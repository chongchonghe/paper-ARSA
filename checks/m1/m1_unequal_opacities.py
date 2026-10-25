"""Verify the unequal-mean-opacity, static-background M1 extension.

All opacity coefficients and density are frozen. The background satisfies
kappa_P B0 = kappa_E E0, but its radiation flux need not vanish. The
source includes every velocity derivative of the paper's O(v/c) four-force
at v0=0. Gas hydrodynamic transport is outside this radiation substep.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp

from m1_stability import matrices, closure
from m1_exchange_bounds import pdars


def source_blocks(u, kP, kE, kF, material_energy, material_momentum):
    C = np.zeros((4, 4))
    C[0, 0] = kP/kE
    C[0, 1:] = (2-kF/kE)*u[1:]/u[0]
    C[1:, 1:] = np.eye(3)+closure(u)
    R = np.diag([kE, kF, kF, kF])  # Units: c_hat*rho=1.
    D = np.diag([material_energy, material_momentum, material_momentum, material_momentum])
    Q = np.block([[D@C, -D], [-R@C, R]])
    transform = np.block([[np.eye(4), D@np.linalg.inv(R)], [-C, np.eye(4)]])
    return C, R, D, Q, transform


def full_source(state, kP, kE, kF, material_energy, material_momentum):
    """Independent nonlinear four-force: state=(B, v/c, E, F/c), c=1.

    The equation is written for gas total-energy perturbations scaled to
    emission perturbations at the static background. Constant balancing
    terms, if used, cancel from its derivative.
    """
    emission, velocity = state[0], state[1:4]
    erad, j = state[4], state[5:8]
    pressure = erad*closure(np.r_[erad, j])
    force_energy = kE*erad-kP*emission+(kF-2*kE)*np.dot(velocity, j)
    force_momentum = (kF*j-(kP*emission+(kF-kE)*erad)*velocity
                      - kF*pressure@velocity)
    return np.r_[material_energy/kE*force_energy,
                 material_momentum/kF*force_momentum,
                 -force_energy, -force_momentum]


def symbolic_check():
    p, e, ff, de, dm = sp.symbols('p e ff de dm', positive=True)
    f, chi = sp.symbols('f chi', real=True)
    C = sp.Matrix([[p/e, (2-ff/e)*f], [0, 1+chi]])
    R, D = sp.diag(e, ff), sp.diag(de, dm)
    Q = (D*C).row_join(-D).col_join((-R*C).row_join(R))
    change = sp.eye(2).row_join(D*R.inv()).col_join((-C).row_join(sp.eye(2)))
    block = sp.zeros(4)
    block[2:, 2:] = R+C*D
    assert (change*Q-block*change).applyfunc(sp.simplify) == sp.zeros(4)
    expected_det = (e+p/e*de)*(ff+(1+chi)*dm)/(e*ff)
    assert sp.simplify(change.det()-expected_det) == 0
    # Differentiate all longitudinal velocity terms independently.
    B, v, E, J = sp.symbols('B v E J', real=True)
    energy = e*E-p*B+(ff-2*e)*v*J
    momentum = ff*J-(p*B+(ff-e)*E+ff*chi*E)*v
    source = sp.Matrix([de/e*energy, dm/ff*momentum, -energy, -momentum])
    jac = source.jacobian([B, v, E, J]).subs({B: e/p, E: 1, J: f, v: 0})
    assert (jac+Q).applyfunc(sp.simplify) == sp.zeros(4)
    return {'source_transform_identity': True,
            'transform_determinant': str(sp.factor(expected_det)),
            'full_static_velocity_derivatives': True}


def checks():
    rng = np.random.default_rng(2404)
    stats = {'states': 0, 'four_force_derivative_checks': 0, 'source_resolvent_checks': 0,
             'amplification_checks': 0, 'maximum_relative_derivative_error': 0.,
             'maximum_transform_residual': 0., 'maximum_resolvent_norm_over_bound': 0.,
             'energy_velocity_coupling_signs': set()}
    ratios = [(.3, 2., .1), (4., .2, 3.), (.5, 1., 2.), (1., 3., 3.), (3., .1, .3)]
    ratios += [tuple(10**rng.uniform(-1.5, 1.5, 3)) for _ in range(10)]
    for beta in [0., .6, .95]:
        for kP, kE, kF in ratios:
            u, a, h, speeds = matrices(beta, rng.normal(size=3))
            de, dm = 10**rng.uniform(-2, 1), 10**rng.uniform(-4, 0)
            C, R, D, Q, transform = source_blocks(u, kP, kE, kF, de, dm)
            B = R+C@D
            expected = np.zeros((8, 8))
            expected[4:, 4:] = B
            residual = float(np.linalg.norm(transform@Q-expected@transform)/(1+np.linalg.norm(transform@Q)))
            assert residual < 1e-12
            stats['maximum_transform_residual'] = max(stats['maximum_transform_residual'], residual)
            # E0=1, and B0 is fixed by radiative energy balance.
            state = np.r_[kE/kP, np.zeros(3), u]
            fd = np.zeros((8, 8))
            for j, direction in enumerate(np.eye(8)):
                eps = 1e-6*max(1., abs(state[j]))
                fd[:, j] = (full_source(state+eps*direction, kP, kE, kF, de, dm)
                            - full_source(state-eps*direction, kP, kE, kF, de, dm))/(2*eps)
            error = float(np.max(abs(fd+Q))/(1+np.max(abs(Q))))
            assert error < 1e-8
            stats['maximum_relative_derivative_error'] = max(stats['maximum_relative_derivative_error'], error)
            stats['four_force_derivative_checks'] += 1
            stats['energy_velocity_coupling_signs'].add(int(np.sign(2-kF/kE)))
            H = np.zeros((8, 8))
            H[:4, :4] = np.eye(4)
            H[4:, 4:] = h
            chol = np.linalg.cholesky(H).T
            ci = np.linalg.inv(chol)
            Qh = chol@Q@ci
            f = np.linalg.norm(u[1:]/u[0])
            chi = (1+3*beta**2)/(3+beta**2)
            aa = kE+kP/kE*de
            bb = kF+(1+chi)*dm
            offdiag = abs(2-kF/kE)*f*dm
            KB = 1+offdiag/(np.sqrt(aa)+np.sqrt(bb))**2
            K = float(np.linalg.cond(transform@ci))*KB
            M = float(np.linalg.svd(Qh, compute_uv=False)[0])
            Cstar = K*M*(1+3*K)/2
            ell = min(kE, kF)/2
            Atime = 3*np.sqrt(3/2)*K*K*ell
            Btime = 9/4*K*K*ell*ell
            for s in np.logspace(-4, 4, 9):
                resolvent = np.linalg.inv(np.eye(8)+s*Qh)
                ratio = float(np.linalg.svd(resolvent, compute_uv=False)[0])/K
                assert ratio <= 1+1e-8
                stats['maximum_resolvent_norm_over_bound'] = max(stats['maximum_resolvent_norm_over_bound'], ratio)
                stats['source_resolvent_checks'] += 1
            for dt in [.01, .1, 1., 10.]:
                for dim in [1, 2, 3]:
                    theta = rng.uniform(-np.pi, np.pi, dim)
                    nu0 = 1/sum(speeds[:dim])
                    wr = nu0*(np.dot(speeds[:dim], 1-np.cos(theta))*np.eye(4)
                              + 1j*np.einsum('d,dij->ij', np.sin(theta), a[:dim]))
                    Wbar = np.zeros((8, 8), dtype=complex)
                    Wbar[4:, 4:] = wr
                    Wbar = chol@Wbar@ci
                    Gbar = pdars(dt*Qh, Wbar)
                    G = pdars(dt*Qh, (1+ell*dt)*Wbar)
                    assert np.linalg.svd(Gbar, compute_uv=False)[0] <= 1+Cstar*dt+1e-8
                    assert np.linalg.svd(G-Gbar, compute_uv=False)[0] <= Atime*dt+Btime*dt*dt+1e-8
                    assert np.log(np.linalg.svd(G, compute_uv=False)[0]) <= (Cstar+Atime)*dt+1e-8
                    # Both original channel conditions imply the min(q_E,q_F) envelope.
                    allowed = min(1+kE*dt/(2+kP/kE*de*dt),
                                  1+kF*dt/(2+(1+chi)*dm*dt))
                    assert allowed <= 1+ell*dt+1e-12
                    stats['amplification_checks'] += 1
            stats['states'] += 1
    stats['energy_velocity_coupling_signs'] = sorted(stats['energy_velocity_coupling_signs'])
    return stats


if __name__ == '__main__':
    result = {'symbolic': symbolic_check(), 'numerical': checks()}
    dest = Path(__file__).resolve().parents[2]/'data/m1/unequal_opacities.json'
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
