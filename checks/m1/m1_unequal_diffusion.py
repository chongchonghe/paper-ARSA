"""Unequal-opacity diffusion theorem: exact source identities and independent checks.

alpha_F is the opacity scale; opacity ratios and material response ratios
are fixed when this scale increases. All velocity couplings at v0=0 are kept.
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp

from m1_stability import matrices, closure
from m1_exchange_bounds import pdars
from m1_unequal_opacities import full_source


def blocks(u, p, ratio, gamma, ue, um):
    C = np.zeros((4, 4))
    C[0, 0] = p
    C[0, 1:] = (2-ratio)*u[1:]/u[0]
    C[1:, 1:] = gamma*np.eye(3)+closure(u)
    rates = np.diag([1/ratio, 1., 1., 1.])
    M = np.diag([ue, um, um, um])
    Q = np.block([[M@rates@C, -M@rates], [-rates@C, rates]])
    return C, rates, M, Q


def constants(p, ratio, gamma, ue, um):
    v = 1+um
    e = 1+ue+ue*abs(2-ratio)*v
    er = 1+p*e+abs(2-ratio)*v
    jr = 1+(gamma+1)*v
    Kphys = max(v, e, er, jr)
    cQphys = (max(ue/ratio, um)+max(1/ratio, 1))*(abs(2-ratio)+2)
    return Kphys, cQphys


def symbolic():
    p, ratio, eta, ue, um, q = sp.symbols('p ratio eta ue um q', positive=True)
    f, chi = sp.symbols('f chi', real=True)
    gamma = 1+(p*eta-1)/ratio
    C = sp.Matrix([[p, (2-ratio)*f], [0, gamma+chi]])
    R, M = sp.diag(1/ratio, 1), sp.diag(ue, um)
    Q = (M*R*C).row_join(-M*R).col_join((-R*C).row_join(R))
    B, v, E, J = sp.symbols('B v E J', real=True)
    forceE = E/ratio-p*B/ratio+(1-2/ratio)*v*J
    forceJ = J-(p*B/ratio+(1-1/ratio)*E+chi*E)*v
    rhs = sp.Matrix([ue*forceE, um*forceJ, -forceE, -forceJ])
    jac = rhs.jacobian([B, v, E, J]).subs({B: eta, E: 1, J: f, v: 0})
    assert (jac+Q).applyfunc(sp.simplify) == sp.zeros(4)
    # Elimination with distinct radiation energy and flux relaxation times.
    Tq = q*R*(sp.eye(2)+q*R).inv()
    A = sp.eye(2)+M*Tq*C
    amap = A.inv()
    bmap = amap*M*Tq
    resolvent = amap.row_join(bmap).col_join(
        (Tq*C*amap).row_join(sp.eye(2)-Tq+Tq*C*bmap))
    assert ((sp.eye(4)+q*Q)*resolvent-sp.eye(4)).applyfunc(sp.simplify) == sp.zeros(4)
    # Fixed isotropic coordinates reduce to the original channel problem.
    C0 = sp.diag(p, gamma+sp.Rational(1, 3))
    Q0 = Q.subs({f: 0, chi: sp.Rational(1, 3)})
    S = sp.diag(C0, sp.eye(2))
    Xi = M*C0
    channelQ = (R*Xi).row_join(-R*Xi).col_join((-R).row_join(R))
    assert (S*Q0-channelQ*S).applyfunc(sp.simplify) == sp.zeros(4)
    return dict(full_nonequilibrium_four_force=True,
                two_relaxation_time_elimination=True,
                fixed_isotropic_channel_identity=True,
                gamma='1+(p*eta-1)/ratio',
                adequate_drag_condition='gamma >= 0 (sufficient, not necessary)')


def numerical():
    rng = np.random.default_rng(20261025)
    stats = dict(states=0, derivative_checks=0, reference_mode_checks=0,
                 diffusion_mode_checks=0, max_derivative_error=0.,
                 max_reference_energy_norm=0., max_resolvent_over_bound=0.,
                 max_perturbation_over_bound=0., max_identity_residual=0.,
                 max_generator_scaled_difference=0.)
    examples = []
    for j in range(30):
        p, ratio = 10**rng.uniform(-1, 1, 2)
        if j < 3:
            ratio = [.25, 2., 5.][j]  # both signs and zero of energy-velocity coupling
        gamma = 1. if j % 3 == 0 else max(0., 1-1/ratio)+10**rng.uniform(-1, .5)
        eta = (1+ratio*(gamma-1))/p
        ue, um = 10**rng.uniform(-3, -.1), 10**rng.uniform(-5, -1)
        direction = rng.normal(size=3); direction /= np.linalg.norm(direction)
        _, A0, _, _ = matrices(0., direction)
        C0 = np.diag([p]+[gamma+1/3]*3)
        M = np.diag([ue]+[um]*3)
        Hr = np.diag([1., 3., 3., 3.])
        H0 = np.block([[C0@Hr@np.linalg.inv(M), np.zeros((4, 4))],
                       [np.zeros((4, 4)), Hr]])
        chol = np.linalg.cholesky(H0).T; ci = np.linalg.inv(chol)
        # Fixed equivalence between physical block-max norm and this Euclidean norm.
        # max(|energy|,||velocity||,|E|,||J||) <= ||x||_2 <= 2 max(...).
        norm_transfer = 2*np.linalg.cond(chol)
        Kphys, cQphys = constants(p, ratio, gamma, ue, um)
        K, cQ, cW = norm_transfer*Kphys, norm_transfer*cQphys, norm_transfer*18
        xiE, xiF = p*ue, (gamma+1/3)*um
        Lstar = 1+1/max(xiE, xiF)
        d0 = 1+2*Lstar
        dactual = d0+cW*np.sqrt(3)*Lstar
        c1 = cQ*(K/4*(1+d0*d0)+K*K/2*d0*dactual)
        c2 = cW*K/2*(K*dactual+d0)
        for alpha in [1., 100., 1e4, 1e6]:
            # Fixed B_f=0.2, tau_*=10 and f<=0.2; increasing stiffness drives f down.
            f = .2/alpha
            beta = 3*f/(2+np.sqrt(4-3*f*f))
            u, A, _, speeds = matrices(beta, direction)
            C, rates, M, Qs = blocks(u, p, ratio, gamma, ue, um)
            u0 = np.r_[1., np.zeros(3)]
            _, _, _, Qs0 = blocks(u0, p, ratio, gamma, ue, um)
            Qh, Q0h = chol@Qs@ci, chol@Qs0@ci
            # Independent nonlinear four-force derivative, outside balance as well.
            state = np.r_[eta, np.zeros(3), u]
            fd = np.empty((8, 8))
            for k, ej in enumerate(np.eye(8)):
                eps = 1e-6*max(1., abs(state[k]))
                plus = full_source(state+eps*ej, p/ratio, 1/ratio, 1., ue/ratio, um)
                minus = full_source(state-eps*ej, p/ratio, 1/ratio, 1., ue/ratio, um)
                fd[:, k] = (plus-minus)/(2*eps)
            err = np.linalg.norm(fd+Qs, np.inf)/(1+np.linalg.norm(Qs, np.inf))
            assert err < 2e-8
            stats['max_derivative_error'] = max(stats['max_derivative_error'], float(err))
            stats['derivative_checks'] += 1
            assert np.linalg.norm(Qh-Q0h, 2) <= cQ*f+1e-12
            for dim in [1, 2, 3]:
                dx = 10/alpha
                # Solve the implicit two-channel timestep inequality in nu.
                def allowed(nu):
                    qF = 10*nu; qE = qF/ratio
                    return min(1+qE/(2+xiE*qE), 1+qF/(2+xiF*qF))
                lo, hi = 0., np.sqrt(3)*Lstar/dim
                for _ in range(70):
                    mid = (lo+hi)/2
                    if dim*mid/np.sqrt(3) <= allowed(mid): lo = mid
                    else: hi = mid
                nu = .95*lo; dt = nu*dx; qF = alpha*dt
                Q, Q0 = qF*Qh, qF*Q0h
                B0 = 2*np.eye(8)+Q0
                chB = np.linalg.cholesky(B0).T
                chBi = np.linalg.inv(chB)
                for theta in [np.zeros(dim), np.full(dim, np.pi), rng.uniform(-np.pi, np.pi, dim)]:
                    Wr = nu*(np.dot(speeds[:dim], 1-np.cos(theta))*np.eye(4)
                              + 1j*np.einsum('d,dij->ij', np.sin(theta), A[:dim]))
                    Wr0 = nu*(np.sum(1-np.cos(theta))/np.sqrt(3)*np.eye(4)
                               + 1j*np.einsum('d,dij->ij', np.sin(theta), A0[:dim]))
                    W = np.zeros((8, 8), complex); W[4:, 4:] = Wr
                    W0 = np.zeros((8, 8), complex); W0[4:, 4:] = Wr0
                    W, W0 = chol@W@ci, chol@W0@ci
                    D, D0 = np.eye(8)-W, np.eye(8)-W0
                    R, R0 = np.linalg.inv(np.eye(8)+Q), np.linalg.inv(np.eye(8)+Q0)
                    H, Href = np.linalg.inv(2*np.eye(8)+Q), np.linalg.inv(2*np.eye(8)+Q0)
                    G, G0 = pdars(Q, W), pdars(Q0, W0)
                    en = np.linalg.norm(chB@G0@chBi, 2)
                    assert en <= 1+1e-8
                    stats['max_reference_energy_norm'] = max(stats['max_reference_energy_norm'], float(en))
                    ratioK = max(np.linalg.norm(R, 2), 2*np.linalg.norm(H, 2))/K
                    assert ratioK <= 1+1e-10
                    stats['max_resolvent_over_bound'] = max(stats['max_resolvent_over_bound'], float(ratioK))
                    difference = (H-Href)@(np.eye(8)+D0@R0@D0)+H@(
                        (D-D0)@R@D+D0@(R-R0)@D+D0@R0@(D-D0))
                    residual = np.linalg.norm(G-G0-difference, 2)/(1+np.linalg.norm(G, 2))
                    assert residual < 1e-10
                    stats['max_identity_residual'] = max(stats['max_identity_residual'], float(residual))
                    bound = c1*qF*f+c2*f*dim*nu
                    diffnorm = np.linalg.norm(G-G0, 2)
                    assert diffnorm <= bound*(1+1e-10)+1e-12
                    assert bound <= dt*.2*(c1+c2*dim/10)*(1+1e-10)
                    stats['max_perturbation_over_bound'] = max(stats['max_perturbation_over_bound'], float(diffnorm/bound))
                    stats['max_generator_scaled_difference'] = max(stats['max_generator_scaled_difference'], float(diffnorm/dt))
                    stats['reference_mode_checks'] += 1
                    stats['diffusion_mode_checks'] += 1
        stats['states'] += 1
        if j < 3:
            examples.append(dict(p=float(p), opacity_ratio=float(ratio), eta=float(eta),
                                 gamma=float(gamma), limiting_CFL=float(Lstar)))
    return dict(statistics=stats, examples=examples,
                scope='Frozen ratios; gamma>=0. Numerical checks supplement the analytic proof.')


if __name__ == '__main__':
    result = dict(symbolic=symbolic(), numerical=numerical())
    path = Path(__file__).resolve().parents[2]/'data/m1/unequal_diffusion.json'
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
