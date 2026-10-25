"""Exact counterexample to dropping K_F in the manuscript's ARSA limiter.

Isotropic, stationary, energy-balanced ideal gas; kappa_P = kappa_E != kappa_F.
No nonlinear solver is run. The repeated update is the frozen linear PD-ARS
scheme on a periodic 2x2x2 grid, initialized with its checkerboard eigenmode.
Run: python3 checks/m1/p1_dropped_momentum_counterexample.py
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

R = sp.Rational
beta, phi0 = R(9, 10), R(1, 10000)
C, gamma = sp.Integer(100), R(5, 3)  # c / isothermal sound speed
h = R(1, 10000)
L = sp.Integer(9000)  # sum_d nu_d / sqrt(3), evaluated at phi = 1
qE, qF = sp.Integer(200000), R(89991, 5)
lambdaF = h / (3 * (gamma - 1) * C**2)
KE, KF = h*qE, lambdaF*qF
energy = L-qE/(2+KE)
approx_momentum = L-qF/2
exact_momentum = L-qF/(2+KF)
phi_approx = min(1, max(phi0, beta/max(energy, approx_momentum)))
phi_exact = min(1, max(phi0, beta/max(energy, exact_momentum)))
assert phi_approx == 1 and approx_momentum == beta
assert energy < 0 and exact_momentum > 1
assert 0 < KF < R(1, 10000)
assert phi0*L == beta  # ordinary 3D RSLA CFL with a 0.9 safety factor


def amplification(phi, q, K):
    # Rescale material momentum so the exchange-equilibrium vector is (1,1).
    Q = sp.Matrix([[K, -K], [-phi*q, phi*q]])
    D = sp.diag(1, 1-2*phi*L)
    return (2*sp.eye(2)+Q).inv()*(sp.eye(2)+D*(sp.eye(2)+Q).inv()*D)


G = amplification(phi_approx, qF, KF)
Ge = amplification(phi_approx, qE, KE)
Gfixed = amplification(phi_exact, qF, KF)
# det(I-G)<0 gives a real eigenvalue >1, since the characteristic polynomial
# is monic and tends to +infinity on the positive real axis.
certificate = sp.factor((sp.eye(2)-G).det())
assert certificate < 0
assert L <= 1+qE/(2+KE)
assert phi_exact*L < 1+phi_exact*qF/(2+KF)


def mpmatrix(A):
    return mp.matrix([[mp.mpf(str(x.p))/mp.mpf(str(x.q)) for x in row]
                      for row in A.tolist()])


def dominant(A):
    t, d = A[0, 0]+A[1, 1], mp.det(A)
    return (t+mp.sqrt(t*t-4*d))/2


estimates = []
for precision in [60, 90]:
    mp.mp.dps = precision
    eigenvalue = dominant(mpmatrix(G))
    fixed_eigenvalue = dominant(mpmatrix(Gfixed))
    energy_eigenvalue = dominant(mpmatrix(Ge))
    assert eigenvalue > 1 and fixed_eigenvalue < 1 and energy_eigenvalue <= 1
    estimates.append(dict(digits=precision, momentum_gain=mp.nstr(eigenvalue, 40),
                          corrected_gain=mp.nstr(fixed_eigenvalue, 40),
                          energy_gain=mp.nstr(energy_eigenvalue, 40)))

# Execute the two stages directly, with HLL diffusion on a periodic 2^3 grid.
# At isotropy the centered physical-flux difference vanishes on this mode.
mp.mp.dps = 90
A = mpmatrix(G)
g = dominant(A)
v = np.array([float(A[0, 1]/(g-A[0, 0])), 1.])
signs = np.array([(-1.)**(bin(j).count('1')) for j in range(8)])
initial = signs[:, None]*v
steps = 10000


def grid_run(phi):
    phi = float(phi)
    Q = np.array([[float(KF), -float(KF)], [-phi*float(qF), phi*float(qF)]])
    inv1 = np.linalg.inv(np.eye(2)+Q)
    inv2 = np.linalg.inv(2*np.eye(2)+Q)
    weight = phi*float(L)/3  # nu_d sigma_d in each of the 3 directions
    def transport(u):
        z = u.copy()
        for bit in [1, 2, 4]:
            # Both periodic neighbors coincide on a two-cell axis.
            z[:, 1] += weight*(u[np.arange(8)^bit, 1]-u[:, 1])
        return z
    u = initial.copy()
    for _ in range(steps):
        stage1 = transport(u)@inv1.T
        u = (u + transport(stage1))@inv2.T
    amplitude = (signs@u)/8
    gain = float(np.linalg.norm(amplitude)/np.linalg.norm(v))
    # Independently compare with the exact rational amplification matrix.
    target = mpmatrix(amplification(R(str(phi)), qF, KF))**steps*mp.matrix(v.tolist())
    target = np.array([float(x) for x in target])
    relative_error = float(np.linalg.norm(amplitude-target)/np.linalg.norm(target))
    assert relative_error < 2e-7
    return dict(steps=steps, amplitude_gain=gain, matrix_power_relative_error=relative_error)


unstable_run = grid_run(phi_approx)
corrected_run = grid_run(phi_exact)
assert unstable_run['amplitude_gain'] > 1.5
assert corrected_run['amplitude_gain'] < 1
result = dict(
    model='Frozen isotropic equilibrium P1, also the full-M1 linearization at f=0',
    parameters={k:str(v) for k,v in dict(beta=beta, phi0=phi0, c_over_sound_speed=C,
        gamma=gamma, heat_capacity_ratio=h, L_at_phi_one=L, qE_at_phi_one=qE,
        qF_at_phi_one=qF, KE=KE, KF=KF, lambdaF_at_phi_one=lambdaF,
        kappa_P_over_kappa_E=1, kappa_E_over_kappa_F=qE/qF).items()},
    limiter=dict(energy_residual=float(energy), approximate_momentum_residual=float(approx_momentum),
        exact_momentum_residual=float(exact_momentum), approximate_phi=float(phi_approx),
        corrected_phi=float(phi_exact)),
    exact_certificate=dict(det_I_minus_G=str(certificate), negative=True),
    high_precision=estimates, grid_run=unstable_run, corrected_grid_run=corrected_run,
    physical_scope=dict(velocity_background=0, reduced_flux=0, energy_balance=True,
        nominal_RSLA_CFL='nu0=0.9/sqrt(3)',
        flux_cell_optical_depth=float(sp.sqrt(3)*qF/L),
        energy_cell_optical_depth=float(sp.sqrt(3)*qE/L),
        gas_CFL_per_direction=float(L/(sp.sqrt(3)*C)),
        limitation='Radiation-only large step; exceeds ordinary explicit gas CFL. This is not a counterexample to a hydrodynamically CFL-limited full run.'))
path = Path(__file__).resolve().parents[2]/'data/m1/dropped_momentum_counterexample.json'
path.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
