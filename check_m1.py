"""Does the stability analysis of Sec. 3.2 (Proposition 1, eq:bound_md) hold for the M1 closure in any
radiation state, or only for the Eddington flux of eq:flux_split?

About a uniform state the HLL speeds of Honghu.jl and QUOKKA, S = -/+ chat max(0.1, sqrt(T_dd)) from the two sides,
coincide, so the linearized HLL flux is the Rusanov flux with the full M1 Jacobian A_d, and the checkerboard
(theta = pi, sin theta = 0) does not see A_d at all. What uses the Eddington Jacobian is Proposition B.1 (the
checkerboard is the least stable mode), whose proof needs spectral radius rho(A_d) = sigma_d. This script
  1. compares rho(A_d) of the Levermore M1 flux with the HLL speed s_d = max(0.1, sqrt(T_dd)) over states;
  2. scans every Fourier mode of the PD-ARS step in an optically thin cell (transport only,
     G = I/2 + (I + Z)^2/2, Z = -nu sum_d [i sin(theta_d) A_d + s_d (1 - cos theta_d) I]) at
       (a) nu = 1/sqrt(D), the standard RSLA at the largest step of Proposition 1, and
       (b) nu = 1/sum_d s_d, the bound eq:bound_md in a thin cell (the ARSA at beta = 1),
     for D = 1, 2, 3 and radiation states up to |F|/(cE) = 0.999: random directions (fixed seed), every 2 deg
     in the x-y plane in 1D, and directions 30-89 deg from x in the x-y and x-z planes in 2D and 3D, where
     rho(A_x)/s_x is largest (part 1).
Writes data/check_m1.txt (the printout). About 3 minutes.

    python check_m1.py
"""
import itertools
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / 'data' / 'check_m1.txt'
rng = np.random.default_rng(1)
lines = []


def say(s):
    print(s)
    lines.append(s)


def chi(f):                                      # Levermore (1984) M1 Eddington factor
    return (3 + 4 * f**2) / (5 + 2 * np.sqrt(4 - 3 * f**2))


def edd(U):
    E, F = U[0], U[1:]
    f = np.linalg.norm(F) / E
    n = F / np.linalg.norm(F)
    c = chi(f)
    return (1 - c) / 2 * np.eye(3) + (3 * c - 1) / 2 * np.outer(n, n)


def flux(U, d):                                  # M1 flux along d, c = 1: (F_d, E T_d.)
    return np.concatenate(([U[1 + d]], U[0] * edd(U)[d]))


def jac(U, d, h=1e-7):
    J = np.zeros((4, 4))
    for k in range(4):
        e = np.zeros(4)
        e[k] = h
        J[:, k] = (flux(U + e, d) - flux(U - e, d)) / (2 * h)
    return J


def state(f, n):
    return np.concatenate(([1.0], max(f, 1e-9) * n))


def speeds(U, D):
    T = edd(U)
    return [max(0.1, np.sqrt(T[d, d])) for d in range(D)]


# 1. spectral radius of the M1 Jacobian against the HLL speed, along x
worst = (0, None)
for f in np.linspace(0, 0.999, 112):
    for ang in np.linspace(0, 90, 46):
        n = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang)), 0])
        U = state(f, n)
        r = np.max(np.abs(np.linalg.eigvals(jac(U, 0)))) / speeds(U, 1)[0]
        if r > worst[0]:
            worst = (r, (f, ang))
say(f'1. max rho(A_x)/max(0.1, sqrt(T_xx)) = {worst[0]:.3f} at f = {worst[1][0]:.3f}, '
    f'angle(F, x) = {worst[1][1]:.0f} deg (1 for the Eddington flux and for M1 at f = 0)')


# 2. every mode in an optically thin cell
def maxrho(U, D, nu, thetas):
    A = [jac(U, d) for d in range(D)]
    s = speeds(U, D)
    w = (0.0, None)
    for t in thetas:
        Z = -nu * sum(1j * np.sin(t[d]) * A[d] + s[d] * (1 - np.cos(t[d])) * np.eye(4) for d in range(D))
        IZ = np.eye(4) + Z
        r = np.max(np.abs(np.linalg.eigvals(0.5 * np.eye(4) + 0.5 * IZ @ IZ)))
        if r > w[0]:
            w = (r, t)
    return w


for D, nth in ((1, 400), (2, 60), (3, 21)):
    grid = list(itertools.product(np.linspace(0, np.pi, nth)[1:], repeat=D))
    rays = [r * k / np.linalg.norm(k) for r in np.geomspace(1e-3, 0.5, 12)
            for k in rng.normal(size=(8, D))]            # long waves in random directions
    thetas = grid + [tuple(t) for t in rays]
    dirs = [np.eye(3)[0], np.ones(3) / np.sqrt(3)] + [v / np.linalg.norm(v) for v in rng.normal(size=(10, 3))]
    angs = np.radians(np.arange(0, 91, 2) if D == 1 else (30, 60, 75, 80, 86, 89))
    dirs += [np.array([np.cos(a), np.sin(a), 0]) for a in angs]
    if D > 1:
        dirs += [np.array([np.cos(a), 0, np.sin(a)]) for a in angs]
    for label, nu_of in (('RSLA, nu = 1/sqrt(D)', lambda s: 1 / np.sqrt(D)),
                         ('bound, nu = 1/sum_d s_d', lambda s: 1 / sum(s))):
        w = (0.0, None, None, None)
        for f in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99, 0.999):
            for n in dirs:
                U = state(f, n)
                nu = nu_of(speeds(U, D))
                r, t = maxrho(U, D, nu, thetas)
                if r > w[0]:
                    w = (r, f, n.round(2), t)
        r, f, n, t = w
        say(f'2. D = {D}, {label}: max rho - 1 = {r - 1:+.1e} (at f = {f}, n = {n}, theta = {np.round(t, 3)})')

OUT.write_text('\n'.join(lines) + '\n')
print('wrote', OUT)
