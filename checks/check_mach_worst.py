"""Moving background (M != 0, 1D, chi = 1/3, C = 1e4): how far is the closed
form nu_max (eq:numax) from the true boundary nu_true of the full 4x4 system
(max over theta), and from the checkerboard-only boundary nu_pi?

  dev_true = nu_true / nu_max - 1   (negative = closed form is unsafe)
  dev_pi   = nu_pi   / nu_max - 1
  cb_least = nu_true / nu_pi - 1    (negative = some theta < pi is less stable than pi)

usage: python moving_bg.py mode   with mode in {scan, de}
"""
import sys
import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'src'))
from pdars import matrices

C = 1e4
TOL = 1e-9
THETA = np.pi * np.concatenate([np.geomspace(1e-4, 4e-3, 10), np.linspace(5e-3, 1.0, 240)])
I4 = np.eye(4)


def growth(thetas, tau, nu, nuc, P, M, gamma, chi=1 / 3):
    phi, dt = nu / nuc, nuc / C
    A, S = matrices(M, phi, C, P, gamma, chi)
    Q = (dt * tau * S)[1:, 1:]
    sig = np.sqrt(chi) * phi * C
    Ar = A[3:, 3:]
    th = np.atleast_1d(thetas)[:, None, None]
    E = np.broadcast_to(I4.astype(complex), (th.shape[0], 4, 4)).copy()
    E[:, 2:, 2:] -= dt * (1j * np.sin(th) * Ar + (1 - np.cos(th)) * sig * np.eye(2))
    G1 = np.linalg.solve(I4 - Q, E)
    G = np.linalg.solve(I4 - 0.5 * Q, 0.5 * I4 + 0.5 * E @ G1)
    return np.abs(np.linalg.eigvals(G)).max(axis=-1)


def nu_bound(thetas, tau, nuc, P, M, gamma, rel=1e-7):
    ok = lambda nu: growth(thetas, tau, nu, nuc, P, M, gamma).max() <= 1 + TOL
    if ok(nuc):
        return nuc
    lo, hi = 0.0, nuc
    while hi - lo > rel * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def nu_formula(tau, nuc, P, gamma, chi=1 / 3):
    den = np.sqrt(chi) - tau / (2 + 4 * (gamma - 1) * P * nuc * tau)
    return np.inf if den <= 0 else 1 / den


def case(tau, nuc, P, M, gamma):
    nf = min(nu_formula(tau, nuc, P, gamma), nuc)
    npi = nu_bound(np.array([np.pi]), tau, nuc, P, M, gamma)
    ntrue = nu_bound(THETA, tau, nuc, P, M, gamma)
    return nf, npi, ntrue


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'scan'
    if mode == 'scan':
        rng = np.random.default_rng(0)
        rows = []
        for n in range(400):
            tau = 10 ** rng.uniform(-3, 6); nuc = 10 ** rng.uniform(1, 3.5); P = 10 ** rng.uniform(-6, 1)
            gamma = rng.uniform(1.1, 5 / 3); M = rng.uniform(0, 3)
            nf, npi, ntrue = case(tau, nuc, P, M, gamma)
            if nf >= nuc * (1 - 1e-6):
                continue           # boundary not active (phi = 1 stable)
            rows.append((tau, nuc, P, M, gamma, nf, npi, ntrue))
        rows = np.array(rows)
        dev_true = rows[:, 7] / rows[:, 5] - 1
        dev_pi = rows[:, 6] / rows[:, 5] - 1
        cb = rows[:, 7] / rows[:, 6] - 1
        print(f'{len(rows)} active cases, M in [0, 3]:')
        print(f'   nu_true/nu_max - 1: min {dev_true.min():+.3e} (unsafe side), max {dev_true.max():+.3e}')
        print(f'   nu_pi/nu_max   - 1: min {dev_pi.min():+.3e}, max {dev_pi.max():+.3e}')
        print(f'   nu_true/nu_pi  - 1: min {cb.min():+.3e} (checkerboard not least stable if < -1e-6), max {cb.max():+.3e}')
        k = np.argmin(dev_true)
        print('   most unsafe case (tau, nuc, P, M, gamma, nu_max, nu_pi, nu_true):', rows[k])
        for Mlo, Mhi in ((0, 0.5), (0.5, 1), (1, 2), (2, 3)):
            sel = (rows[:, 3] >= Mlo) & (rows[:, 3] < Mhi)
            if sel.any():
                print(f'   M in [{Mlo},{Mhi}): n={sel.sum()}  min dev_true {dev_true[sel].min():+.2e}  max |dev_true| {np.abs(dev_true[sel]).max():.2e}'
                      f'  min cb {cb[sel].min():+.2e}')
    else:
        Mmax = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0

        def obj(x):
            tau, nuc, P = 10 ** x[0], 10 ** x[1], 10 ** x[2]; gamma, M = x[3], x[4]
            nf = min(nu_formula(tau, nuc, P, gamma), nuc)
            if nf >= nuc * (1 - 1e-6):
                return 0.0
            ntrue = nu_bound(THETA, tau, nuc, P, M, gamma, rel=1e-6)
            return ntrue / nf - 1          # minimize -> most unsafe

        bounds = [(-3, 6), (1, 3.5), (-6, 1), (1.1, 5 / 3), (0, Mmax)]
        res = differential_evolution(obj, bounds, seed=1, maxiter=25, popsize=12, tol=1e-12, polish=False)
        tau, nuc, P = 10 ** res.x[0], 10 ** res.x[1], 10 ** res.x[2]; gamma, M = res.x[3], res.x[4]
        nf, npi, ntrue = case(tau, nuc, P, M, gamma)
        print(f'DE most-unsafe search, M <= {Mmax}: tau={tau:.4g} nuc={nuc:.4g} P={P:.4g} gamma={gamma:.4f} M={M:.4f}')
        print(f'   nu_max={nf:.8g} nu_pi={npi:.8g} nu_true={ntrue:.8g}  nu_true/nu_max-1={ntrue/nf-1:+.3e}  nu_true/nu_pi-1={ntrue/npi-1:+.3e}')
        th = THETA
        g = growth(th, tau, ntrue * (1 + 1e-3), nuc, P, M, gamma)
        print(f'   at nu = 1.001 nu_true the fastest-growing mode is theta/pi = {th[np.argmax(g)]/np.pi:.4f}')
