"""Is the checkerboard bound, with the actual M1 HLL speeds, stable in an absorbing cell with anisotropic radiation?

The model of the absorption counterexample of m1_stability.py: gray absorption into a fixed reservoir (the gas holds
the heat capacity, K = 0), the full M1 flux Jacobian A_d of the Levermore closure about a uniform state with reduced
flux f and direction n, HLL diffusion with the actual speeds sigma_d = max(0.1, sqrt(T_dd)), and the PD-ARS update
G = (I + D^2/(1 + q))/(2 + q), D = I - W(theta), on cubic cells. The limiter sets the step from the checkerboard bound,
nu = beta / (S - tau/2) with beta = 0.9, tau the cell optical depth and q = nu tau; where S - tau/2 <= 0 it sets no
bound, and steps up to nu = 1e4 are tested. A case fails if some Fourier mode has |g| > 1 + 1e-9. Three rules for S:
  actual      S = sum_d sigma_d over the active directions (the rule of the paper, eq:phi_sb);
  max_iso     S = sum_d max(sigma_d, 1/sqrt3);
  lambda_max  S = sum_d of the largest |eigenvalue| of A_d.
Three geometries:
  1D, with the flux at any angle to the grid axis; symmetry forbids a transverse flux in a 1D simulation, but the
      counterexamples of m1_stability.py and m1_thick_cells.py are of this kind;
  2D, with the flux in the grid plane, as in a 2D simulation;
  3D, with the flux in any direction (one octant, by symmetry).
Writes data/m1/absorbing_scan.json: the number of cases and, per geometry and rule, every failure with its exact
parameters (b is the closure coordinate, f = 4b/(3 + b^2); n the flux direction). About 5 minutes.

    python3 checks/m1/absorbing_scan.py
"""
import json
from pathlib import Path
import numpy as np
from m1_stability import matrices

BETA, NU_CAP, TOL = 0.9, 1e4, 1e-9
ISO = 1 / np.sqrt(3)
OUT = Path(__file__).resolve().parents[2] / 'data' / 'm1' / 'absorbing_scan.json'


def phases(dim):
    """Fourier phases per direction, the zero mode excluded."""
    t = {1: np.pi * np.r_[np.geomspace(1e-3, 0.02, 15), np.linspace(0.02, 1, 200)],
         2: np.pi * np.r_[0, np.geomspace(0.003, 0.05, 5), np.linspace(0.08, 1, 24)],
         3: np.pi * np.r_[0, 0.01, 0.04, np.linspace(0.1, 1, 12)]}[dim]
    th = np.array(np.meshgrid(*[t] * dim)).reshape(dim, -1)
    return th[:, np.any(th > 0, axis=0)]


def rho_max(A, sig, nu, q, th):
    """Largest |g| over the phases th (one column per mode)."""
    re = sum(sig[d] * (1 - np.cos(th[d])) for d in range(len(sig)))[:, None, None] * np.eye(4)
    im = sum(np.sin(th[d])[:, None, None] * A[d] for d in range(len(sig)))
    D = np.eye(4) - nu * (re + 1j * im)
    G = (np.eye(4) + D @ D / (1 + q)) / (2 + q)
    return np.abs(np.linalg.eigvals(G)).max()


def directions(dim):
    if dim == 1:                                     # cosine with the grid axis
        return [[mu, np.sqrt(1 - mu * mu), 0.0] for mu in np.linspace(0, 1, 9)]
    if dim == 2:                                     # in the grid plane
        return [[np.cos(p), np.sin(p), 0.0] for p in np.radians(np.linspace(0, 45, 7))]
    return [[np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)]
            for t in np.radians([0, 20, 45, 70, 90]) for p in np.radians([0, 22.5, 45])]


def scan(dim, bs, taus):
    th = phases(dim)
    out = {'cases': 0, 'failures': {'actual': [], 'max_iso': [], 'lambda_max': []}}
    for b in bs:                                     # f = 4b/(3 + b^2)
        for n in directions(dim):
            u, a, h, speeds = matrices(b, n)
            A, sig = a[:dim], speeds[:dim]
            lam = [np.abs(np.linalg.eigvals(A[d])).max() for d in range(dim)]
            for tau in taus:
                out['cases'] += 1
                for rule, S in (('actual', sum(sig)), ('max_iso', sum(max(s, ISO) for s in sig)),
                                ('lambda_max', sum(lam))):
                    den = S - tau / 2
                    for nu in ([BETA / den] if den > 0 else list(np.geomspace(1, NU_CAP, 9))):
                        r = rho_max(A, sig, nu, nu * tau, th)
                        if r > 1 + TOL:
                            out['failures'][rule].append({'f': 4 * b / (3 + b * b), 'b': float(b), 'n': [float(x) for x in n],
                                                          'tau': float(tau), 'nu': float(nu), 'g': float(r)})
                            break
    return out


if __name__ == '__main__':
    result = {'1D': scan(1, np.r_[0.0, np.linspace(0.05, 0.95, 10), 0.98], np.geomspace(1e-3, 1e3, 25)),
              '2D': scan(2, np.r_[0.0, np.linspace(0.1, 0.9, 6), 0.95, 0.98], np.geomspace(1e-3, 1e3, 19)),
              '3D': scan(3, np.r_[0.0, np.linspace(0.2, 0.9, 4), 0.95, 0.98], np.geomspace(1e-3, 1e3, 13))}
    OUT.write_text(json.dumps(result, indent=1) + '\n')
    for dim, r in result.items():
        print(dim, 'cases:', r['cases'], ' failures:', {k: len(v) for k, v in r['failures'].items()})
    print('wrote', OUT.relative_to(OUT.parents[2]))
