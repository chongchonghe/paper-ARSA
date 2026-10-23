"""Is the checkerboard theta = pi always the mode that sets the stability
boundary of the PD-ARS update? (Appendix of the paper.)

For each parameter case:
  nu_pi   boundary of the theta = pi mode alone (bisection on the 5x5 G)
  test 1  stable for all theta on a fine grid at every nu in [1e-3, 1-d] * nu_pi
  test 2  at nu = (1+d) nu_pi the largest |eigenvalue| over theta occurs at pi
Cases where phi = 1 is stable at theta = pi (nu_pi = nuc) only need test 1.
The trivially neutral density mode (eigenvalue exactly 1) is excluded by
taking the eigenvalues of the (v, T, E, F) block.
Writes data/theta_scan.csv (Table B.3) and data/theta_curves.csv (Fig. B.10).
"""
import itertools, csv
import numpy as np
from pdars import matrices
from paper_style import DATA

C = 1e4
D = 1e-3                       # relative bracket around nu_pi
TOL = 1e-9                     # stability tolerance on the spectral radius
THETA = np.pi * np.concatenate([np.geomspace(1e-4, 4e-3, 20), np.linspace(5e-3, 1.0, 400)])
I4 = np.eye(4)


def nu_formula(tau, Pi, gamma, chi):
    """Closed-form bound, eq:numax (same as check_formula.nu_formula)."""
    den = np.sqrt(chi) - tau / (2 + 4 * (gamma - 1) * Pi * tau)
    return np.inf if den <= 0 else 1 / den


def growth(thetas, tau, nu, nuc, P, M, gamma, chi):
    """Spectral radius of G over the (v, T, E, F) block, one value per theta."""
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


def nu_pi(tau, nuc, P, M, gamma, chi):
    ok = lambda nu: growth(np.pi, tau, nu, nuc, P, M, gamma, chi)[0] <= 1 + TOL
    if ok(nuc):
        return nuc
    lo, hi = 0.0, nuc
    while hi - lo > 1e-9 * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def run_case(tau, nuc, P, M, gamma, chi):
    npi = nu_pi(tau, nuc, P, M, gamma, chi)
    below = npi * np.geomspace(1e-3, 1 - D, 30)
    gbelow = max(growth(THETA, tau, nu, nuc, P, M, gamma, chi).max() for nu in below)
    active = npi < nuc
    if active:
        g = growth(THETA, tau, npi * (1 + D), nuc, P, M, gamma, chi)
        theta_star = THETA[np.argmax(g)] / np.pi
    else:
        theta_star = np.nan
    f = min(nu_formula(tau, P * nuc, gamma, chi), nuc)
    return dict(M=M, chi=chi, gamma=gamma, nuc=nuc, P=P, tau=tau, nu_pi=npi, active=int(active),
                max_g_below=gbelow, theta_star=theta_star, nu_formula=f)


def scan():
    rows = []
    for M, chi, gamma, nuc, P, tau in itertools.product(
            (0.0, 0.5, 3.0), (1 / 3, 1.0), (5 / 3, 1.4), (30., 300., 3000.),
            (1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1, 10), (1e-3, 0.1, 0.5, 1, 2, 5, 20, 1e3, 1e6)):
        if M != 0 and chi != 1 / 3:
            continue                      # matrices() has chi != 1/3 for M = 0 only
        rows.append(run_case(tau, nuc, P, M, gamma, chi))
    with open(DATA / 'theta_scan.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)
    return rows


def nu_crit(theta, tau, nuc, P, gamma=5 / 3, chi=1 / 3):
    """Largest nu <= nuc for which the single mode theta is stable (scan, then bisection)."""
    ok = lambda nu: growth(theta, tau, nu, nuc, P, 0.0, gamma, chi)[0] <= 1 + TOL
    grid = np.geomspace(1e-2, nuc, 200)
    bad = [nu for nu in grid if not ok(nu)]
    if not bad:
        return nuc
    hi = bad[0]; lo = grid[grid < hi][-1] if (grid < hi).any() else 0.0
    while hi - lo > 1e-6 * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def curves():
    """Stability limit nu_crit(theta) of each Fourier mode, for the figure.
    Cells are labelled by Lambda = 4 (gamma-1) P nuc, so that K = Lambda tau."""
    out = []
    nuc, gamma = 1e4, 5 / 3
    th = np.pi * np.concatenate([np.geomspace(1e-3, 0.05, 15), np.linspace(0.06, 1.0, 60)])
    for tau, Lam in ((0.1, 0.0), (1.0, 0.0), (1.0, 1.5), (10.0, 10.0)):
        P = max(Lam / (4 * (gamma - 1) * nuc), 1e-14)
        npi = nu_crit(np.pi, tau, nuc, P)
        out += [dict(tau=tau, Lam=Lam, nuc=nuc, nu_pi=npi, theta=t / np.pi, nu_crit=nu_crit(t, tau, nuc, P))
                for t in th]
    with open(DATA / 'theta_curves.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=out[0].keys())
        w.writeheader(); w.writerows(out)


if __name__ == '__main__':
    rows = scan()
    curves()
    act = [r for r in rows if r['active']]
    print(f"cases: {len(rows)}, boundary active (nu_pi < nuc): {len(act)}")
    print("max over cases of max_theta rho below nu_pi - 1:", max(r['max_g_below'] for r in rows) - 1)
    print("active cases with theta* = pi:", sum(r['theta_star'] == 1.0 for r in act), "/", len(act))
    print("group                  cases active  max(rho-1) below  theta*=pi  max|nu_pi/nu_f-1|")
    for M, chi in ((0.0, 1 / 3), (0.0, 1.0), (0.5, 1 / 3), (3.0, 1 / 3)):
        g = [r for r in rows if r['M'] == M and r['chi'] == chi]
        a = [r for r in g if r['active']]
        print(f"M={M:<3} chi={chi:.3f}   {len(g):5d} {len(a):6d}  {max(r['max_g_below'] for r in g) - 1:16.1e}"
              f"  {sum(r['theta_star'] == 1.0 for r in a):4d}/{len(a):<4d}"
              f"  {max(abs(r['nu_pi'] / r['nu_formula'] - 1) for r in g):.1e}")
