"""Linear waves of the RHD system with the derived stability-bound phi.

A wave of optical depth tau = calL/k (per radian) resolved by N cells has cell
optical depth tau_c = 2 pi tau / N.  The step is the RSLA step,
c dt/dx = nuc = C_r/phi0 (1D), unless a reduced step is requested.
phi = min(1, max(phi0, beta*phimax)), phimax = 1/(nuc*(sig - tau_c/(2+K))),
K = 4 (gamma-1) P nuc tau_c, sig = 1/sqrt(3) (P1)."""
import numpy as np
from rhd_waves import eigensystem, classify_modes

SIG = 1 / np.sqrt(3)

def phi_sb(tau, P, phi0, gamma=5/3, N=32, Cr=0.4, beta=0.9, nuc=None):
    nuc = Cr / phi0 if nuc is None else nuc
    tc = 2 * np.pi * tau / N
    K = 4 * (gamma - 1) * P * nuc * tc
    den = SIG - tc / (2 + K)
    if den <= 0:
        return 1.0
    return min(1.0, max(phi0, beta / (nuc * den)))

def modes(tau, M, C, P, gamma, phi):
    lam, vec = eigensystem(tau, M, C, P, gamma, phi)
    return classify_modes(lam, vec)

if __name__ == '__main__':
    C, M, gamma, phi0 = 1e4, 0.0, 5/3, 3e-3
    taus = np.logspace(-3, 7, 11)
    for P in (1e-4, 1e-2, 1e-1, 1.0):
        print(f"--- P = {P}")
        for tau in taus:
            ex = modes(tau, M, C, P, gamma, 1.0)
            rs = modes(tau, M, C, P, gamma, phi0)
            ph = phi_sb(tau, P, phi0)
            ir = modes(tau, M, C, P, gamma, ph)
            err = lambda m: max(abs(m['acoustic'][1] - ex['acoustic'][1]) / abs(ex['acoustic'][1]),
                                abs(m['entropy'] - ex['entropy']) / abs(ex['entropy']))
            print(f"tau={tau:8.1e} phi={ph:8.4g}  err RSLA={err(rs):9.3g}  IRSLA={err(ir):9.3g}")


def tracked_modes(taus, M, C, P, gamma, phif):
    """Eigenvalues along a dense tau grid, each mode followed by continuity from
    the low-tau end. There the modes are classified by speed relative to the flow:
    the two radiation modes are the fastest, |Re(Omega) - M| ~ phi C sigma >> 1, the
    entropy mode is the slowest, and the acoustic pair is in between. (Eigenvector
    weights are ambiguous there when radiation dominates the heat capacity.)
    Returns dict name -> array over taus for rad+, rad-, ac+, ac-, ent."""
    from rhd_waves import eigensystem as es
    lam, _ = es(taus[0], M, C, P, gamma, phif(taus[0]))
    lam = lam[np.argsort(np.abs(lam.real - M))]      # slowest first
    ent, ac, rad = lam[0], np.sort_complex(lam[1:3]), np.sort_complex(lam[3:5])
    first = np.array([rad[1], rad[0], ac[1], ac[0], ent])
    return _track(taus, first, M, C, P, gamma, phif)


def _track(taus, prev, M, C, P, gamma, phif):
    import itertools
    from rhd_waves import eigensystem as es
    names = ['rad+', 'rad-', 'ac+', 'ac-', 'ent']
    out = [prev]
    perms = list(itertools.permutations(range(5)))
    for t in taus[1:]:
        lam, _ = es(t, M, C, P, gamma, phif(t))
        scale = np.maximum(np.abs(prev), 1e-12)
        best = min(perms, key=lambda p: sum(abs(lam[p[k]] - prev[k]) / scale[k] for k in range(5)))
        prev = np.array([lam[best[k]] for k in range(5)])
        out.append(prev)
    out = np.array(out)
    return {n: out[:, k] for k, n in enumerate(names)}
