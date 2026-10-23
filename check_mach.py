import numpy as np
from pdars import amplification, matrices, THETA
from check_formula import nu_formula

C = 1e4
def cont_growth(theta, tau_c, nu, nuc, P, M, gamma=5/3):
    phi = nu / nuc; dt = nuc / C
    A, S = matrices(M, phi, C, P, gamma)
    om = theta * np.linalg.eigvals(A + 1j * tau_c * S)
    return max(1.0, np.max(np.exp(om.imag * dt)))

def stable(tau, nu, nuc, P, M, hydro):
    for th in THETA:
        G = amplification(th, tau, nu, nuc, P, C=C, M=M, hydro=hydro)
        if np.max(np.abs(np.linalg.eigvals(G))) > cont_growth(th, tau, nu, nuc, P, M) * (1 + 1e-9):
            return False
    return True

def numax(tau, nuc, P, M, hydro):
    lo, hi = 1e-3, nuc
    if stable(tau, hi, nuc, P, M, hydro): return hi
    while hi - lo > 1e-4 * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if stable(tau, mid, nuc, P, M, hydro) else (lo, mid)
    return lo

worst = 0
for M in (0.0, 0.5, 3.0):
    for hydro in (False, True):
        for P in (1e-4, 1e-2, 1.0):
            for tau in (0.01, 0.5, 1.0, 3.0, 100.0):
                nuc = 300.
                n = numax(tau, nuc, P, M, hydro); f = min(nu_formula(tau, P * nuc), nuc)
                worst = max(worst, abs(n - f) / f)
                if abs(n - f) / f > 3e-3:
                    print(f"M={M} hydro={hydro} P={P} tau={tau}: num={n:.5g} formula={f:.5g}")
print("worst:", worst)
