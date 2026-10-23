"""Compare the closed-form PD-ARS bound with the full 5x5 von Neumann boundary."""
import numpy as np
from pdars import nu_max

def nu_formula(tau, Pi, gamma=5/3, chi=1/3):
    sig = np.sqrt(chi)
    Lam = 4 * (gamma - 1) * Pi
    den = sig - tau / (2 + Lam * tau)
    return np.inf if den <= 0 else 1 / den

worst = 0
for chi in (1/3, 1.0):
    for gamma in (5/3, 1.4):
        for nuc in (30., 300., 3000.):
            for P in (1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1, 10):
                for tau in (1e-3, 0.1, 0.5, 1, 2, 5, 20, 1e3, 1e6):
                    n = nu_max(tau, nuc, P, chi=chi, gamma=gamma, rel=1e-4)
                    f = min(nu_formula(tau, P * nuc, gamma, chi), nuc)
                    err = abs(n - f) / f
                    worst = max(worst, err)
                    if err > 2e-3:
                        print(f"chi={chi:.3f} g={gamma} nuc={nuc} P={P} tau={tau}: num={n:.5g} formula={f:.5g}")
print("worst relative deviation:", worst)
