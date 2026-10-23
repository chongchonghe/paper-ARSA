"""Check the closed-form checkerboard boundaries of the alternative schemes of Sec. 3.2
(forward Euler + implicit exchange, SSP-RK2 with the exchange once at the end, Lax-Wendroff
+ implicit exchange) against the numerically determined boundary of the full 5x5 amplification
matrix, maximized over theta. Reports the largest relative difference per scheme."""
import numpy as np
from pdars import nu_max

TH = np.pi * np.linspace(0.02, 1.0, 50)


def closed(scheme, tau, nuc, P, s, g1):
    K = 4 * g1 * P * nuc * tau
    if scheme in ("pdars", "euler"):      # nu [s - tau/(2+K)] <= 1
        br = s - tau / (2 + K)
        nu = 1 / br if br > 0 else np.inf
    elif scheme == "rk2end":              # a <= 2
        nu = 1 / s
    elif scheme == "lw":                  # nu^2 s^2 - nu tau/(2+K) <= 1
        b = tau / (2 + K)
        nu = (b + np.sqrt(b**2 + 4 * s**2)) / (2 * s**2)
    return min(nu, nuc)


if __name__ == "__main__":
    rng = np.random.default_rng(5)
    worst = {k: 0.0 for k in ("pdars", "euler", "rk2end", "lw")}
    for _ in range(40):
        s = rng.choice([1 / np.sqrt(3), 1.0]); g1 = 2 / 3
        nuc = 10 ** rng.uniform(1.5, 3); P = 10 ** rng.uniform(-6, -1); tau = 10 ** rng.uniform(-2, 3)
        for sch in worst:
            num = nu_max(tau, nuc, P, scheme=sch, thetas=TH, chi=s**2)
            cf = closed(sch, tau, nuc, P, s, g1)
            worst[sch] = max(worst[sch], abs(num / cf - 1))
    print("max relative difference, numerical vs closed form:", worst)
