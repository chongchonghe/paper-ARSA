"""Stability of the IRSLA with the wave-speed correction of the energy dissipation
(Sec. 3.3), over every wavelength, in a uniform medium.

The limiter replaces sigma by eps*sigma, eps = min(1, 1/tau_c), in the energy
bracket. That is the checkerboard bound of the corrected scheme, but with the energy
dissipation lowered an intermediate wavelength (theta ~ pi/2) can grow where the
radiation's inertia is not negligible, K_F = (1 + chi) P nuc tau_c / C^2 (it vanishes
as C grows). Honghu.jl therefore counts the correction only where
2 E kappa_F dt / c = 2 P nuc tau_c / C^2 <= KF_MAX and bounds the cell as if
uncorrected elsewhere. For each threshold this reports the cases in which the phi the
limiter picks (beta = 0.9) is unstable at some theta, and the cases in which the
correction raises phi above the uncorrected limiter's value."""
import itertools
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))   # pdars.py
from pdars import growth

G1 = 2 / 3
TH = np.pi * np.linspace(0.01, 1.0, 200)
GRID = list(itertools.product((1 / 3, 1.0), np.logspace(0, 6, 13), (30, 300, 3000),
                              (1e-2, 1.0, 10.0, 30.0), (1e3, 1e4, 1e5)))


def phi_limiter(tau, nuc, P, chi, eps):
    K = 4 * G1 * P * nuc * tau
    b = eps * np.sqrt(chi) - tau / (2 + K)
    return 1.0 if b <= 0 else min(1.0, 0.9 / (nuc * b))


def scan(kf_max):
    bad, gain = [], []
    for chi, tau, nuc, P, C in GRID:
        phi_off = phi_limiter(tau, nuc, P, chi, 1.0)
        if kf_max is None:          # the uncorrected scheme
            phi, wc = phi_off, False
        else:
            KF = 2 * P * nuc * tau / C**2
            phi = phi_limiter(tau, nuc, P, chi, min(1, 1 / tau)) if KF <= kf_max else phi_off
            wc = True
        g = max(growth(t, tau, phi * nuc, nuc, P, chi=chi, C=C, ap=wc) for t in TH)
        if g > 1 + 1e-9:
            bad.append((chi, tau, nuc, P, C, phi, g))
        if phi > 1.01 * phi_off:
            gain.append(phi / phi_off)
    return bad, gain


if __name__ == "__main__":
    print(f"{len(GRID)} cases: chi = 1/3, 1; tau_c = 1 ... 1e6; nuc = 30 ... 3000; "
          "P = 0.01 ... 30; C = 1e3 ... 1e5")
    for kf in (None, np.inf, 10.0, 1.0, 0.3, 0.1):
        bad, gain = scan(kf)
        label = "no correction" if kf is None else f"K_F <= {kf:g}"
        extra = (f" (median factor {np.median(gain):.3g})"
                 if gain else "")
        print(f"{label:>14}: unstable {len(bad):3d}; correction raises phi in "
              f"{len(gain):3d}{extra}")
