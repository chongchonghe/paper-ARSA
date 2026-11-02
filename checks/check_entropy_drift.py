"""Checks the drift of thermal structures under the RSLA (paper Sec. 2.2: a speed of order (h/phi) v, and a drift of
order D/phi times the size of the structure before it diffuses away, whatever h), from the linear waves of the RHD
equations on a moving background (src/waves.py, Mach M = 5, C = 1e4, gamma = 5/3). For the entropy mode, which
carries a temperature structure with the flow:
  drift    (M - Re Omega)/M, the speed at which it falls behind the flow, as a fraction of the flow speed; the exact
           equations give zero, the RSLA (h/gamma)(1/phi - 1) to first order in h/phi;
  life     the drift accumulated before the mode decays, (M - Re Omega)/(-Im Omega) radians, against 3 D (1 - phi)/phi
           with D = M tau/C the dynamic-diffusion number of the mode.
tau is the optical depth per radian, from static (D << phi) to dynamic diffusion (D >> 1). Prints the ratios
drift/(h (1 - phi)/(gamma phi)) and life/(3 D (1 - phi)/phi), which are 1 to first order in h/phi.

    python checks/check_entropy_drift.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))   # waves.py
from waves import modes

C, M, GAMMA = 1e4, 5.0, 5 / 3
worst_exact = 0.0
print(f"M = {M}, C = {C:g}, gamma = {GAMMA:.4g}; h = 4 (gamma - 1) P at T = rho = 1")
print(f"{'phi':>6s} {'h/phi':>7s} {'tau':>7s} {'D':>8s} {'drift':>10s} {'drift/eq':>9s} {'life':>9s} {'life/eq':>8s}")
for phi in (0.01, 0.1):
    for x in (0.001, 0.01, 0.1):
        P = x * phi / (4 * (GAMMA - 1))
        for tau in (1e1, 1e2, 1e3, 1e4, 1e6):
            ex = modes(tau, M, C, P, GAMMA, 1.0)['entropy']
            rs = modes(tau, M, C, P, GAMMA, phi)['entropy']
            worst_exact = max(worst_exact, abs(ex.real - M) / M)
            drift = (M - rs.real) / M
            life = (M - rs.real) / (-rs.imag)
            D = M * tau / C
            print(f"{phi:6g} {x:7g} {tau:7.0e} {D:8.3g} {drift:10.3e} {drift / (x * (1 - phi) / GAMMA):9.4f}"
                  f" {life:9.3g} {life / (3 * D * (1 - phi) / phi):8.4f}")
print(f"exact equations: largest |Re Omega - M|/M = {worst_exact:.1e}")
