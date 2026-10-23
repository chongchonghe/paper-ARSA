"""Sec. 2.4: the high-tau numbers of Fig. 2 (fig:waves_RSLA) from the parameters of Sec. 2.3, checked
against the dispersion relation (waves.tracked_modes, as in fig_waves_rsla.py), at tau = 1e4.
  h = 4 (gamma - 1) P (eq:h_ratio) and h/phi;
  acoustic speed sqrt(gamma_phi), gamma_phi = 1 + (gamma - 1)/(1 + h/phi) (eq:gamma_phi): gas and radiation
      compressed together, with the apparent heat-capacity ratio h/phi;
  slow mode -Im Omega = 1/Pe_phi at M = 1 (eq:slow_rate, eq:Pe_phi), with h -> h/gamma at Pe_phi >> 1;
  acoustic decay about (gamma_phi - 1)/2 times the slow-mode rate (sound damped by heat conduction).

    python check_waves_heat_capacity.py
"""
import numpy as np
from waves import tracked_modes

C, M, gamma, P, phi = 1e4, 0.0, 5 / 3, 5e-4, 1e-3     # the parameters of fig_waves_rsla.py
taus = np.logspace(-6, 6, 1201)
i = int(np.argmin(abs(taus - 1e4)))
tau = taus[i]
h = 4 * (gamma - 1) * P
print(f'h = {h:.4g}, h/phi = {h / phi:.4g}, h/gamma = {h / gamma:.3g}; tau = {tau:.3g}')
for name, f in (('exact', 1.0), ('RSLA', phi)):
    m = tracked_modes(taus, M, C, P, gamma, lambda t: f)
    x = h / f
    g_phi = 1 + (gamma - 1) / (1 + x)
    pe = 3 * tau / C * (1 / h + 1 / f)                  # eq:Pe_phi at M = 1
    slow_num = min(-m[k][i].imag for k in ('ent', 'rad+', 'rad-'))   # the RSLA carries it in a radiation branch
    slow = C * h / (3 * tau * (gamma + x))              # eq:slow_rate with h/gamma
    ac = (g_phi - 1) / 2 * slow
    print(f'{name:5s}: Pe_phi = {pe:.3g}, 1/Pe_phi = {1 / pe:.3g}; acoustic speed {m["ac+"][i].real:.4f}, '
          f'sqrt(gamma_phi) = {np.sqrt(g_phi):.4f}; slow decay {slow_num:.4g}, closed form {slow:.4g}; '
          f'acoustic decay {-m["ac+"][i].imag:.4g}, (gamma_phi - 1)/2 x slow = {ac:.4g}')
x = h / phi
print(f'exact/RSLA: slow decay (gamma + h/phi)/(gamma + h) = {(gamma + x) / (gamma + h):.3f}; '
      f'(gamma_1 - 1)/(gamma_phi - 1) = (1 + h/phi)/(1 + h) = {(1 + x) / (1 + h):.3f}; '
      f'product {(gamma + x) / (gamma + h) * (1 + x) / (1 + h):.2f}')
