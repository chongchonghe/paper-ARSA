"""The middle layers of the convective envelopes of the Sun and an RGB star, halfway in radius between the base of
the envelope and the surface: the stellar rows of Table 1 (tab:media; data/media_parameters.md, Sec. C), with
the base of the envelope of an AGB star. Fully ionized gas plus radiation, mu = 0.6.
  Sun: Model S (Christensen-Dalsgaard et al. 1996, Science 272, 1286; file cptrho.l5bi.d.15c) at r = 0.85 R_sun,
      halfway between the base of the convection zone (0.713 R_sun) and the surface.
  RGB star (1 M_sun, 30 L_sun, T_eff = 4600 K, so R = 8.6 R_sun; M_r = 0.3 M_sun, the helium core): from the
      base of the envelope (T = 2e6 K, rho = 0.03 g/cm^3), at the radius r_b where a point-mass envelope with
      nabla_ad = 0.4 has this T, T_b = nabla_ad mu m_H G M_r (1/r_b - 1/R)/k, outward along the adiabat
      ln(T^(3/2)/rho) + h/2 = const (the entropy of gas plus radiation, in units of k/(mu m_H)), with
      d(1/r) = dP/(rho G M_r), to r = (r_b + R)/2. With h << 1 the adiabat gives h ~ T^3/rho ~ rho.
  AGB star: h along the adiabat from the base of the envelope (T = 7e7 K, rho = 3 g/cm^3): with h >> 1 the
      radiation carries most of the entropy and h hardly falls, so Table 1 keeps the base; the point-mass
      envelope does not apply, since the envelope holds most of the mass.
In each layer: the pressure scale height H_p = P r^2/(rho G M_r), the convective velocity v = (L/(4 pi r^2 rho))^(1/3),
the opacity kappa = 0.34 + K rho T^(-7/2) (electron scattering plus a Kramers term, K set so that kappa = 20 at the
base of the solar convection zone, where the radiation carries the luminosity), L = kappa rho H_p, D = L v/c and
h = 4 a_r T^3/(rho c_v), c_v = 2e8 erg/g/K as in Table 1. Writes no file; run_all.sh saves the printout in
data/convective_layers.txt.

    python calc_convective_layers.py
"""
import numpy as np

A, C, G, KB, MH = 7.5657e-15, 2.99792458e10, 6.674e-8, 1.380649e-16, 1.6726e-24
MSUN, LSUN, RSUN = 1.989e33, 3.828e33, 6.957e10
MU, CV, KES = 0.6, 2.0e8, 0.34
KRAM = (20 - KES) / (0.189 * 2.2e6**-3.5)               # kappa = 20 at the base of the solar convection zone


def h_of(T, rho):
    return 4 * A * T**3 / (rho * CV)


def pressure(T, rho):
    return rho * KB * T / (MU * MH) + A * T**4 / 3


def entropy(T, rho):                                    # in units of k/(mu m_H), gas plus radiation
    return np.log(T**1.5 / rho) + 4 * A * T**3 * MU * MH / (3 * rho * KB)


def rho_on_adiabat(T, s):                               # entropy falls with rho: bisection in log rho
    lo, hi = -30.0, 5.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if entropy(T, 10**mid) > s else (lo, mid)
    return 10**((lo + hi) / 2)


def outward(T_b, rho_b, M, r_b, r_end, n=40000):        # base (T_b, rho_b, r_b) outward along the adiabat to r_end
    s, inv_r = entropy(T_b, rho_b), 1 / r_b
    T_prev, rho_prev = T_b, rho_b
    for T in np.geomspace(T_b, 1e3, n)[1:]:
        rho = rho_on_adiabat(T, s)
        dP = pressure(T, rho) - pressure(T_prev, rho_prev)
        inv_r += dP / (0.5 * (rho + rho_prev) * G * M)     # dP/dr = -rho G M/r^2, i.e. d(1/r) = dP/(rho G M)
        T_prev, rho_prev = T, rho
        if 1 / inv_r >= r_end:
            return T, rho, 1 / inv_r
    raise ValueError('the envelope ends before r_end')


def layer(name, T, rho, M, L, r, P=None):
    P = pressure(T, rho) if P is None else P
    Hp = P * r**2 / (rho * G * M)
    v = (L / (4 * np.pi * r**2 * rho))**(1 / 3)
    kap = KES + KRAM * rho * T**-3.5
    tau = kap * rho * Hp
    print(f'{name:34s} r = {r / RSUN:6.3g} R_sun  T = {T:8.3g} K  rho = {rho:8.3g} g/cm^3  P_rad/P_gas = '
          f'{A * T**4 / 3 / (rho * KB * T / (MU * MH)):7.2g}  kappa = {kap:6.3g} cm^2/g  H_p = {Hp:8.3g} cm  '
          f'v = {v / 1e5:6.3g} km/s  L = {tau:8.3g}  D = {tau * v / C:8.3g}  h = {h_of(T, rho):8.3g}')


print(f'middle layers; kappa = {KES} + {KRAM:.3g} rho T^-3.5; mu = {MU}, c_v = {CV:g}')
# Model S at r = 0.85 R_sun (cptrho.l5bi.d.15c): T, rho, P
layer('Sun (Model S, r = 0.85 R_sun)', 9.62e5, 0.0538, 1.0 * MSUN, LSUN, 0.85 * RSUN, P=7.04e12)
# RGB star
M, L, T_b, rho_b = 0.3 * MSUN, 30 * LSUN, 2e6, 3e-2
R = np.sqrt(30) * (5772 / 4600)**2 * RSUN
r_b = 1 / (T_b * KB / (0.4 * MU * MH * G * M) + 1 / R)
print(f'RGB star: R = {R / RSUN:.3g} R_sun, base T = {T_b:g} K, rho = {rho_b:g} g/cm^3 at r_b = {r_b / RSUN:.3g} R_sun, '
      f'h = {h_of(T_b, rho_b):.3g}')
T, rho, r = outward(T_b, rho_b, M, r_b, (r_b + R) / 2)
layer('RGB star (halfway, r = (r_b + R)/2)', T, rho, M, L, r)
# AGB star: h along the adiabat from the base
s = entropy(7e7, 3.0)
print('AGB star, h along the adiabat from the base: ' +
      ', '.join(f'{T:.0e} K: {h_of(T, rho_on_adiabat(T, s)):.2g}' for T in (7e7, 1e7, 1e6, 1e5)))
