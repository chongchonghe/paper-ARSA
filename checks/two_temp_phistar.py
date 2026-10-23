"""Estimate of phi_* (eq:phi_star) for matter at two temperatures (sssec:two_temp).

In a thick cell the damping term alpha tau_E/(2 + K) tends to alpha tau_E/K, so the
results of Sec. 3.3 hold with h -> h/alpha, h of eq:h_2T:
    h/alpha = 4 a_r T_d^3 (psi/alpha) / (rho c_v),   phi_* = (h/alpha) s / (1 + h/alpha),
with s = nu0 sum_d sigma_d = 1 (the RSLA at its largest stable step, no beta, as tab:media).

psi and alpha need a model of T_d(T, E). For this estimate only: dust without heat
capacity in balance, rho c (kappa_P a_r T_d^4 - kappa_E E) = Gamma, with the gas-dust rate
Gamma = K_gd n^2 sqrt(T) (T - T_d) (Keto & Field 2005 form; K_gd = 2.5e-34 cgs, n = rho/(mu m_H),
the value QUOKKA uses). Then
    psi   = Gamma_T / (Gamma_d + B),   alpha = Gamma_d / (Gamma_d + B),   B = 4 a_r T_d^3 rho kappa_P c,
    psi/alpha = Gamma_T/Gamma_d = 3/2 - T_d/(2T),
with kappa_P = kappa_E = 2e-4 T_d^2 cm^2/g (Bell & Lin 1994, ice grains, per gas mass).
phi_* depends only on psi/alpha; alpha sets how thick a cell must be (tau_c -> alpha tau_c)."""
import numpy as np

AR, KB, MH, C = 7.5657e-15, 1.380649e-16, 1.6726e-24, 2.99792458e10
KGD, PHI0 = 2.5e-34, 3e-3

# medium, T_gas [K], T_d [K], n [cm^-3], mu
ROWS = [
    ('Warm neutral medium', 8000, 20, 1, 1.27),
    ('Cold neutral medium', 80, 20, 50, 1.27),
    ('HII region', 1e4, 50, 1e3, 0.61),
    ('Molecular cloud', 15, 10, 1e3, 2.33),
    ('Dense core', 10, 10, 1e5, 2.33),
    ('Prestellar core, centre', 8, 8, 1e6, 2.33),
    ('Protostellar envelope', 100, 100, 1e8, 2.33),
    ('First hydrostatic core', 100, 100, 1e12, 2.33),
    ('Protoplanetary disk, 1 au', 150, 150, 2.6e14, 2.33),
]

print('%-26s %8s %6s %9s %6s %9s %9s %9s %9s' % ('medium', 'T', 'T_d', 'n', 'psi/a', 'h/alpha',
                                                'phi_*', 'phimax', 'alpha'))
for name, T, Td, n, mu in ROWS:
    rho = n * mu * MH
    rcv = 1.5 * n * KB
    kP = 2e-4 * Td**2
    Gd = KGD * n**2 * np.sqrt(T)                         # -dGamma/dT_d
    GT = Gd * (1.5 - Td / (2 * T))                       # dGamma/dT
    B = 4 * AR * Td**3 * rho * kP * C
    alpha, psi = Gd / (Gd + B), GT / (Gd + B)
    heff = 4 * AR * Td**3 * (psi / alpha) / rcv          # h/alpha
    phis = heff / (1 + heff)
    den = 1 - PHI0 / heff
    phim = 1.0 if den <= 0 else min(1.0, PHI0 / den)
    print('%-26s %8.3g %6.3g %9.2g %6.3f %9.2g %9.2g %9.2g %9.2g' % (name, T, Td, n, psi / alpha, heff,
                                                                   phis, phim, alpha))
