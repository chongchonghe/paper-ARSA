"""Table 1 (tab:media): order-of-magnitude parameters of astrophysical systems and the gain of the ARSA.
Also used by calc_rsla_accuracy.py (Fig. 1), which imports the systems and the functions.

Inputs (rho, T, c_v, kappa, tau, v) are order-of-magnitude estimates, good to a factor of 2-3;
the star-formation and disk rows, and the convective envelopes of stars, come from
data/media_parameters.md. Only optically thick systems are listed. This script recomputes from them
    h      = 4 a_r T^3 / (rho c_v),        the radiation-to-gas heat-capacity ratio (paper eq:h_ratio),
    tau v/c,                                the dynamic-diffusion indicator (dynamic where >~ 1),
    phi_*  = h / (1 + h),                   eq:phi_star at the largest stable RSLA step (nu0 sum sigma = 1),
    phi/phi0                                the limiter eq:phi_sb at phi0 = 1e-2 in a cell of optical depth
                                            tau/100 (100 cells across the system), gray, isotropic radiation
                                            in 3D, RSLA step nu0 = beta/sqrt(3), beta = 0.9.
The printout also gives, for Sec. 3.4, the cell optical depth tau_c = tau/100, phi_* of the limiter
(with beta: h nu0 sum sigma/(1 + beta h)) and tau_thresh of eq:tau_one with beta (inf where phi_* > phi0).
kappa includes electron scattering in the ionized rows.

    python make_table_media.py      # writes tables/tab_media.tex
"""
from pathlib import Path
import numpy as np

AR, C = 7.5657e-15, 2.99792458e10                     # cgs
PHI0, BETA = 1e-2, 0.9
SIGSUM, NU0 = np.sqrt(3.0), 0.9 / np.sqrt(3.0)
CELLS = 100                                           # cells across the system for the phi/phi0 column
OUT = Path(__file__).resolve().parents[1] / 'tables' / 'tab_media.tex'

# specific heat at constant volume [erg/g/K], by composition
CV = {'ion': 2.0e8, 'warm': 8e7, 'cold': 5e7, 'kilonova': 6e6}

# name (LaTeX), rho [g/cm^3], T [K], c_v key, kappa [cm^2/g] (None: n/a), tau across the system (None: << 1),
# v [cm/s], v as printed [km/s or c]
# Rows commented out were removed from Table 1 by the authors (2026-10-21; the optically thin rows and the
# solar core, whose v ~ 0 leaves D undefined, on 2026-11-03).
SYSTEMS = [
    # ('Hot Jupiter atmosphere', 1e-7, 1.5e3, 'warm', 0.01, 1, 3e5, '3'),
    # star formation and disks: data/media_parameters.md (tau = kappa rho L; v in the laboratory frame,
    # Keplerian around 1 Msun in the disk rows); c_v as above, so h is lower than there by 1.4-2.3
    ('Massive protostellar envelope', 1e-14, 100, 'cold', 1, 15, 3e5, '3'),
    ('First hydrostatic core', 5e-12, 100, 'cold', 1, 375, 1e5, '1'),
    # ('PPD atmosphere, 1 au', 1e-14, 200, 'warm', 1, 5e-3, 3e6, '30'),
    ('PPD midplane, 0.1 au', 1e-7, 1000, 'warm', 1, 3e3, 9.4e6, '94'),
    ('PPD midplane, 1 au', 1e-9, 280, 'warm', 1, 500, 3e6, '30'),
    # ('PPD outer region, 100 au', 1e-15, 20, 'cold', 1, 0.10, 3e5, '3'),
    ('AGN dusty torus', 1e-18, 300, 'warm', 5, 15, 2e7, '200'),
    # ('Solar core', 150, 1.5e7, 'sun', 1, 1e12, 0.0, r'$\approx 0$'),
    # ('SN IIP ejecta, 50 d', 1e-12, 8e3, 'ion', 0.2, 3e2, 4e8, r'$4\times10^{3}$'),
    # ('Star-forming clump', 3e-19, 20, 'cold', 1, 1, 5e5, '5'),
    ('Massive-star envelope', 1e-9, 2e5, 'ion', 1, 1e2, 1e6, '10'),
    # convective envelopes of stars: data/media_parameters.md, Sec. C (tau = kappa rho H_p, v from the convective
    # flux; tau and v to an order of magnitude). The middle layers: data/convective_layers.txt
    ('Sun, middle of convection zone', 0.054, 9.6e5, 'ion', 100, 1.9e10, 1.2e4, '0.1'),
    ('RGB star, middle of envelope', 8e-4, 1.8e5, 'ion', 500, 2.8e10, 4.8e4, '0.5'),
    ('AGB star, base of envelope', 3, 7e7, 'ion', 0.34, 4.5e9, 6e5, '6'),
    # ('Cold molecular cloud', 2e-22, 15, 'cold', 0.1, 7e-5, 1e5, '1'),
    # ('Radiation-mediated shock', 7e-11, 1.7e5, 'ion', 0.34, 340, 5e8, r'$5\times10^{3}$'),
    # ('AGN disk, 20 $r_g$', 1e-9, 5e5, 'ion', 0.34, 3e2, 0.004 * C, r'$0.004c$'),
    # ('Kilonova ejecta, 1 d', 2e-13, 1e4, 'kilonova', 10, 9e2, 0.2 * C, r'$0.2c$'),
    ('ULIRG nucleus', 2.5e-19, 80, 'warm', 1, 32, 1e7, '100'),    # tau = kappa rho L: 160 at kappa = 5
    # ('Super-Eddington inflow', 5e-12, 2.6e5, 'ion', 0.34, 27, 0.14 * C, r'$0.14c$'),
    # ('Diffuse ISM', 1e-24, 1e4, 'ion', None, None, 1e6, '10'),    # warm medium, n ~ 0.5 cm^-3
]


def heat_ratio(rho, T, cv):
    return 4 * AR * T**3 / (rho * CV[cv])


def phi(tau_c, h, phi0=PHI0):
    """eq:phi_sb, gray: the energy condition binds; phi = 1 where its bracket is not positive."""
    nuc = NU0 / phi0
    tau_c, h = np.broadcast_arrays(np.asarray(tau_c, float), np.asarray(h, float))
    K = h * nuc * tau_c
    D = SIGSUM - tau_c / (2 + K)
    with np.errstate(divide='ignore'):
        p = np.where(D > 0, BETA / (nuc * np.where(D > 0, D, 1.0)), 1.0)
    return np.minimum(1.0, np.maximum(phi0, p))


def rows():
    """The computed columns, sorted by h."""
    out = []
    for name, rho, T, cv, kappa, tau, v, vtex in SYSTEMS:
        h = heat_ratio(rho, T, cv)
        tvc = None if tau is None or v == 0 else tau * v / C
        gain = 1.0 if tau is None else float(phi(tau / CELLS, h)) / PHI0
        out.append(dict(name=name, rho=rho, T=T, kappa=kappa, tau=tau, v=v, vtex=vtex, h=h, tvc=tvc,
                        phis=h / (1 + h), gain=gain))
    return sorted(out, key=lambda r: r['h'])


def tex(x, sig=1, lo=1e-2, hi=1e3):
    """Number in LaTeX with sig significant figures: plain for lo <= x < hi, else a x 10^b."""
    if lo <= x < hi:
        return '$%g$' % float('%.*g' % (sig, x))
    m, e = ('%.*e' % (sig - 1, x)).split('e')
    m, e = float(m), int(e)
    if m == 10:
        m, e = 1, e + 1
    return '$10^{%d}$' % e if m == 1 else r'$%g\times10^{%d}$' % (m, e)


def main():
    lines = [r'\begin{table}', r'    \centering', r'    \small', r'    \setlength{\tabcolsep}{2.5pt}',
             r'    \begin{tabular}{lcccccccc}', r'    \hline',
             r'    System & $\rho$ & $T$ & $\kappa$ & $v$ & '
             r'$\mathcal{L}$ & $\mathcal{D}$ & $\hat{h}$ & $\phi/\phi_0$ \\',
             r'    & [g cm$^{-3}$] & [K] & [cm$^2$ g$^{-1}$] & [km s$^{-1}$] \\',
             r'    \hline']
    for r in rows():
        kap = '--' if r['kappa'] is None else tex(r['kappa'], 2, hi=1e2)
        tau = r'$\ll 1$' if r['tau'] is None else tex(r['tau'], 2, hi=1e3)
        tvc = '--' if r['v'] == 0 else r'$\ll 1$' if r['tvc'] is None else tex(r['tvc'], hi=1e2)
        g = r['gain']
        gain = '$%.0f$' % g if g > 10 else '$%.1f$' % g
        lines.append('    %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\'
                     % (r['name'], tex(r['rho'], 2), tex(r['T'], 2, hi=1e4), kap, r['vtex'], tau, tvc, tex(r['h'], 2),
                        gain))
        h = r['h']
        tau_c = float('nan') if r['tau'] is None else r['tau'] / CELLS
        ps = h * NU0 * SIGSUM / (1 + BETA * h)                  # phi_* of the limiter
        thr = (2 * (SIGSUM - BETA * PHI0 / NU0) * PHI0 / ((1 + BETA * h) * (PHI0 - ps)) if ps < PHI0
               else float('inf'))                               # eq:tau_one with beta
        print('%-26s h = %8.2g  tau v/c = %8s  phi_* = %6.3g  phi/phi0 (tau/%d) = %6.3g  tau_c = %7.3g  '
              'phi_*(beta) = %7.3g  tau_thresh = %5.3g'
              % (r['name'], h, '-' if r['tvc'] is None else '%.2g' % r['tvc'], r['phis'], CELLS, g, tau_c, ps, thr))
    lines += [r'    \hline', r'    \end{tabular}',
              r'    \caption{Representative parameters of optically thick astrophysical systems, good to a factor of '
              r'2--3 (for the stars, $\mathcal{L}$ and $v$ to an order of magnitude), and the '
              r'gain of the ARSA, sorted by the heat-capacity ratio $\hat{h}$. Columns are, from left to right, the system '
              r'being described (here PPD = protoplanetary disk, AGN = active galactic nucleus, ULIRG = ultraluminous '
              r'infrared galaxy, RGB = red giant branch, AGB = asymptotic giant branch), its characteristic density $\rho$, temperature $T$, specific opacity $\kappa$, '
              r'velocity $v$, optical depth across the system $\mathcal{L}$, dynamic diffusion number '
              r'$\mathcal{D} = \mathcal{M}\mathcal{L}/\mathcal{C} = \mathcal{L}v/c$, ratio of radiation to gas heat capacity '
              r'$\hat{h} = 4a_rT^3/(\rho c_v)$, and the reduction factor that the limiter of the ARSA sets at '
              r'$\phi_0 = 10^{-2}$ in a cell of optical depth $\mathcal{L}/100$, i.e., with 100 cells across the system, '
              r'for isotropic radiation in three dimensions at the RSLA step $\nu_0 = 0.9/\sqrt{3}$, $\phi/\phi_0$. '
              r'We adopt $c_v = 2\times10^8$~erg~g$^{-1}$~K$^{-1}$ for ionized hydrogen and helium, $8\times10^7$ and '
              r'$5\times10^7$ for warm and cold molecular gas. The stellar rows are halfway in radius between the base '
              r'of the convective envelope and the surface, in the present Sun \citep{ChristensenDalsgaard96a} and '
              r'in a $1\,M_\odot$ star on the lower RGB, and at the base of the envelope of an AGB star with hot '
              r'bottom burning \citep{Karakas14b, Doherty17a}; their $\mathcal{L}$ is across a pressure scale height.',
              r'    \label{tab:media}', r'    }', r'\end{table}']
    OUT.write_text('\n'.join(lines) + '\n')
    print('wrote', OUT.relative_to(OUT.parents[1]))


if __name__ == '__main__':
    main()
