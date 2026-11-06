"""Fig. of the paper (fig:expanding_sphere): the interior of the opaque expanding sphere of Roth & Kasen (2015), with
absorption, from data/expanding_sphere.csv (written by sim_expanding_sphere.jl). Top: the box mean of T_gas against
the outer radius r_out [cm] for the RSLA (blue, filled) and the ARSA (orange, open) at each phi0 (markers), with
the exact solution (black line). Bottom: the relative error of the box mean against the exact solution. Writes figures/expanding_sphere.pdf.

    python make_expanding_sphere.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, C1, C2, read, col, save

MARKERS = {3e-4: 'o', 3e-3: 's'}
R0 = 1e13                    # r_out at t = 0 [cm]


@plt.rc_context(STYLE)
def expanding_sphere():
    d = read('expanding_sphere.csv')
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(0.55 * TEXTWIDTH, 3.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1]))
    ex = [r for r in d if r['method'] == 'rsla' and float(r['phi0']) == 3e-4]
    a0.loglog(R0 * col(ex, 'r_out_R0'), col(ex, 'T_exact_K'), '-', color='k', lw=1, label='exact', zorder=1)
    for phi0, m in MARKERS.items():
        for meth, c, lab in (('rsla', C1, 'RSLA'), ('arsa', C2, 'ARSA')):
            rows = [r for r in d if r['method'] == meth and float(r['phi0']) == phi0]
            r, T, Te = R0 * col(rows, 'r_out_R0'), col(rows, 'Tg_mean_K'), col(rows, 'T_exact_K')
            kw = dict(color=c, marker=m, ms=3.5, lw=0.8, mfc='none' if meth == 'arsa' else c, zorder=2)
            a0.loglog(r, T, ls='none', label=rf'{lab}, $\phi_0 = {phi0:g}$', **kw)
            a1.loglog(r[1:], np.abs(T[1:] / Te[1:] - 1), **kw)
    a0.set_ylabel(r'$T$ [K]')
    a1.set_ylabel(r'$|T - T_{\rm exact}|/T_{\rm exact}$')
    a1.set_xlabel(r'$r_{\rm out}$ [cm]')
    a0.legend(ncol=1, loc='lower left', handlelength=2.2)
    save(fig, 'expanding_sphere')


if __name__ == '__main__':
    expanding_sphere()
