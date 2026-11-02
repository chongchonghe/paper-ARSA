"""Figs. 7 and 8 of the paper (fig:shocks_Tgas, fig:shocks_Trad): the radiative shocks of Melon Fuksman et
al. (2021), from data/mf_shocks.csv and data/mf_shocks_profiles.csv (written by sim_mf_shocks.jl). Writes
figures/shocks_Tgas.pdf and figures/shocks_Trad.pdf.

    python make_mf_shocks.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from paper_style import STYLE, TEXTWIDTH, C1, C2, read, col, save


@plt.rc_context(STYLE)
def mf_shocks():
    """Two figures, T_gas and T_rad, of the radiative shocks of Melon Fuksman et al. (2021) on
    64 cells (data/mf_shocks*.csv, from sim_mf_shocks.jl). Left column: test 1 (subcritical). Right
    column: the supercritical shock with 10 x their mass opacity, at phi0 = 0.03 (test 3, solid)
    and at the phi0 at which the ARSA matches the error of the RSLA at 0.03 (test 2, dashed).
    RSLA blue, ARSA orange.
    Rows: the temperature, its relative error against chat = c, and phi."""
    p = read('mf_shocks_profiles.csv')
    sm = {int(r['test']): r for r in read('mf_shocks.csv')}
    columns = ((r'Subcritical, $\kappa\rho = 3.1\times10^{-10}$ cm$^{-1}$', 4, (1,)),
               (r'Supercritical, $\kappa\rho = 3.1\times10^{-9}\,\rho/\rho_0$ cm$^{-1}$', 3, (3, 2)))
    for pre, name, fname in (('T_', r'T_{\rm gas}', 'shocks_Tgas'), ('Tr_', r'T_{\rm rad}', 'shocks_Trad')):
        fig, axes = plt.subplots(3, 2, figsize=(0.85 * TEXTWIDTH, 5.0), sharex='col',
                                 gridspec_kw=dict(height_ratios=[2, 1.2, 1]))
        for j, (title, xmax, tests) in enumerate(columns):
            a0, a1, a2 = axes[:, j]
            for t, ls in zip(tests, ('-', '--')):
                rows = [r for r in p if int(r['test']) == t]
                x = col(rows, 'x') / 1e10
                phi0 = float(sm[t]['phi0'])
                ref = col(rows, pre + 'ref')
                if ls == '-':
                    a0.plot(x, ref, '-', color='0.75', lw=4, label=r'$\hat{c} = c$', zorder=1)
                for k, c, lab in (('rsla', C1, 'RSLA'), ('arsa', C2, 'ARSA')):
                    a0.plot(x, col(rows, pre + k), ls, color=c, lw=1.2, label=rf'{lab}, $\phi_0 = {phi0:.2g}$', zorder=2)
                    a1.semilogy(x, np.abs(col(rows, pre + k) - ref) / ref, ls, color=c, lw=1.2)
                a2.semilogy(x, col(rows, 'phi_arsa'), ls, color=C2, lw=1.2)
                a2.axhline(phi0, color=C1, lw=0.8, ls=ls)
            a0.set_title(title)
            hl = sorted(zip(*a0.get_legend_handles_labels()),        # reference, RSLA, then ARSA
                        key=lambda q: ('RSLA' in q[1]) + 2 * ('ARSA' in q[1]))
            a0.legend(*zip(*hl), loc='upper right', handlelength=2.5)
            a1.set_ylim(1e-5, 3)
            a2.set_ylim(5e-4, 2)
            a2.set_xlim(0, xmax)
            a2.set_xlabel(r'$x$ [$10^{10}$ cm]')
            for a in (a0, a1, a2):
                a.grid(True, color='0.92', lw=0.5)
            if j == 0:
                a0.set_ylabel(rf'${name}$ [K]')
                a1.set_ylabel('relative difference')
                a2.set_ylabel(r'$\phi$')
        save(fig, fname)


mf_shocks()
