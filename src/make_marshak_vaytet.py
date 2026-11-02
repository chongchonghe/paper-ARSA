"""Fig. 10 of the paper (fig:marshak_vaytet): the multigroup Marshak wave, from
data/marshak_vaytet_profiles.csv, _scores.csv and _groups.csv (written by sim_marshak_vaytet.jl). Writes
figures/marshak_vaytet.pdf.

    python make_marshak_vaytet.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, C1, C2, read, col, save


def vaytet():
    """Multigroup Marshak wave. Top: gas temperature at phi0 = 0.01 (profiles.csv) and the error against
    phi0 (scores.csv). Bottom, one column per group (profiles.csv, groups.csv): the radiation
    (brightness) temperature of the group, with its opacity and cell optical depth, and below it phi.
    RSLA blue, ARSA orange."""
    from matplotlib.ticker import FixedLocator, NullLocator
    gp = read('marshak_vaytet_profiles.csv'); s = read('marshak_vaytet_scores.csv'); info = read('marshak_vaytet_groups.csv')

    def tex(v):
        if 0.01 <= v < 1000:
            return f'{v:.3g}'
        m, e = f'{v:.2e}'.split('e')
        return rf'{float(m):.3g}\times10^{{{int(e)}}}'

    plt.rcParams.update(STYLE)
    fig = plt.figure(figsize=(TEXTWIDTH, 4.9))
    sf_top, sf_bot = fig.subfigures(2, 1, height_ratios=[1, 1.35])
    a0, a1 = sf_top.subplots(1, 2)
    runs = {k: [r for r in gp if r['run'] == k] for k in ('ref', 'rsla', 'arsa')}
    x = col(runs['ref'], 'x')
    a0.plot(x, col(runs['ref'], 'T_gas'), '-', color='0.75', lw=4, label=r'$\hat{c} = c$')
    a0.plot(x, col(runs['rsla'], 'T_gas'), '-', color=C1, lw=1.3, label='RSLA')
    a0.plot(x, col(runs['arsa'], 'T_gas'), '-', color=C2, lw=1.3, label='ARSA')
    a0.set_xlabel(r'$x$ [cm]'); a0.set_ylabel(r'$T_{\rm gas}$ [K]'); a0.set_xlim(0, 20); a0.set_xticks(range(0, 21, 5))
    a0.set_title(r'$\phi_0 = 0.01$')
    a0.legend()
    ph = col(s, 'phi0')
    a1.loglog(ph, 100 * col(s, 'L1_rsla'), 's-', color=C1, mfc='none', lw=1.3, ms=5, label='RSLA')
    a1.loglog(ph, 100 * col(s, 'L1_arsa'), 'o-', color=C2, mfc='none', lw=1.3, ms=5, label='ARSA')
    a1.set_xlabel(r'$\phi_0$'); a1.set_ylabel(r'$L_1$ error in $T_{\rm gas}$ [\%]')
    a1.xaxis.set_major_locator(FixedLocator(ph)); a1.xaxis.set_minor_locator(NullLocator())
    a1.set_xticklabels([f'{v:g}' for v in ph])
    a1.yaxis.set_major_locator(FixedLocator([3, 5, 10, 20, 30])); a1.yaxis.set_minor_locator(NullLocator())
    a1.set_yticklabels(['3', '5', '10', '20', '30'])
    a1.legend()
    for a in (a0, a1):
        a.grid(True, which='major', color='0.92', lw=0.5)

    b = sf_bot.subplots(2, 4, sharex=True, sharey='row', height_ratios=[2, 1])
    for g in range(1, 5):
        b0, b1 = b[0, g - 1], b[1, g - 1]
        b0.plot(x, col(runs['ref'], f'T_g{g}'), '-', color='0.75', lw=4)
        b0.plot(x, col(runs['rsla'], f'T_g{g}'), '-', color=C1, lw=1.2)
        b0.plot(x, col(runs['arsa'], f'T_g{g}'), '-', color=C2, lw=1.2)
        b1.semilogy(x, col(runs['arsa'], f'phi_g{g}'), '-', color=C2, lw=1.2)
        b1.axhline(0.01, color=C1, lw=0.8)
        q = info[g - 1]
        b0.set_title(f'group {g}')
        b0.text(0.96, 0.95, rf'$\kappa = {tex(float(q["kappa"]))}$ cm$^2$ g$^{{-1}}$' + '\n'
                + rf'$\tau_{{\rm c}} = {tex(float(q["tau_c"]))}$',
                transform=b0.transAxes, ha='right', va='top', fontsize=7.5, linespacing=1.3)
        b1.set_ylim(5e-3, 2); b1.set_xlim(0, 20)
        b1.set_xlabel(r'$x$ [cm]')
        for a in (b0, b1):
            a.grid(True, color='0.92', lw=0.5)
    b[0, 0].set_ylabel(r'$T_{{\rm rad},g}$ [K]')
    b[1, 0].set_ylabel(r'$\phi_g$')
    save(fig, 'marshak_vaytet')


vaytet()
