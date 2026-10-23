"""Fig. 6 of the paper (fig:diffusion_pulse): the thick pulse at rest and advected, from
data/diffusion_pulse_profiles.csv (written by sim_diffusion_pulse.jl). Writes figures/diffusion_pulse.pdf.

    python make_diffusion_pulse.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, C1, C2, read, col, save


def diffusion_pulse():
    """The pulse of sim_diffusion_pulse.jl (T_lo = 1e6 K, phi0 = 0.01): left at rest
    (static diffusion), right advected at v0 = 3e7 cm/s (dynamic diffusion), against x - v0 t wrapped
    into the periodic box, so that the exact pulse stays at the centre on the initial profile. Rows: T,
    its relative error against chat = c, and phi. RSLA blue, ARSA orange."""
    p = [r for r in read('diffusion_pulse_profiles.csv') if float(r['T_lo']) == 1e6 and float(r['phi0']) == 0.01]
    L = 1024.0
    columns = ((0.0, r'At rest, $t = 4.8\times10^{-4}$ s', r'$x$ [cm]'),
               (3e7, r'$v = 3\times10^{7}$ cm s$^{-1}$, $t = 4.8\times10^{-5}$ s', r'$x - vt$ [cm]'))
    plt.rcParams.update(STYLE)
    fig, axes = plt.subplots(3, 2, figsize=(0.85 * TEXTWIDTH, 5.2), sharex='col',
                             gridspec_kw=dict(height_ratios=[2, 1.2, 1]))
    for j, (v0, title, xlabel) in enumerate(columns):
        rows = [r for r in p if float(r['v0']) == v0]
        x0, Ti = col(rows, 'x'), col(rows, 'T_initial')
        xc = (x0 - v0 * float(rows[0]['t']) + L / 2) % L - L / 2      # comoving position
        rows = [rows[i] for i in np.argsort(xc)]
        x = np.sort(xc)
        ref = col(rows, 'T_ref')
        a0, a1, a2 = axes[:, j]
        a0.plot(x0, Ti / 1e6, ':', color='0.5', lw=1, label='initial')
        a0.plot(x, ref / 1e6, '-', color='0.75', lw=4, label=r'$\hat{c} = c$', zorder=1)
        for k, c, lab in (('rsla', C1, 'RSLA'), ('arsa', C2, 'ARSA')):
            a0.plot(x, col(rows, 'T_' + k) / 1e6, '-', color=c, lw=1.2, label=rf'{lab}, $\phi_0 = 0.01$', zorder=2)
            a1.semilogy(x, np.abs(col(rows, 'T_' + k) - ref) / ref, '-', color=c, lw=1.2)
        a2.semilogy(x, col(rows, 'phi_arsa'), '-', color=C2, lw=1.2)
        a2.axhline(0.01, color=C1, lw=1)
        a0.set_title(title)
        a1.set_ylim(1e-6, 2)
        a2.set_ylim(3e-3, 3)
        a2.set_xlim(-200, 200)
        a2.set_xlabel(xlabel)
        for a in (a0, a1, a2):
            a.grid(True, color='0.92', lw=0.5)
        if j == 0:
            a0.set_ylabel(r'$T$ [$10^6$ K]')
            a1.set_ylabel('relative error')
            a2.set_ylabel(r'$\phi$')
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='outside lower center', ncol=4)
    save(fig, 'diffusion_pulse')


diffusion_pulse()
