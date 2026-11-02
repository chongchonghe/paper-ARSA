"""Fig. 1 of the paper (fig:accuracy): when the RSLA is accurate and where the ARSA can help, with the systems
of Table 1, from data/rsla_accuracy_systems.csv and data/rsla_accuracy_curves.csv (written by
calc_rsla_accuracy.py). Writes figures/rsla_accuracy.pdf.

    python make_rsla_accuracy.py
"""
from matplotlib.lines import Line2D
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, read, col, save
INK, MUTED = '#1f1f1f', '#555555'
RED, RED_INK = '#f8d6d1', '#8e1b1b'           # RSLA inaccurate, the ARSA cannot help
AMBER, AMBER_INK = '#fde5ad', '#8a4a00'       # the ARSA can help
GREEN, GREEN_INK = '#d8eed3', '#1f5e2a'       # RSLA accurate
PHI0, PE_T = 1e-3, 0.1                        # as in calc_rsla_accuracy.py; PHI0 is the phi of the RSLA
H_LO, H_HI = 0.1 * PHI0, 0.1
DRIFT_D = 0.1 * PHI0                          # as in calc_rsla_accuracy.py: phi >> D
XLIM, YLIM = (1e-11, 1e8), (1e-7, 1e5)     # x to 1e8: room for the column of labels; y in h (P = 3h/8 up to 4e4)


def bold(text):
    return '\n'.join(rf'\textbf{{{line}}}' for line in text.split('\n'))


c = read('rsla_accuracy_curves.csv')
hc, xc = col(c, 'h'), col(c, 'tau_v_c_Pe')

plt.rcParams.update(STYLE)
fig, ax = plt.subplots(figsize=(0.85 * TEXTWIDTH, 4.3))
for lo, hi, color in ((YLIM[0], H_LO, GREEN), (H_LO, H_HI, AMBER), (H_HI, YLIM[1], RED)):
    m = (hc >= lo) & (hc <= hi)
    hh = np.concatenate(([lo], hc[m], [hi]))
    xx = np.interp(np.log(hh), np.log(hc), xc)
    if color == GREEN:                        # the RSLA is accurate only where advected structures do not drift
        ax.fill_betweenx(hh, xx, DRIFT_D, color=GREEN, lw=0)
        ax.fill_betweenx(hh, DRIFT_D, XLIM[1], color=AMBER, lw=0)
    else:
        ax.fill_betweenx(hh, xx, XLIM[1], color=color, lw=0)
ax.plot([DRIFT_D, DRIFT_D], [YLIM[0], H_LO], color=INK, lw=1.0, ls='--')
ax.plot(xc, hc, color=INK, lw=1.6)
pe_label = ax.text(1.5e-9, 1.5e-7, r'$\mathrm{Pe}_\phi = 0.1$', color=INK, ha='left', va='bottom',
                   rotation_mode='anchor')

# the regions
ax.text(2e-11, 50, bold('radiation is faster than\nthe flow; temperature\ntracks the flow; accuracy\ndepends on the problem\n'
        '($\\mathrm{Pe}_\\phi < 0.1$)'), color=INK, fontsize=8, ha='left', va='center',
        bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=MUTED, lw=0.8))
ax.text(0.3, 2.0, bold('RSLA inaccurate:\ndiffusion too slow;\nARSA cannot help'), color=RED_INK,
        ha='center', va='center')
ax.text(XLIM[1] / 1.5, 1.2e-3, bold('ARSA can help:\nwith $\\phi_0$ above $\\hat{h}$,\n$\\phi = 1$ in thick cells'),
        color=AMBER_INK, ha='right', va='center')
ax.text(DRIFT_D / 1.5, 3e-7, bold('RSLA\naccurate'), color=GREEN_INK, ha='right', va='center')
ax.text(DRIFT_D * 2, 2.2e-7, bold('advected structures drift: $\\mathcal{D} > 0.1\\,\\phi$'), color=AMBER_INK,
        ha='left', va='center', fontsize=7)

LABEL = {   # label text and position beside the point: (text, x, y, ha)
    'PPD midplane, 1 au': ('PPD midplane, 1 au', 0.1, 1.1e-5, 'left'),
    'PPD midplane, 0.1 au': ('PPD midplane, 0.1 au', 1.8, 3.4e-6, 'left'),
    'Hot Jupiter atmosphere': ('hot Jupiter atmosphere', 2e-5, 3.5e-5, 'left'),
    'First hydrostatic core': ('first hydrostatic core', 2.5e-3, 1.2e-4, 'left'),
    'PPD outer region, 100 au': ('PPD, 100 au', 5e-7, 4.8e-3, 'right'),
    'Massive protostellar envelope': ('massive protostellar envelope', 3e-4, 0.061, 'left'),
    'PPD atmosphere, 1 au': ('PPD atmosphere, 1 au', 2.5e-7, 0.3, 'right'),
    'Cold molecular cloud': ('cold molecular cloud', 5e-10, 1e4, 'left'),
}
COLUMN = {   # labels in a column at the right, with leader lines, bottom to top (in the order of their h)
    'Star-forming clump': 'star-forming clump', 'ULIRG nucleus': 'ULIRG nucleus',
    'Massive-star envelope': 'massive-star envelope', 'AGN dusty torus': 'AGN dusty torus'}
systems = read('rsla_accuracy_systems.csv')
ycol = dict(zip(COLUMN, np.geomspace(20, 2e4, len(COLUMN))))
xcol = 5e2
for r in systems:
    name, x, y = r['system'], float(r['tau_v_c']), float(r['h'])
    if name == 'Solar core':          # v ~ 0
        ax.annotate('', xy=(XLIM[0], y), xytext=(XLIM[0] * 30, y), arrowprops=dict(arrowstyle='->', color=INK, lw=1.4))
        ax.text(XLIM[0] * 1.5, y * 1.6, 'solar core ($v \\approx 0$):\ndiffusion-driven', color=INK, fontsize=8, ha='left', va='bottom')
        continue
    if name == 'Diffuse ISM':          # L << 1, P ~ 6e13: above the plot
        ax.annotate('', xy=(3e-11, YLIM[1]), xytext=(3e-11, YLIM[1] / 30), arrowprops=dict(arrowstyle='->', color=INK, lw=1.4))
        ax.text(6e-11, YLIM[1] / 2.5, 'diffuse ISM', color=INK, fontsize=8, ha='left', va='center')
        continue
    thin = float(r['tau']) <= 1                 # optically thin: outside the diffusion argument
    ax.plot([x], [y], 'o', ms=6, color=INK, mfc='none' if thin else INK, mew=1.3)
    if name in LABEL:
        txt, xx, yy, ha = LABEL[name]
        ax.text(xx, yy, txt, ha=ha, va='center', fontsize=8, color=INK)
    else:
        yy = ycol[name]
        ax.plot([x * 1.6, xcol / 1.3], [y, yy], color=MUTED, lw=0.6)
        ax.text(xcol, yy, COLUMN[name], ha='left', va='center', fontsize=8, color=INK)

ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=6, color=INK, label='optically thick'),
                   Line2D([], [], marker='o', ls='', ms=6, color=INK, mfc='none', mew=1.3,
                          label='optically thin')],
          loc='lower left', bbox_to_anchor=(0.0, 0.17), handletextpad=0.3)
ax.set(xscale='log', yscale='log', xlim=XLIM, ylim=YLIM,
       xlabel=r'dynamic-diffusion number $\mathcal{D} = \mathcal{M}\mathcal{L}/\mathcal{C}$',
       ylabel=r'radiation-to-gas pressure ratio $\mathcal{P}$')
ax.set_xticks(10.0 ** np.arange(-10, 9, 2))
# the data and regions are in h; the axis is P = 3 h / 8, i.e. h = 4 (gamma - 1) P with gamma = 5/3
P_PER_H = 3 / 8
ax.set_yticks(10.0 ** np.arange(-6, 6, 2) / P_PER_H, labels=[rf'$10^{{{n}}}$' for n in range(-6, 6, 2)])
ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(10.0 ** np.arange(-7, 5) / P_PER_H))
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.canvas.draw()                  # the curve has slope 1 below h ~ 1e-5: align the label with it on the page
w, h = ax.get_window_extent().size
pe_label.set_rotation(np.degrees(np.arctan((h / np.log10(YLIM[1] / YLIM[0])) / (w / np.log10(XLIM[1] / XLIM[0])))))

save(fig, 'rsla_accuracy')
