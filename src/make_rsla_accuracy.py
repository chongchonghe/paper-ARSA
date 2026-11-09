"""Fig. 1 of the paper (fig:accuracy): when the RSLA is accurate and where the ARSA can help, with the optically
thick systems of Table 1, from data/rsla_accuracy_systems.csv (written by calc_rsla_accuracy.py). Writes
figures/rsla_accuracy.pdf.

    python make_rsla_accuracy.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, read, save
INK, MUTED = '#1f1f1f', '#555555'
RED, RED_INK = '#f8d6d1', '#8e1b1b'           # neither the RSLA nor the ARSA is accurate
AMBER, AMBER_INK = '#fde5ad', '#8a4a00'       # the ARSA can help
GREEN, GREEN_INK = '#d8eed3', '#1f5e2a'       # RSLA accurate
PHI0 = 1e-2                                   # as in calc_rsla_accuracy.py; PHI0 is the phi of the RSLA
H_LO, H_HI = 0.1 * PHI0, 0.1
DRIFT_D = 0.1 * PHI0                          # as in calc_rsla_accuracy.py: phi >> D
XLIM, YLIM = (1e-8, 1e7), (1e-7, 1e5)         # y in h (P = 3h/8 up to 4e4)


def bold(text):
    return '\n'.join(rf'\textbf{{{line}}}' for line in text.split('\n'))


plt.rcParams.update(STYLE)
fig, ax = plt.subplots(figsize=(0.85 * TEXTWIDTH, 4.3))
ax.fill_between(XLIM, H_HI, YLIM[1], color=RED, lw=0)
ax.fill_between(XLIM, YLIM[0], H_HI, color=AMBER, lw=0)
ax.fill_between([XLIM[0], DRIFT_D], YLIM[0], H_LO, color=GREEN, lw=0)
ax.plot([DRIFT_D, DRIFT_D], [YLIM[0], H_LO], color=INK, lw=1.0, ls='--')
ax.plot([XLIM[0], DRIFT_D], [H_LO, H_LO], color=INK, lw=1.0, ls='--')

# the regions
ax.text(2e-8, 2.0, bold('neither RSLA nor ARSA accurate:\nradiation holds much of the heat'), color=RED_INK,
        ha='left', va='center')
ax.text(2e-8, 1.5e-2, bold('RSLA inaccurate;\nARSA accurate in thick cells\n($\\tau_{\\rm c} > 3.5$ in 3D) if $\\phi_0 > \\hat{h}$'),
        color=AMBER_INK, ha='left', va='center')
ax.text(DRIFT_D / 1.5, 1e-5, bold('RSLA\naccurate'), color=GREEN_INK, ha='right', va='center')
ax.text(DRIFT_D * 2, 2.2e-7, bold('advected structures drift: $\\mathcal{D} > 0.1\\,\\phi_0$'), color=AMBER_INK,
        ha='left', va='center', fontsize=7)

LABEL = {   # label text and position beside the point: (text, x, y, ha)
    'PPD midplane, 1 au': ('PPD midplane, 1 au', 0.05, 2.4e-5, 'center'),
    'PPD midplane, 0.1 au': ('PPD midplane, 0.1 au', 1.8, 3.8e-6, 'left'),
    'First hydrostatic core': ('first hydrostatic core', 2.2e-3, 1.2e-4, 'left'),
    'Massive protostellar envelope': ('massive protostellar\nenvelope', 2.8e-4, 0.061, 'left'),
    'Sun, middle of convection zone': ('Sun, middle of\nconvection zone', 7.4e3, 0.012, 'center'),
    'RGB star, middle of envelope': ('RGB star, middle\nof envelope', 4.4e4, 2.2e-4, 'center'),
    'AGB star, base of envelope': ('AGB star, base\nof envelope', 9e4, 1.2, 'center'),
    'ULIRG nucleus': ('ULIRG nucleus', 0.02, 7.7e2, 'left'),
    'Massive-star envelope': ('massive-star envelope', 2e-3, 1.2e3, 'right'),
    'AGN dusty torus': ('AGN dusty torus', 0.02, 1e4, 'left'),
}
for r in read('rsla_accuracy_systems.csv'):
    name, x, y = r['system'], float(r['tau_v_c']), float(r['h'])
    ax.plot([x], [y], 'o', ms=6, color=INK)
    txt, xx, yy, ha = LABEL[name]
    ax.text(xx, yy, txt, ha=ha, va='center', fontsize=8, color=INK, linespacing=1.0)

ax.set(xscale='log', yscale='log', xlim=XLIM, ylim=YLIM,
       xlabel=r'dynamic-diffusion number $\mathcal{D} = \mathcal{M}\mathcal{L}/\mathcal{C}$',
       ylabel=r'radiation-to-gas pressure ratio $\mathcal{P}$')
ax.set_xticks(10.0 ** np.arange(-8, 8, 2))
# the data and regions are in h; the axis is P = 3 h / 8, i.e. h = 4 (gamma - 1) P with gamma = 5/3
P_PER_H = 3 / 8
ax.set_yticks(10.0 ** np.arange(-6, 6, 2) / P_PER_H, labels=[rf'$10^{{{n}}}$' for n in range(-6, 6, 2)])
ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(10.0 ** np.arange(-7, 5) / P_PER_H))
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

save(fig, 'rsla_accuracy')
