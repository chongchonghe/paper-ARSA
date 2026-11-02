"""Fig. B.10 of the paper (fig:theta): stability limit of each Fourier mode, normalized by that of the
checkerboard, from data/theta_curves.csv (written by calc_theta_scan.py). Writes figures/theta_scan.pdf."""
import pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import STYLE, TEXTWIDTH, DATA, save

LINES = [('#0072B2', '-'), ('#D55E00', '--'), ('#009E73', '-.'), ('0.15', ':')]   # Okabe-Ito, with dashes

plt.rcParams.update(STYLE)
d = pd.read_csv(DATA / 'theta_curves.csv')
fig, ax = plt.subplots(figsize=(0.6 * TEXTWIDTH, 2.9))
for (key, g), (col, ls) in zip(d.groupby(['tau', 'Lam'], sort=False), LINES):
    tau, Lam = key        # Lam = 4 (gamma - 1) P nu_c = h nu_c
    ax.semilogy(g.theta, g.nu_crit / g.nu_pi, ls=ls, lw=1.5, color=col,
                label=rf'$\tau_{{\rm c}} = {tau:g}$, $\hat{{h}}\nu_c = {Lam:g}$')
ax.axhline(1, color='0.45', lw=0.8, ls=':')
ax.set_xlim(0, 1)
ax.set_ylim(0.7, 1.5e5)
ax.set_xlabel(r'phase change per cell $\theta/\pi$')
ax.set_ylabel(r'$\nu_{\rm crit}(\theta)/\nu_\pi$')
ax.grid(True, which='major', color='0.9', lw=0.5)
ax.legend(loc='upper right', handlelength=2.6)
save(fig, 'theta_scan')
