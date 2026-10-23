"""Figure 2 (fig:waves_RSLA): dimensionless linear wave frequencies of the RHD equations and
of the RSLA, as a function of the optical depth per radian tau, for C = 1e4, M = 0, gamma = 5/3.
Two versions: the paper's, with P = 5e-4 and phi = 1e-3 (the bottom-left case of Fig. 4), written
to figures/waves_RSLA.pdf; and one with the parameters of the original figure, P = 0.1 and
phi = 3e-3, written to figures/waves_RSLA_P0.1.pdf for comparison with it.
Top: Re(Omega) of all five modes (symmetric-log axis); bottom: -Im(Omega). Modes are followed by
eigenvalue continuity from the low-tau end (waves.tracked_modes), for the exact equations (solid)
and the RSLA (dashed). Gray bands: tau > phi C (light) and tau > C (dark). Style of the original
figure: one column, Computer Modern fonts."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from waves import tracked_modes
from paper_style import STYLE

C, M, gamma = 1e4, 0.0, 5 / 3
VERSIONS = (('waves_RSLA', 5e-4, 1e-3),          # the paper's: the bottom-left case of fig:waves_IRSLA
            ('waves_RSLA_P0.1', 0.1, 3e-3))      # the parameters of the original figure
taus = np.logspace(-6, 6, 1201)
COL = dict(rad='C0', ac='C1', ent='C2')
KIND = {'rad+': 'rad', 'rad-': 'rad', 'ac+': 'ac', 'ac-': 'ac', 'ent': 'ent'}
plt.rcParams.update({**STYLE, 'figure.constrained_layout.use': False})   # keep the two panels joined (hspace=0)


def make(name, P, phi):
    ex = tracked_modes(taus, M, C, P, gamma, lambda t: 1.0)
    rs = tracked_modes(taus, M, C, P, gamma, lambda t: phi)

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(3.5, 4.5), sharex=True, gridspec_kw=dict(hspace=0))
    for ax in (a0, a1):
        ax.axvspan(phi * C, taus[-1], color='0.86', lw=0)
        ax.axvspan(C, taus[-1], color='0.72', lw=0)
    for key, kind in KIND.items():
        for modes, ls in ((ex, '-'), (rs, '--')):
            a0.plot(taus, modes[key].real, ls, color=COL[kind], lw=1.5)
            a1.plot(taus, -modes[key].imag, ls, color=COL[kind], lw=1.5)

    a0.set_xscale('log')
    a0.set_yscale('symlog', linthresh=3, linscale=1)
    a0.set_ylim(-1e4, 1e4)
    a0.set_yticks([-1e3, -1e2, -1e1, 0, 1e1, 1e2, 1e3, 1e4])
    a0.set_yticklabels([r'$-10^3$', r'$-10^2$', r'$-10^1$', r'$0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
    a0.set_ylabel(r'$\mathrm{Re}(\hat\Omega)$')
    a1.set_yscale('log')
    a1.set_ylim(1e-4, 1e4)
    a1.set_yticks([1e-4, 1e-2, 1, 1e2])            # no label at the join with the panel above
    a1.set_xlim(taus[0], taus[-1])
    a1.set_xticks([1e-6, 1e-3, 1, 1e3, 1e6])
    a1.set_xlabel(r'$\tau$')
    a1.set_ylabel(r'$-\mathrm{Im}(\hat\Omega)$')
    a1.text((phi * C)**0.15 * C**0.85, 3e2, r'$\tau > \phi\mathcal{C}$', rotation=-90, ha='center', va='center')
    a1.text(np.sqrt(C * taus[-1]), 3e2, r'$\tau > \mathcal{C}$', rotation=-90, ha='center', va='center')

    h1 = [Line2D([], [], color=COL[k], lw=1.5, label=l) for k, l in (('rad', 'Radiation'), ('ac', 'Acoustic'), ('ent', 'Entropy'))]
    h2 = [Line2D([], [], color='k', lw=1.5, ls=ls, label=l) for ls, l in (('-', 'Exact'), ('--', 'RSLA'))]
    leg = a0.legend(handles=h1, loc='upper right', frameon=True, framealpha=0.9)
    a0.add_artist(leg)
    a0.legend(handles=h2, loc='lower right', frameon=True, framealpha=0.9)
    fig.savefig('figures/%s.pdf' % name, bbox_inches='tight')
    fig.savefig('figures/' + name + '.png', dpi=150, bbox_inches='tight')
    print('wrote', name, 'P =', P, 'phi =', phi)

    # numbers quoted in the text, at tau = 1e4: the slow (diffusive) decay, which the exact
    # equations carry in the entropy mode and the RSLA in one radiation branch, and the
    # acoustic phase speed and decay rate
    i = int(np.argmin(abs(taus - 1e4)))
    print('  tau = %.3g: exact entropy decay %.3g; RSLA radiation decays %.3g, %.3g; acoustic speed '
          'exact %.3f, RSLA %.3f; acoustic decay exact %.3g, RSLA %.3g'
          % (taus[i], -ex['ent'][i].imag, -rs['rad+'][i].imag, -rs['rad-'][i].imag,
             ex['ac+'][i].real, rs['ac+'][i].real, -ex['ac+'][i].imag, -rs['ac+'][i].imag))


for v in VERSIONS:
    make(*v)
