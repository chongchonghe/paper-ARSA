"""Fig. 5 of the paper (fig:linear_waves_code): phase speed and decay rates of linear waves measured in
Honghu.jl (data/linear_waves_rhd.csv, from Honghu.jl's examples/linear_waves_rhd.jl) against the dispersion
relation, for the exact equations, the RSLA and the ARSA. Writes figures/linear_waves_code.pdf.

    python make_linear_waves.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from paper_style import STYLE, TEXTWIDTH, C1, C2, REF, read, col, save


def linear_waves():
    from waves import tracked_modes
    rows = read('linear_waves_rhd.csv')
    C, P, gamma, phi0, N, Cr, beta = 1e3, 1e-2, 5 / 3, 0.03, 64, 0.9 * np.sqrt(3), 0.9
    sig = 1 / np.sqrt(3)
    def phi_code(tau):
        dx = 1 / N; chi = tau * 2 * np.pi
        cdt = Cr * dx / phi0
        hrat = 4 * P / (1 / (gamma - 1))
        K = hrat * cdt * chi
        D = sig / dx - chi / (2 + K)
        return 1.0 if D <= 0 else min(1.0, max(phi0, beta / (cdt * D)))
    taus = np.logspace(-1.3, 3.3, 181)
    th = {'exact': tracked_modes(taus, 0.0, C, P, gamma, lambda t: 1.0),
          'rsla': tracked_modes(taus, 0.0, C, P, gamma, lambda t: phi0),
          'irsla': tracked_modes(taus, 0.0, C, P, gamma, phi_code)}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH, 2.45))
        sty = {'exact': ('-', REF, 'o'), 'rsla': ('--', C1, 's'), 'irsla': (':', C2, '^')}
        for s, (ls, c, mk) in sty.items():
            ac = th[s]['ac+']; en = th[s]['ent']
            axes[0].semilogx(taus, ac.real, ls, color=c, lw=1.4)
            axes[1].loglog(taus, -ac.imag, ls, color=c, lw=1.4)
            axes[2].loglog(taus, -en.imag, ls, color=c, lw=1.4)
            sa = [r for r in rows if r['scheme'] == s and r['mode'] == 'acoustic']
            se = [r for r in rows if r['scheme'] == s and r['mode'] == 'entropy']
            axes[0].semilogx(col(sa, 'tau'), col(sa, 'Re_meas'), mk, color=c, mfc='none', ms=5, mew=1)
            axes[1].loglog(col(sa, 'tau'), -col(sa, 'Im_meas'), mk, color=c, mfc='none', ms=5, mew=1)
            axes[2].loglog(col(se, 'tau'), -col(se, 'Im_meas'), mk, color=c, mfc='none', ms=5, mew=1)
        for a, title, ylab in zip(axes, ('acoustic phase speed', 'acoustic decay rate', 'entropy decay rate'),
                                  (r'$\mathrm{Re}(\hat\Omega)$', r'$-\mathrm{Im}(\hat\Omega)$',
                                   r'$-\mathrm{Im}(\hat\Omega)$')):
            a.set_title(title); a.set_ylabel(ylab)
            a.set_xlabel(r'$\tau$'); a.grid(True, which='major', color='0.92', lw=0.5)
        h = [Line2D([], [], color=c, ls=ls, lw=1.4, marker=mk, mfc='none', ms=5, mew=1, label=l)
             for (ls, c, mk), l in zip(sty.values(), ('exact ($\\hat{c} = c$)', 'RSLA', 'ARSA'))]
        fig.legend(handles=h, loc='outside lower center', ncol=3, handlelength=3)
        save(fig, 'linear_waves_code')


linear_waves()
