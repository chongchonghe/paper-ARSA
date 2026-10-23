"""Figures of the numerical tests (paper Sec. 5), from the CSV files written by
Honghu.jl's examples/irsla_paper.jl and examples/linear_waves_rhd.jl."""
import csv, sys, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from pdars import growth
from paper_style import STYLE, TEXTWIDTH

FIG = 'figures/'
C1, C2, C3, C4, C5 = '#1f4e9c', '#ef8a2c', '#1baf7a', '#eda100', '#e87ba4'   # C1 dark, C2 light: apart in grayscale
REF = '0.35'

def read(name):
    with open('data/' + name) as f:
        rows = list(csv.DictReader(f))
    return rows

def col(rows, k, typ=float):
    return np.array([typ(r[k]) for r in rows])

def save(fig, name):
    fig.savefig(FIG + name + '.pdf', bbox_inches='tight')
    fig.savefig(FIG + name + '.png', dpi=110, bbox_inches='tight')
    print('wrote', name)

# --------------------------------------------------------------- checkerboard
def energy_block_g(phi, cdt, sdx_over_sig, chi, P, gamma=5 / 3):
    """Largest |eigenvalue| of the PD-ARS amplification of the checkerboard's
    energy block (T, E): the mode an E checkerboard excites."""
    sig = 1 / np.sqrt(3)
    a = 2 * sig * phi * cdt * sdx_over_sig          # transport: E -> (1 - a) E
    q = phi * cdt * chi                              # E damping per step
    lam = 4 * (gamma - 1) * P / phi                  # heat-capacity ratio / phi
    I = np.eye(2)
    Q = q * np.array([[-lam, lam / 4], [4, -1]])
    Z = np.diag([0.0, -a])
    E = I + Z
    G1 = np.linalg.solve(I - Q, E)
    G = np.linalg.solve(I - Q / 2, I / 2 + E @ G1 / 2)
    ev = np.linalg.eigvals(G)
    # the neutral thermal mode (|g| -> 1 as P -> 0) is not excited by an E
    # checkerboard at small P; report the radiation mode, the one that moves
    return max(abs(ev[np.argmax(np.abs(ev - 1))]), abs(ev[np.argmin(np.abs(ev - 1))])) \
        if P > 1e-3 else abs(ev[np.argmax(np.abs(ev - 1))])

@plt.rc_context(STYLE)
def checkerboard():
    rows = read('checkerboard.csv')
    fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH, 2.4), sharey=True)
    cases = [(0.3, 1e-6), (1.0, 1e-6), (1.0, 1e-1)]
    DIMS = ((1, C1, 'o', 6), (2, C2, 's', 4.5), (3, C3, '^', 3.5))   # nested symbols: coincident points stay visible
    for ax, (tau, P) in zip(axes, cases):
        for D, colr, mk, ms in DIMS:
            sel = [r for r in rows if int(r['D']) == D and float(r['tau']) == tau and
                   abs(float(r['Pcal']) - P) < 1e-12 and r['limiter'] == 'fixed']
            fac = col(sel, 'fac'); g = col(sel, 'g')
            ax.plot(fac, g, mk, ms=ms, mfc='none', mec=colr, mew=1)
            # theory: the checkerboard in D dimensions is the 1D theta = pi mode with
            # tau -> tau/D and the CFL numbers summed over directions
            n, phi0 = 16, 1e-3
            dtc = 0.4 / (phi0 * D * n)            # c dt, with dx = 1/n
            nuc = dtc * D * n
            chi = tau * n
            K = 4 * P * (2 / 3) * dtc * chi
            phib = 1 / (dtc * (D * n / np.sqrt(3) - chi / (2 + K)))
            ff = np.linspace(0.8, 1.25, 120)
            gt = [energy_block_g(f * phib, dtc, D * n, chi, P) for f in ff]
            ax.plot(ff, gt, '-', color=colr, lw=1.2)
            # the earlier limiter, placed at its own phi / phi_b
            e = [r for r in rows if int(r['D']) == D and float(r['tau']) == tau and
                 abs(float(r['Pcal']) - P) < 1e-12 and r['limiter'] == 'earlier']
            if e and float(e[0]['g']) > 1.0 + 1e-6:
                ax.plot([], [])
        ax.axhline(1, color='0.5', lw=0.8, ls=':')
        ax.axvline(1, color='0.5', lw=0.8, ls=':')
        ax.set_title(rf'$\tau = {tau:g}$, $\mathcal{{P}} = 10^{{{int(np.log10(P))}}}$')
        ax.set_xlabel(r'$\phi/\phi_{\max}$')
        ax.set_ylim(0.75, 1.3)
        ax.grid(True, color='0.92', lw=0.5)
    axes[0].set_ylabel('growth factor per sub-step')
    h = [Line2D([], [], color=c, lw=1.2, marker=mk, ms=ms, mfc='none', mew=1) for D, c, mk, ms in DIMS]
    h += [Line2D([], [], color='0.3', lw=1.2),
          tuple(Line2D([], [], color='0.3', ls='none', marker=mk, ms=ms, mfc='none', mew=1) for D, c, mk, ms in DIMS)]
    fig.legend(h, [f'{D}D' for D, *_ in DIMS] + ['von Neumann theory', 'Honghu.jl'], loc='outside lower center',
               ncol=5, handler_map={tuple: HandlerTuple(ndivide=None)})
    save(fig, 'checkerboard')

# ------------------------------------------------------ linear waves in code
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

# ------------------------------------------------------- radiative shocks (MF21)
@plt.rc_context(STYLE)
def mf_shocks():
    """Two figures, T_gas and T_rad, of the radiative shocks of Melon Fuksman et al. (2021) on
    64 cells (data/mf_shocks*.csv, from mf_shocks.jl). Left column: test 1 (subcritical). Right
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


if __name__ == '__main__':
    which = sys.argv[1:] or ['checkerboard', 'linear_waves', 'mf_shocks']
    for w in which:
        globals()[w]()
