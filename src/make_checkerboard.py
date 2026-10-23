"""Fig. 3 of the paper (fig:checkerboard): growth factor of a checkerboard in 1D, 2D and 3D against
phi/phi_max, measured in Honghu.jl (data/checkerboard.csv, from Honghu.jl's examples/irsla_paper.jl, section
checkerboard) and from the von Neumann analysis of the PD-ARS step. Writes figures/checkerboard.pdf.

    python make_checkerboard.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from paper_style import STYLE, TEXTWIDTH, C1, C2, C3, read, col, save


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


checkerboard()
