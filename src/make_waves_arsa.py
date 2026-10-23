"""Figure 4 (fig:waves_IRSLA): linear waves under the ARSA at the RSLA time step, phi0 = 1e-3,
in the style of Fig. 2 (make_waves_rsla.py: Computer Modern fonts, the same colors, panels without
gaps, gray bands tau > phi0 C and tau > C). Four cases (2 x 2): P = 1e-4 at M = 0 and 5 (h = 0.27 phi0,
phi_* < phi0: thick waves reach phi = 1); P = 5e-4, M = 0 (h = 1.3 phi0, phi_* just above phi0:
phi saturates at 6 phi0; the case of Fig. 2); P = 10, M = 0 (phi_* ~ 1 >> phi0: phi = phi0).
Each case has three stacked panels: phase speed relative to the flow, Re(Omega) - M, of all five modes on a
symmetric-log axis (linear within +-1), broken at +-15 to +-3e3, decay rate -Im(Omega)
of the forward modes, and the phi of eq:phi_waves. Modes of the reduced equations are matched to the exact ones by
eigenvector (see matched); for tau >= 1e3 the eigenvalues are computed at 40 digits (hp_eig),
since the matter decay rates there are many orders of magnitude below the radiation ones."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import itertools
import mpmath as mp
from waves import phi_sb, tracked_modes
import rhd_waves
from paper_style import STYLE, TEXTWIDTH, save
PERMS = list(itertools.permutations(range(5)))

mp.mp.dps = 40


def hp_eig(tau, M, C, P, gamma, phi):
    """Eigenvalues and eigenvectors of (A + iS), rhd_waves.eigensystem, at 40 digits for tau >= 1e3."""
    if tau < 1e3:
        return np.linalg.eig(rhd_waves.advection_matrix(M, phi, C, gamma)
                             + 1j * rhd_waves.source_matrix(M, phi, C, P, gamma, tau))
    f = mp.mpf
    tau, M, C, P, phi, g1 = f(tau), f(M), f(C), f(P), f(phi), f(gamma) - 1
    D = 3 * C**2 - M**2
    A = mp.matrix([[M, 1, 0, 0, 0], [1, M, 1, 0, 0], [0, g1, M, 0, 0],
                   [0, 0, 0, 0, phi * C], [0, 0, 0, phi * C / 3, 0]])
    S = mp.matrix([[0, 0, 0, 0, 0],
                   [0, -4 * C * P / D, -4 * P * M / C, -P * M / (3 * C), P],
                   [0, 0, -4 * (C**2 - M**2) * P * g1 / C, (3 * C**2 + M**2) * P * g1 / (3 * C), -2 * P * M * g1],
                   [0, 4 * phi * C * M / D, 4 * phi * C, -phi * C, phi * M],
                   [0, 4 * phi * C**2 / D, 4 * phi * M, phi * M / 3, -phi * C]])
    lam, vec = mp.eig(A + 1j * tau * S)
    return (np.array([complex(x) for x in lam]),
            np.array([[complex(vec[i, j]) for j in range(5)] for i in range(5)]))


rhd_waves.eigensystem = hp_eig       # also used by waves.tracked_modes
eigensystem = hp_eig

C, gamma, phi0, N, Cr = 1e4, 5 / 3, 1e-3, 32, 0.9 * np.sqrt(3)   # Cr = nu0: 0.9 x the largest stable RSLA step for sigma = 1/sqrt(3)
COL = dict(radiation='C0', acoustic='C1', entropy='C2')
plt.rcParams.update({**STYLE, 'figure.constrained_layout.use': False})   # keep the stacked panels joined (hspace=0)
taus = np.logspace(-6, 6, 361)
sub = slice(0, None, 10)
CASES = ((1e-4, 0.0), (1e-4, 5.0), (5e-4, 0.0), (10.0, 0.0))     # (P, M)


def exact_modes(M, P):
    """Exact (phi = 1) modes tracked by eigenvalue continuity, with eigenvectors."""
    tm = tracked_modes(taus, M, C, P, gamma, lambda t: 1.0)
    keys = ['rad+', 'rad-', 'ac+', 'ac-', 'ent']
    vecs = []
    for i, t in enumerate(taus):
        lam, vec = eigensystem(t, M, C, P, gamma, 1.0)
        idx = [int(np.argmin(np.abs(lam - tm[k][i]))) for k in keys]
        vecs.append(vec[:, idx])
    return tm, np.array(vecs)


def matched(M, P, phif, vex, w=(10.0, 0.2, 0.2, 0.2, 0.2)):
    """Modes of the reduced system, each matched at every tau to the exact mode whose
    eigenvector it overlaps most (so a mode is compared with the same physical mode,
    even where the reduced and exact eigenvalues follow different branches). A small
    continuity term, w[k] per decade of change from the previous tau, keeps each mode on
    one branch: where the reduced radiation pair splits into two purely damped modes,
    both overlap the exact radiation mode similarly, and the match would otherwise flip
    between them. The plotted forward radiation mode (first key) therefore gets a strong weight;
    all other modes a weak one, so that their match stays set by the eigenvectors."""
    keys = ['rad+', 'rad-', 'ac+', 'ac-', 'ent']
    out = {k: np.zeros(len(taus), complex) for k in keys}
    prev = None
    for i, t in enumerate(taus):
        lam, vec = eigensystem(t, M, C, P, gamma, phif(t))
        a = vex[i] / np.linalg.norm(vex[i], axis=0)
        b = vec / np.linalg.norm(vec, axis=0)
        ov = np.abs(a.conj().T @ b)                       # ov[k, j]: exact k vs reduced j
        if prev is None:
            jump = np.zeros((5, 5))
        else:
            with np.errstate(divide='ignore'):
                jump = np.abs(np.log10(np.maximum(np.abs(lam)[None, :], 1e-300)
                                       / np.maximum(np.abs(prev)[:, None], 1e-300)))
        best = max(PERMS, key=lambda p: sum(ov[k, p[k]] - w[k] * jump[k, p[k]] for k in range(5)))
        prev = np.array([lam[best[k]] for k in range(5)])
        for k, key in enumerate(keys):
            out[key][i] = lam[best[k]]
    return out


def pick(tm):
    return dict(rad=tm['rad+'], ac=tm['ac+'], ent=tm['ent'])


def speed(z, key):
    """Phase speed. The radiation mode is shown only where it propagates (underdamped,
    |Im| < |Re|); in thick waves it is a strongly damped diffusive mode with no
    meaningful phase speed. The matter modes are shown everywhere."""
    re = np.array(z.real, float)
    if key.startswith('rad'):
        re[np.abs(z.imag) >= np.abs(z.real)] = np.nan
    return re


def broken(y):
    """Break a line where the matched mode switches branch (jump by > 10x between
    neighbouring taus), so that no vertical connector is drawn."""
    y = np.array(y, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        jump = np.abs(np.diff(np.log10(np.abs(y)))) > 1
    y[1:][jump] = np.nan
    return y


def ptex(P):
    if 1 <= P < 1e3:
        return '%g' % P
    m, e = ('%.0e' % P).split('e')
    return (r'10^{%d}' % int(e)) if m == '1' else (r'%s\times10^{%d}' % (m, int(e)))


fig = plt.figure(figsize=(0.75 * TEXTWIDTH, 7.0))
sfs = fig.subfigures(2, 2, hspace=0, wspace=0).ravel()
for jc, (sf, (P, M)) in enumerate(zip(sfs, CASES)):
    phif = lambda t: phi_sb(t, P, phi0, gamma, N, Cr)
    tex, vex = exact_modes(M, P)
    rs_all, ir_all = matched(M, P, lambda t: phi0, vex), matched(M, P, phif, vex)
    ex, rs, ir = pick(tex), pick(rs_all), pick(ir_all)
    phi = np.array([phif(t) for t in taus])
    # the phase-speed panel is symmetric-log (linear within +-1) with the y axis broken at
    # +-15 to +-3e3: exact forward radiation above, exact backward radiation below, the rest between
    outer = sf.add_gridspec(2, 1, height_ratios=[1.6, 1.55], hspace=0, left=0.22, right=0.97, top=0.94, bottom=0.11)
    top = outer[0].subgridspec(3, 1, height_ratios=[0.25, 1.1, 0.25], hspace=0.07)
    bot = outer[1].subgridspec(2, 1, height_ratios=[1, 0.55], hspace=0)
    a2 = sf.add_subplot(bot[1])
    a1 = sf.add_subplot(bot[0], sharex=a2)
    a0d = sf.add_subplot(top[2], sharex=a2)
    a0 = sf.add_subplot(top[1], sharex=a2)
    a0u = sf.add_subplot(top[0], sharex=a2)
    for ax in (a0u, a0, a0d, a1):
        ax.tick_params(labelbottom=False)
    for ax in (a0u, a0, a0d, a1, a2):
        ax.axvspan(phi0 * C, taus[-1], color='0.86', lw=0)
        ax.axvspan(C, taus[-1], color='0.72', lw=0)
    for key in ('rad+', 'rad-', 'ac+', 'ac-', 'ent'):   # phase speed relative to the flow, all five modes
        c = COL[{'r': 'radiation', 'a': 'acoustic', 'e': 'entropy'}[key[0]]]
        for ax in (a0u, a0, a0d):
            ax.plot(taus, speed(tex[key], key) - M, '-', color=c, lw=1.5)
            ax.plot(taus, speed(rs_all[key], key) - M, '--', color=c, lw=1.5)
            ax.plot(taus[sub], (speed(ir_all[key], key) - M)[sub], 'o', ms=3.5, mfc='none', mec=c, mew=1)
    for key, name in (('rad', 'radiation'), ('ac', 'acoustic'), ('ent', 'entropy')):
        c = COL[name]
        a1.plot(taus, -ex[key].imag, '-', color=c, lw=1.5)
        a1.plot(taus, broken(-rs[key].imag), '--', color=c, lw=1.5)
        a1.plot(taus[sub], -ir[key].imag[sub], 'o', ms=3.5, mfc='none', mec=c, mew=1)
    a2.plot(taus, phi, '-', color='k', lw=1.5)
    a2.axhline(phi0, color='k', ls=':', lw=1)
    g = phi[-1] / phi0
    a2.text(0.03, 0.78, r'$\phi/\phi_0 \in [1, %s]$' % ('%.0f' % g if g >= 10 else '%.2g' % g),
            transform=a2.transAxes, ha='left', va='center')
    a2.set_xscale('log')
    for ax in (a1, a2):
        ax.set_yscale('log')
    for ax in (a0u, a0, a0d):
        ax.set_yscale('symlog', linthresh=1, linscale=1)
        ax.minorticks_off()
    a0u.set_ylim(3e3, 1e4); a0.set_ylim(-15, 15); a0d.set_ylim(-1e4, -3e3)
    a0u.set_yticks([1e4]); a0u.set_yticklabels([r'$10^4$'])
    a0.set_yticks([-10, -1, 0, 1, 10]); a0.set_yticklabels([r'$-10$', r'$-1$', r'$0$', r'$1$', r'$10$'])
    a0d.set_yticks([-1e4]); a0d.set_yticklabels([r'$-10^4$'], va='bottom')
    a0u.set_yticklabels([r'$10^4$'], va='top')
    a0u.spines['bottom'].set_visible(False); a0.spines['top'].set_visible(False)
    a0.spines['bottom'].set_visible(False); a0d.spines['top'].set_visible(False)
    a0u.tick_params(bottom=False); a0.tick_params(top=False, bottom=False); a0d.tick_params(top=False)
    kw = dict(color='k', clip_on=False, lw=1)        # break marks
    for x in (0, 1):
        a0u.plot([x - 0.015, x + 0.015], [-0.07, 0.07], transform=a0u.transAxes, **kw)
        a0.plot([x - 0.015, x + 0.015], [0.985, 1.015], transform=a0.transAxes, **kw)
        a0.plot([x - 0.015, x + 0.015], [-0.015, 0.015], transform=a0.transAxes, **kw)
        a0d.plot([x - 0.015, x + 0.015], [0.93, 1.07], transform=a0d.transAxes, **kw)
    dmin = min(np.nanmin(-ex[k].imag) for k in ('ac', 'ent'))
    a1.set_ylim(10 ** np.floor(np.log10(dmin)), 1e6)
    a2.set_ylim(3e-4, 3)
    a2.set_yticks([1e-3, 1e-1])                    # no label at the join with the panel above
    a2.set_xlim(taus[0], taus[-1])
    a2.set_xticks([1e-6, 1e-3, 1, 1e3, 1e6])
    a2.set_xlabel(r'$\tau$')
    a0u.set_title(r'$\mathcal{P} = %s$, $\mathcal{M} = %g$' % (ptex(P), M))
    if jc % 2 == 0:
        a0.set_ylabel(r'$\mathrm{Re}(\hat\Omega) - \mathcal{M}$')
        a1.set_ylabel(r'$-\mathrm{Im}(\hat\Omega)$')
        a2.set_ylabel(r'$\phi$')
    if jc == 0:
        for x, txt in ((np.sqrt(phi0 * C * C), r'$\tau > \phi_0\mathcal{C}$'), (np.sqrt(C * taus[-1]), r'$\tau > \mathcal{C}$')):
            a0.text(x, 4, txt, ha='center', va='center', fontsize=7.5)   # in the bands, nothing between 1 and 10 here
        h1 = [Line2D([], [], color=COL[n], lw=1.5, label=l)
              for n, l in (('radiation', 'Radiation'), ('acoustic', 'Acoustic'), ('entropy', 'Entropy'))]
        h2 = [Line2D([], [], color='k', lw=1.5, ls='-', label='Exact'),
              Line2D([], [], color='k', lw=1.5, ls='--', label='RSLA'),
              Line2D([], [], color='k', ls='none', marker='o', ms=3.5, mfc='none', label='ARSA')]

fig.legend(handles=h1 + h2, loc='lower center', ncol=6, bbox_to_anchor=(0.5, 0.995), handlelength=1.6, columnspacing=0.9)
save(fig, 'waves_ARSA')
