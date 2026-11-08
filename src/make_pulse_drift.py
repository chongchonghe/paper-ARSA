"""Fig. of the paper (fig:pulse_drift): the advected pulse after WIDTHS widths traveled, one row per
phi0/P at the peak (3 and 300), from data/pulse_drift.json (written by sim_pulse_drift.jl), and below them, set
apart, the weak pulse (amplitude 0.1, shown divided by it) at rest at phi0/P = 3 at the end of its run, from data/pulse_rest.json (sim_pulse_rest.jl), over its
whole box. (Until 2026-11-08 also phi0/P = 30 and 3000, and the pulse at D = 10 phi0 in the bottom row.)
Each drift row gives the widths traveled at its end. The x-axis is the real
(unwrapped) position, broken into windows of the same width, 2 HALF pulse widths, so that the pulse has the same
width on the page in every window: the start (the initial pulse, dashed), and at the end the exact solution
(phi = 1), the ARSA and the RSLA, in one window when the RSLA lags by less than 2 HALF w (the exact pulse HALF w
from its right edge) and in two otherwise.
Writes figures/pulse_drift.pdf.

    python make_pulse_drift.py
"""
import json
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from paper_style import STYLE, TEXTWIDTH, C1, C2, DATA, save

WIDTHS = 256                 # the widths traveled shown in the rows of pulse_drift.json
HALF = 7                     # half width of a window [w]
SCALE = 0.8                  # the figure is SCALE times the text width, so that its text prints larger
AX_H, GAP = 1.30 * SCALE, 0.23 * SCALE   # height of a panel and gap between rows [in]; twice the gap above the last row
TOP, BOTTOM = 0.37 * SCALE, 0.53 * SCALE # margins for the legend and the axis label [in]
LEFT = 0.49                  # margin for the y label [in]

d = json.load(open(DATA / 'pulse_drift.json'))
w, L, xs = d['w'], d['L'], np.array(d['x'])
dr = json.load(open(DATA / 'pulse_rest.json'))
rows = [(r, r['widths'].index(WIDTHS), rf'$\phi_0/\mathcal{{P}} = {r["phi0_P"]}$') for r in d['runs']]
GRID = [0, 1]                # grid row of each drift row; grid row 2 is the gap above the row of the pulse at rest


def place(T, centre):        # the box profile at the image of each cell nearest the real position centre
    x = xs + L * np.round((centre - xs) / L)
    i = np.argsort(x)
    return x[i], np.array(T)[i]


def breaks(left, right):     # the marks of a broken axis between two neighbouring windows
    left.spines['right'].set_visible(False); right.spines['left'].set_visible(False)
    right.tick_params(axis='y', left=False, labelleft=False)
    kw = dict(color='k', lw=0.6, clip_on=False)
    for ax, x0 in ((left, 1), (right, 0)):
        for y0 in (0, 1):
            ax.plot([x0 - 0.012, x0 + 0.012], [y0 - 0.04, y0 + 0.04], transform=ax.transAxes, **kw)


def lagtext(v):
    return f'{v:.0f}' if v >= 10 else f'{v:.1f}' if v >= 1 else f'{v:.2f}'


with plt.rc_context(STYLE):
    H = TOP + BOTTOM + 3 * AX_H + 3 * GAP
    fig = plt.figure(figsize=(SCALE * TEXTWIDTH, H))
    fig.set_layout_engine('none')
    ratios = [1, 1, 1e-4, 1]                         # grid row 2: no height, so the gap there is 2 GAP
    gs = fig.add_gridspec(4, 3, height_ratios=ratios, left=LEFT / (SCALE * TEXTWIDTH), right=0.99, bottom=BOTTOM / H, top=1 - TOP / H,
                          wspace=0.06, hspace=GAP / (AX_H * sum(ratios) / 4))
    for (r, K, label), i in zip(rows, GRID):
        T_lo = r['T_lo']
        y = lambda T: (np.array(T) - T_lo) / T_lo
        xr, xl = r['x_exact'][K], r['x_rsla'][K]
        merged = xr - xl < 2 * HALF * w
        start = fig.add_subplot(gs[i, 0])
        ends = [fig.add_subplot(gs[i, 1:])] if merged else [fig.add_subplot(gs[i, 1]), fig.add_subplot(gs[i, 2])]
        centres = [None] if merged else [xl, xr]       # None: the exact pulse HALF w from the right edge
        axes = [start] + ends
        fig.canvas.draw()
        cm_per_in = 2 * HALF * w / (start.get_position().width * fig.get_figwidth())
        for ax, c in zip(axes, [0.0] + centres):
            span = ax.get_position().width * fig.get_figwidth() * cm_per_in
            c = xr + HALF * w - span / 2 if c is None else c
            ax.set_xlim(c - span / 2, c + span / 2)
            ax.set_ylim(-0.08, 1.12)
            ax.set_xticks([t for t in np.arange(100 * np.ceil((c - span / 2) / 100), c + span / 2, 100)
                           if abs(t - c) < span / 2 - 20])        # none within 20 cm of the breaks
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'${v:,.0f}$'.replace(',', r'\,')))
        for a, b in zip(axes[:-1], axes[1:]):
            breaks(a, b)
        start.plot(*place(y(r['T_exact'][0]), 0.0), color='k', ls='--', lw=1.0)
        start.set_ylabel(r'$(T - T_{\rm lo})/T_{\rm lo}$')
        start.set_yticks([0, 0.5, 1])
        for ax in ends:
            ax.plot(*place(y(r['T_exact'][K]), xr), color='k', lw=1.8)
            ax.plot(*place(y(r['T_arsa'][K]), r['x_arsa'][K]), color=C2, lw=0.9)
            ax.plot(*place(y(r['T_rsla'][K]), xl), color=C1, lw=1.2)
        lag = r['lag_rsla_w'][K]
        ends[0].text(0.03, 0.93, label + f', after {r["widths"][K]} widths', transform=ends[0].transAxes, va='top')
        ends[0].text(0.03, 0.78, rf'RSLA lag ${lagtext(lag)}\,w$', transform=ends[0].transAxes, va='top', color=C1)
    # the pulse at rest: the initial pulse, and at the end the exact solution, the ARSA and the RSLA, over the box
    T_lo, xr_, amp = dr['T_lo'], np.array(dr['x']), dr['amp']
    y = lambda T: (np.array(T) - T_lo) / (amp * T_lo)
    start, end = fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1:])
    for ax in (start, end):
        ax.set_xlim(-dr['L'] / 2, dr['L'] / 2)
        ax.set_ylim(-0.08, 1.12)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'${v:,.0f}$'.replace(',', r'\,')))
    start.set_xlim(-HALF * w, HALF * w); start.set_xticks([-100, 0, 100])
    end.set_xticks(np.arange(-300, 301, 100))
    end.tick_params(axis='y', labelleft=False)
    start.plot(xr_, y(dr['T_exact'][0]), color='k', ls='--', lw=1.0)
    start.set_ylabel(rf'$(T - T_{{\rm lo}})/({amp:g}\,T_{{\rm lo}})$')
    start.set_yticks([0, 0.5, 1])
    end.plot(xr_, y(dr['T_exact'][-1]), color='k', lw=1.8)
    end.plot(xr_, y(dr['T_arsa'][-1]), color=C2, lw=0.9)
    end.plot(xr_, y(dr['T_rsla'][-1]), color=C1, lw=1.2)
    end.text(0.03, 0.93, rf'weak pulse at rest, $\phi_0/\mathcal{{P}} = {dr["phi0_P"]}$, $t = {dr["t"][-1] * 1e3:.1f}$ ms',
             transform=end.transAxes, va='top')
    ex = np.array(dr['T_exact'][-1]) - T_lo
    l1 = lambda T: np.sum(np.abs(np.array(T) - T_lo - ex)) / np.sum(np.abs(ex))
    end.text(0.03, 0.78, rf'$L_1$ error of $T - T_{{\rm lo}}$: RSLA ${100 * l1(dr["T_rsla"][-1]):.0f}\%$',
             transform=end.transAxes, va='top', color=C1)
    end.text(0.03, 0.63, rf'ARSA ${100 * l1(dr["T_arsa"][-1]):.2f}\%$', transform=end.transAxes, va='top', color=C2)
    fig.text(0.5, 0.008, rf'real position $x$ [cm] (pulse width $w = {w:g}$ cm)',
             ha='center')
    fig.legend(handles=[Line2D([], [], color='k', ls='--', lw=1.0), Line2D([], [], color='k', lw=1.8),
                        Line2D([], [], color=C2, lw=0.9), Line2D([], [], color=C1, lw=1.2)],
               labels=['initial', 'exact', 'ARSA', 'RSLA'],
               loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.0))
    save(fig, 'pulse_drift')
