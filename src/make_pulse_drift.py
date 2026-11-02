"""Fig. 6 of the paper (fig:pulse_drift): the advected pulse after WIDTHS widths traveled, one row per
phi0/P at the peak (3 to 3000), from data/pulse_drift.json (written by sim_pulse_drift.jl), and below them, set
apart, the pulse at D = 10 phi0 after WIDTHS_LONG widths, from data/pulse_drift_long.json (sim_pulse_drift_long.jl).
Each row gives the widths traveled at its end. The x-axis is the real
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
WIDTHS_LONG = 2048           # and in the row of pulse_drift_long.json
HALF = 7                     # half width of a window [w]
AX_H, GAP = 1.30, 0.23       # height of a panel and gap between rows [in]; twice the gap above the last row
TOP, BOTTOM = 0.37, 0.53     # margins for the legend and the axis label [in]

d = json.load(open(DATA / 'pulse_drift.json'))
w, L, xs = d['w'], d['L'], np.array(d['x'])
dl = json.load(open(DATA / 'pulse_drift_long.json'))
assert dl['w'] == w and dl['L'] == L and len(dl['x']) == len(xs)
rows = [(r, r['widths'].index(WIDTHS), rf'$\phi_0/\mathcal{{P}} = {r["phi0_P"]}$') for r in d['runs']]
rl = dl['runs'][0]
rows.append((rl, rl['widths'].index(WIDTHS_LONG),
             rf'$\phi_0/\mathcal{{P}} = {rl["phi0_P"]}$' + '\n' + rf'$\mathcal{{D}} = {dl["D"] / dl["phi0"]:.0f}\,\phi_0$'))
GRID = [0, 1, 2, 3, 5]       # grid row of each row; grid row 4 is the gap above the last row


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
    H = TOP + BOTTOM + 5 * AX_H + 5 * GAP
    fig = plt.figure(figsize=(TEXTWIDTH, H))
    fig.set_layout_engine('none')
    ratios = [1, 1, 1, 1, 1e-4, 1]                   # grid row 4: no height, so the gap there is 2 GAP
    gs = fig.add_gridspec(6, 3, height_ratios=ratios, left=0.075, right=0.99, bottom=BOTTOM / H, top=1 - TOP / H,
                          wspace=0.06, hspace=GAP / (AX_H * sum(ratios) / 6))
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
        start.text(0.04, 0.93, label, transform=start.transAxes, va='top')
        for ax in ends:
            ax.plot(*place(y(r['T_exact'][K]), xr), color='k', lw=1.8)
            ax.plot(*place(y(r['T_arsa'][K]), r['x_arsa'][K]), color=C2, lw=0.9)
            ax.plot(*place(y(r['T_rsla'][K]), xl), color=C1, lw=1.2)
        lag = r['lag_rsla_w'][K]
        ends[0].text(0.03, 0.93, f'after {r["widths"][K]} widths', transform=ends[0].transAxes, va='top')
        ends[0].text(0.03, 0.78, rf'RSLA lag ${lagtext(lag)}\,w$', transform=ends[0].transAxes, va='top', color=C1)
    fig.text(0.5, 0.008, rf'real position $x$ [cm] (pulse width $w = {w:g}$ cm)',
             ha='center')
    fig.legend(handles=[Line2D([], [], color='k', ls='--', lw=1.0), Line2D([], [], color='k', lw=1.8),
                        Line2D([], [], color=C2, lw=0.9), Line2D([], [], color=C1, lw=1.2)],
               labels=['initial', 'exact', 'ARSA', 'RSLA'],
               loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.0))
    save(fig, 'pulse_drift')
