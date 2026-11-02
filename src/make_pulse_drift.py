"""Fig. 7 of the paper (fig:pulse_drift): the advected pulse after K crossings of the periodic box, one row per
phi0/P at the peak (1 to 1000), from data/pulse_drift.json (written by sim_pulse_drift.jl). The x-axis is the real
(unwrapped) position, broken into windows of the same width, 2 HALF pulse widths, so that the pulse has the same
width on the page in every window: the start (the initial pulse, dashed), and at the end the reference (chat = c),
the ARSA and the RSLA, in one window when the RSLA lags by less than 2 HALF w and in two otherwise. Where the
ARSA at phi0 = 0.01 saturates (phi0/P = 1), the row shows the ARSA at the higher phi0 of the run arsa_hi.
Writes figures/pulse_drift.pdf.

    python make_pulse_drift.py
"""
import json
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter
from paper_style import STYLE, TEXTWIDTH, C1, C2, DATA, save

K = 12                       # the crossing shown
HALF = 6                     # half width of a window [w]
ROW_H = 1.8                  # height of a row [in]

d = json.load(open(DATA / 'pulse_drift.json'))
w, L, xs = d['w'], d['L'], np.array(d['x'])
runs = d['runs']


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
    fig = plt.figure(figsize=(TEXTWIDTH, len(runs) * ROW_H + 0.55))
    fig.set_layout_engine('none')
    gs = fig.add_gridspec(len(runs), 3, left=0.075, right=0.99, bottom=0.065, top=0.955, wspace=0.06, hspace=0.42)
    for i, r in enumerate(runs):
        T_lo = r['T_lo']
        y = lambda T: (np.array(T) - T_lo) / T_lo
        arsa = 'arsa_hi' if 'x_arsa_hi' in r else 'arsa'
        xr, xl = r['x_ref'][K], r['x_rsla'][K]
        merged = xr - xl < 2 * HALF * w
        start = fig.add_subplot(gs[i, 0])
        ends = [fig.add_subplot(gs[i, 1:])] if merged else [fig.add_subplot(gs[i, 1]), fig.add_subplot(gs[i, 2])]
        centres = [(xr + xl) / 2] if merged else [xl, xr]
        axes = [start] + ends
        fig.canvas.draw()
        cm_per_in = 2 * HALF * w / (start.get_position().width * fig.get_figwidth())
        for ax, c in zip(axes, [0.0] + centres):
            span = ax.get_position().width * fig.get_figwidth() * cm_per_in
            ax.set_xlim(c - span / 2, c + span / 2)
            ax.set_ylim(-0.08, 1.12)
            ax.xaxis.set_major_locator(MultipleLocator(100))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'${v:,.0f}$'.replace(',', r'\,')))
        for a, b in zip(axes[:-1], axes[1:]):
            breaks(a, b)
        start.plot(*place(y(r['T_ref'][0]), 0.0), color='k', ls='--', lw=1.0)
        start.set_ylabel(r'$(T - T_{\rm lo})/T_{\rm lo}$')
        start.set_yticks([0, 0.5, 1])
        start.text(0.04, 0.93, rf'$\phi_0/\mathcal{{P}} = {r["phi0_P"]}$', transform=start.transAxes, va='top')
        for ax in ends:
            ax.plot(*place(y(r['T_ref'][K]), xr), color='k', lw=1.8)
            ax.plot(*place(y(r[f'T_{arsa}'][K]), r[f'x_{arsa}'][K]), color=C2, lw=0.9)
            ax.plot(*place(y(r['T_rsla'][K]), xl), color=C1, lw=1.2)
        lag = r['lag_rsla_w'][K]
        ends[0].text(0.03 if not merged else 0.03, 0.93, rf'RSLA lag ${lagtext(lag)}\,w$',
                     transform=ends[0].transAxes, va='top', color=C1)
        if arsa == 'arsa_hi':
            ends[-1].text(0.97, 0.93, rf'ARSA at $\phi_0 = {r["phi0_arsa_hi"]:g}$', transform=ends[-1].transAxes,
                          va='top', ha='right', color=C2)
    fig.text(0.5, 0.008, rf'real position $x$ [cm], after ${K * L / w:.0f}$ widths traveled ($w = {w:g}$ cm)',
             ha='center')
    fig.legend(handles=[Line2D([], [], color='k', ls='--', lw=1.0), Line2D([], [], color='k', lw=1.8),
                        Line2D([], [], color=C2, lw=0.9), Line2D([], [], color=C1, lw=1.2)],
               labels=['initial', r'reference ($\hat{c} = c$)', 'ARSA', 'RSLA'],
               loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.0))
    save(fig, 'pulse_drift')
