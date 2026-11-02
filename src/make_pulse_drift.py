"""Fig. 7 of the paper (fig:pulse_drift): the error of the advected pulse against the number of crossings of the box, for
the RSLA and the ARSA at phi0 = 0.01, with h/phi0 at the peak of the pulse from 1e-3 to 0.3, from data/pulse_drift.csv
(written by sim_pulse_drift.jl). Writes figures/pulse_drift.pdf.

    python make_pulse_drift.py
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from paper_style import STYLE, TEXTWIDTH, C2, read, col, save

rows = read('pulse_drift.csv')
cases = sorted(set(col(rows, 'h_phi0_peak')))
blues = plt.get_cmap('Blues')(np.linspace(0.35, 1.0, len(cases)))     # RSLA: light to dark as h/phi0 grows

with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(0.6 * TEXTWIDTH, 2.7))
    for x, color in zip(cases, blues):
        sel = [r for r in rows if float(r['h_phi0_peak']) == x]
        k = col(sel, 'crossing')
        ax.semilogy(k, 100 * col(sel, 'L1_T_rsla'), 'o-', color=color, ms=3, lw=1.3,
                    label=rf'${x:g}$' + (' (amber)' if x > 0.1 else ''))
        ax.semilogy(k, 100 * col(sel, 'L1_T_arsa'), '-', color=C2, lw=0.9, alpha=0.8)
    ax.set(xlabel='crossings of the box', ylabel=r'$L_1$ error relative to the pulse [\%]', xlim=(0.7, 8.3))
    ax.grid(True, which='major', color='0.92', lw=0.5)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color=C2, lw=0.9)); labels.append('ARSA, every $\\hat{h}/\\phi_0$')
    ax.legend(handles, labels, title=r'RSLA, $\hat{h}/\phi_0$ at the peak', loc='center left',
              bbox_to_anchor=(1.02, 0.5), handlelength=2.2)
    save(fig, 'pulse_drift')
