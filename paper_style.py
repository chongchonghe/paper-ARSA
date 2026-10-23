"""Matplotlib style of the paper's figures: Computer Modern through LaTeX, as in the text, with 9 pt labels
and 8 pt ticks and legends. Make each figure at the width it has in the paper (a fraction of TEXTWIDTH),
so that LaTeX prints it at its natural size and these font sizes hold.

    plt.rcParams.update(STYLE)   or   with plt.rc_context(STYLE): ...
"""
TEXTWIDTH = 500.484 / 72.27          # \textwidth of the paper (from its main.log) [in]

STYLE = {
    'text.usetex': True, 'font.family': 'serif', 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8, 'legend.frameon': False,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'figure.constrained_layout.use': True, 'savefig.pad_inches': 0.02,
}
