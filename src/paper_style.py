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

# ------------------------------------------------- paths and helpers of the figure scripts
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA, FIGURES, TABLES = ROOT / 'data', ROOT / 'figures', ROOT / 'tables'
C1, C2, C3 = '#1f4e9c', '#ef8a2c', '#1baf7a'    # RSLA C1 (blue), ARSA C2 (orange): apart in grayscale; C3 green
REF = '0.35'                                    # the reference runs (chat = c)


def read(name):
    """The rows of data/<name>, a CSV file with a header, as dicts of strings."""
    with open(DATA / name) as f:
        return list(csv.DictReader(f))


def col(rows, k, typ=float):
    return np.array([typ(r[k]) for r in rows])


def save(fig, name):
    """Write figures/<name>.pdf."""
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f'{name}.pdf', bbox_inches='tight')
    print(f'wrote figures/{name}.pdf')
