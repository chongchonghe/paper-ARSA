"""Checks how sharply the reduction factor of the ARSA varies between neighbouring cells in the tests (paper Sec. 4:
every run stays stable although phi jumps between cells). Reads the phi of the ARSA at the final time from
data/mf_shocks_profiles.csv (column phi_arsa; tests 1-3 of sim_mf_shocks.jl) and data/marshak_vaytet_profiles.csv
(columns phi_g1-phi_g4, run arsa; phi0 = 0.01). Prints, for each, the range of phi and the largest ratio of phi
between neighbouring cells.

    python checks/check_phi_jumps.py
"""
import csv
from pathlib import Path
import numpy as np

DATA = Path(__file__).resolve().parents[1] / 'data'


def jumps(p):
    p = np.asarray(p)
    return p.min(), p.max(), np.max(np.maximum(p[1:] / p[:-1], p[:-1] / p[1:]))


rows = list(csv.DictReader(open(DATA / 'mf_shocks_profiles.csv')))
phi0 = {r['test']: r['phi0'] for r in csv.DictReader(open(DATA / 'mf_shocks.csv'))}
for t in sorted({r['test'] for r in rows}):
    lo, hi, j = jumps([float(r['phi_arsa']) for r in rows if r['test'] == t])
    print(f'shock test {t} (phi0 = {float(phi0[t]):.3g}): phi in [{lo:.3g}, {hi:.3g}], '
          f'largest ratio between neighbouring cells {j:.0f}')
rows = [r for r in csv.DictReader(open(DATA / 'marshak_vaytet_profiles.csv')) if r['run'] == 'arsa']
for g in range(1, 5):
    lo, hi, j = jumps([float(r[f'phi_g{g}']) for r in rows])
    print(f'Marshak wave, group {g} (phi0 = 0.01): phi in [{lo:.3g}, {hi:.3g}], '
          f'largest ratio between neighbouring cells {j:.0f}')
