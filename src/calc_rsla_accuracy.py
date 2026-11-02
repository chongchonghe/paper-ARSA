"""When the RSLA is accurate and where the ARSA can help: the systems of Table 1 (tab:media) in the plane of
the dynamic-diffusion number tau v/c and the heat-capacity ratio h, for the smallest reduction factor of the
RSLA considered, phi = 1e-3 (PHI0; the worst case: any larger phi is more accurate). Curves and regions:
  Pe_phi = 0.1, where Pe_phi = 3 (tau v/c) (1/h + 1/phi) (eq:Pe_phi) is the Peclet number of the heat transport
      under the RSLA: the thermal diffusion time of a structure (eq:T_diffusion: t_phi = 3 tau l (1/h + 1/phi)/c)
      over the time the flow takes to cross it; left of the curve even the slowed diffusion is ten times faster
      than the flow, and the temperature tracks the flow, for any phi >= 1e-3; whether the RSLA is accurate
      there depends on the problem (heating faster than the flow, or evolution set by diffusion itself);
  right of it, by the criterion phi >> max(P, D) of eq:rsla_valid, with h for P and the same factor 10 for both:
      h < 0.1 phi = 1e-4 and tau v/c < 0.1 phi: the RSLA is accurate: its slowdown (1 + h/phi)/(1 + h)
          (eq:slowdown), about h/phi, is below 10%, and structures carried by the flow do not drift (Sec. 2.2);
      h < 0.1 phi and tau v/c > 0.1 phi: such structures drift; the ARSA can help, as below;
      1e-4 < h < 0.1: the ARSA can help: run with phi0 above h (up to 0.1), it runs thick cells at phi = 1;
      h > 0.1: the RSLA is inaccurate for any useful phi (by eq:slowdown, 10% accuracy needs phi > 0.48 at h = 0.1).
Optically thin systems (tau <= 1, flagged by the column tau) lie outside the diffusion argument: there the
RSLA delays radiation by the light-crossing time l/c_hat instead.
Writes data/rsla_accuracy_systems.csv and data/rsla_accuracy_curves.csv.

    python calc_rsla_accuracy.py
"""
import numpy as np
from paper_style import DATA
from make_table_media import rows                          # the systems of Table 1

PHI0 = 1e-3                                           # the smallest reduction factor phi considered
PE_T = 0.1                                            # threshold of Pe_phi
H_LO, H_HI = 0.1 * PHI0, 0.1                          # edges of the region where the ARSA can help
DRIFT_D = 0.1 * PHI0                                  # tau v/c above which advected structures drift (phi >> D)

h = np.logspace(-8, 8, 400)
with open(DATA / 'rsla_accuracy_curves.csv', 'w') as f:
    f.write('h,tau_v_c_Pe\n')                         # Pe_phi = PE_T: tau v/c = PE_T / (3 (1/h + 1/phi))
    for hh in h:
        f.write(f'{hh:.6g},{PE_T / (3 * (1 / hh + 1 / PHI0)):.6g}\n')

print(f'phi = {PHI0}; Pe_phi = {PE_T} at tau v/c = {PE_T}/(3 (1/h + 1/phi)): {PE_T / 3:.3g} h for h << phi, '
      f'{PE_T * PHI0 / 3:.3g} for h >> phi; the ARSA can help for {H_LO:g} < h < {H_HI:g}, and for h < {H_LO:g} '
      f'where tau v/c > 0.1 phi = {DRIFT_D:.3g}')
print(f'{"system":30s} {"tau v/c":>9s} {"h":>9s} {"Pe_phi":>9s} {"t_phi/t_1":>10s}  region')
with open(DATA / 'rsla_accuracy_systems.csv', 'w') as f:
    f.write('system,tau,tau_v_c,h,Pe_phi,slowdown,region\n')
    for r in rows():
        hh, tvc = r['h'], r['tvc']
        slow = (1 + hh / PHI0) / (1 + hh)
        if tvc is None:                                # the solar core (v ~ 0) and the IGM (tau << 1)
            tvc, pp = float('nan'), float('nan')
        else:
            pp = 3 * tvc * (1 / hh + 1 / PHI0)
        region = ('no flow / thin' if np.isnan(pp) else 'fast diffusion' if pp < PE_T else
                  'RSLA accurate' if hh < H_LO and tvc < DRIFT_D else 'ARSA can help' if hh < H_HI else 'RSLA inaccurate')
        tau = float('nan') if r['tau'] is None else r['tau']
        f.write(f'"{r["name"]}",{tau:.4g},{tvc:.4g},{hh:.4g},{pp:.4g},{slow:.4g},{region}\n')
        print(f'{r["name"]:30s} {tvc:9.2g} {hh:9.2g} {pp:9.2g} {slow:10.3g}  {region}')
print('wrote data/rsla_accuracy_curves.csv, data/rsla_accuracy_systems.csv')
