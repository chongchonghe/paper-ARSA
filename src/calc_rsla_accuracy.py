"""When the RSLA is accurate and where the ARSA can help: the optically thick systems of Table 1 (tab:media) in
the plane of the dynamic-diffusion number D = tau v/c and the heat-capacity ratio h, for the smallest reduction
factor of the RSLA considered, phi = 1e-3 (PHI0; the worst case: any larger phi is more accurate). The regions
follow the criterion phi >> max(P, D) of eq:rsla_valid, with h for P and the same factor 10 for both:
  h < 0.1 phi = 1e-4 and D < 0.1 phi: the RSLA is accurate: its slowdown (1 + h/phi)/(1 + h) (eq:slowdown),
      about h/phi, is below 10%, and structures carried by the flow do not drift (Sec. 2.2);
  h < 0.1 phi and D > 0.1 phi: such structures drift; the ARSA can help, as below;
  1e-4 < h < 0.1, at any D: the ARSA can help: run with phi0 above h (up to 0.1), it runs thick cells at phi = 1;
  h > 0.1, at any D: neither is accurate for any useful phi (by eq:slowdown, 10% accuracy needs phi > 0.48 at
      h = 0.1).
Writes data/rsla_accuracy_systems.csv.

    python calc_rsla_accuracy.py
"""
from paper_style import DATA
from make_table_media import rows                          # the systems of Table 1

PHI0 = 1e-3                                           # the smallest reduction factor phi considered
H_LO, H_HI = 0.1 * PHI0, 0.1                          # edges of the region where the ARSA can help
DRIFT_D = 0.1 * PHI0                                  # D above which advected structures drift (phi >> D)

print(f'phi = {PHI0}; the RSLA is accurate for h < {H_LO:g} and D < {DRIFT_D:g}; the ARSA can help for '
      f'{H_LO:g} < h < {H_HI:g} at any D, and for h < {H_LO:g} where D > {DRIFT_D:g}; neither for h > {H_HI:g}')
print(f'{"system":30s} {"D":>9s} {"h":>9s} {"t_phi/t_1":>10s}  region')
with open(DATA / 'rsla_accuracy_systems.csv', 'w') as f:
    f.write('system,tau,tau_v_c,h,slowdown,region\n')
    for r in rows():
        hh, tvc = r['h'], r['tvc']
        slow = (1 + hh / PHI0) / (1 + hh)
        region = ('RSLA accurate' if hh < H_LO and tvc < DRIFT_D else 'ARSA can help' if hh < H_HI else
                  'neither accurate')
        f.write(f'"{r["name"]}",{r["tau"]:.4g},{tvc:.4g},{hh:.4g},{slow:.4g},{region}\n')
        print(f'{r["name"]:30s} {tvc:9.2g} {hh:9.2g} {slow:10.3g}  {region}')
print('wrote data/rsla_accuracy_systems.csv')
