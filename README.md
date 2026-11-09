# Scripts and data for "An Adaptive Reduced Speed of Light Approximation"

This repository makes every figure, Table 1 and the key numbers of C.-C. He, M. R. Krumholz and B. D. Wibking, "An Adaptive Reduced Speed of Light Approximation" (submitted to the Journal of Computational Physics). The simulations use [Honghu.jl](https://github.com/chongchonghe/Honghu.jl) v0.1.0 ([doi:10.5281/zenodo.23569082](https://doi.org/10.5281/zenodo.23569082)), pinned in `Manifest.toml`.

## Quick start

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # Julia 1.12; installs Honghu.jl v0.1.0 and its dependencies
python3 -m pip install -r requirements.txt            # Python 3.13
./run_all.sh figures    # the figures and Table 1, from the data in this repository: about 20 s
./run_all.sh            # everything, the simulations first: about 25 min with 14 threads
```

`PY=<python> ./run_all.sh` picks another Python. The figures set matplotlib's `text.usetex`, so they need a LaTeX installation. Every script runs from any directory.

## Layout

```
run_all.sh                        the single entry point
src/sim_*.jl                      simulations with Honghu.jl; write data/
src/calc_*.py                     calculations (linear theory, scans); write data/
src/make_*.py                     figures and Table 1, from data/ (or linear theory); write figures/, tables/
src/pdars.py rhd_waves.py waves.py paper_style.py    modules the scripts import
checks/                           scripts that check statements in the text; no output files, except checks/m1/
lean/                             Lean 4 proofs of the stability results (B. Wibking)
data/  figures/  tables/          generated, in git
copy_to_paper.sh                  for the authors: copies figures/ and tables/ into the paper's source
```

## Which script makes what

| Paper | Output | Script | Data (step) |
|---|---|---|---|
| Table 1 | `tables/tab_media.tex` | `make_table_media.py` | estimates in the script; `media_parameters.md`, `convective_layers.txt` |
| Fig. 1 | `rsla_accuracy.pdf` | `make_rsla_accuracy.py` | `rsla_accuracy_systems.csv` (`rsla_accuracy`) |
| Fig. 2 | `waves_RSLA.pdf` | `make_waves_rsla.py` | linear theory |
| Fig. 3 | `waves_ARSA.pdf` | `make_waves_arsa.py` | linear theory |
| Fig. 4 | `linear_waves_code.pdf` | `make_linear_waves.py` | `linear_waves_rhd.csv` (`linear_waves_rhd`) |
| Fig. 5 | `pulse_drift.pdf` | `make_pulse_drift.py` | `pulse_drift.json` (`pulse_drift`), `pulse_rest.json` (`pulse_rest`, bottom row) |
| Fig. 6 | `expanding_sphere.pdf` | `make_expanding_sphere.py` | `expanding_sphere.csv` (`expanding_sphere`) |
| Figs. 7, 8 | `shocks_Tgas.pdf`, `shocks_Trad.pdf` | `make_mf_shocks.py` | `mf_shocks*.csv` (`mf_shocks`) |
| Fig. 9 | `marshak_vaytet.pdf` | `make_marshak_vaytet.py` | `marshak_vaytet_*.csv` (`marshak_vaytet`) |
| Fig. A.10, Table A.2 | `theta_scan.pdf`; Table A.2 typed in | `make_theta_scan.py` | `theta_curves.csv`, `theta_scan.csv` (`theta_scan`) |

Figures are in `figures/`. Two figures are made but not in the current paper: `checkerboard.pdf` (step `checkerboard`) and `waves_RSLA_P0.1.pdf` (Fig. 2 at $\mathcal{P} = 0.1$, $\phi = 3\times10^{-3}$).

## The steps

`./run_all.sh` runs the steps in parallel, each in its own process, then `floor_energy`, then the figures; `./run_all.sh <step> ...` runs only the steps named. Each step writes its printout to `data/<step>.txt`, after a header with the date, the commit of this repository (`-dirty`: uncommitted changes in `src/` or `Manifest.toml`) and the Honghu.jl version. The data in git were made with the commit that their header names. Times are from one `./run_all.sh` on an 18-core Apple laptop.

| Step | Script | Writes (in `data/`) | Time |
|---|---|---|---|
| `linear_waves_rhd` | Honghu.jl `examples/linear_waves_rhd.jl` | `linear_waves_rhd.csv` | 1 min |
| `pulse_drift` | `sim_pulse_drift.jl` (worker processes, one thread each) | `pulse_drift.json` | 3 min |
| `pulse_rest` | `sim_pulse_rest.jl` (two workers) | `pulse_rest.json` | 12 min |
| `expanding_sphere` | `sim_expanding_sphere.jl` (four workers; `THREADS=<n>`) | `expanding_sphere.csv` | 5 min |
| `mf_shocks` | `sim_mf_shocks.jl` | `mf_shocks*.csv` | 10 min |
| `marshak_vaytet` | `sim_marshak_vaytet.jl` | `marshak_vaytet_*.csv` | 3 min |
| `rsla_accuracy`, `convective_layers`, `theta_scan`, `m1_stability` | `calc_*.py` | `rsla_accuracy_systems.csv`, `theta_scan.csv`, `theta_curves.csv`; printouts | 1 s–2 min |
| `m1_checks` | `checks/m1/*.py` | `m1/*.json` | 5 min |
| `table_media`, `waves_rsla` | `make_table_media.py`, `make_waves_rsla.py` | printouts | 5 s |
| `checkerboard`, `arsa_phi0`, `diffusion_pulse`, `floor_energy`, `pulse_drift_substep`, `waves_heat_capacity` | supporting runs, not quoted in the current paper | see each script's header | 10 s–25 min |

## Data files

All files are in `data/`. Units are cgs, temperatures in K, errors are fractions (the paper gives percent). Each `.txt` is the printout of its step.

- **`linear_waves_rhd.csv`**: `tau` (optical depth per radian), `scheme` (`exact`, `rsla`, `irsla` = ARSA), `mode`, `phi`, measured and theoretical `Re`, `Im` of the frequency.
- **`pulse_drift.json`**: the advected pulse, $\phi_0 = 0.01$, $\mathcal{D} = 14.4$. Top level: `phi0`, `D`, `w`, `L`, `v0`, `t_out`, `x`. `runs` (one per $\phi_0/\mathcal{P}$ at the peak, 3 and 300): `T_lo`, `P_peak`, `h_peak`, `arsa_phi_min`, `widths` (0–256), the centroid lags `lag_rsla_w`, `lag_arsa_w` [widths], and the profiles `T_exact`, `T_rsla`, `T_arsa`.
- **`pulse_rest.json`**: the weak pulse at rest: `phi0`, `phi0_P`, `amp`, `kappa`, `w`, `L`, `T_lo`, `t_w` ($3\kappa\rho_0 w^2/c$), `t`, `x`, and `T_exact`, `T_rsla`, `T_arsa` per output.
- **`expanding_sphere.csv`**: box means at 13 radii: `phi0`, `method`, `r_out_R0`, `t_s`, `Tg_mean_K`, `Tg_min_K`, `Tg_max_K`, `Tr_mean_K`, `phi_min`, `phi_max`, `T_exact_K` (eq:sphere_adiabat).
- **`mf_shocks.csv`**: one row per test (1 subcritical; 3 supercritical at $\phi_0 = 0.03$; 2 supercritical at the equal-accuracy $\phi_0$): `phi0`, `h_max`, `tau_c_min`, `tau_c_max`, `L1_rsla`, `L1_arsa`, `L1_arsa_fluxform`, `cells_phi1`. `_profiles.csv`: $T$, $T_{\rm rad}$ and `phi_arsa` at $t_{\rm end}$. `_match.csv`: the equal-accuracy $\phi_0$ and `speedup_substeps`.
- **`marshak_vaytet_scores.csv`**: `phi0`, `L1_rsla`, `L1_arsa`, and the RSLA of equal accuracy (`phi0_rsla_matched`, `speedup`). `_profiles.csv`: $T$, $T_{{\rm rad},g}$, $\phi_g$ at $\phi_0 = 0.01$. `_groups.csv`: `kappa`, `tau_c` per group. `_budget.csv`: share [%] of each group in the flux and the heating. `_attribution.csv`: $\hat c$ reduced in some groups only.
- **`rsla_accuracy_systems.csv`**: per system of Table 1: `tau`, `tau_v_c` ($\mathcal{D}$), `h`, `slowdown` (eq:slowdown), `region`.
- **`theta_scan.csv`**: per case, `nu_pi`, `active`, `max_g_below`, `theta_star`, `nu_formula` (eq:numax) and the parameters. `theta_curves.csv`: $\nu_{\rm crit}(\theta)$ for Fig. A.10.
- **`media_parameters.md`**: derivation and sources of the star-formation, disk and stellar rows of Table 1 (input).

## The tests

All runs except the linear waves are 3D: the pulse, shocks and Marshak wave as slabs of $N\times2\times2$ cubic cells, periodic in $y$ and $z$; the sphere as a $32^3$ box. The RSLA and the ARSA run at $\nu_0 = 0.9/\sqrt{3}$ (Honghu.jl $C_r = 3\nu_0 = 0.9\sqrt{3}$), the ARSA with $\beta = 0.9$; the $\hat c = c$ references at $\nu = 0.4$ ($C_r = 1.2$). The 1D linear waves use the $P_1$ closure at $\nu_0 = 0.9\sqrt{3}$.

- **Advected pulse** (Sec. 5.2, Fig. 5): a Gaussian pulse, $w = 24$ cm, in pressure balance, $\rho = 1.2$ g cm$^{-3}$, $\kappa = 500$ cm$^2$ g$^{-1}$, $512$ cells over 1024 cm, moving at $3\times10^7$ cm/s for 256 widths; $T_{\rm lo}$ set by bisection so that $\phi_0/\mathcal{P} = 3$ and 300 at the peak. Reference: the exact solution, i.e. the enthalpy diffusion of the pulse at rest in column mass, shifted by $vt$. Position: the centroid of $T - T_{\rm lo}$.
- **Weak pulse at rest** (Sec. 5.2, Fig. 5 bottom): amplitude 0.1, $\phi_0/\mathcal{P} = 3$ at the peak, $\kappa = 30$ cm$^2$ g$^{-1}$, 384 cells over 768 cm, to $1000\,t_w = 2.1$ ms; exact reference as above, without the shift.
- **Expanding sphere** (Sec. 5.3, Fig. 6): the interior of the sphere of Roth & Kasen (2015), with absorbing opacity and $v_{\rm out} = 10$ km/s, in homologous flow; exact reference (the adiabat of gas plus radiation). Caveat: in a uniform flow the test checks the heat capacity of the reduced radiation and the work on it, not transport across gradients.
- **Radiative shocks** (Sec. 5.4, Figs. 7, 8): the shocks of Melon Fuksman et al. (2021) on 64 cells; the supercritical one with a constant mass opacity ten times theirs.
- **Marshak wave** (Sec. 5.5, Fig. 9): the multigroup wave of Vaytet et al. (2011), four groups, 64 cells. Equal accuracy: a bisection in $\phi_0$ for the RSLA with the error of the ARSA; the speedup is the ratio of sub-steps.
- **Fig. 1**: the systems of Table 1 in the plane of $\mathcal{D}$ and $\mathcal{P} = 3\hat{h}/8$ ($\gamma = 5/3$) at $\phi_0 = 10^{-2}$; regions from eq:rsla_valid with a factor of 10. Caveat: positions inherit the uncertainties of Table 1 (a factor of 2–3).

## Key numbers and their sources

| Claim in the paper | Source |
|---|---|
| Pulse (Sec. 5.2): RSLA lags 13 and 0.28 widths after 256 widths traveled; ARSA within $0.04\,w$ | `pulse_drift.txt` (lines `256 w traveled`) |
| Weak pulse at rest: $L_1$ RSLA $15\%$, ARSA $0.064\%$; RSLA peak $19\%$ too high; $t = 8.9\,t_{\rm diff}$; $\tau_{\rm c} = 65$–$72$ | `pulse_rest.txt` (line `t =  1000 t_w`, and its header lines) |
| Sphere (Sec. 5.3): RSLA $+270\%$ and $+6.2\%$, ARSA $0.35\%$ and $0.067\%$; $\phi_* = 2.6\times10^{-4}$, $\tau_{\rm thresh} = 23$ and $3.8$; $\tau_{\rm c}$, $\mathcal{D} = 63$ | `expanding_sphere.txt` |
| Shocks (Sec. 5.4): $L_1$ $0.78\%$/$0.62\%$ (subcritical), $14\%$/$3.8\%$ (supercritical); 53 of 64 cells can reach $\phi = 1$ | `mf_shocks.csv` (`L1_rsla`, `L1_arsa`); `mf_shocks.txt` (lines `isotropic`) |
| Shocks: equal accuracy at $\phi_0 = 0.011$, speedup 2.6 | `mf_shocks_match.csv` |
| Sec. 4: the flux form is worse ($52\%$); $\phi$ jumps by up to 65 times between cells and stays stable | `mf_shocks.csv` (`L1_arsa_fluxform`); `mf_shocks_profiles.csv` (`phi_arsa`), printed by `checks/check_phi_jumps.py` |
| Marshak (Sec. 5.5): RSLA $29\%$, $13\%$, $5.3\%$; ARSA $25\%$, $9.5\%$, $3.3\%$; speedups 1.4, 1.5, 1.6 | `marshak_vaytet_scores.csv` |
| Marshak: group 4 carries $96\%$ of the flux; $\hat c$ reduced in group 4 only $9.3\%$, in groups 1–3 only $5.9\%$ | `marshak_vaytet_budget.csv`, `marshak_vaytet_attribution.csv` |
| Table 1 and Fig. 1: $\hat{h}$, $\mathcal{D}$, $\phi/\phi_0$ of each system | `table_media.txt`, `rsla_accuracy.txt` |
| Sec. 2.3: RSLA phase speed 1.13 at high $\tau$ | `waves_rsla.txt` |
| Sec. 3.1.7: the stable limit at most $0.6\%$ below the bound (moving background, up to $6\times10^{-3}$ at $\mathcal{M} = 3$); matter transport shifts it by up to $2.4\%$ | `checks/check_mach_worst.py`; `checks/check_mach.py` (`hydro=True`) |
| App. A.2: no gap in 2D and 3D; in 1D, flux at $75^\circ$ or more, 5 of 1375 cases up to $90\%$, 2 of 2700 unstable at the limiter step | `m1/absorbing_scan.json` (rules `gap` and `actual`) |
| App. A.2: thin anisotropic M1 cells stable in 1–3D up to $\lvert\vec F\rvert/cE = 0.999$; Jacobian speed up to 1.7 times $T_{dd}^{1/2}$ | `m1_stability.txt` |
| Table A.2 | `theta_scan.csv` |
| Sec. 3.1.4: dropping $K_F$ shifts the condition by $0.40$ at $\nu\tau_F = 1.8\times10^4$ | `m1/dropped_momentum_counterexample.json` |
| Sec. 3.1.9: wave-speed correction unstable in 110 of 936 cases | `checks/check_wavespeed.py` |
| Sec. 3.3: multigroup checkerboard stable in 4000 random cases | `checks/check_multigroup.py` |
| App. C: no step condition suffices for $f > 2\sqrt{3}/5$ in thick cells | `m1/finite_flux_obstruction.json`; Lean: `lean/RSLA/M1Full.lean`, `M1UnequalFull.lean` |

## Checks

The scripts in `checks/` check statements in the text; each says in its header what it checks, and only `checks/m1/` writes files (`data/m1/*.json`, by the step `m1_checks`). Run them with `python3 checks/<script>` or `julia --project=. -t 4 checks/<script>`. They cover the closed-form bound against the full von Neumann analysis (`check_formula.py`, `check_mach*.py`, `check_md.py`, `check_checkerboard_*.py`, `sym_energy.py`), the alternative schemes and variants (`check_schemes.py`, `check_multigroup.py`, `check_two_temp.py`, `check_wavespeed.py`), the drift of Sec. 2.2 (`check_entropy_drift.py`), checkerboard runs in Honghu.jl (`checkerboard.jl`, `cb_*.jl`), and the M1 analysis of App. C by B. Wibking (`checks/m1/`).

## Lean proofs

`lean/` holds B. Wibking's Lean 4 formalization of Propositions 1, A.1 and B.1, their lemmas, and the full-M1 theorems of App. C; `lean/README.md` says what is proved and under which hypotheses. With [elan](https://github.com/leanprover/elan) installed, `lean/scripts/check.sh` builds it (Lean 4.19.0, about 4 minutes) and audits every theorem for extra axioms. It is not part of `./run_all.sh`.

## Citation and license

If you use these scripts or data, please cite the paper and Honghu.jl; `CITATION.cff` has both. The code is under the MIT license (`LICENSE`).
