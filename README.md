# Scripts and data for "An Adaptive Reduced Speed of Light Approximation"

This repository makes every figure, table and quoted number of C.-C. He, M. R. Krumholz and B. D. Wibking, "An Adaptive Reduced Speed of Light Approximation" (submitted to the Journal of Computational Physics). The simulations use [Honghu.jl](https://github.com/chongchonghe/Honghu.jl) v0.1.0 ([doi:10.5281/zenodo.23569082](https://doi.org/10.5281/zenodo.23569082)), pinned in `Manifest.toml`.

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
README.md  LICENSE  CITATION.cff  RELEASING.md
requirements.txt                  the Python packages, with versions
Project.toml  Manifest.toml       the Julia environment (Honghu.jl v0.1.0)
run_all.sh                        the single entry point
src/
  sim_*.jl                        simulations with Honghu.jl; write data/
  calc_*.py                       calculations (linear theory, scans); write data/
  make_*.py                       figures and Table 1, from data/ (or linear theory); write figures/, tables/
  pdars.py  rhd_waves.py  waves.py  paper_style.py     modules the scripts import
checks/                           scripts that check statements in the text; no output files, except checks/m1/
  m1/                            checks of the M1 analysis of App. C (B. Wibking); write data/m1/
lean/                             Lean 4 proofs of the stability results (B. Wibking)
data/                             all data and printouts (generated, in git)
figures/                          every figure of the paper, as PDF (generated, in git)
tables/                           Table 1 (generated, in git)
copy_to_paper.sh                  for the authors: copies figures/ and tables/ into the paper's source
```

## Which script makes what

| Paper | Output | Script | Data |
|---|---|---|---|
| Fig. 1 | `figures/rsla_accuracy.pdf` | `make_rsla_accuracy.py` | `rsla_accuracy_systems.csv` (`calc_rsla_accuracy.py`) |
| Fig. 2 | `figures/waves_RSLA.pdf` | `make_waves_rsla.py` | none (linear theory) |
| Fig. 3 | `figures/checkerboard.pdf` | `make_checkerboard.py` | `checkerboard.csv` (Honghu.jl `examples/irsla_paper.jl`, section `checkerboard`) |
| Fig. 4 | `figures/waves_ARSA.pdf` | `make_waves_arsa.py` | none (linear theory) |
| Fig. 5 | `figures/linear_waves_code.pdf` | `make_linear_waves.py` | `linear_waves_rhd.csv` (Honghu.jl `examples/linear_waves_rhd.jl`) |
| Fig. 6 | `figures/pulse_drift.pdf` | `make_pulse_drift.py` | `pulse_drift.json` (`sim_pulse_drift.jl`), `pulse_rest.json` (`sim_pulse_rest.jl`, bottom row) |
| Sec. 5.3 figure (fig:expanding_sphere) | `figures/expanding_sphere.pdf` | `make_expanding_sphere.py` | `expanding_sphere.csv` (`sim_expanding_sphere.jl`) |
| Figs. 7, 8 | `figures/shocks_Tgas.pdf`, `figures/shocks_Trad.pdf` | `make_mf_shocks.py` | `mf_shocks.csv`, `mf_shocks_profiles.csv` (`sim_mf_shocks.jl`) |
| Fig. 9 | `figures/marshak_vaytet.pdf` | `make_marshak_vaytet.py` | `marshak_vaytet_*.csv` (`sim_marshak_vaytet.jl`) |
| Fig. B.10 | `figures/theta_scan.pdf` | `make_theta_scan.py` | `theta_curves.csv` (`calc_theta_scan.py`) |
| Table 1 | `tables/tab_media.tex` | `make_table_media.py` | the estimates typed into the script; the star-formation and disk rows and the stellar rows from `data/media_parameters.md` (the middle layers from `convective_layers.txt`) |
| Table 2 | typed into the paper | none | see [Table 2 and numbers in the text](#table-2-and-numbers-in-the-text) |
| Table B.3 | typed into the paper | none | `theta_scan.csv` (`calc_theta_scan.py`) |

Scripts are in `src/`, data in `data/`. `make_waves_rsla.py` also writes `figures/waves_RSLA_P0.1.pdf`, Fig. 2 at the parameters of an earlier version of the figure ($\mathcal{P} = 0.1$, $\phi = 3\times10^{-3}$), which the paper does not show.

## The steps

`./run_all.sh` runs the steps below in parallel, each in its own process, then `floor_energy`, which uses the reference solutions of three of them, then the figures. `./run_all.sh <step> ...` runs only the steps named. Each step writes its printout, after a header, to `data/<step>.txt`. The header gives the date, the commit of this repository (`-dirty`: uncommitted changes in `src/` or `Manifest.toml`), and the Honghu.jl version and tree hash. The data in git were made with the commit that their header names.

| Step | Script | Writes (in `data/`) | Time |
|---|---|---|---|
| `checkerboard` | Honghu.jl `examples/irsla_paper.jl`, section `checkerboard` | `checkerboard.csv` | 10 s |
| `linear_waves_rhd` | Honghu.jl `examples/linear_waves_rhd.jl` | `linear_waves_rhd.csv` | 1 min |
| `mf_shocks` | `sim_mf_shocks.jl` | `mf_shocks*.csv` | 10 min |
| `arsa_phi0` | `sim_arsa_phi0.jl` | `arsa_phi0.txt` (the table is the printout) | 2 min |
| `diffusion_pulse` | `sim_diffusion_pulse.jl` | `diffusion_pulse_profiles.csv` | 25 min |
| `marshak_vaytet` | `sim_marshak_vaytet.jl` | `marshak_vaytet_*.csv` | 3 min |
| `floor_energy` | `sim_floor_energy.jl` | `floor_energy.csv` | 1 min |
| `rsla_accuracy` | `calc_rsla_accuracy.py` | `rsla_accuracy_systems.csv` | 1 s |
| `convective_layers` | `calc_convective_layers.py` | `convective_layers.txt` (the printout: the middle layers of the stellar rows of Table 1) | 1 s |
| `theta_scan` | `calc_theta_scan.py` | `theta_scan.csv`, `theta_curves.csv` | 1 min |
| `m1_stability` | `calc_m1_stability.py` | `m1_stability.txt` (the printout) | 2 min |
| `waves_heat_capacity` | `calc_waves_heat_capacity.py` | `waves_heat_capacity.txt` (the printout) | 1 s |
| `m1_checks` | `checks/m1/*.py` | `m1/*.json` | 5 min |
| `pulse_rest` | `src/sim_pulse_rest.jl` (two single-threaded worker processes) | `pulse_rest.json` | 12 min |
| `pulse_drift_substep` | `checks/check_pulse_drift_substep.jl` (after `pulse_drift`, whose exact solution it reads) | `pulse_drift_substep.txt` (the printout) | 3 min |
| `pulse_drift` | `sim_pulse_drift.jl` (eight runs on eight single-threaded worker processes) | `pulse_drift.json` | 3 min |
| `expanding_sphere` | `sim_expanding_sphere.jl` (four worker processes; `THREADS=<n>` threads each, default 1) | `expanding_sphere.csv` | 5 min |
| `table_media`, `waves_rsla` | `make_table_media.py`, `make_waves_rsla.py` | `table_media.txt`, `waves_rsla.txt` (printouts quoted in the text) | 5 s |

The times are from one `./run_all.sh` on an Apple M-series laptop with 18 cores, all steps at once. Each step uses one thread, except the two Honghu.jl examples, which use four, and `pulse_drift`, which uses eight worker processes (Honghu.jl spreads its grid loops over the threads of a process, so several runs in one multi-threaded process slow each other down).

## Data files

All files are in `data/`. Units are cgs, temperatures are in K, and errors are fractions (the paper gives them in percent). Each `.txt` file is the printout of its step, with the derived parameters and the errors of every run.

- **`checkerboard.csv`**: growth factor per sub-step of a checkerboard. `D` (dimensions), `tau` (cell optical depth), `Pcal` ($\mathcal{P}$), `limiter` (`fixed` $\phi$, or the limiter of earlier drafts), `fac` ($\phi/\phi_{\max}$), `phi`, `g` (growth factor).
- **`linear_waves_rhd.csv`**: complex frequency of linear waves measured in Honghu.jl and from the dispersion relation. `tau` (optical depth per radian), `scheme` (`exact`, `rsla`, `irsla` for the ARSA), `mode` (`acoustic`, `entropy`), `phi`, `Re_meas`, `Im_meas`, `Re_theory`, `Im_theory`.
- **`mf_shocks.csv`**: one row per shock test. `test`, `u` [cm/s], `t_end` [s], `phi0`, `h_max` ($\hat{h}$), `P_max` ($\mathcal{P}$), `tau_c_min`, `tau_c_max`, `L1_rsla`, `L1_arsa`, `L1_Trad_rsla`, `L1_Trad_arsa`, `L1_arsa_fluxform`, `cells_phi1` (cells at $\phi = 1$). `mf_shocks_profiles.csv`: `test`, `x` [cm], $T$ and $T_{\rm rad}$ of the reference, the RSLA and the ARSA, and `phi_arsa`, at $t_{\rm end}$. `mf_shocks_phi0.csv`: test 3 at other `phi0`. `mf_shocks_match.csv`: the $\phi_0$ at which the ARSA matches the error of the RSLA, with the sub-steps and the speedup.
- **`arsa_phi0.txt`**: for the supercritical shock, one row per $\phi_0$: $\phi_*$, the threshold of eq:tau_one, the cells predicted and found at $\phi = 1$, the range of $\phi$, the sub-steps, and the errors of the RSLA and the ARSA in $T$ and $T_{\rm rad}$ [%]. The header lines explain the columns.
- **`diffusion_pulse_profiles.csv`**: $T$ at the final time of every pulse run. `T_lo` [K], `v0` [cm/s], `t` [s], `phi0`, `x` [cm], `T_initial`, `T_ref` ($\hat c = c$), `T_rsla`, `T_arsa`, `phi_arsa`. Each (`T_lo`, `v0`, `phi0`) is one set of runs; `sim_floor_energy.jl` uses the references `T_ref` at `T_lo = 1e6`.
- **`marshak_vaytet_scores.csv`**: `phi0`, `L1_rsla`, `L1_arsa`, `sub_rsla`, `sub_arsa`, and the RSLA of equal accuracy: `phi0_rsla_matched`, `L1_rsla_matched`, `sub_rsla_matched`, `speedup` (`sub_rsla_matched / sub_arsa`). `_profiles.csv`: `run` (`ref`, `rsla`, `arsa`), `x` [cm], `T_gas`, `T_g1`–`T_g4`, `phi_g1`–`phi_g4`, at $\phi_0 = 0.01$. `_groups.csv`: `group`, `nu_lo`, `nu_hi` [Hz], `kappa` [cm$^2$ g$^{-1}$], `tau_c`. `_budget.csv`: in the reference, the share [%] of each group in the flux at four positions at $t_{\rm end}$, and in the heating of the gas behind and ahead of the front ($T_{\rm gas} = 650$ K). `_attribution.csv`: `reduced_in`, `phi_thin` ($\phi$ of group 4), `phi_thick` ($\phi$ of groups 1–3), `L1`, `substeps`.
- **`floor_energy.csv`**: one row per case. `case`, `L1_conservative` and `L1_floor` (the ARSA with Honghu.jl's conservative floor and with the floor everywhere), `events_below_floor`, `substeps`, `floor_energy_fraction` (the energy the floor adds, as a fraction of a reference energy; see below).
- **`expanding_sphere.csv`**: the box means of the expanding-sphere runs at 13 outer radii. `phi0`, `method` (`rsla`, `arsa`), `r_out_R0` ($r_{\rm out}/r_{\rm out,i}$), `t_s` [s], `Tg_mean_K`, `Tg_min_K`, `Tg_max_K` (gas temperature over the box), `Tr_mean_K`, `phi_min`, `phi_max` (range of $\phi$ in the box over the run so far), `T_exact_K` (the adiabat of gas plus radiation, eq:sphere_adiabat).
- **`rsla_accuracy_systems.csv`**: one row per system of Table 1: `system`, `tau` (optical depth across the system), `tau_v_c` ($\mathcal{D}$), `h`, `slowdown` ($t_\phi/t_1$, eq:slowdown), `region`.
- **`theta_scan.csv`**: for each parameter case, `nu_pi` (the stability boundary of the checkerboard, $\theta = \pi$), `active` (whether it is below `nuc`), `max_g_below` (largest growth factor over $\theta$ below `nu_pi`), `theta_star` (the $\theta/\pi$ of the least stable mode just above it) and `nu_formula` (the closed form eq:numax), with the parameters `M`, `chi`, `gamma`, `nuc`, `P`, `tau`. `theta_curves.csv`: $\nu_{\rm crit}(\theta)$ for the four cases of Fig. B.10.
- **`pulse_rest.json`**: the weak pulse at rest. Top level: `phi0`, `phi0_P`, `amp` (0.1), `kappa`, `w`, `L` [cm], `T_lo` [K], `t_w` [s] ($3\kappa\rho_0 w^2/c$), `t` [s] (the outputs, 0 to $1000\,t_w$), `x` [cm] (the cell centres), and `T_exact`, `T_rsla`, `T_arsa` [K], one profile per output.
- **`pulse_drift.json`**: the advected pulse over 256 pulse widths traveled, with output every 32 widths. Top level: `phi0`, `D` ($\kappa\rho w v/c$), `w`, `L` [cm], `v0` [cm/s], `t_out` [s] (the time between outputs), `x` [cm] (the cell centres). `runs`: one entry per $\phi_0/\mathcal{P}$ at the peak, with `phi0_P`, `T_lo` [K], `P_peak`, `h_peak`, `mach`, `arsa_phi_min` (the smallest $\phi$ of the ARSA), `widths` (widths traveled, 0 to 256), and at every output `x_exact`, `x_rsla`, `x_arsa` (real position of the centroid [cm]), `lag_rsla_w`, `lag_arsa_w` (lag behind the exact solution, in pulse widths), and the profiles on `x`: `T_exact`, `T_rsla`, `T_arsa`, `Tr_rsla`, `Tr_arsa` [K].
- **`media_parameters.md`**: the derivation of the star-formation and disk systems of Table 1, and (Sec. C) of the bases of the convective envelopes of the Sun, an RGB star and an AGB star, with their sources (input, not generated).

## The tests

The shock, pulse and Marshak runs are 3D slabs ($N\times2\times2$ cubic cells, periodic in $y$ and $z$). The RSLA and the ARSA run at $\nu_0 = 0.9/\sqrt{3}$ per direction ($C_r = 0.9\sqrt{3}$), the 3D RSLA step for isotropic radiation, and the $\hat c = c$ references at $\nu_0 = 0.4$ per direction.

- **Pulse** (`sim_diffusion_pulse.jl`, Sec. 5.3: the pulse at rest and the clipping of $\lvert\vec{F}\rvert$; also the set-up of `sim_pulse_drift.jl` and the references of `sim_floor_energy.jl`): the advecting pulse of He et al. (2024a), Honghu.jl's `radhydro_pulse_dyn`, with a background ten times colder, so that the gas holds the heat capacity. A Gaussian pulse $T = T_{\rm lo}[1 + e^{-x^2/2w^2}]$, $w = 24$ cm, in pressure balance, $\rho = 1.2$ g cm$^{-3}$ in the background; a periodic box of 1024 cm with $128\times2\times2$ cells, $\kappa = 500$ cm$^2$ g$^{-1}$ ($\tau_{\rm c} = 2400$–$4800$), $\mu = 2.33$, $\gamma = 5/3$, Levermore closure, $v/c$ terms on. Moving runs ($v_0 = 3\times10^7$ cm/s, to $t = 4.8\times10^{-5}$ s, 1.4 box crossings) at $T_{\rm lo} = 10^5$ and $10^6$ K and $\phi_0 = v_0/c$, $2v_0/c$ and $0.01$; a run at rest ($v_0 = 0$, to $t = 4.8\times10^{-4}$ s) at $T_{\rm lo} = 10^6$ K and $\phi_0 = 0.01$. Error: $\sum|T - T_{\rm ref}| / \sum(T_{\rm ref} - T_{\rm lo})$. The end of `diffusion_pulse.txt` also has the ARSA with $\lvert\vec{F}\rvert$ clipped to $cE$ before the exchange update, which Sec. 4 quotes. Caveat: the centroid lag in the printout is meaningless when the RSLA pulse falls more than half a box behind (the runs at $\phi_0 = v_0/c$ and $2v_0/c$); their L1 errors are valid.
- **Marshak wave** (`sim_marshak_vaytet.jl`, Fig. 9, Sec. 5.5): the non-equilibrium multigroup Marshak wave of Vaytet et al. (2011), Honghu.jl's `marshak_vaytet` with a transverse-cell option (with one transverse cell it reproduces `marshak_vaytet` exactly). $64\times2\times2$ cells over 20 cm, 4 log-spaced groups from $6\times10^{10}$ to $6\times10^{14}$ Hz, $\kappa = 2000\,(\nu/4\times10^{13}\,{\rm Hz})^{-2}$ cm$^2$ g$^{-1}$, static medium, to $t = 1.36\times10^{-7}$ s. Equal accuracy: for each $\phi_0$, a 14-step bisection in $\log\phi_0$ finds the larger $\phi_0$ at which the RSLA has the error of the ARSA; the speedup is the ratio of sub-steps. The heating of the gas by group $g$ is $c\rho\kappa_g(E_g - f_g a_r T^4)$, summed over the steps of the reference. The attribution runs take the step of the reference and fix $\phi$ per group (group 4 is the one with $\rho\kappa\Delta x < 1$). $T_{{\rm rad},g}$ is the brightness temperature of the group, $a_r T^4 f_g(T) = E_g$, with $f_g$ the Planck fraction of the group. Caveat: the regular steps in $T_{{\rm rad},4}$, about 4 cells apart, appear in the $\hat c = c$ reference too; they come from the transport of the free-streaming group, not from the ARSA, and their cause has not been checked.
- **Drift of the advected pulse** (`sim_pulse_drift.jl`, Fig. 6, Sec. 5.3): the moving pulse above on $512\times2\times2$ cells for 256 widths traveled (6 crossings of the box), with output every 32 widths, at $\phi_0 = 0.01$, with $T_{\rm lo}$ chosen by bisection so that $\phi_0/\mathcal{P}$ at the peak is 3 and 300 (until 2026-11-08 also 30 and 3000) ($\rho$, $\kappa$, $w$, $v$ fixed, so $\kappa\rho w v/c = 14.4$; $\hat h < \phi_0$ everywhere, so the ARSA keeps $\phi = 1$). The reference is the exact solution ($\phi = 1$): by Galilean invariance, the pulse at rest, shifted by $vt$. At rest the pulse stays in pressure balance and in the diffusion limit, so in the column mass $m$ ($dm = \rho\,dx$) the enthalpy per mass $h = c_pT + \tfrac43 a_rT^4/\rho$ obeys $\partial_t h = (c/3\kappa)\,\partial_m^2(a_rT^4)$; the script solves this on 8192 cells (converged to $2\times10^{-4}$ in $T$). The position of a pulse is the centroid of $T - T_{\rm lo}$ (a circular mean in the box, then the mean within 8 widths of it, which leaves out the small ripples far from the pulse), unwrapped from four samples per output; that of the exact solution is $vt$.
- **Pulse at rest** (`sim_pulse_rest.jl`, bottom row of the pulse figure, ssec:test_thick): a weak pulse, $T = T_{\rm lo}[1 + 0.1\,e^{-x^2/(2w^2)}]$, with $w = 24$ cm, $\phi_0/\mathcal{P} = 3$ at the peak ($\hat{h} = 0.61$--$0.89\,\phi_0$ over the whole pulse) and $v = 0$ ($\mathcal{D} = 0$), $\kappa = 30$ cm$^2$ g$^{-1}$ ($\tau_{\rm c} = 65$--$72$, above $\tau_{\rm thresh}$), in a box of 768 cm on $384\times2\times2$ cells, to $1000\,t_w = 2.1$ ms, $t_w = 3\kappa\rho_0 w^2/c$. The RSLA slows the diffusion of the whole pulse by $(1 + \hat{h}/\phi_0)/(1 + \hat{h})$, up to 1.9; there is no drift. The reference is the exact enthalpy diffusion of `sim_pulse_drift.jl`, without the shift. (Until 2026-11-08 the full pulse, amplitude 1, in a box of 384 cm to $3000\,t_w$: RSLA $14.5\%$, ARSA $0.07\%$.)
- **Expanding sphere** (`sim_expanding_sphere.jl`, Sec. 5.3): the interior of the opaque expanding sphere of Roth & Kasen (2015, Sec. 6.3), with their $\rho$, $T$, $\mu$, $\gamma$, $\kappa$ and expansion factor 11, but absorbing opacity ($\kappa_P = \kappa_E = \kappa_F$) and $v_{\rm out} = 10$ km/s instead of $10^4$ km/s. A $32^3$ box $[-10^{13}, 10^{13}]^3$ cm in homologous flow $\vec{v} = \vec{x}/(t + t_0)$, ghost cells copied from the nearest interior cell with the homologous velocity at their own position; Levermore closure, $v/c$ terms on; RSLA and ARSA at $\phi_0 = 3\times10^{-4}$ and $3\times10^{-3}$ at the 3D RSLA step. The reference is exact. Caveat: in a uniform flow the test checks the heat capacity of the reduced radiation and the work on it, not transport across gradients.
- **Radiative shocks** (`sim_mf_shocks.jl`, Figs. 7 and 8, Sec. 5.4): the shocks of Melon Fuksman et al. (2021) on 64 cells. Test 1 is the subcritical shock; tests 2 and 3 the supercritical shock with 10 times their mass opacity, at $\phi_0 = 0.03$ (test 3) and at the $\phi_0$ at which the ARSA matches the error of the RSLA at 0.03 (test 2). `sim_arsa_phi0.jl` sweeps $\phi_0$ on the same set-up.
- **Flooring $E$ before the exchange update** (`sim_floor_energy.jl`, Sec. 4): the set-ups of the pulse (advected at $\phi_0 = 0.01$, $v_0/c$ and $2v_0/c$, and at rest at $0.01$, all at $T_{\rm lo} = 10^6$ K), the Marshak wave ($\phi_0 = 0.03$, $0.01$, $0.003$) and the supercritical shock (its two $\phi_0$), included from their scripts. Honghu.jl's `floor_rad_energy!` is replaced by a version that either keeps Honghu.jl's conservative rule (floor $E$ only where $\rho e + \sum_g E_g c/\hat c_g \leq 0$) or raises $E$ to the floor wherever it is at or below it, as an earlier algorithm did. The floor changes the energy the exchange update conserves by $\sum (E_{\rm floor} - E)\,c/\hat c\,\Delta V$; this is given as a fraction of the thermal energy of the pulse, of the heat the Marshak slab gains by $t_{\rm end}$, or of the kinetic energy that flows into the shock over the run, $\tfrac{1}{2}\rho_0 u^3 t\,A$. The errors use the stored $\hat c = c$ references, so no reference is rerun.
- **Fig. 1** (`calc_rsla_accuracy.py`): the systems of Table 1 in the plane of the dynamic-diffusion number $\mathcal{D} = \mathcal{M}\mathcal{L}/\mathcal{C} = \rho\kappa\ell v/c$ and the pressure ratio $\mathcal{P} = 3\hat{h}/8$ ($\gamma = 5/3$), for $\phi = 10^{-2}$, the smallest reduction factor of the RSLA considered. Only optically thick systems are shown. The regions follow from the criterion $\phi \gg \max(\mathcal{P}, \mathcal{D})$ of eq:rsla_valid, with $\hat{h}$ for $\mathcal{P}$ and the same factor 10 for both: the RSLA is accurate for $\hat{h} < 0.1\,\phi$ and $\mathcal{D} < 0.1\,\phi$, where its slowdown $(1 + \hat{h}/\phi)/(1 + \hat{h})$ (eq:slowdown) is below 10% and structures carried by the flow do not drift (Sec. 2.2); the ARSA can help for $0.1\,\phi < \hat{h} < 0.1$ at any $\mathcal{D}$, and for $\hat{h} < 0.1\,\phi$ with $\mathcal{D} > 0.1\,\phi$ (right of the dashed line); neither is accurate above $\hat{h} = 0.1$, at any $\mathcal{D}$. For the star-formation and disk rows, $\hat{h}$ is recomputed with the specific heats of Table 1, since `media_parameters.md` takes $c_v = k/(\mu m_{\rm H})$, too small for H$_2$–He gas. Caveat: the positions inherit the uncertainties of Table 1 (a factor of 2–3; for the stellar rows, an order of magnitude in $\mathcal{D}$).

## Table 2 and numbers in the text

These values are typed into the paper by hand. All files are in `data/`.

| Where in the paper | Source |
|---|---|
| Table 2, subcritical shock ($\phi_0 = 10^{-3}$, CSV test 1) and supercritical shock ($\phi_0 = 0.03$ and the equal-accuracy $0.011$, CSV tests 3 and 2) | `mf_shocks.csv`, columns `L1_rsla`, `L1_arsa`, `L1_arsa_fluxform` |
| Sec. 5.4, errors of the supercritical shock at $\phi_0 = 0.1$ to $10^{-3}$ (the sweep) | `arsa_phi0.txt`, columns `rsla_T%`, `arsa_T%` |
| Sec. 5.4, dimensionless numbers ($\hat{h}$, $\mathcal{P}$, $\phi_*$, $\tau_{\rm c}$, the threshold of eq:tau_one, cells at $\phi = 1$) | `mf_shocks.txt` |
| Sec. 4, flux form ($52\%$, supercritical shock at $\phi_0 = 0.03$) | `mf_shocks.csv`, column `L1_arsa_fluxform`; the ARSA with $E$ floored before the exchange update, not quoted in the paper, is at the end of `mf_shocks.txt` |
| Sec. 5.3, the pulse at rest (errors $1.6\%$ and $0.0059\%$) and with $\lvert\vec{F}\rvert$ clipped before the exchange update ($4.6\%$), also the realizability paragraph of Sec. 4 | `diffusion_pulse.txt` |
| Table 2, Marshak wave; Sec. 5.5, speedups at equal accuracy (1.4, 1.5, 1.6) | `marshak_vaytet_scores.csv` |
| Sec. 5.5, the per-group $\phi$ and equilibrium, the share of group 4 in the flux, the errors with $\hat c$ reduced in some groups only | `marshak_vaytet_profiles.csv`, `marshak_vaytet.txt`, `marshak_vaytet_budget.csv`, `marshak_vaytet_attribution.csv` |
| Sec. 5.4, equal accuracy on the supercritical shock: the ARSA at $\phi_0 = 0.011$ matches the RSLA at $0.03$ (speedup 2.6) | `mf_shocks_match.csv` |
| Sec. 4, energy a floor before the exchange update would add ($4\%$ and $15\%$ on the shock, about $10^{-4}$ or less on the pulse, none on the Marshak wave) | `floor_energy.csv`, column `floor_energy_fraction` |
| Sec. 2.3, numbers read from Fig. 2 at $\tau = 10^4$ | `waves_rsla.txt` |
| Sec. 2.3, slow-mode rates at $\tau = 10^4$ ($2.66\times10^{-4}$ exact, $1.48\times10^{-4}$ RSLA, factor 1.8) and $\hat{h}/\gamma = 8.0\times10^{-4}$ | `waves_rsla.txt` (`exact entropy decay`, second `RSLA radiation decays`); eq:slow_rate with $\hat{h}/\gamma$ from the parameters of Fig. 2 ($\mathcal{P} = 5\times10^{-4}$, $\gamma = 5/3$, $\phi = 10^{-3}$, $\mathcal{C} = 10^4$) gives the same values |
| Sec. 2.3, $\hat{h} = 1.3\times10^{-3}$, $\hat{h}/\phi = 1.3$, $\sqrt{\gamma_\phi} = 1.13$ (eq:gamma_phi), the damping factors 2.3 and 1.8 and their product 4, and slow-mode rate $= 1/\mathrm{Pe}_\phi$ at $\mathcal{M} = 1$ | `waves_heat_capacity.txt` (closed forms against the dispersion relation at $\tau = 10^4$) |
| Appendix B.2, last paragraph: the M1 Jacobian's largest characteristic speed up to 1.7 times $T_{dd}^{1/2}$; every mode of a thin cell stable at $\nu_0 = 1/\sqrt{D}$ and at $\nu\sum_d\sigma_d = 1$ in 1--3D for $\lvert\vec{F}\rvert/(cE) \leq 0.999$ | `m1_stability.txt` |
| Sec. 3.4, last paragraph, and Sec. 6.2 ($\tau_{\rm c} = 0.15$, massive protostellar envelope), at $\phi_0 = 10^{-2}$: $\tau_{\rm c} = 5$ and 30 and $\tau_{\rm thresh} \approx 3.4$ (PPD midplanes), $\tau_{\rm c} = 3.75$ and $\tau_{\rm thresh} = 3.5$ (first hydrostatic core), $\phi_* = 2\times10^{-3}$ and $10^{-3}$ and $\tau_{\rm thresh} \approx 4$ (Sun and RGB star, middle of the convective envelope) | `table_media.txt`, columns `tau_c`, `phi_*(beta)`, `tau_thresh`, `phi/phi0` |
| Sec. 3.2.4, why $K_F$ is kept: at $\nu\tau_F = 1.8\times10^4$ and $K_F = 9\times10^{-5}$ dropping it moves the momentum condition by $0.40$, against a margin of $0.1$ | `m1/dropped_momentum_counterexample.json` (`checks/m1/p1_dropped_momentum_counterexample.py`): `qF_at_phi_one`, `KF`, `approximate_momentum_residual` (0.9) and `exact_momentum_residual` (1.305); the shift is their difference |
| Sec. 3.1.6, after Proposition 1, the example: in an absorbing cell with $\lvert\vec{F}\rvert/(cE) = 0.94$ at $80^\circ$ to the axis of a 1D grid, the mode $\theta = 3\pi/4$ grows by $1.5\%$ per step at 0.999 of the checkerboard bound, and the stable limit lies $0.8\%$ below that bound | `m1/checks.json` (`checks/m1/m1_stability.py`), `absorption_counterexample`: `reduced_flux`, `normal_cosine` (1/6), `amplification_modulus`, `stable_limit_over_checkerboard_bound` |
| App. B.2, absorbing anisotropic cells: every mode stable in all 1197 (2D) and 1365 (3D) cases; 2 of 2700 unstable in 1D, with reduced fluxes above 0.97 at $75^\circ$ and $83^\circ$ to the axis | `m1/absorbing_scan.json` (`checks/m1/absorbing_scan.py`), rule `actual` |
| Sec. 3.1.6 after Proposition 1, and App. B.2: at the bound itself ($\beta = 1$) no gap in 2D and 3D, none in 1D with the flux within $68^\circ$ of the axis; in 1D a gap in 5 of 1375 cases, all at $75^\circ$ or more with $\lvert\vec{F}\rvert/(cE) > 0.9$, up to $90\%$ | `m1/absorbing_scan.json`, `gap`: `cases_with_bound`, `positive` (`n[0]` is the cosine with the axis), `max_gap_within_68deg` |
| App. C, the obstruction for strongly directed radiation: $f > 2\sqrt3/5$; the example $f = 80/91$, $\hat{h} = \phi = 1$, isothermal sound speed $c/61$ | `m1/finite_flux_obstruction.json` (`checks/m1/m1_finite_flux_obstruction.py`, with $\mu_E = 1$, $\mu_M = 10^{-4}$ and $\eta = 1$ in the script): `flux_threshold`, `particular_discriminant` $< 0$; $c/a_{\rm iso} = (3\mu_E/8\mu_M)^{1/2} = 61$ for an ideal gas with $\gamma = 5/3$. The theorem of App. C is proved in `lean/RSLA/M1Full.lean` and `M1UnequalFull.lean` |
| Sec. 5.3 (expanding sphere, ssec:test_sphere): $\hat{h} = 2.6\times10^{-4}$ at the start; $\tau_{\rm c} = 1.2\times10^5$ and 89; $\mathcal{D} = 63$; $\hat{h}/\phi_0 = 0.85$ and $0.085$; $\phi_* = 2.6\times10^{-4}$ and $\tau_{\rm thresh} = 23$ and $3.8$ (as Table 3: $\beta = 1$, $\nu_0 = 1/\sqrt{3}$); at $r_{\rm out} = 11\,r_{\rm out,i}$ the RSLA $+270\%$ and $+6.2\%$ (306 K and 88 K against 83 K), the ARSA $0.35\%$ and $0.067\%$, $\phi = 1$ everywhere; Table 3 row: the same errors and RSLA/ARSA $= 780$ and $93$; $\hat{c} = 9\,v_{\rm out}$ at $\phi_0 = 3\times10^{-4}$ (from the parameters) | `expanding_sphere.txt` |
| Table 3, Speedup column: supercritical shock $2.6$; Marshak wave $1.6$, $1.5$, $1.4$ (at $\phi_0 = 0.03$, $0.01$, $0.003$); expanding sphere and advecting pulse estimated as $1/\phi_0$ ($3300$, $330$, $100$); subcritical shock $\approx 1$ (RSLA/ARSA $= 1.2$) | `mf_shocks_match.csv` (`speedup_substeps`, with `phi0_arsa` = 0.0114); `marshak_vaytet_scores.csv` (`speedup`); $1/\phi_0$ from the parameters |
| Table 3, advecting-pulse rows ($\phi_0/\mathcal{P} = 3$ and 300): $\hat{h}$ (largest), $\tau_{\rm c} = 600$--$1200$, $\tau_{\rm thresh}$ ($\beta = 1$), and the error, the lag after 256 widths traveled in pulse widths, of the RSLA and the ARSA and their ratio | `pulse_drift.txt` (the lines `Table 3:`) |
| ssec:test_thick, the weak pulse at rest: $\hat{h} = 0.61$--$0.89\,\phi_0$, $\tau_{\rm c} = 65$--$72$, $\tau_{\rm thresh} \leq 28.5$ ($\beta = 1$), $\phi = 1$ in the ARSA; at $1000\,t_w = 2.1$ ms the $L_1$ error of $T - T_{\rm lo}$ RSLA $15\%$, ARSA $0.064\%$; peaks $0.375$ (RSLA) and $0.316$ (exact) times the amplitude; $t = 8.9\,t_{\rm diff}$, $t_{\rm diff} = 3\kappa\rho_0 w^2/(\hat{h}c) = t_w/\hat{h}$ at the peak ($1000\,t_w$ times $\hat{h} = 8.9\times10^{-3}$) | `pulse_rest.txt` (the lines `t = 1000 t_w` and `Table 3:`) |
| Sec. 5.3, drift of the advected pulse: the centroid lags of the RSLA and the ARSA after 256 widths traveled; $\mathcal{D} = 14$ | `pulse_drift.txt` (the lines `256 w traveled`); `pulse_drift.json` (`lag_rsla_w`, `lag_arsa_w`); $\mathcal{D} = \kappa\rho w v/c = 14.4$ on the first line of `pulse_drift.txt` ($\mathcal{D}/\phi_0 = 1.4\times10^3$ in the figure); the Mach numbers (3.8 to 12) on its `==` lines |
| Fig. 6 caption, the top-row ARSA deviation comes from the radiation substeps: $c\,\Delta t/\Delta x \approx 50$ at $\phi = 1$; $L_1$ after 32 widths traveled $4.4\%$, $0.45\%$ with substeps ten times shorter, $4.8\%$ with a seven times shorter hydrodynamic step | `pulse_drift_substep.txt` |
| Sec. 6.2, PPD midplanes and first hydrostatic core: $\mathcal{D}$ ($0.13\,\phi_0$ for the core, with $\phi_0 = 10^{-2}$), $t_\phi/t_1$ (1.01 for the core); stellar convective envelopes: $\mathcal{D} = 8\times10^3$ (Sun) and $4\times10^4$ (RGB), $\hat{h} = 2.5\times10^{-3}$ and $1.1\times10^{-3}$, $\hat{h} \approx 17$ (AGB, base) | `rsla_accuracy.txt`, columns `D`, `h`, `t_phi/t_1` |
| Sec. 6.2, $\hat{h}$ of the AGB envelope stays above unity out to $T = 10^5$ K | `convective_layers.txt`, last line |
| Table B.3 | `theta_scan.csv` |
| Table 2, columns $\tau_{\rm c}$, $\hat{h}$, $\phi_*$ (largest $\phi_*$, isotropic radiation in 3D): shocks | `mf_shocks.txt`: `tau_c in [...]`, `h max`, `isotropic largest phi_*` |
| Table 2, column $\tau_{\rm thresh}$ (eq:tau_one, isotropic radiation), shocks | supercritical at $\phi_0 = 0.03$: `arsa_phi0.txt`, columns `thr_min_a`, `thr_max_a` (cells with an infinite threshold: 64 minus `nstar_a`); subcritical and the matched $\phi_0$: `mf_shocks.txt`, lines `isotropic` |
| same, Marshak wave | `marshak_vaytet.txt`: line `reference at t_end` ($\hat{h}$, $\phi_*$) and the lines `eq:tau_one threshold` below it ($\tau_{\rm thresh}$); `marshak_vaytet_groups.csv`, column `tau_c` |

## Checks

The scripts in `checks/` check statements in the text; each says in its header what it checks, and none writes a file. Run them with `python3 checks/<script>` or `julia --project=. -t 4 checks/<script>`.

- The stability bound: the closed form eq:numax against the full von Neumann analysis (`check_formula.py`), at nonzero Mach numbers (`check_mach.py`, `check_mach_worst.py`), and in 2D and 3D (`check_md.py`, `check_checkerboard_md.py`); the checkerboard as the least stable mode (`check_checkerboard_norm.py`); the energy block in closed form (`sym_energy.py`).
- Variants: the alternative schemes of Sec. 3.2 (`check_schemes.py`), several groups sharing one gas temperature (`check_multigroup.py`), matter at two temperatures (`check_two_temp.py`, `two_temp_phistar.py`), and a wave-speed correction of the energy dissipation (`check_wavespeed.py`).
- The drift of Sec. 2.2: the entropy mode of the linear waves on a moving background falls behind the flow at $(\hat{h}/\gamma)(1/\phi - 1)$ times its speed, and by $3\mathcal{D}(1 - \phi)/\phi$ radians before it decays, to first order in $\hat{h}/\phi$ (`check_entropy_drift.py`).
- Checkerboard runs in Honghu.jl against the von Neumann prediction (`checkerboard.jl`), with the earlier and the current limiter, the standard RSLA at large $\mathcal{P}$, the flux checkerboard, long runs and other reconstructions (`cb_*.jl`).

- The M1 analysis of App. C, by B. Wibking (`checks/m1/`, run by the step `m1_checks`): the M1 flux Jacobian and its HLL symmetrizer (`m1_stability.py`, which also has the absorbing-cell counterexample), the exchange with M1 (`m1_exchange_bounds.py`, `m1_unequal_opacities.py`), the diffusion-regime theorem with equal and unequal opacity means (`m1_thick_cells.py`, `m1_unequal_diffusion.py`, `m1_diffusion_limits.py`), the obstruction for strongly anisotropic radiation (`m1_finite_flux_obstruction.py`) and the counterexample to dropping $K_F$ from the limiter (`p1_dropped_momentum_counterexample.py`). Each writes its results to `data/m1/<name>.json`; `data/m1/checks.json` is the output of `m1_stability.py`, and `data/m1/lean_validation.json`, written by hand, records the Lean build. `absorbing_scan.py` tests the checkerboard bound with the actual HLL speeds in absorbing cells with anisotropic radiation, in 1D, 2D and 3D (`data/m1/absorbing_scan.json`).

## Lean proofs

`lean/` holds B. Wibking's Lean 4 formalization of Propositions 1, A.1 and B.1, Lemmas A.1–A.5 and B.1–B.3, and the two full-M1 theorems of App. C; `lean/README.md`, `lean/FORMALIZATION.md` and `lean/M1_FORMALIZATION.md` say what is proved and under which hypotheses. With [elan](https://github.com/leanprover/elan) installed, fetch the pinned mathlib cache as `lean/README.md` describes, then `lean/scripts/check.sh` builds the project (Lean 4.19.0; about 4 minutes) and audits every theorem for axioms beyond the standard three. It is not part of `./run_all.sh`.

## Citation and license

If you use these scripts or data, please cite the paper and Honghu.jl; `CITATION.cff` has both. The code is under the MIT license (`LICENSE`).
