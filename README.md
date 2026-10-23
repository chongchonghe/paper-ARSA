# Scripts and data for "An Adaptive Reduced Speed of Light Approximation"

This repository makes every figure, table and quoted number of C.-C. He, M. R. Krumholz and B. D. Wibking, "An Adaptive Reduced Speed of Light Approximation" (submitted to the Journal of Computational Physics). The simulations use [Honghu.jl](https://github.com/chongchonghe/Honghu.jl) v0.1.0 ([doi:10.5281/zenodo.23569082](https://doi.org/10.5281/zenodo.23569082)), pinned in `Manifest.toml`.

## Quick start

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # Julia 1.12; installs Honghu.jl v0.1.0 and its dependencies
python3 -m pip install -r requirements.txt            # Python 3.13
./run_all.sh figures    # the figures and Table 1, from the data in this repository: about 20 s
./run_all.sh            # everything, the simulations first: about 25 min with 14 threads
```

`PY=<python> ./run_all.sh` picks another Python. The figures set matplotlib's `text.usetex`, so they need a LaTeX installation.

The outputs go to fixed places:

- `figures/`: every figure of the paper, as PDF, with PNG previews of most of them.
- `tables/tab_media.tex`: Table 1.
- `data/` and the `fig_*/` and `floor_energy/` folders: the simulation data, the printouts the paper quotes numbers from (`metadata.txt` and `data/*.txt`), and the figures made from each folder.
- `logs/`: the printout of each run of `./run_all.sh` (not in git).

Each printout starts with a header: the date, the commit of this repository (`-dirty`: uncommitted changes to the script), and the Honghu.jl version and tree hash. The data in git were made with the commit that their header names.

## Which script makes what

| Paper | Output | Script | Input |
|---|---|---|---|
| Fig. 1 | `figures/rsla_accuracy.pdf` | `fig_rsla_accuracy/plot.py` | `fig_rsla_accuracy/*.csv` (`run.py`; see its README) |
| Fig. 2 | `figures/waves_RSLA.pdf` | `fig_waves_rsla.py` | none (linear theory) |
| Fig. 3 | `figures/checkerboard.pdf` | `fig_tests.py checkerboard` | `data/checkerboard.csv` (Honghu.jl `examples/irsla_paper.jl`, section `checkerboard`) |
| Fig. 4 | `figures/waves_ARSA.pdf` | `fig_waves.py` | none (linear theory) |
| Fig. 5 | `figures/linear_waves_code.pdf` | `fig_tests.py linear_waves` | `data/linear_waves_rhd.csv` (Honghu.jl `examples/linear_waves_rhd.jl`) |
| Fig. 6 | `figures/diffusion_pulse.pdf` | `fig_diffusion_pulse/plot.py` | `fig_diffusion_pulse/profiles.csv` (`run.jl`; see its README) |
| Figs. 7, 8 | `figures/shocks_Tgas.pdf`, `figures/shocks_Trad.pdf` | `fig_tests.py mf_shocks` | `data/mf_shocks.csv`, `data/mf_shocks_profiles.csv` (`mf_shocks.jl`) |
| Fig. 9 | `figures/marshak_vaytet.pdf` | `fig_marshak_vaytet/plot.py` | `fig_marshak_vaytet/*.csv` (`run.jl`; see its README) |
| Table 1 | `tables/tab_media.tex` | `table_media.py` | the estimates typed into the script |
| Table 2 | typed into the paper | none | see [Table 2 and numbers in the text](#table-2-and-numbers-in-the-text) |
| Table B.3 | typed into the paper | none | `data/theta_scan.csv` (`check_theta.py`) |
| Fig. B.10 | `figures/theta_scan.pdf` | `fig_theta.py` | `data/theta_curves.csv` (`check_theta.py`) |

`fig_waves_rsla.py` also writes `figures/waves_RSLA_P0.1.pdf`, Fig. 2 at the parameters of an earlier version of the figure ($\mathcal{P} = 0.1$, $\phi = 3\times10^{-3}$), which the paper does not show.

## The simulations

The shock, pulse and Marshak runs are 3D slabs ($N\times2\times2$ cubic cells, periodic in $y$ and $z$). `./run_all.sh` starts these runs in parallel, each in its own process, and runs `floor_energy/run.sh` after them, because it uses their reference solutions. To run one alone, from the top folder:

| Run | Writes | Time |
|---|---|---|
| `julia --project=. -t 4 $EXAMPLES/irsla_paper.jl --no-plot --outdir=data checkerboard` | `data/checkerboard.csv` | 10 s |
| `julia --project=. -t 4 $EXAMPLES/linear_waves_rhd.jl --no-plot --outdir=data` | `data/linear_waves_rhd.csv` | 1 min |
| `./run_mf_shocks.sh` | `data/mf_shocks*.csv`, `data/mf_shocks.txt` | 10 min |
| `julia --project=. -t 1 arsa_phi0.jl` | `data/arsa_phi0.txt` | 2 min |
| `fig_diffusion_pulse/run.sh` | `fig_diffusion_pulse/profiles.csv`, `metadata.txt` | 25 min |
| `fig_marshak_vaytet/run.sh` | `fig_marshak_vaytet/*.csv`, `metadata.txt` | 3 min |
| `fig_rsla_accuracy/run.sh` | `fig_rsla_accuracy/*.csv`, `metadata.txt` | 1 s |
| `floor_energy/run.sh` (after the shocks, the pulse and the Marshak wave) | `floor_energy/results.csv`, `metadata.txt` | 1 min |
| `python3 check_theta.py` | `data/theta_scan.csv`, `data/theta_curves.csv` | 1 min |
| `python3 check_m1.py` | `data/check_m1.txt` | 2 min |

Here `EXAMPLES=$(julia --project=. -e 'using Honghu; print(pkgdir(Honghu, "examples"))')`, the examples folder of the installed Honghu.jl. The times are from one `./run_all.sh` on an Apple M-series laptop with 18 cores, all runs at once; each run uses one thread, the first two four. `./run_all.sh` saves the printouts of the two Honghu.jl examples in `data/checkerboard.txt` and `data/linear_waves_rhd.txt`.

## Table 2 and numbers in the text

These values are typed into the paper by hand. Errors in the CSV files are fractions; the paper gives them in percent.

| Where in the paper | Source |
|---|---|
| Table 2, subcritical shock ($\phi_0 = 10^{-3}$, CSV test 1) and supercritical shock ($\phi_0 = 0.03$ and the equal-accuracy $0.011$, CSV tests 3 and 2) | `data/mf_shocks.csv`, columns `L1_rsla`, `L1_arsa`, `L1_arsa_fluxform` |
| Sec. 5.4, errors of the supercritical shock at $\phi_0 = 0.1$ to $10^{-3}$ (the sweep) | `data/arsa_phi0.txt`, columns `rsla_T%`, `arsa_T%` |
| Sec. 5.4, dimensionless numbers ($\hat{h}$, $\mathcal{P}$, $\phi_*$, $\tau_{\rm c}$, the threshold of eq:tau_one, cells at $\phi = 1$) | `data/mf_shocks.txt`, the printout of `mf_shocks.jl` |
| Sec. 4, flux form ($52\%$, supercritical shock at $\phi_0 = 0.03$) | `data/mf_shocks.csv`, column `L1_arsa_fluxform`; the ARSA with $E$ floored before the exchange update, not quoted in the paper, is at the end of `data/mf_shocks.txt` |
| Table 2, pulse at rest and advected; Sec. 5.3 and the realizability paragraph of Sec. 4 | `fig_diffusion_pulse/metadata.txt` |
| Table 2, Marshak wave; Sec. 5.5, speedups at equal accuracy (1.4, 1.5, 1.6) | `fig_marshak_vaytet/scores.csv` |
| Sec. 5.4, equal accuracy on the supercritical shock: the ARSA at $\phi_0 = 0.011$ matches the RSLA at $0.03$ (speedup 2.6) | `data/mf_shocks_match.csv` (`mf_shocks.jl`) |
| Sec. 4, energy a floor before the exchange update would add ($4\%$ and $15\%$ on the shock, about $10^{-4}$ or less on the pulse, none on the Marshak wave) | `floor_energy/results.csv`, column `floor_energy_fraction` |
| Sec. 2.4, numbers read from Fig. 2 at $\tau = 10^4$ | `data/fig_waves_rsla.txt` |
| Sec. 2.4, slow-mode rates at $\tau = 10^4$ ($2.66\times10^{-4}$ exact, $1.48\times10^{-4}$ RSLA, factor 1.8) and $\hat{h}/\gamma = 8.0\times10^{-4}$ | `data/fig_waves_rsla.txt` (`exact entropy decay`, second `RSLA radiation decays`); eq:slow_rate with $\hat{h}/\gamma$ from the parameters of Fig. 2 ($\mathcal{P} = 5\times10^{-4}$, $\gamma = 5/3$, $\phi = 10^{-3}$, $\mathcal{C} = 10^4$) gives the same values |
| Sec. 2.4, $\hat{h} = 1.3\times10^{-3}$, $\hat{h}/\phi = 1.3$, $\sqrt{\gamma_\phi} = 1.13$ (eq:gamma_phi), the damping factors 2.3 and 1.8 and their product 4, and slow-mode rate $= 1/\mathrm{Pe}_\phi$ at $\mathcal{M} = 1$ | `data/check_waves_heat_capacity.txt` (closed forms against the dispersion relation at $\tau = 10^4$) |
| Appendix B.2, last paragraph: the M1 Jacobian's largest characteristic speed up to 1.7 times $T_{dd}^{1/2}$; every mode of a thin cell stable at $\nu_0 = 1/\sqrt{D}$ and at $\nu\sum_d\sigma_d = 1$ in 1--3D for $\lvert\vec{F}\rvert/(cE) \leq 0.999$ | `data/check_m1.txt` |
| Sec. 3.4, last paragraph, and Sec. 6.2 ($\tau_{\rm c} = 0.15$, massive protostellar envelope): $\tau_{\rm c} = 5$ and 30 and $\tau_{\rm thresh} \approx 3.5$ (PPD midplanes), $\tau_{\rm c} = 3.75$, $\tau_{\rm thresh} = 3.9$ and $\phi = 32\,\phi_0$ (first hydrostatic core), $\phi_* = 4\times10^{-3}$ and $1.3\,\phi_0$ (solar core) | `data/table_media.txt`, columns `tau_c`, `phi_*(beta)`, `tau_thresh`, `phi/phi0` |
| Table B.3 | `data/theta_scan.csv` |
| Table 2, columns $\tau_{\rm c}$, $\hat{h}$, $\phi_*$ (largest $\phi_*$, isotropic radiation in 3D): shocks | `data/mf_shocks.txt`: `tau_c in [...]`, `h max`, `isotropic largest phi_*` |
| Table 2, column $\tau_{\rm thresh}$ (eq:tau_one, isotropic radiation), shocks | supercritical at $\phi_0 = 0.03$: `data/arsa_phi0.txt`, columns `thr_min_a`, `thr_max_a` (cells with an infinite threshold: 64 minus `nstar_a`); subcritical and the matched $\phi_0$: `data/mf_shocks.txt`, lines `isotropic` |
| same, pulse (at rest and advected) | `fig_diffusion_pulse/metadata.txt`, block `T_lo = 1e+06 K`: `h in [...]`, `tau_c in [...]`, `largest phi_*`, and for $\tau_{\rm thresh}$ the line `phi0 = 0.01`: `finite thresholds in [...]` |
| same, Marshak wave | `fig_marshak_vaytet/metadata.txt`: line `reference at t_end` ($\hat{h}$, $\phi_*$) and the lines `eq:tau_one threshold` below it ($\tau_{\rm thresh}$); `groups.csv`, column `tau_c` |

## Other files

- `pdars.py` (the PD-ARS amplification matrix), `rhd_waves.py` (the linear waves of the radiation-hydrodynamics equations) and `waves.py` are modules that the figure scripts import; `paper_style.py` is the matplotlib style of the figures; `common.sh` holds the shell functions of the run scripts.
- `arsa_phi0.jl` writes `data/arsa_phi0.txt`, a table of the derived variables of the ARSA ($\phi_*$, the threshold of eq:tau_one, the cells predicted and found at $\phi = 1$, the range of $\phi$, the sub-steps) and the errors of the RSLA and the ARSA as a function of $\phi_0$, for the two shocks. It uses the set-up of `mf_shocks.jl`.
- `check_*.py`, `sym_energy.py`, `two_temp_phistar.py`, `checkerboard.jl` and `cb_*.jl` check statements in the text. Except for `check_theta.py`, `check_m1.py` and `check_waves_heat_capacity.py`, which `./run_all.sh` runs, they make no figure, table or file; run them with `python3 <script>` or `julia --project=. -t 4 <script>`. Each says in its header what it checks.
- `copy_to_paper.sh` copies `figures/*.pdf` and `tables/tab_media.tex` into the paper's LaTeX source (`../writing-latex-v2` by default), for the authors.

## Citation and license

If you use these scripts or data, please cite the paper and Honghu.jl; `CITATION.cff` has both. The code is under the MIT license (`LICENSE`).
