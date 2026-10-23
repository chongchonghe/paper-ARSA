# Fig. 9: multigroup Marshak wave

Paper: Fig. 9 (`fig:marshak_vaytet`), §5.5 (`ssec:test_marshak`), and the three Marshak rows of Table 2 (`tab:tests`). §5.5 also quotes the speedups at equal accuracy (`scores.csv`), the per-group $\phi$ and equilibrium (`profiles.csv`, `metadata.txt`), the share of group 4 in the flux (`budget.csv`) and the errors with $\hat c$ reduced in some groups only (`attribution.csv`).

## Reproduce

```bash
$PY plot.py       # figure only, from the data in this folder; about 2 s
./run.sh         # simulations: rewrites the CSV files and metadata.txt; about 6 min
```

Run `plot.py` with the Python of `../requirements.txt`. `run.sh` needs Julia and the environment of `../Project.toml` (Honghu.jl v0.1.0; see `../README.md`). Both scripts work from any directory.

## Files

| File | Made by | Content |
|---|---|---|
| `run.jl` | | the simulations and the equal-accuracy bisection |
| `run.sh` | | runs `run.jl` and saves its printout in `metadata.txt` |
| `scores.csv` | `run.jl` | errors and sub-steps against $\phi_0$, and the RSLA $\phi_0$ of equal accuracy |
| `profiles.csv` | `run.jl` | profiles at $\phi_0 = 0.01$ |
| `groups.csv` | `run.jl` | the four groups: edges, $\kappa$, $\tau_{\rm c}$ |
| `budget.csv` | `run.jl` | reference: share of each group in the flux and in the heating of the gas |
| `attribution.csv` | `run.jl` | error with $\hat c = 0.01c$ in group 4 only, groups 1–3 only, all groups |
| `metadata.txt` | `run.sh` | header with date and versions; the same numbers as text, and $\phi$ of each group |
| `plot.py` | | reads the three CSV files, writes the figure |
| `marshak_vaytet.pdf`, `.png` | `plot.py` | the figure; `plot.py` also copies the PDF to `../figures/` |

- `scores.csv` columns: `phi0`, `L1_rsla`, `L1_arsa` (fractions), `sub_rsla`, `sub_arsa`, `phi0_rsla_matched`, `L1_rsla_matched`, `sub_rsla_matched`, `speedup` (`sub_rsla_matched / sub_arsa`).
- `profiles.csv` columns: `run` (`ref`, `rsla`, `arsa`), `x` [cm], `T_gas` [K], `T_g1`–`T_g4` [K], `phi_g1`–`phi_g4`.
- `groups.csv` columns: `group`, `nu_lo`, `nu_hi` [Hz], `kappa` [cm$^2$ g$^{-1}$], `tau_c`.
- `budget.csv` columns: `quantity` (flux at a position at $t_{\rm end}$, or gas heating over the run behind or ahead of the front at $T_{\rm gas} = 650$ K), `x` [cm], `share_g1`–`share_g4` [%].
- `attribution.csv` columns: `reduced_in`, `phi_thin` ($\phi$ of group 4), `phi_thick` ($\phi$ of groups 1–3), `L1` (fraction), `substeps`.

## The test

The non-equilibrium multigroup Marshak wave of Vaytet et al. (2011), Honghu.jl's `marshak_vaytet`:

- Setup: a 3D slab of $64\times2\times2$ cubic cells over 20 cm, periodic in $y$ and $z$, 4 log-spaced groups from $6\times10^{10}$ to $6\times10^{14}$ Hz, $\kappa = 2000\,(\nu/4\times10^{13}\,{\rm Hz})^{-2}$ cm$^2$ g$^{-1}$, static medium, to $t = 1.36\times10^{-7}$ s.
- Steps: the RSLA and the ARSA at $\nu_0 = 0.9/\sqrt{3}$ per direction ($C_r = 0.9\sqrt{3}$); the $\hat c = c$ reference at $\nu_0 = 0.4$ per direction. The set-up is a copy of Honghu.jl's `marshak_vaytet` with a transverse-cell option; with one transverse cell it reproduces `marshak_vaytet` exactly.
- **Equal accuracy:** for each $\phi_0$, a 14-step bisection in $\log\phi_0$ finds the larger $\phi_0$ at which the RSLA has the error of the ARSA. The speedup is the ratio of sub-steps.
- **Budget and attribution:** the heating of the gas by group $g$ is $c\rho\kappa_g(E_g - f_g a_r T^4)$, summed over the steps of the reference. The attribution runs use the step of the reference, so that the time step adds no error; they fix $\phi$ per group with the limiter `GroupPhi` defined in `run.jl` (group 4 is the group with $\rho\kappa\Delta x < 1$).
- **$T_{{\rm rad},g}$:** the brightness temperature of the group, the $T$ with $a_r T^4 f_g(T) = E_g$, where $f_g$ is the Planck fraction of the group.

**Caveat:** the regular steps in $T_{{\rm rad},4}$, about 4 cells apart, appear in the $\hat c = c$ reference too. They come from the transport of the free-streaming group, not from the ARSA; their cause has not been checked.

These runs replace the sections `vaytet` and `equal_accuracy` (Marshak rows) of `irsla_paper.jl` for this figure. The results agree with `data/vaytet_scores.csv` and `data/equal_accuracy.csv` to all printed digits.
