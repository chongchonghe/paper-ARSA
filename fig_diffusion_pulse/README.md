# Fig. 6: thick pulse at rest and advected

Paper: Fig. 6 (`fig:diffusion_pulse`), §5.3 (`ssec:test_thick`), the two pulse rows of Table 2 (`tab:tests`), and the realizability paragraph of §4.

## Reproduce

```bash
$PY plot.py       # figure only, from the data in this folder; about 2 s
./run.sh         # simulations: rewrites profiles.csv and metadata.txt; about 30 min
```

Run `plot.py` with the Python of `../requirements.txt`. `run.sh` needs Julia and the environment of `../Project.toml` (Honghu.jl v0.1.0; see `../README.md`). Both scripts work from any directory.

## Files

| File | Made by | Content |
|---|---|---|
| `run.jl` | | the simulations |
| `run.sh` | | runs `run.jl` and saves its printout in `metadata.txt` |
| `profiles.csv` | `run.jl` | $T$ at the final time for every run (columns below) |
| `metadata.txt` | `run.sh` | header with date and versions; derived variables and errors of every run |
| `plot.py` | | reads `profiles.csv`, writes the figure |
| `diffusion_pulse.pdf`, `.png` | `plot.py` | the figure; `plot.py` also copies the PDF to `../figures/` |

`profiles.csv` columns: `T_lo` [K], `v0` [cm/s], `t` [s] (final time), `phi0`, `x` [cm], `T_initial`, `T_ref` ($\hat c = c$), `T_rsla`, `T_arsa` [K], `phi_arsa`. Each (`T_lo`, `v0`, `phi0`) is one set of runs. The figure uses `T_lo = 1e6` and `phi0 = 0.01`.

## The test

The advecting pulse of He et al. (2024a), Honghu.jl's `radhydro_pulse_dyn`, with a background ten times colder, so that the gas holds the heat capacity:

- Gaussian pulse $T = T_{\rm lo}[1 + e^{-x^2/2w^2}]$, $w = 24$ cm, in gas plus radiation pressure balance, $\rho = 1.2$ g cm$^{-3}$ in the background.
- 3D slab: periodic box of 1024 cm with $128\times2\times2$ cubic cells, the same in every column, $\kappa = 500$ cm$^2$ g$^{-1}$ ($\tau_{\rm c} = 2400$–$4800$), $\mu = 2.33$, $\gamma = 5/3$, Levermore closure, $v/c$ terms on.
- Steps: the RSLA and the ARSA at $\nu_0 = 0.9/\sqrt{3}$ per direction ($C_r = 0.9\sqrt{3}$), the 3D RSLA step for isotropic radiation; the references at $\nu_0 = 0.4$ per direction.
- Two kinds of run:
  - **Moving** ($v_0 = 3\times10^7$ cm/s, to $t = 4.8\times10^{-5}$ s, 1.4 box crossings): at $T_{\rm lo} = 10^5$ and $10^6$ K, and at $\phi_0 = v_0/c$, $2v_0/c$ and $0.01$.
  - **At rest** ($v_0 = 0$, to $t = 4.8\times10^{-4}$ s): at $T_{\rm lo} = 10^6$ K and $\phi_0 = 0.01$.
- **Error:** relative to the pulse, $\sum|T - T_{\rm ref}| / \sum(T_{\rm ref} - T_{\rm lo})$.
- **Flux-limited variant:** the end of `metadata.txt` also has the ARSA with $|F|$ clipped to $cE$ before the exchange update, which §4 quotes.

**Caveat:** a centroid lag in `metadata.txt` is meaningless when the RSLA pulse falls more than half a box behind (the runs at $\phi_0 = v_0/c$ and $2v_0/c$). The L1 errors of those runs are valid.
