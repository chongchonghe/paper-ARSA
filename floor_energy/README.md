# Flooring the radiation energy before the exchange update

Paper: §4 (`sec:implementation`), the paragraph on the transport stage leaving $E < 0$: how much energy a floor before the exchange update would add on the supercritical shock (§5.4), the pulse (§5.3) and the Marshak wave (§5.5). No figure.

## Reproduce

```bash
./run.sh         # rewrites results.csv and metadata.txt; about 1 min
```

`run.sh` needs Julia and the environment of `../Project.toml` (Honghu.jl v0.1.0; see `../README.md`). It runs from any directory.

## Files

| File | Made by | Content |
|---|---|---|
| `run.jl` | | the runs |
| `run.sh` | | runs `run.jl` and saves its printout in `metadata.txt` |
| `results.csv` | `run.jl` | one row per case (columns below) |
| `metadata.txt` | `run.sh` | header with date and versions; the same numbers as text |

`results.csv` columns: `case`, `L1_conservative` and `L1_floor` (fractions; error of $T_{\rm gas}$ of the ARSA with Honghu.jl's conservative rule and with the floor everywhere), `events_below_floor` (cell, group and call with $E$ at or below the floor before the exchange update, in the conservative run; cells left at the floor by the post-update repair count too), `substeps`, `floor_energy_fraction` (energy the floor adds, below).

## The computation

- **Set-ups:** the 3D slabs of `fig_diffusion_pulse/` (advected at $\phi_0 = 0.01$, $v_0/c$ and $2v_0/c$, and at rest at $\phi_0 = 0.01$, all at $T_{\rm lo} = 10^6$ K), `fig_marshak_vaytet/` ($\phi_0 = 0.03$, $0.01$, $0.003$) and `mf_shocks.jl` (the supercritical shock at its two $\phi_0$). `run.jl` includes their set-up code, so it uses the same grids and steps.
- **Two rules:** Honghu.jl's `floor_rad_energy!` is replaced by a version with a switch:
  - conservative: Honghu.jl's rule, floor $E$ only where $\rho e + \sum_g E_g c/\hat c_g \leq 0$;
  - floor: raise $E$ to the floor wherever it is at or below it, as the earlier algorithm did.
- **Added energy:** the floor changes the energy that the exchange update conserves by $\sum (E_{\rm floor} - E)\,c/\hat c\,\Delta V$, summed over the floored cells and calls. Because $E_{\rm floor} = a_r \times 10^{-4}$ is negligible, this is the negative radiation energy that the transport stage leaves, weighted by $1/\phi$.
- **Normalization:** the added energy is given as a fraction of the thermal energy of the pulse, $\sum \rho c_v (T_{\rm initial} - T_{\rm lo})\,\Delta V$; of the heat the Marshak slab gains by $t_{\rm end}$ in the reference; or of the kinetic energy that flows into the shock over the run, $\tfrac{1}{2}\rho_0 u^3 t\,A$.
- **References:** the errors use the stored $\hat c = c$ references of the three folders, so no reference is rerun.
