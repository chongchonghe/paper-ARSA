# Fig. 1: where the RSLA is accurate and where the ARSA can help, with the systems of Table 1

Paper: Fig. 1 (`fig:accuracy`), §2.3 (`ssec:heat_capacity`), with eq:Pe_phi; the systems are those of Table 1 (`tab:media`, §2.3); also meant for talks and posters. `plot.py` copies the PDF to `../figures/rsla_accuracy.pdf`. The paper quotes no number from this folder.

## Reproduce

```bash
$PY plot.py       # figure only, from the data in this folder; about 2 s
./run.sh         # the curve and the systems: rewrites systems.csv, curves.csv and metadata.txt; about 2 s
```

Run both with the Python of `../requirements.txt`; `run.sh` uses `$PY` (default `python3`). Both scripts work from any directory. `run.py` imports the systems from `../table_media.py` (Table 1).

## Files

| File | Made by | Content |
|---|---|---|
| `parameters.md` | given | the star-formation and disk systems (S1–S6) that `../table_media.py` uses |
| `run.py` | | the curve and the position and region of each system |
| `run.sh` | | runs `run.py` and saves its printout in `metadata.txt` |
| `curves.csv` | `run.py` | the curve $\mathrm{Pe}_\phi = 0.1$: `h`, `tau_v_c_Pe` |
| `systems.csv` | `run.py` | one row per system: `system`, `tau` (optical depth across the system; `nan` for the diffuse ISM, $\tau \ll 1$), `tau_v_c`, `h`, `Pe_phi`, `slowdown` ($t_\phi/t_1$, eq:slowdown), `region`; `nan` for the solar core ($v \approx 0$) and the IGM ($\tau \ll 1$) |
| `metadata.txt` | `run.sh` | header with date and version; the same numbers as text |
| `plot.py` | | reads `curves.csv` and `systems.csv`; writes the figure |
| `rsla_accuracy.pdf`, `.png` | `plot.py` | the figure |

## The diagram

Axes: the dynamic-diffusion number $\mathcal{M}\mathcal{L}/\mathcal{C} = \rho\kappa\ell v/c$ of each system and its pressure ratio $\mathcal{P} = 3\hat{h}/8$, i.e. the heat-capacity ratio $\hat{h}$ of Table 1 with $\gamma = 5/3$. The curve and the regions are computed in $\hat{h}$ and drawn for the smallest reduction factor of the RSLA considered, $\phi = 10^{-3}$ (`PHI0` in the scripts), the worst case: a larger $\phi$ is more accurate everywhere.

- **$\mathrm{Pe}_\phi = 0.1$:** by eq:T_diffusion, a temperature structure of size $\ell$ and optical depth $\tau$ diffuses away under the RSLA in $t_\phi = 3\tau\ell\,(1/\hat{h} + 1/\phi)/c$; the flow crosses it in $\ell/v$. Their ratio is the Péclet number of the heat transport under the RSLA, $\mathrm{Pe}_\phi = 3(\tau v/c)(1/\hat{h} + 1/\phi)$ (eq:Pe_phi), and $\mathrm{Pe}_\phi = 0.1$ is the curve $\tau v/c = 0.1/[3(1/\hat{h} + 1/\phi)]$: $\tau v/c \approx \hat{h}/30$ for $\hat{h} \ll \phi$ and $3.3\times10^{-5}$ for $\hat{h} \gg \phi$. Left of it (white), even the slowed diffusion is ten times faster than the flow, and the temperature tracks the flow for any $\phi \geq 10^{-3}$. The error of the RSLA then does not matter if the flow is the fastest driver of the temperature; transient heating faster than the slowed diffusion, or evolution driven by diffusion itself (gas at rest), still suffers the slowdown of eq:slowdown in optically thick gas.
- **Right of it**, by the slowdown $t_\phi/t_1 = (1 + \hat{h}/\phi)/(1 + \hat{h})$ of eq:slowdown:
  - green, $\hat{h} < 0.1\,\phi = 10^{-4}$: the error of the RSLA, about $\hat{h}/\phi$, is below 10% for any $\phi \geq 10^{-3}$;
  - amber, $10^{-4} < \hat{h} < 0.1$: the ARSA can help. Run with a reduction factor $\phi_0$ above $\hat{h}$ (up to $0.1$), it runs thick cells at $\phi = 1$; at $\phi_0$ below $\hat{h}$, $\phi$ saturates at $\phi_0/[1 - \phi_0/(\beta \hat{h})]$, which cuts the error by the fraction $\phi_0/(\beta \hat{h})$;
  - red, $\hat{h} > 0.1$: the RSLA is inaccurate for any useful $\phi$ (eq:rsla_accurate with $\epsilon = 0.1$ needs $\phi > 0.48$ at $\hat{h} = 0.1$), and the ARSA cannot help.

## The systems

From Table 1: `../table_media.py` holds them: order-of-magnitude estimates and, for the star-formation and disk rows, the values of `parameters.md`. For those rows the geometry of `parameters.md` ($\tau = \kappa\rho L$, $v$ in the laboratory frame, Keplerian in the disk) is used as given, but $\hat{h}$ is recomputed with the specific heats of Table 1 ($8\times10^7$ erg g$^{-1}$ K$^{-1}$ above 100 K, $5\times10^7$ at and below), since `parameters.md` takes $c_v = k/(\mu m_{\rm H})$, too small for H$_2$–He gas; its $\hat{h}$ is higher by a factor of 1.4–2.3. The PPD midplanes lie in the green region with either specific heat ($\hat{h} = 4\times10^{-6}$ and $8\times10^{-6}$ here, $9\times10^{-6}$ and $2\times10^{-5}$ in `parameters.md`).

## Caveats

- $\mathrm{Pe}_\phi$ compares the diffusion time with the time in which the flow crosses a structure of optical depth $\tau$, taken as the whole system; smaller structures have smaller $\tau$ and smaller $\mathrm{Pe}_\phi$.
- Optically thin systems ($\tau \leq 1$, open symbols) lie outside the diffusion argument behind eq:T_diffusion and $\mathrm{Pe}_\phi$: there the RSLA delays radiation by the light-crossing time $\ell/\hat c$, and the relevant check is $v \ll \hat c$ and $\ell/\hat c$ short compared with the time scales of the problem. All systems of Table 1 left of the curve are of this kind (cold molecular cloud, PPD atmosphere, PPD at 100 au, star-forming clump), as is the hot Jupiter atmosphere ($\tau = 1$, green). Suggested by B. Wibking, 2026-10-16.
- "Thick cells": the amber region marks where the ARSA reaches $\phi = 1$ in cells thicker than $\tau_{\rm thresh} \approx 3.5$; systems in it with thin cells (the hot Jupiter atmosphere) do not gain.
- The solar core has $v \approx 0$: no flow, so $\mathrm{Pe}_\phi \to 0$ on a dynamical time; the relevant time is then the thermal time of the star.
- The positions inherit the uncertainties of Table 1 (a factor of 2–3) and its caveats (electron scattering in $\kappa$, relativistic speeds).
