## B. Star formation and accretion disks

The characteristic velocity \(v\) is now defined as the **fluid velocity relative to the laboratory frame**, including orbital motion. For the protoplanetary disks, I adopt Keplerian velocities around a \(1\,M_\odot\) star.

| ID | Astrophysical system | \(\rho\) | \(v\) | \(\kappa\) | \(L\) | \(T\) | \(\mu\) | \(n\) | \(\tau\) | \(\beta\tau\) | \(E_r\) | \(\eta\) | \(h\) | Regime |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| S1 | Massive protostellar envelope | 1e−14 | 3 | 1 | 1.5e15 | 100 D | 2.33 | 2.6e9 | 15 | 1.5e−4 | 7.57e−7 | 1 | 0.085 | Static diffusion |
| S2 | First hydrostatic core | 5e−12 | 1 | 1 | 7.5e13 | 100 G/D | 2.33 | 1.3e12 | 375 | 0.00125 | 7.57e−7 | 1 | 1.7e−4 | Static diffusion |
| S3 | PPD atmosphere, 1 AU | 1e−14 | **30** | 1 | 5e11 | 200 D | 2.33 | 2.6e9 | 0.005 | **5.0e−7** | 4e−5 | 3.3 | 0.68 | Streaming |
| S4a | PPD midplane, 0.1 AU | 1e−7 | **94** | 1 | 3e10 | 1000 G/D | 2.33 | 2.6e16 | 3000 | **0.94** | 7.57e−3 | 1 | 8.5e−6 | Transition |
| S4 | PPD midplane, 1 AU | 1e−9 | **30** | 1 | 5e11 | 280 G/D | 2.33 | 2.6e14 | 500 | **0.050** | 4.65e−5 | 1 | 1.9e−5 | Static diffusion |
| S5 | PPD outer region, 100 AU | 1e−15 | **3** | 1 | 1e14 | 20 D | 2.33 | 2.6e8 | 0.10 | **1.0e−6** | 4e−10 | 0.33 | 0.0068 | Streaming |
| S6 | AGN dusty torus | 1e−18 | **200** | 5 | 3e18 | 300 D | 2.33 | 2.6e5 | 15 | **0.010** | 6.13e−5 | 1 | 2.3e4 | Static diffusion |

**Units:**
- \(\rho\): g cm\(^{-3}\)
- \(v\): km s\(^{-1}\)
- \(\kappa\): cm\(^2\) g\(^{-1}\)
- \(L\): cm
- \(T\): K
- \(n\): cm\(^{-3}\)
- \(E_r\): erg cm\(^{-3}\)

**Velocity assumptions:**

- **S1–S2:** Representative bulk infall velocities; rotation is not explicitly modeled.
- **S3–S5, S4a:** Keplerian orbital velocities, \(v=\sqrt{GM_\star/r}\), with \(M_\star=1\,M_\odot\).
- **S6:** Representative orbital velocity of approximately 200 km/s at 1 pc around a \(10^7\,M_\odot\) SMBH.

All other physical parameters remain unchanged.

The key change is that the PPD midplane at 1 AU now has

\[
\boxed{\beta\tau=0.05}
\]

rather than \(1.7\times10^{-4}\).

It remains in the static diffusion regime, but is considerably closer to the dynamic diffusion boundary.

**Note:** This classification uses the laboratory-frame velocity as requested. For a rotating flow, the significance of radiation advection for the energy evolution still depends on the direction of radiation transport and the spatial structure of the flow.


---

## B. Star formation and accretion disks (updated)

Added **S4a: PPD midplane at 0.1 AU**, assuming a \(1\,M_\odot\) central star. All velocities are measured in the laboratory frame, including orbital motion.

| ID      | Astrophysical system          | \(\rho\) |  \(v\) | \(\kappa\) |    \(L\) |         \(T\) |  \(\mu\) |      \(n\) |  \(\tau\) | \(\beta\tau\) |     \(E_r\) | \(\eta\) | \(h\)      | Regime         |
|---------|-------------------------------|---------:|-------:|-----------:|---------:|--------------:|---------:|-----------:|----------:|--------------:|------------:|---------:|------------|----------------|
| S1      | Massive protostellar envelope |    1e−14 |      3 |          1 |   1.5e15 |         100 D |     2.33 |      2.6e9 |        15 |        1.5e−4 |     7.57e−7 |        1 | 0.085      | Static         |
| S2      | First hydrostatic core        |    5e−12 |      1 |          1 |   7.5e13 |       100 G/D |     2.33 |     1.3e12 |       375 |       0.00125 |     7.57e−7 |        1 | 1.7e−4     | Static         |
| S3      | PPD atmosphere, 1 AU          |    1e−14 |     30 |          1 |     5e11 |         200 D |     2.33 |      2.6e9 |     0.005 |        5.0e−7 |        4e−5 |      3.3 | 0.68       | Streaming      |
| S4      | PPD midplane, 1 AU            |     1e−9 |     30 |          1 |     5e11 |       280 G/D |     2.33 |     2.6e14 |       500 |         0.050 |     4.65e−5 |        1 | 1.9e−5     | Static         |
| **S4a** | **PPD midplane, 0.1 AU**      | **1e−7** | **94** |      **1** | **3e10** | **1,000 G/D** | **2.33** | **2.6e16** | **3,000** |      **0.94** | **7.57e−3** |    **1** | **8.5e−6** | **Transition** |
| S5      | PPD outer region, 100 AU      |    1e−15 |      3 |          1 |     1e14 |          20 D |     2.33 |      2.6e8 |      0.10 |        1.0e−6 |       4e−10 |     0.33 | 0.0068     | Streaming      |
| S6      | AGN dusty torus               |    1e−18 |    200 |          5 |     3e18 |         300 D |     2.33 |      2.6e5 |        15 |         0.010 |     6.13e−5 |        1 | 2.3e4      | Static         |

Units: \(\rho\) in g/cm³; \(v\) in km/s; \(\kappa\) in cm²/g; \(L\) in cm; \(T\) in K; \(n\) in cm⁻³; \(E_r\) in erg/cm³.

### Assumptions for the new 0.1 AU row

For a solar-mass star,

\[
v_K=\sqrt{\frac{GM_\odot}{r}}
\simeq94.2\ \mathrm{km\,s^{-1}}.
\]

At \(T=1000\) K, the approximate hydrostatic scale height is

\[
H=\frac{c_s}{\Omega_K}
\simeq3.0\times10^{10}\ \mathrm{cm}.
\]

Adopting a representative midplane density of \(10^{-7}\ \mathrm{g\,cm^{-3}}\) and \(\kappa_R=1\ \mathrm{cm^2\,g^{-1}}\),

\[
\tau=\kappa_R\rho H=3000,
\]

giving

\[
\boxed{\beta\tau=\frac{94.2}{299792}\times3000\simeq0.94.}
\]

This places the illustrative inner disk **near the static–dynamic diffusion transition**.

The assumed \(1000\) K temperature is plausible for a hot inner disk, although dust sublimation makes the actual opacity strongly temperature- and composition-dependent. Detailed inner-disk models show substantial variation in the thermal structure near the sublimation region. [OUP Academic](https://academic.oup.com/mnras/article/420/2/1541/984786?utm_source=chatgpt.com)

The result is sensitive to opacity and density: increasing either by a factor of ten would give \(\beta\tau\simeq9.4\), formally placing the gas in dynamic diffusion. This makes the inner PPD a potentially interesting transition-regime test for RHD algorithms.


**Added 2026-10-15:** row S4a (PPD midplane, 0.1 AU), given in the conversation of that day.
