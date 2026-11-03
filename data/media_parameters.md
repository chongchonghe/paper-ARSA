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


---

## C. Bases of the convective envelopes of stars (added 2026-11-03)

Three rows of Table 1, all fully ionized ($c_v = 2\times10^8$ erg g$^{-1}$ K$^{-1}$, $\mu = 0.6$), at the base of the convective envelope, where the radiation alone just carries the luminosity $L_\star$ along the adiabatic gradient ($\nabla_{\rm rad} = \nabla_{\rm ad} \approx 0.4$). This fixes the opacity for given $T$, total pressure $P$ and enclosed mass $M_r$:

$$\kappa = \frac{16\pi a c G M_r T^4 \nabla_{\rm ad}}{3 L_\star P}.$$

The size is the pressure scale height $H_p = P r^2/(\rho G M_r)$, so $\mathcal{L} = \kappa\rho H_p$, and the velocity is set by the convective flux, $v \approx [L_\star/(4\pi r^2\rho)]^{1/3}$. $\mathcal{L}$ and $v$ are good to an order of magnitude only, mainly through the radius $r$ of the base; $\mathcal{D} \propto r^{4/3}$. Every row has $\mathcal{D} \gg \phi_0$, so this does not move any of them between regions of Fig. 1.

| System | $T$ [K] | $\rho$ [g cm$^{-3}$] | $M_r$ [$M_\odot$] | $L_\star$ [$L_\odot$] | $r$ [$R_\odot$] | $\kappa$ [cm$^2$ g$^{-1}$] | $H_p$ [cm] | $\mathcal{L}$ | $v$ [km s$^{-1}$] | $P_{\rm rad}/P_{\rm gas}$ |
|---|---|---|---|---|---|---|---|---|---|---|
| Sun | $2.2\times10^6$ | 0.2 | 0.98 | 1 | 0.713 | 20 | $5.7\times10^9$ | $2.3\times10^{10}$ | 0.09 | $10^{-3}$ |
| RGB, $1\,M_\odot$, lower branch | $2\times10^6$ | $3\times10^{-2}$ | 0.3 | 30 | 0.5 | 1.0 | $8.4\times10^9$ | $2.5\times10^8$ | 0.6 | $5\times10^{-3}$ |
| AGB with hot bottom burning | $7\times10^7$ | 3 | 1.0 | $5\times10^4$ | 0.064 | 0.34 | $4.5\times10^9$ | $4.5\times10^9$ | 6 | 2 |

- **Sun:** $T$, $\rho$ and $r$ at the base of the convection zone from the standard solar model of Bahcall, Pinsonneault & Basu (2001, ApJ 555, 990), $T \approx 2.2\times10^6$ K and $\rho \approx 0.19$ g cm$^{-3}$ at $0.71\,R_\odot$. The relation above gives $\kappa = 20$, as OPAL gives there; electron scattering alone ($\kappa = 0.34$) would let the radiation carry 60 $L_\odot$.
- **RGB:** $T$ and $\rho$ in the ranges $(1$–$5)\times10^6$ K and $10^{-3}$–$10^{-1}$ g cm$^{-3}$ at the base of the envelope (given by the user, 2026-11-03); $\rho = 3\times10^{-2}$ chosen by the user. With $M_r$ equal to the helium core mass and $L_\star = 30\,L_\odot$, the relation gives $\kappa = 1$, the lower end of the range $1$–$10$ given by the user. $r = 0.5\,R_\odot$ is assumed.
- **AGB:** hot bottom burning needs the base of the envelope above $5\times10^7$ K, up to $10^8$ K in the most massive AGB stars (Karakas & Lattanzio 2014, PASA 31, e030; Doherty et al. 2017, PASA 34, e056, who give $\rho \sim 10$ g cm$^{-3}$ near $10^8$ K for super-AGB stars). The gas is fully ionized, so $\kappa$ is close to electron scattering. A density $10^{-3}$–$10^{-1}$ g cm$^{-3}$ at this temperature would make $P_{\rm rad}/P_{\rm gas} \sim 600$ and $H_p$ larger than the star; $\rho = 3$ gives $P_{\rm rad}/P_{\rm gas} = 2$ and, with $M_r = 1\,M_\odot$ and $L_\star = 5\times10^4\,L_\odot$, $\kappa = 0.28$ from the relation, consistent with electron scattering. $r$ is set where $H_p = r$. The low-mass and the hot-bottom-burning AGB rows given by the user both have $\hat{h} \gg 0.1$; they are merged into this one row.

### C.1 Middle layers (2026-11-03): the Sun and the RGB star in Table 1

Table 1 lists the Sun and the RGB star halfway in radius between the base of the convective envelope and the surface, and keeps the AGB star at the base. `src/calc_convective_layers.py` (printout: `data/convective_layers.txt`) computes the layers:

- **Sun:** Model S (Christensen-Dalsgaard et al. 1996, Science 272, 1286) at $r = 0.85\,R_\odot$: $T = 9.6\times10^5$ K, $\rho = 0.054$ g cm$^{-3}$.
- **RGB star** ($1\,M_\odot$, $30\,L_\odot$, $T_{\rm eff} = 4600$ K, so $R = 8.6\,R_\odot$): from the base row above, placed at $r_b = 0.76\,R_\odot$ where a point-mass envelope has $T = 2\times10^6$ K, outward along the adiabat of gas plus radiation, $\ln(T^{3/2}/\rho) + \hat{h}/2 = \mathrm{const}$, with hydrostatic balance, to $r = (r_b + R)/2 = 4.7\,R_\odot$: $T = 1.8\times10^5$ K, $\rho = 8\times10^{-4}$ g cm$^{-3}$. The base radius $0.5\,R_\odot$ assumed in the table above is too small: from it the envelope cools to $10^3$ K before reaching the middle.
- **Why $\hat{h}$ falls outward:** along the adiabat, while $\hat{h} \ll 1$, $T^{3/2}/\rho$ is constant, so $\hat{h} \propto T^3/\rho \propto \rho$. Where $\hat{h} \gg 1$ the radiation carries most of the entropy and $\hat{h}$ hardly falls: from the base of the AGB envelope, $\hat{h} = 17$, 12, 6.3 and 1.7 at $T = 7\times10^7$, $10^7$, $10^6$ and $10^5$ K. The AGB star stays in the red region at any layer, and a point-mass envelope does not describe it (the envelope holds most of the mass), so Table 1 keeps its base.
- **Opacity:** electron scattering plus a Kramers term, $\kappa = 0.34 + K\rho T^{-7/2}$, with $K$ set by $\kappa = 20$ at the base of the solar convection zone: $\kappa \approx 100$ (Sun) and 500 (RGB star), uncertain by a few times (the Kramers law misses the iron bump near $2\times10^5$ K). The RGB layer lies near the He II ionization zone, so its $c_v$ and $\hat{h}$ are uncertain by a similar factor.
- $H_p = P r^2/(\rho G M_r)$ and $v = [L_\star/(4\pi r^2\rho)]^{1/3}$ as at the base: $\mathcal{L} = 1.9\times10^{10}$ and $2.8\times10^{10}$, $v = 0.12$ and 0.5 km s$^{-1}$, $\mathcal{D} = 7\times10^3$ and $4\times10^4$, $\hat{h} = 2.5\times10^{-3}$ and $1.1\times10^{-3}$.
