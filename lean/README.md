# Lean formalization of the ARSA paper

Formalization of all 11 named mathematical results in **An Adaptive Reduced Speed
of Light Approximation**, from `Improved_RSLA (4).pdf` (23 October 2026).
The proofs use Lean 4.19.0 and mathlib v4.19.0, pinned by `lean-toolchain` and
`lake-manifest.json`.

Written by B. Wibking in the paper's LaTeX project (commits `07a488f` to `fae6805`,
23-24 October 2026) and moved here on 25 October 2026, with his Python checks
(now `../checks/m1/`) and their results (now `../data/m1/`). The companion notes
`simplified_proofs.tex`, `m1_stability_extension.tex` and
`m1_finite_flux_extension.tex` stay in the paper's LaTeX project.

The project also includes the full M1 diffusion-scaling finite-time extension,
for both equal and unequal opacity means.
See [M1_FORMALIZATION.md](M1_FORMALIZATION.md) for its statement, proof map,
and explicit algebraic model-identification hypotheses.

The project contains proofs, not just theorem statements. The build includes a
transitive axiom audit of every theorem in the `RSLA` namespace. Only the standard
foundational axioms `propext`, `Classical.choice`, and `Quot.sound` are allowed;
placeholders and additional axioms fail the audit. Undeclared implicit parameters
are disabled.

## Build and verify

In this prepared workspace:

```sh
./scripts/check.sh
```

The script uses `lake` on PATH, or the project-local `.tools/lake` installed during
formalization. It builds the full library and reruns the axiom audit.

On another machine with [elan](https://github.com/leanprover/elan) installed:

```sh
lake exe cache get
./scripts/check.sh
```

If the cache executable encounters the macOS linker error observed during this
setup, use the interpreted downloader instead:

```sh
lake env lean --run .lake/packages/mathlib/Cache/Main.lean get
```

Import the complete development with `import RSLA`. Generated dependencies,
toolchains, and temporary PDF extracts are ignored by `.gitignore`.

## Paper-to-proof map

PDF page numbers below count the highlights page; printed page numbers are one lower.
All declaration names below are prefixed with `RSLA.`.

| Paper result | PDF page | Lean declaration | Representation |
| --- | --- | --- | --- |
| Proposition 1 | 23 | `cfl_sufficient`, `cfl_iff`, `m1_sharp_witness` | Sharp CFL bound from nonnegative signal speeds and unit trace; a normalized streaming M1 witness attains equality. |
| Lemma A.1 | 47 | `matrix_jury` | Arbitrary real 2×2 matrix; all complex characteristic roots have norm at most one. |
| Lemma A.2 | 47 | `determinant_split` | Arbitrary finite matrix size over any field, with invertible B. |
| Lemma A.3 | 47 | `EnergyBlock.condition_i` | First Jury condition holds exactly up to the boundary. |
| Lemma A.4 | 48 | `EnergyBlock.condition_ii` | Second Jury determinant is strictly positive for every real a. |
| Lemma A.5 | 48 | `EnergyBlock.condition_iii` | The amplification determinant is strictly between zero and one throughout the stable interval. |
| Proposition A.1 | 46 | `EnergyBlock.proposition_A1` | Exact complex spectral stability of the actual 2×2 amplification matrix, using q and λ. |
| Lemma B.1 | 50 | `two_stage_contract`, `Channels.full_contract` | General inner-product proof and a concrete complex-channel proof with an explicit energy-loss term. |
| Lemma B.2 | 50–51 | `AllModes.all_modes_stage_iff` | One-stage contraction for every Fourier phase iff every channel satisfies its checkerboard bound. |
| Proposition B.1 | 49 | `Checkerboard.all_modes_stability_iff`, `AllModes.common_weighted_norm` | Both directions of all-mode complex spectral stability, plus a common phase-independent weighted length. |
| Lemma B.3 | 51–52 | `Strict.lemma_B3` | Strict damping of interior 1D modes, with the positive-step assumption made explicit. |

The energy-block matrix is defined in [RSLA/EnergyBlock.lean](RSLA/EnergyBlock.lean).
Its determinant formulas are derived from matrix multiplication and inversion;
they are not assumed. The complex Fourier symbol is explicitly defined in
[RSLA/StaticTransport.lean](RSLA/StaticTransport.lean), including the energy/flux
couplings in every direction. The necessity proof embeds a complex eigenvector
of an individual checkerboard block into the full state.

## Explicit hypotheses and the strict-damping edge case

Appendix A uses `q > 0`, `K > 0`, and `a ≥ 0`, where `K = λ q`. The named
`proposition_A1` theorem also exposes the paper's original `λ > 0` parameter.
The boundary is `2*(2+q+K)/(2+K)`. `EnergyBlock.courant_iff` additionally verifies
the rearrangement into Equation 62.

Appendix B uses strictly positive `q_k` and `K_k`, nonnegative directional
transport weights `w_d = a_d/2 = ν_d σ_d`, and the static, fixed, diagonal-closure
Fourier symbol. The variables are already rescaled as in the paper. Each exchange
channel's gas variable has been scaled so that its equilibrium is `u=v`. Its
squared B-energy is

```text
(q/K)|u|² + |v|² + q|u-v|².
```

The checkerboard transport parameter is `a = 2∑_d w_d = 2∑_d ν_d σ_d`;
in one dimension this is `a = 2νσ`. One explicit transport stage multiplies the
checkerboard radiation perturbation by `1-a`.

The common squared full-step energy is the sum over channels of

```text
2(q/K)|u|² + 2|v|² + q|u-v|².
```

[RSLA/WeightedNorm.lean](RSLA/WeightedNorm.lean) identifies the square root of this
energy with a Euclidean norm of a linear embedding, and proves positive
definiteness, absolute homogeneity, and the triangle inequality.

The fixed-point argument proves Lemma B.3 with **a > 0**, in addition to positive
channel couplings and `0 < θ < π`. The displayed statement in the paper does not
explicitly state this positive-step condition, although its proof uses it.
`Strict.zero_step_neutral` verifies an exchange-only counterexample in the
independent-parameter model: at `a=0`, `q=K=1`, a nonzero exchange equilibrium has
eigenvalue one at every phase. Thus the formalization makes this implicit
condition explicit. The density perturbation is excluded from Appendix B's
state, as in the paper; its neutral eigenvalue would also preclude strict damping.

## Scope of the result

These are formal proofs about the specified linear amplification model. The
formalization does not derive the frozen linear model from the nonlinear RHD
PDEs, verify a simulation implementation, or prove the empirical accuracy and
speed claims. The trace-based CFL result concerns Equation 66; it does not extend
the all-mode theorem to anisotropic M1 flux Jacobians or moving backgrounds.
No additional hypotheses about the desired stability conclusion are assumed in
the assembled all-mode theorems.

[FORMALIZATION.md](FORMALIZATION.md) records the source fingerprint and verification
process. [RSLA/AxiomAudit.lean](RSLA/AxiomAudit.lean) implements the audit.

## Informal proof note

`simplified_proofs.tex` (in the paper's LaTeX project) gives alternative informal proofs
of Lemmas B.1, A.5, and B.3, including the positive-transport clarification and a
conditional multidimensional consequence. It defines the state, exchange and
transport matrices, stage and full-step operators, and weighted norms before
using them. It distinguishes the von Neumann eigenvalue condition from uniform
power boundedness and foregrounds the common norm proved for static P1.

The note also gives a positive-opacity Jordan-block example for a singular
frozen Eddington tensor, relates it to Levermore's streaming endpoint, and
records the implementation and Git history of Quokka's HLL speed floor.
It also explains how direct stellar irradiation can remain directional in an
optically thick absorption cell, including the distinction between transient
and continuously driven illumination.
This additional section is checked symbolically and numerically, not in Lean.

## Run the Jordan-block example (script not yet in this repository)

With NumPy, SymPy, and Matplotlib available:

```sh
python3 examples/jordan_block.py
```

The script runs the full two-stage update on an 8 by 32 periodic grid for 1000
steps, with strictly positive exchange coefficients. It verifies the exact
Jordan chain and compares the complete grid field and Fourier amplitude against
it. With no transverse speed floor, the radiation-energy amplitude grows by a
factor of 2001. With a speed floor of 0.1 times the reduced light speed, this
particular mode has transient growth followed by decay. The model is a linear
frozen-tensor system using the limiting Rusanov flux, not a Quokka simulation.

Results (not yet in this repository): `output/jordan/results.json`, `growth.csv`
and `growth.png`.
