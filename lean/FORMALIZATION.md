# Formalization evidence

## Source and reproducibility

- Source: `Improved_RSLA (4).pdf`, 61 PDF pages, dated 23 October 2026.
- Title: *An Adaptive Reduced Speed of Light Approximation*.
- Initial source snapshot SHA-256: `a3b0ad1780748687684fd631ae0018c94df18cccc66b73da00fd683b91e26e44`.
- The source PDF was not edited by this work. Its on-disk fingerprint changed
  during the session; a comparison against a later snapshot found identical
  normalized extracted text on all 61 pages, with highlights visible in the
  later PDF. The fingerprint above identifies the initial snapshot.
- Lean: 4.19.0, commit `6caaee842e94`, arm64-apple-darwin23.6.0 release.
- mathlib: v4.19.0, commit `c44e0c8ee63ca166450922a373c7409c5d26b00b`.
- All transitive package revisions are recorded in `lake-manifest.json`.

## Completed proof structure

The inventory contains 11 named results. The complete paper-to-declaration map is
in README.md. The development includes:

1. The closed-disk Jury criterion, including repeated real and nonreal roots, and
   its characteristic-determinant formulation for arbitrary real 2×2 matrices.
2. The noncommutative determinant factorization in arbitrary finite dimension.
3. The actual energy-block amplification matrix, all three determinant identities,
   their signs, and both directions of the exact spectral stability boundary.
4. The trace/Cauchy–Schwarz CFL bound, its converse, and a normalized M1 streaming
   witness attaining equality.
5. The explicit complex multidimensional Fourier symbol and its squared-norm
   estimate; the gas-variable square completion; and both directions of the
   all-mode one-stage contraction criterion.
6. The concrete two-stage update and its exact energy-loss identity, yielding a
   common weighted length for every phase.
7. Both directions of the all-mode spectral stability theorem, with checkerboard
   necessity proved by embedding a channel eigenvector in the full state.
8. A Euclidean embedding of the common weighted length, proving that it satisfies
   positive definiteness, absolute homogeneity, and the triangle inequality.
9. Strict damping in one dimension, proved through the absence of stage
   eigenvalues +1 and -1 and the strictness term in the full-step energy identity.

## Source clarification

Lemma B.3 has an explicit positive transport assumption `a > 0` in Lean. A checked
zero-transport neutral example documents the need for this condition when q, K,
and a are treated as independent parameters. Positive channel couplings and the
exclusion of the density mode are also explicit. See README.md for the model
boundary and coordinate scaling.

## Validation

Each module was compiled as it was developed. The final verification command is:

```sh
./scripts/check.sh
```

This builds every theorem through `RSLA.lean` with `autoImplicit = false` and then
runs a fresh transitive dependency audit of every theorem in namespace `RSLA`.
The audit rejects all axioms except `propext`, `Classical.choice`, and `Quot.sound`.
The original-paper audit covered 148 RSLA theorems. The subsequent full M1
diffusion-scaling extension is documented in M1_FORMALIZATION.md.
The checked source contains no proof placeholders, additional axioms, unsafe
proof code, or native-decision shortcuts.

Intermediate compiler diagnostics and extracted/rendered PDF pages are in the
ignored `tmp/` directory. Symbolic-algebra checks were used to guide the proof
construction; identities in the Lean development are proved independently by Lean.
The later singular-tensor Jordan-block section of the informal note is separate:
its identities are checked with exact SymPy arithmetic, and its periodic-grid
experiment is executed by `examples/jordan_block.py` (not yet in this repository). It is not part of the Lean
coverage or the 148-theorem audit.

## Environment notes

The toolchain was installed under `.tools/elan` without changing the user's shell
configuration. The generated mathlib cache executable hit a macOS linker error;
running `Cache/Main.lean` through the Lean interpreter downloaded the cache.
Transient disk exhaustion required stopping dependency rebuilds and restoring the
complete downloaded cache. The final proofs use the pinned, restored libraries.

## Full M1 extension (24 October 2026)

`RSLA.M1Full.finite_time` and `RSLA.M1UnequalFull.finite_time` assemble
the equal- and unequal-opacity diffusion-scaling results from the
explicit M1 source and HLL symbol. They derive every norm certificate and
the isotropic reference power bound. See M1_FORMALIZATION.md for its exact
algebraic coordinate/inverse hypotheses, module map, and scope boundaries.
The extended axiom audit covers 336 RSLA theorems.
