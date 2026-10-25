# Full M1 diffusion-scaling theorem in Lean

The main declarations are `RSLA.M1Full.finite_time` and
`RSLA.M1UnequalFull.finite_time`. They prove the finite-time estimates
in Theorems 1 and 2 of
`m1_stability_extension.tex` for the specified frozen linearized model.
The source and physical-flux matrices are explicit definitions; no
amplification bound, reference power bound, source-resolvent norm bound,
or M1 flux Lipschitz estimate is an assumption of this declaration.

The theorem uses Lean 4.19.0 and the pinned mathlib v4.19.0. Run
`./scripts/check.sh` to build it and audit the transitive axioms of all
RSLA theorems. In the Overleaf checkout this project is under `lean/`.

## Mathematical statement

With q = dt chat alpha, alpha = rho kappa, and reduced flux f, assume:

- a static frozen background, equal gray opacities, and the stated full
  first-order velocity-perturbation couplings;
- positive reference channel ratios, fixed material response ratios,
  eta = B0/E0 >= 0, and a Levermore closure coordinate 0 <= b < 1;
- alpha f <= B and alpha dx_i >= tau > 0;
- nonnegative directional weights nu_i <= dt chat/dx_i;
- sum_i nu_i/sqrt(3) <= 1 + q/(2 + xi_max q).

Then the powers of the full PD-ARS amplification operator obey

    ||G^n|| <= K0 exp(C chat B n dt).

The displayed Lean constants contain no opacity, grid spacing, timestep,
Fourier phase, or beam direction. They depend on fixed material data,
tau, dimension, and a fixed change of coordinates. Inactive directions
can have nu_i = 0, so transverse moment components need not be removed.
The enlarged limiting CFL is 1 + 1/xi_max. This is a finite-time result,
not a claim that ||G|| <= 1 or that the exponential constant is sharp.

## Unequal opacity means

The original static P1 theorem already accepts distinct channel rates.
The new M1 theorem uses alpha_F = rho kappa_F, fixed positive ratios
p = kappa_P/kappa_E and s = kappa_F/kappa_E, and

    gamma = 1 + (p eta - 1)/s >= 0,
    q_F = dt chat alpha_F,    q_E = q_F/s,
    xi_E = p mu_E,           xi_F = (gamma + 1/3) mu_M.

Here mu_E and mu_M are normalized to the energy and momentum radiation
relaxation rates separately. Energy balance means p eta = 1, so gamma = 1.
The sufficient gamma >= 0 condition also covers some unbalanced states;
positive opacity means alone do not imply it.

Replace the equal-opacity CFL by

    L0 <= min(1 + q_E/(2 + xi_E q_E), 1 + q_F/(2 + xi_F q_F)),

and impose alpha_F f <= B and alpha_F dx_i >= tau. The finite-time
estimate is uniform in the common opacity scale. Its constants may depend
on the fixed opacity ratios; no uniform limit in zero or infinite ratios
is claimed. The thick-limit enlargement is 1 + 1/max(xi_E,xi_F).

The explicit source uses

    C = [[p, (2-s) f n^T], [0, gamma I + T]],
    Rbar = diag(1/s, 1, 1, 1),
    Q = q_F [[M Rbar C, -M Rbar], [-Rbar C, Rbar]].

All first-order velocity couplings at zero background velocity are retained.
In the Lean declarations `rateE = 1/s`, `ell = 2-s`, and `gamma` is the
coefficient above. The fixed reference coordinates use
C0 = diag(p, gamma+1/3, ...), with channel ratios xi_E and xi_F.
The Lean power-bound constant conservatively uses Lstar = 1+1/xi_E;
the note's sharper constant using max(xi_E,xi_F) follows by taking the
stricter of the two channel caps. The timestep assumption is identical.
The algebraic identification boundary below applies to both declarations.

## Explicit algebraic inputs and interpretation

The formal theorem takes the actual source inverses as continuous linear
maps satisfying their left/right inverse identities. It also takes a
fixed continuous linear equivalence S identifying the isotropic reference
with the already verified P1 channels. These are algebraic identification
hypotheses, not stability hypotheses. The isotropic source and transport
identities are named `hQ₀` and `hD₀`; `hr₀action` identifies the explicit
channel formula for the reference source solve.

For Cartesian moment components choose

    C0 = diag(1, eta + 1/3, ...),
    P  = diag(1, sqrt(3), ...),
    u  = P C0 m,   v = P r,
    xi_none = mu_E,   xi_some = (eta + 1/3) mu_M.

Apply this transformation to real and imaginary amplitudes and combine
them into complex u and v. Then S is independent of opacity, beam direction,
phase, and timestep. The reference norm is

    sum_k (|u_k|^2/xi_k + |v_k|^2).

The theorem quantifies over any S satisfying the identification equations;
it does not construct this particular Cartesian equivalence inside Lean.
Likewise it verifies stability of the defined linear operators, not that
Quokka's C++ implementation computes those operators. The note documents
the linearization, and independent finite differences check the physical
and two-state HLL flux Jacobians. These boundaries should not be confused
with a fully verified nonlinear solver or a Lean verification of C++.

## Proof map

| Module | Kernel-checked content |
| --- | --- |
| `M1Coefficients` | Levermore coefficient bounds; positive Eddington tensor; full physical-flux Jacobian isotropic-limit bound; HLL speed floor estimate. |
| `M1Source` | Exact source elimination, retaining velocity couplings; uniform full and half source-resolvent bounds; source perturbation bound. |
| `M1Fourier` | All-phase HLL symbol perturbation estimate, with the physical Jacobian and HLL diffusion treated separately. |
| `M1Coordinates` | Real/imaginary duplication, zero material-block embedding, Fourier rotation, and fixed-coordinate norm transfer. |
| `M1Model` | Explicit physical source and Fourier symbols, combining the concrete M1 estimates. |
| `M1Reference` | Original ARSA CFL implies the actual P1 all-phase stage and full-step energy inequalities. |
| `M1ReferenceNorm` | Genuine opacity-independent Hilbert norm; exact channel-solve and full-update identities; all-phase reference power bound. |
| `M1Energy` | First-step estimate removing the opacity-dependent energy-norm comparison; uniform reference power bound. |
| `M1Diffusion` | Noncommutative power perturbation, exact inverse/PD-ARS difference identities, and diffusion scaling. |
| `M1Theorem` | Abstract operator assembly and explicit-CFL assembly with primitive norm certificates. |
| `M1Full` | Concrete equal-opacity M1 assembly deriving all norm certificates from the model formulas. |
| `M1UnequalReference` | Distinct relaxation rates in the all-phase P1 energy and fixed-norm power estimates. |
| `M1UnequalSource` | Explicit unequal-opacity source, two-rate elimination, uniform implicit-solve bounds, and source perturbation. |
| `M1UnequalFull` | Concrete unequal-opacity assembly deriving every norm certificate. |

The reference power proof is not a spectral-radius argument. It controls
operator powers directly and therefore does not overlook Jordan blocks.
The perturbation proof does not assume that source and transport commute,
that the M1 operator is normal, or that the anisotropic source shares the
transport symmetrizer.

## Scope not covered

This formalization covers the diffusion-scaling finite-time theorem.
The focused note also contains growing-mode and checkerboard counterexamples.
They have analytic and symbolic/numerical validation and are not claimed as
Lean theorems. Earlier versions additionally discussed compatible norms,
absorption-only bounds, and less uniform unequal-opacity bounds; their supporting scripts and
results remain in the repository history and validation package.
Neither changing frozen coefficients
between steps nor nonlinear Quokka streaming-transient experiments are
proved or executed here. The note specifies the latter experiments to run.

## Independent validation of unequal means

`../checks/m1/m1_unequal_diffusion.py` checks the full nonequilibrium four-force
Jacobian, separate-rate implicit elimination, and fixed isotropic channel
identification symbolically. It also checks 120 finite-difference source
Jacobians and 1080 cases each of the reference energy bound and full-M1
perturbation estimate, across 30 fixed-ratio states and opacity scales
from 1 to 10^6. These finite numerical samples supplement the proof;
they are not a proof over all Fourier phases. Recorded results are in
`../data/m1/unequal_diffusion.json`.
