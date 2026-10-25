import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-! Lemma B.1 in a real inner product space. Complex vector spaces are covered
by restriction of scalars to ℝ, so the quadratic form is the real Hermitian form.
No square root of an operator or spectral theorem is needed for this proof. -/
namespace RSLA
open scoped RealInnerProductSpace

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The B-weighted squared energy. -/
def energy (B : E →ₗ[ℝ] E) (x : E) : ℝ := ⟪x, B x⟫

lemma energy_cross_bound (B : E →ₗ[ℝ] E)
    (hsym : ∀ x y, ⟪x, B y⟫ = ⟪y, B x⟫)
    (hpos : ∀ x, 0 ≤ energy B x) (x y : E) :
    2 * ⟪x, B y⟫ ≤ energy B x + energy B y := by
  have h := hpos (x-y)
  simp only [energy, map_sub, inner_sub_left, inner_sub_right] at h ⊢
  rw [hsym y x] at h
  linarith

/-- Lemma B.1. `T = B⁻¹D`; the equation for G is precisely
`(I+B)G = I + DB⁻¹D = I + B T²`. -/
theorem two_stage_contract (B : E →ₗ[ℝ] E) (T G : E → E)
    (hsym : ∀ x y, ⟪x, B y⟫ = ⟪y, B x⟫)
    (hpos : ∀ x, 0 ≤ energy B x)
    (hstage : ∀ x, energy B (T x) ≤ energy B x)
    (hupdate : ∀ x, G x + B (G x) = x + B (T (T x))) (x : E) :
    ‖G x‖^2 + energy B (G x) ≤ ‖x‖^2 + energy B x := by
  have htwice := (hstage (T x)).trans (hstage x)
  have hcross := energy_cross_bound B hsym hpos (G x) (T (T x))
  have hord : 2 * ⟪G x, x⟫ ≤ ‖G x‖^2 + ‖x‖^2 := by
    have hh : 0 ≤ ⟪G x - x, G x - x⟫ := real_inner_self_nonneg
    simp only [inner_sub_left, inner_sub_right, real_inner_self_eq_norm_sq] at hh
    linarith [real_inner_comm x (G x)]
  have heq := congrArg (fun y => ⟪G x, y⟫) (hupdate x)
  simp only [inner_add_right, real_inner_self_eq_norm_sq] at heq
  change ‖G x‖^2 + energy B (G x) = ⟪G x, x⟫ + ⟪G x, B (T (T x))⟫ at heq
  linarith

/-- The common weighted length in Lemma B.1 is also nonincreasing. -/
theorem two_stage_weighted_norm (B : E →ₗ[ℝ] E) (T G : E → E)
    (hsym : ∀ x y, ⟪x, B y⟫ = ⟪y, B x⟫)
    (hpos : ∀ x, 0 ≤ energy B x)
    (hstage : ∀ x, energy B (T x) ≤ energy B x)
    (hupdate : ∀ x, G x + B (G x) = x + B (T (T x))) (x : E) :
    Real.sqrt (‖G x‖^2 + energy B (G x)) ≤ Real.sqrt (‖x‖^2 + energy B x) :=
  Real.sqrt_le_sqrt (two_stage_contract B T G hsym hpos hstage hupdate x)

end RSLA
