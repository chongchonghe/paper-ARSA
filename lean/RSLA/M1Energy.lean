import RSLA.TwoStage
import Mathlib.Analysis.NormedSpace.OperatorNorm.Basic

/-! Remove the opacity-dependent comparison constant from the reference
PD-ARS energy estimate by bounding its first full step separately. -/
noncomputable section
namespace RSLA.M1Energy
open scoped RealInnerProductSpace
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

def sourceEnergy (Q : E → E) (x : E) : ℝ := 2*‖x‖^2 + ⟪x, Q x⟫

theorem resolvent_contract (Q R : E → E)
    (hpos : ∀ x, 0 ≤ ⟪x, Q x⟫) (hres : ∀ x, R x+Q (R x)=x) (x : E) :
    ‖R x‖ ≤ ‖x‖ := by
  have hi : ‖R x‖^2 ≤ ⟪R x, x⟫ := by
    conv_rhs => arg 2; rw [← hres x]
    rw [inner_add_right, real_inner_self_eq_norm_sq]
    linarith [hpos (R x)]
  have hh := hi.trans (real_inner_le_norm (R x) x)
  by_cases hz : ‖R x‖=0
  · rw [hz]; exact norm_nonneg _
  · have hp : 0 < ‖R x‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hz)
    apply (mul_le_mul_left hp).mp
    simpa [pow_two] using hh

theorem half_resolvent_bound (Q A : E →L[ℝ] E)
    (hpos : ∀ x, 0 ≤ ⟪x,Q x⟫) (hres : ∀ x, (2:ℝ) • A x+Q (A x)=x) :
    ‖A‖ ≤ 1/2 := by
  apply ContinuousLinearMap.opNorm_le_bound _ (by norm_num)
  intro x
  have he : 2*‖A x‖^2 ≤ ⟪A x,x⟫ := by
    conv_rhs => arg 2; rw [← hres x]
    rw [inner_add_right, inner_smul_right, real_inner_self_eq_norm_sq]
    linarith [hpos (A x)]
  have hh := he.trans (real_inner_le_norm (A x) x)
  by_cases hz : ‖A x‖=0
  · rw [hz]; positivity
  · have hp : 0 < ‖A x‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hz)
    have hm : ‖A x‖*(2*‖A x‖) ≤ ‖A x‖*‖x‖ := by nlinarith
    have hh2 := (mul_le_mul_left hp).mp hm
    linarith

/-- Lift the one-stage P1 energy inequality through the already checked
two-stage identity; this is the precise reference hypothesis used below. -/
theorem stage_implies_energy_contract (Q R D G : E →L[ℝ] E)
    (hsym : ∀ x y, ⟪x, Q y⟫ = ⟪y, Q x⟫)
    (hpos : ∀ x, 0 ≤ ⟪x, Q x⟫) (hres : ∀ x, R x+Q (R x)=x)
    (hupdate : ∀ x, (2 : ℝ) • G x + Q (G x) = x + D (R (D x)))
    (hstage : ∀ x, ‖R (D x)‖^2+⟪R (D x), Q (R (D x))⟫ ≤
      ‖x‖^2+⟪x, Q x⟫) (x : E) : sourceEnergy Q (G x) ≤ sourceEnergy Q x := by
  let B : E →ₗ[ℝ] E := (1+Q).toLinearMap
  have hb (y : E) : B y = y+Q y := rfl
  have he (y : E) : RSLA.energy B y = ‖y‖^2+⟪y, Q y⟫ := by
    simp only [RSLA.energy, hb, inner_add_right, real_inner_self_eq_norm_sq]
  have hbsym : ∀ y z, ⟪y, B z⟫ = ⟪z, B y⟫ := by
    intro y z
    simp only [hb, inner_add_right]
    rw [real_inner_comm y z, hsym]
  have hbpos : ∀ y, 0 ≤ RSLA.energy B y := by
    intro y; rw [he]; exact add_nonneg (sq_nonneg _) (hpos _)
  have hbupdate : ∀ y, G y+B (G y) = y+B (R (D (R (D y)))) := by
    intro y
    simp only [hb]
    rw [hres]
    have hh := hupdate y
    simpa only [two_smul, add_assoc] using hh
  have hh := RSLA.two_stage_contract B (fun y => R (D y)) G hbsym hbpos
    (by intro y; simpa only [he] using hstage y) hbupdate x
  simp only [he, sourceEnergy] at hh ⊢
  linarith

theorem first_step_energy (Q R D G : E → E) (d : ℝ)
    (hd : 0 ≤ d) (hpos : ∀ x, 0 ≤ ⟪x, Q x⟫)
    (hR : ∀ x, ‖R x‖ ≤ ‖x‖) (hD : ∀ x, ‖D x‖ ≤ d*‖x‖)
    (hupdate : ∀ x, (2 : ℝ) • G x + Q (G x) = x + D (R (D x))) (x : E) :
    sourceEnergy Q (G x) ≤ ((1+d^2)^2/2)*‖x‖^2 := by
  have hF : ‖x+D (R (D x))‖ ≤ (1+d^2)*‖x‖ := by
    calc
      _ ≤ ‖x‖+‖D (R (D x))‖ := norm_add_le _ _
      _ ≤ ‖x‖+d*‖R (D x)‖ := add_le_add_left (hD _) _
      _ ≤ ‖x‖+d*‖D x‖ := by gcongr; exact hR _
      _ ≤ ‖x‖+d*(d*‖x‖) := by gcongr; exact hD _
      _ = _ := by ring
  have heq : sourceEnergy Q (G x) = ⟪G x, x+D (R (D x))⟫ := by
    rw [← hupdate x]
    simp [sourceEnergy, inner_add_right, inner_smul_right, real_inner_self_eq_norm_sq]
  have hupper : sourceEnergy Q (G x) ≤ ‖G x‖*((1+d^2)*‖x‖) := by
    rw [heq]
    exact (real_inner_le_norm _ _).trans (mul_le_mul_of_nonneg_left hF (norm_nonneg _))
  have hlower : 2*‖G x‖^2 ≤ sourceEnergy Q (G x) := by
    dsimp [sourceEnergy]; linarith [hpos (G x)]
  have hhalf : 2*‖G x‖ ≤ (1+d^2)*‖x‖ := by
    by_cases hz : ‖G x‖=0
    · simp only [hz, mul_zero]; positivity
    · have hp : 0 < ‖G x‖ := lt_of_le_of_ne (norm_nonneg _) (Ne.symm hz)
      apply (mul_le_mul_left hp).mp
      nlinarith
  nlinarith [mul_le_mul_of_nonneg_right hhalf (by positivity : 0 ≤ (1+d^2)*‖x‖)]

theorem iterated_energy (Q G : E →L[ℝ] E)
    (hcontract : ∀ x, sourceEnergy Q (G x) ≤ sourceEnergy Q x) :
    ∀ (n : ℕ) (x : E), sourceEnergy Q ((G^n) x) ≤ sourceEnergy Q x := by
  intro n
  induction n with
  | zero => intro x; simp
  | succ n ih =>
    intro x
    rw [pow_succ', ContinuousLinearMap.mul_apply]
    exact (hcontract _).trans (ih x)

/-- Uniform in Q: the resulting power bound contains no opacity-dependent
comparison of the energy norm to the base norm. -/
theorem uniform_power_bound (Q R D G : E →L[ℝ] E) (d : ℝ)
    (hd : 1 ≤ d) (hpos : ∀ x, 0 ≤ ⟪x, Q x⟫)
    (hR : ∀ x, ‖R x‖ ≤ ‖x‖) (hD : ∀ x, ‖D x‖ ≤ d*‖x‖)
    (hupdate : ∀ x, (2 : ℝ) • G x + Q (G x) = x + D (R (D x)))
    (hcontract : ∀ x, sourceEnergy Q (G x) ≤ sourceEnergy Q x) (n : ℕ) :
    ‖G^n‖ ≤ (1+d^2)/2 := by
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  cases n with
  | zero =>
    simp only [pow_zero, ContinuousLinearMap.one_apply]
    have hsq : 1 ≤ d^2 := by nlinarith
    nlinarith [norm_nonneg x]
  | succ n =>
    have hh := iterated_energy Q G hcontract n (G x)
    rw [← ContinuousLinearMap.mul_apply, ← pow_succ] at hh
    have hf := first_step_energy Q R D G d (by linarith) hpos hR hD hupdate x
    have hp := hpos ((G^(n+1)) x)
    dsimp [sourceEnergy] at hh hf
    have hsq : ‖(G^(n+1)) x‖^2 ≤ (((1+d^2)/2)*‖x‖)^2 := by nlinarith
    exact (sq_le_sq₀ (norm_nonneg _) (by positivity)).mp hsq

/-- The reference power estimate follows from positivity, the actual implicit
solve/update identities, and the one-stage energy inequality. -/
theorem reference_power_bound (Q R D G : E →L[ℝ] E) (d : ℝ)
    (hd : 1 ≤ d) (hsym : ∀ x y, ⟪x, Q y⟫ = ⟪y, Q x⟫)
    (hpos : ∀ x, 0 ≤ ⟪x, Q x⟫) (hres : ∀ x, R x+Q (R x)=x)
    (hD : ∀ x, ‖D x‖ ≤ d*‖x‖)
    (hupdate : ∀ x, (2 : ℝ) • G x + Q (G x) = x + D (R (D x)))
    (hstage : ∀ x, ‖R (D x)‖^2+⟪R (D x), Q (R (D x))⟫ ≤
      ‖x‖^2+⟪x, Q x⟫) (n : ℕ) : ‖G^n‖ ≤ (1+d^2)/2 :=
  uniform_power_bound Q R D G d hd hpos (resolvent_contract Q R hpos hres) hD hupdate
    (stage_implies_energy_contract Q R D G hsym hpos hres hupdate hstage) n

theorem functional_uniform_power (Q R D G : E → E) (d : ℝ)
    (hd : 1 ≤ d) (hpos : ∀ x, 0 ≤ ⟪x,Q x⟫)
    (hR : ∀ x, ‖R x‖ ≤ ‖x‖) (hD : ∀ x, ‖D x‖ ≤ d*‖x‖)
    (hupdate : ∀ x, (2:ℝ) • G x+Q (G x)=x+D (R (D x)))
    (hcontract : ∀ x, sourceEnergy Q (G x) ≤ sourceEnergy Q x)
    (n : ℕ) (x : E) : ‖G^[n] x‖ ≤ ((1+d^2)/2)*‖x‖ := by
  have hi : ∀ (k : ℕ) (y : E), sourceEnergy Q (G^[k] y) ≤ sourceEnergy Q y := by
    intro k
    induction k with
    | zero => intro y; simp
    | succ k ih =>
      intro y
      rw [Function.iterate_succ_apply']
      exact (hcontract _).trans (ih y)
  cases n with
  | zero =>
    simp only [Function.iterate_zero_apply]
    have hd2 : 1 ≤ d^2 := by nlinarith
    nlinarith [norm_nonneg x]
  | succ k =>
    have hh := hi k (G x)
    have hf := first_step_energy Q R D G d (by linarith) hpos hR hD hupdate x
    have hp := hpos (G^[k] (G x))
    dsimp [sourceEnergy] at hh hf
    have hs : ‖G^[k] (G x)‖^2 ≤ (((1+d^2)/2)*‖x‖)^2 := by nlinarith
    simpa only [Function.iterate_succ_apply] using
      (sq_le_sq₀ (norm_nonneg _) (by positivity)).mp hs

end RSLA.M1Energy
