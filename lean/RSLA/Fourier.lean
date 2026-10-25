import RSLA.EnergyBlock
import Mathlib.Algebra.Order.Chebyshev
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-! Scalar and channelwise ingredients of Appendix B. -/
noncomputable section
namespace RSLA.Fourier
open scoped BigOperators

lemma phase_product {z w : ℝ} (hz : 0 ≤ z) (hz' : z ≤ 2)
    (hw : 0 ≤ w) (hw' : w ≤ 2) : z*w ≤ z+w := by
  nlinarith [mul_nonneg hz (show 0 ≤ 2-w by linarith),
    mul_nonneg hw (show 0 ≤ 2-z by linarith)]

/-- The finite-sum estimate in B.8, with w_d = a_d/2. -/
theorem fourier_scalar_bound {ι : Type*} [Fintype ι] [DecidableEq ι]
    (w z : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (hz : ∀ i, 0 ≤ z i) (hz' : ∀ i, z i ≤ 2) :
    (∑ i, w i*z i)^2 + ∑ i, (w i)^2*z i*(2-z i) ≤
      (2*∑ i, w i) * (∑ i, w i*z i) := by
  have hpair : ∀ i j, w i*w j*z i*z j +
      (if i=j then (w i)^2*z i*(2-z i) else 0) ≤ w i*w j*(z i+z j) := by
    intro i j
    by_cases hij : i=j
    · subst j; simp only [ite_true]; nlinarith
    · simp only [hij, ite_false, add_zero]
      have h := mul_le_mul_of_nonneg_left (phase_product (hz i) (hz' i) (hz j) (hz' j))
        (mul_nonneg (hw i) (hw j))
      nlinarith
  have hh := Finset.sum_le_sum (s := Finset.univ) (fun i _ =>
    Finset.sum_le_sum (s := Finset.univ) (fun j _ => hpair i j))
  have hl : (∑ i, ∑ j, (w i*w j*z i*z j +
      if i=j then (w i)^2*z i*(2-z i) else 0)) =
      (∑ i, w i*z i)^2 + ∑ i, (w i)^2*z i*(2-z i) := by
    simp only [Finset.sum_add_distrib, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
    congr 1
    rw [pow_two, Finset.sum_mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    apply Finset.sum_congr rfl
    intro j _
    ring
  have hr : (∑ i, ∑ j, w i*w j*(z i+z j)) =
      (2*∑ i, w i) * (∑ i, w i*z i) := by
    simp_rw [mul_add, Finset.sum_add_distrib]
    have h₁ : (∑ i, ∑ j, w i*w j*z i) = (∑ i, w i*z i) * (∑ i, w i) := by
      rw [Finset.sum_mul_sum]
      congr 1; ext i
      congr 1; ext j
      ring
    have h₂ : (∑ i, ∑ j, w i*w j*z j) = (∑ i, w i) * (∑ i, w i*z i) := by
      rw [Finset.sum_mul_sum]
      congr 1; ext i
      congr 1; ext j
      ring
    rw [h₁, h₂]; ring
  rwa [hl, hr] at hh

/-- B.8 for actual Fourier phases, including arbitrary anisotropic grid weights. -/
theorem fourier_trig_bound {ι : Type*} [Fintype ι] [DecidableEq ι]
    (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i) :
    (∑ i, w i*(1-Real.cos (θ i)))^2 + ∑ i, (w i)^2*(Real.sin (θ i))^2 ≤
      (2*∑ i, w i) * (∑ i, w i*(1-Real.cos (θ i))) := by
  have h := fourier_scalar_bound w (fun i => 1-Real.cos (θ i)) hw
    (fun i => by dsimp; linarith [Real.cos_le_one (θ i)])
    (fun i => by dsimp; linarith [Real.neg_one_le_cos (θ i)])
  have hi : ∀ i, (w i)^2*(Real.sin (θ i))^2 =
      (w i)^2*(1-Real.cos (θ i))*(2-(1-Real.cos (θ i))) := by
    intro i
    have hh : (Real.sin (θ i))^2 = 1 - (Real.cos (θ i))^2 := by
      nlinarith [Real.sin_sq_add_cos_sq (θ i)]
    rw [hh]
    ring
  simp_rw [hi]
  exact h

/-- B-weighted energy of one real channel, in the equivalent scaling where the
exchange equilibrium is u=v. Complex channels are the sum of this form on real
and imaginary parts. -/
def channelEnergy (q K u v : ℝ) : ℝ := q/K*u^2 + v^2 + q*(u-v)^2

def stageGas (q K u r : ℝ) : ℝ := ((1+q)*u+K*r)/(1+q+K)
def stageRad (q K u r : ℝ) : ℝ := (q*u+(1+K)*r)/(1+q+K)

/-- Exact completion of the gas-variable square giving B.7.
Here w is the radiation transport increment, so r=v-w. -/
theorem channel_schur_identity {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v w : ℝ) :
    channelEnergy q K u v -
      channelEnergy q K (stageGas q K u (v-w)) (stageRad q K u (v-w)) =
    q*(2+q+K)/(1+q+K) * (u-v+w/(2+q+K))^2 +
      2*v*w - (2+K)/(2+q+K)*w^2 := by
  unfold channelEnergy stageGas stageRad
  field_simp
  <;> ring

/-- The eliminated coefficient is exactly 2/a_* as in B.7. -/
theorem schur_coefficient {q K : ℝ} (hq : 0 < q) (hK : 0 < K) :
    (2+K)/(2+q+K) = 2 / EnergyBlock.boundary q K := by
  unfold EnergyBlock.boundary
  field_simp
  <;> ring

/-- At the checkerboard, the remaining channel form is nonnegative exactly
up to the channel boundary. The zero-step case is included. -/
theorem checkerboard_schur_iff {q K a : ℝ} (hq : 0 < q) (hK : 0 < K) (ha : 0 ≤ a) :
    0 ≤ 2*a - (2+K)/(2+q+K)*a^2 ↔ a ≤ EnergyBlock.boundary q K := by
  have hp : 0 < 2+q+K := by positivity
  have hpK : 0 < 2+K := by positivity
  have heq : 2*a - (2+K)/(2+q+K)*a^2 = a*(2*(2+q+K)-(2+K)*a)/(2+q+K) := by
    field_simp <;> ring
  rw [heq, le_div_iff₀ hp, EnergyBlock.boundary, le_div_iff₀ hpK]
  simp only [zero_mul]
  constructor
  · intro h
    by_cases hzero : a=0
    · subst a; nlinarith
    · have := nonneg_of_mul_nonneg_right h (lt_of_le_of_ne ha (Ne.symm hzero))
      nlinarith
  · intro h
    exact mul_nonneg ha (by nlinarith)

/-- The two positive determinant factors in B.10 for an interior phase.
The explicit assumption a>0 corrects the zero-step omission in Lemma B.3. -/
theorem interior_determinants_pos {a z K KF aE aF : ℝ}
    (ha : 0 < a) (hz : 0 < z) (hz' : z < 2) (hK : 0 < K) (hKF : 0 < KF)
    (hE : a ≤ aE) (hF : a ≤ aF) :
    0 < a^2/2*K*KF*z ∧
    0 < (2+K)*(2+KF)*((aE-a*z/2)*(aF-a*z/2)+a^2/4*z*(2-z)) := by
  have hza : a*z/2 < a := by nlinarith [mul_pos ha (show 0 < 2-z by linarith)]
  have h₁ : 0 < aE-a*z/2 := by linarith
  have h₂ : 0 < aF-a*z/2 := by linarith
  have hz2 : 0 < 2-z := by linarith
  constructor <;> positivity

end RSLA.Fourier
