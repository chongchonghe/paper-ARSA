import RSLA.M1Diffusion
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-! All-phase HLL symbol perturbation: the physical-flux Jacobians vary with
M1, and the HLL speed variations enter their own diffusion term. -/
noncomputable section
namespace RSLA.M1Fourier
open scoped BigOperators
variable {A : Type*} [NormedRing A] [NormedAlgebra ℝ A]
variable {ι : Type*} [Fintype ι]

def symbol (ν σ θ : ι → ℝ) (J P : A) (a : ι → A) : A :=
  ∑ i, ν i • ((σ i*(1-Real.cos (θ i))) • P + Real.sin (θ i) • (J*a i))

lemma real_norm_smul (s : ℝ) (x : A) : ‖s • x‖=|s| * ‖x‖ := by
  simpa only [Real.norm_eq_abs] using norm_smul s x

theorem symbol_difference (ν σ σ₀ θ : ι → ℝ) (J P : A) (a a₀ : ι → A) :
    symbol ν σ θ J P a-symbol ν σ₀ θ J P a₀ =
      ∑ i, ν i • (((σ i-σ₀ i)*(1-Real.cos (θ i))) • P+
        Real.sin (θ i) • (J*(a i-a₀ i))) := by
  unfold symbol
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  simp only [mul_sub]
  module

theorem symbol_difference_bound (ν σ σ₀ θ : ι → ℝ) (J P : A) (a a₀ : ι → A)
    (f cσ cA : ℝ) (hν : ∀ i, 0 ≤ ν i) (hf : 0 ≤ f) (hcσ : 0 ≤ cσ) (hcA : 0 ≤ cA)
    (hJ : ‖J‖ ≤ 1) (hP : ‖P‖ ≤ 1)
    (hσ : ∀ i, |σ i-σ₀ i| ≤ cσ*f) (ha : ∀ i, ‖a i-a₀ i‖ ≤ cA*f) :
    ‖symbol ν σ θ J P a-symbol ν σ₀ θ J P a₀‖ ≤ (2*cσ+cA)*f*(∑ i, ν i) := by
  rw [symbol_difference]
  calc
    _ ≤ ∑ i, ‖ν i • (((σ i-σ₀ i)*(1-Real.cos (θ i))) • P+
        Real.sin (θ i) • (J*(a i-a₀ i)))‖ := norm_sum_le _ _
    _ ≤ ∑ i, ν i*((2*cσ+cA)*f) := by
      apply Finset.sum_le_sum
      intro i _
      rw [real_norm_smul, abs_of_nonneg (hν i)]
      apply mul_le_mul_of_nonneg_left _ (hν i)
      have hz : |1-Real.cos (θ i)| ≤ 2 := by
        rw [abs_of_nonneg (by linarith [Real.cos_le_one (θ i)])]
        linarith [Real.neg_one_le_cos (θ i)]
      have hs : |Real.sin (θ i)| ≤ 1 := Real.abs_sin_le_one _
      calc
        _ ≤ ‖((σ i-σ₀ i)*(1-Real.cos (θ i))) • P‖+
            ‖Real.sin (θ i) • (J*(a i-a₀ i))‖ := norm_add_le _ _
        _ = |σ i-σ₀ i| * |1-Real.cos (θ i)| * ‖P‖+
            |Real.sin (θ i)| * ‖J*(a i-a₀ i)‖ := by
          rw [real_norm_smul, real_norm_smul, abs_mul]
        _ ≤ (cσ*f)*2*1+1*(1*(cA*f)) := by
          gcongr
          · exact hσ i
          · exact (norm_mul_le _ _).trans (by gcongr; exact ha i)
        _ = (2*cσ+cA)*f := by ring
    _ = (2*cσ+cA)*f*(∑ i, ν i) := by rw [← Finset.sum_mul]; ring

end RSLA.M1Fourier
