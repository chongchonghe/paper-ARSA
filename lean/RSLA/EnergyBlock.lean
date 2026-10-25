import RSLA.Jury

/-! Appendix A. `K = λ q`; positivity of q and K is equivalent to the
paper's positivity of q and λ. Matrices are exactly (59)--(61). -/
namespace RSLA.EnergyBlock
open Matrix

noncomputable def B (q K : ℝ) : Mat2 := !![1 + K, -K / 4; -4 * q, 1 + q]
noncomputable def D (a : ℝ) : Mat2 := !![1, 0; 0, 1 - a]
noncomputable def G (q K a : ℝ) : Mat2 :=
  (1 + B q K)⁻¹ * (1 + D a * (B q K)⁻¹ * D a)
noncomputable def boundary (q K : ℝ) : ℝ := 2 * (2 + q + K) / (2 + K)
noncomputable def detValue (q K a : ℝ) : ℝ :=
  ((2 + K) * (1 - a) ^ 2 + 2 + 2 * q + K) / (2 * (1 + q + K) * (2 + q + K))
noncomputable def minusValue (q K a : ℝ) : ℝ :=
  K * a * (2 * (2 + q + K) - (2 + K) * a) / (2 * (1 + q + K) * (2 + q + K))
noncomputable def plusNumerator (q K a : ℝ) : ℝ :=
  (2 + K)^2 * a^2 - 2 * (K^2 + K*q + 4*K + 4) * a + 4 * (2+q+K)^2

lemma det_B (q K : ℝ) : (B q K).det = 1 + q + K := by
  simp [B, det_fin_two]; ring
lemma det_one_add_B (q K : ℝ) : (1 + B q K).det = 2 * (2 + q + K) := by
  simp [B, det_fin_two, add_apply, one_apply]; ring

lemma G_entries (q K a : ℝ) : G q K a =
    (2 * (2 + q + K))⁻¹ • !![2 + q, K / 4; 4 * q, 2 + K] *
      (1 + D a * ((1 + q + K)⁻¹ • !![1 + q, K / 4; 4 * q, 1 + K]) * D a) := by
  unfold G
  rw [Matrix.inv_def (1 + B q K), Matrix.inv_def (B q K),
    det_B, det_one_add_B]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [B, D, adjugate_fin_two, Ring.inverse_eq_inv, Matrix.add_apply, Matrix.one_apply,
      Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two] <;> ring <;> simp

/-- Equation A.11, derived from the actual amplification matrix. -/
theorem det_G {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    (G q K a).det = detValue q K a := by
  rw [G_entries]
  simp [D, detValue, det_fin_two, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
    Matrix.add_apply, Matrix.smul_apply, Matrix.one_apply]
  field_simp
  <;> ring

/-- Equation A.8 in the variables q and K. -/
theorem det_one_sub_G {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    (1 - G q K a).det = minusValue q K a := by
  have h₁ : 1 + q + K ≠ 0 := ne_of_gt (by positivity)
  have h₂ : 2 + q + K ≠ 0 := ne_of_gt (by positivity)
  rw [G_entries]
  simp [D, minusValue, det_fin_two, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
    Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply]
  field_simp
  <;> ring

/-- Equation A.9, including its normalization by the positive denominator. -/
theorem det_one_add_G {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    (1 + G q K a).det = plusNumerator q K a / (2 * (1+q+K) * (2+q+K)) := by
  have h₁ : 1 + q + K ≠ 0 := ne_of_gt (by positivity)
  have h₂ : 2 + q + K ≠ 0 := ne_of_gt (by positivity)
  rw [G_entries]
  simp [D, plusNumerator, det_fin_two, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
    Matrix.add_apply, Matrix.smul_apply, Matrix.one_apply]
  field_simp
  <;> ring

lemma boundary_pos {q K : ℝ} (hq : 0 < q) (hK : 0 < K) : 0 < boundary q K := by
  unfold boundary; positivity
lemma two_le_boundary {q K : ℝ} (hq : 0 < q) (hK : 0 < K) : 2 ≤ boundary q K := by
  rw [boundary, le_div_iff₀ (by positivity : 0 < 2 + K)]
  linarith

/-- Lemma A.3. -/
theorem condition_i {q K a : ℝ} (hq : 0 < q) (hK : 0 < K) (ha : 0 ≤ a) :
    0 ≤ (1 - G q K a).det ↔ a ≤ boundary q K := by
  rw [det_one_sub_G hq hK, minusValue,
    le_div_iff₀ (by positivity : 0 < 2 * (1+q+K) * (2+q+K))]
  simp only [zero_mul]
  rw [boundary, le_div_iff₀ (by positivity : 0 < 2 + K)]
  constructor
  · intro h
    by_cases haz : a = 0
    · subst a; nlinarith
    · have hKa : 0 < K*a := mul_pos hK (lt_of_le_of_ne ha (Ne.symm haz))
      have := nonneg_of_mul_nonneg_right h hKa
      nlinarith
  · intro h
    exact mul_nonneg (mul_nonneg hK.le ha) (by nlinarith)

lemma plusNumerator_pos {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    0 < plusNumerator q K a := by
  have hdisc : 0 < (K^2+K*q+4*K+4*q+4) * (3*K^2+3*K*q+12*K+4*q+12) := by positivity
  have hid : (2+K)^2 * plusNumerator q K a =
      ((2+K)^2*a - (K^2+K*q+4*K+4))^2 +
      (K^2+K*q+4*K+4*q+4) * (3*K^2+3*K*q+12*K+4*q+12) := by
    unfold plusNumerator; ring
  have hp : 0 < (2+K)^2 * plusNumerator q K a := by
    rw [hid]; positivity
  exact pos_of_mul_pos_right hp (sq_nonneg _)

/-- Lemma A.4, for every real transport parameter. -/
theorem condition_ii {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    0 < (1 + G q K a).det := by
  rw [det_one_add_G hq hK]
  exact div_pos (plusNumerator_pos hq hK a) (by positivity)

lemma detValue_pos {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    0 < detValue q K a := by unfold detValue; positivity

lemma detValue_boundary {q K : ℝ} (hq : 0 < q) (hK : 0 < K) :
    detValue q K (boundary q K) = (2+2*q+K) / ((2+K)*(1+q+K)) := by
  unfold detValue boundary
  field_simp
  <;> ring

/-- Lemma A.5. The determinant is strictly between zero and one, even at both endpoints. -/
theorem condition_iii {q K a : ℝ} (hq : 0 < q) (hK : 0 < K)
    (ha : 0 ≤ a) (hab : a ≤ boundary q K) :
    0 < (G q K a).det ∧ (G q K a).det < 1 := by
  rw [det_G hq hK]
  refine ⟨detValue_pos hq hK a, ?_⟩
  have hb := two_le_boundary hq hK
  have hs : (1-a)^2 ≤ (1-boundary q K)^2 := by
    have := mul_nonneg (show 0 ≤ boundary q K - a by linarith)
      (show 0 ≤ boundary q K + a - 2 by linarith)
    nlinarith
  have hle : detValue q K a ≤ detValue q K (boundary q K) := by
    unfold detValue
    apply div_le_div_of_nonneg_right _ (by positivity)
    nlinarith
  apply lt_of_le_of_lt hle
  rw [detValue_boundary hq hK, div_lt_one (by positivity : 0 < (2+K)*(1+q+K))]
  have : 0 < K*(2+q+K) := by positivity
  nlinarith

/-- Proposition A.1: the exact complex spectral stability boundary of (61). -/
theorem stability_iff {q K a : ℝ} (hq : 0 < q) (hK : 0 < K) (ha : 0 ≤ a) :
    MatrixStable (G q K a) ↔ a ≤ boundary q K := by
  rw [matrix_jury]
  constructor
  · intro h; exact (condition_i hq hK ha).mp h.1
  · intro hab
    obtain ⟨hl, hu⟩ := condition_iii hq hK ha hab
    exact ⟨(condition_i hq hK ha).mpr hab, (condition_ii hq hK a).le,
      abs_le.mpr ⟨by linarith, hu.le⟩⟩

/-- Proposition A.1 using the parameters printed in the paper. -/
theorem proposition_A1 {q lam a : ℝ} (hq : 0 < q) (hlam : 0 < lam) (ha : 0 ≤ a) :
    MatrixStable (G q (lam*q) a) ↔ a ≤ 2*(2+q+lam*q)/(2+lam*q) :=
  stability_iff hq (mul_pos hlam hq) ha

/-- Equation 62: solving the boundary for the local Courant number. -/
theorem courant_iff {ν σ τ K : ℝ} (hK : 0 ≤ K) :
    2*ν*σ ≤ boundary (ν*τ) K ↔ ν*(σ - τ/(2+K)) ≤ 1 := by
  have hp : 0 < 2 + K := by positivity
  have heq : ν*(σ - τ/(2+K)) = (ν*(σ*(2+K)-τ))/(2+K) := by
    field_simp <;> ring
  rw [boundary, le_div_iff₀ hp, heq, div_le_iff₀ hp]
  constructor <;> intro h <;> nlinarith

end RSLA.EnergyBlock
