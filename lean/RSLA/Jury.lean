import Mathlib.Analysis.Complex.Polynomial.Basic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Tactic

/-! The closed-unit-disk Jury criterion, including repeated and nonreal roots. -/
namespace RSLA

/-- Roots of the real monic characteristic quadratic lie in the closed unit disk. -/
def QuadraticStable (t d : ℝ) : Prop :=
  ∀ z : ℂ, z ^ 2 - (t : ℂ) * z + (d : ℂ) = 0 → ‖z‖ ≤ 1

lemma norm_le_one_iff (z : ℂ) : ‖z‖ ≤ 1 ↔ z.re ^ 2 + z.im ^ 2 ≤ 1 := by
  have h := Complex.normSq_eq_norm_sq z
  rw [Complex.normSq_apply] at h
  constructor <;> intro hz
  · nlinarith [norm_nonneg z]
  · nlinarith [norm_nonneg z]

lemma quadratic_parts {t d : ℝ} {z : ℂ}
    (hz : z ^ 2 - (t : ℂ) * z + (d : ℂ) = 0) :
    z.re ^ 2 - z.im ^ 2 - t * z.re + d = 0 ∧ z.im * (2 * z.re - t) = 0 := by
  have hr := congrArg Complex.re hz
  have hi := congrArg Complex.im hz
  simp [pow_two, Complex.mul_re, Complex.mul_im] at hr hi
  constructor <;> nlinarith

lemma real_root_bound {t d x : ℝ} (h₁ : 0 ≤ 1 - t + d)
    (h₂ : 0 ≤ 1 + t + d) (hd : d ≤ 1) (hx : x ^ 2 - t * x + d = 0) :
    -1 ≤ x ∧ x ≤ 1 := by
  have ht : -2 ≤ t ∧ t ≤ 2 := by constructor <;> linarith
  constructor
  · by_contra! h
    have hp : 0 < (-1 - x) * (-x + 1 + t) := mul_pos (by linarith) (by linarith)
    nlinarith
  · by_contra! h
    have hp : 0 < (x - 1) * (x + 1 - t) := mul_pos (by linarith) (by linarith)
    nlinarith

/-- Lemma A.1, expressed in trace and determinant. -/
theorem jury (t d : ℝ) : QuadraticStable t d ↔
    0 ≤ 1 - t + d ∧ 0 ≤ 1 + t + d ∧ |d| ≤ 1 := by
  constructor
  · intro hs
    obtain ⟨s, hsqrt⟩ := IsAlgClosed.exists_pow_nat_eq
      ((t : ℂ) ^ 2 - 4 * (d : ℂ)) (n := 2) (by norm_num)
    let z : ℂ := ((t : ℂ) + s) / 2
    have hz : z ^ 2 - (t : ℂ) * z + (d : ℂ) = 0 := by
      dsimp [z]
      linear_combination hsqrt / 4
    let w : ℂ := (t : ℂ) - z
    have hw : w ^ 2 - (t : ℂ) * w + (d : ℂ) = 0 := by
      dsimp [w]
      linear_combination hz
    have hzn := hs z hz
    have hwn := hs w hw
    have hzsq := (norm_le_one_iff z).mp hzn
    have hwsq := (norm_le_one_iff w).mp hwn
    have hp : z * w = (d : ℂ) := by dsimp [w]; linear_combination -hz
    have hdn : |d| ≤ 1 := by
      calc
        |d| = ‖(d : ℂ)‖ := (Complex.norm_real d).symm
        _ = ‖z‖ * ‖w‖ := by rw [← hp, norm_mul]
        _ ≤ 1 := by nlinarith [norm_nonneg z, norm_nonneg w]
    obtain ⟨hr, hi⟩ := quadratic_parts hz
    have hwr : w.re = t - z.re := by simp [w]
    have hwi : w.im = -z.im := by simp [w]
    rw [hwr, hwi] at hwsq
    rcases mul_eq_zero.mp hi with hy | ht
    · have hbounds : -1 ≤ z.re ∧ z.re ≤ 1 ∧ -1 ≤ t - z.re ∧ t - z.re ≤ 1 := by
        rw [hy] at hzsq hwsq
        constructor
        · nlinarith [sq_nonneg (z.re + 1)]
        constructor
        · nlinarith [sq_nonneg (z.re - 1)]
        constructor <;> nlinarith [sq_nonneg (t - z.re + 1), sq_nonneg (t - z.re - 1)]
      refine ⟨?_, ?_, hdn⟩
      · have := mul_nonneg (show 0 ≤ 1 - z.re by linarith)
            (show 0 ≤ 1 - (t - z.re) by linarith)
        rw [hy] at hr
        nlinarith
      · have := mul_nonneg (show 0 ≤ 1 + z.re by linarith)
            (show 0 ≤ 1 + (t - z.re) by linarith)
        rw [hy] at hr
        nlinarith
    · refine ⟨?_, ?_, hdn⟩ <;>
        nlinarith [sq_nonneg (z.re - 1), sq_nonneg (z.re + 1), sq_nonneg z.im]
  · rintro ⟨h₁, h₂, hd⟩ z hz
    apply (norm_le_one_iff z).mpr
    obtain ⟨hr, hi⟩ := quadratic_parts hz
    have hd' := (abs_le.mp hd).2
    rcases mul_eq_zero.mp hi with hy | ht
    · rw [hy] at hr ⊢
      obtain ⟨hl, hu⟩ := real_root_bound h₁ h₂ hd' (by nlinarith [hr])
      nlinarith [mul_nonneg (by linarith : 0 ≤ 1 - z.re) (by linarith : 0 ≤ 1 + z.re)]
    · nlinarith

abbrev Mat2 := Matrix (Fin 2) (Fin 2) ℝ

def trace2 (G : Mat2) : ℝ := G 0 0 + G 1 1

/-- The usual characteristic-determinant definition of complex spectral stability. -/
def MatrixStable (G : Mat2) : Prop :=
  ∀ z : ℂ, (z • (1 : Matrix (Fin 2) (Fin 2) ℂ) - G.map (algebraMap ℝ ℂ)).det = 0 → ‖z‖ ≤ 1

lemma characteristic_two (G : Mat2) (z : ℂ) :
    (z • (1 : Matrix (Fin 2) (Fin 2) ℂ) - G.map (algebraMap ℝ ℂ)).det =
      z ^ 2 - (trace2 G : ℂ) * z + (G.det : ℂ) := by
  simp [Matrix.det_fin_two, trace2, Matrix.sub_apply, Matrix.smul_apply,
    Matrix.one_apply, Matrix.map_apply]
  ring

/-- Lemma A.1 for arbitrary real 2 by 2 matrices. -/
theorem matrix_jury (G : Mat2) : MatrixStable G ↔
    0 ≤ (1 - G).det ∧ 0 ≤ (1 + G).det ∧ |G.det| ≤ 1 := by
  have hm : (1 - G).det = 1 - trace2 G + G.det := by
    simp [Matrix.det_fin_two, trace2, Matrix.sub_apply, Matrix.one_apply]; ring
  have hp : (1 + G).det = 1 + trace2 G + G.det := by
    simp [Matrix.det_fin_two, trace2, Matrix.add_apply, Matrix.one_apply]; ring
  unfold MatrixStable
  simp only [characteristic_two, hm, hp]
  exact jury (trace2 G) G.det

/-- Lemma A.2, over any field and in arbitrary finite dimension. -/
theorem determinant_split {n 𝕜 : Type*} [Fintype n] [DecidableEq n] [Field 𝕜]
    (B D : Matrix n n 𝕜) (hB : IsUnit B.det) :
    (D * B⁻¹ * D - B).det =
      (-1 : 𝕜) ^ Fintype.card n * (B - D).det * (B + D).det / B.det := by
  have hf : (D - B) * B⁻¹ * (D + B) = D * B⁻¹ * D - B := by
    noncomm_ring [Matrix.mul_nonsing_inv B hB, Matrix.nonsing_inv_mul B hB]
  rw [← hf, Matrix.det_mul, Matrix.det_mul, show D - B = -(B - D) by abel,
    Matrix.det_neg, Matrix.det_nonsing_inv, Ring.inverse_eq_inv, add_comm D B]
  ring

end RSLA
