import RSLA.Jury
import Mathlib.Algebra.Order.Chebyshev

/-! Proposition 1: the sharp uniform CFL bound from the unit trace of the closure. -/
namespace RSLA
open scoped BigOperators

/-- Cauchy--Schwarz applied to directional signal speeds. -/
theorem signal_sum_le {D : ℕ} (σ : Fin D → ℝ)
    (htrace : ∑ i, σ i ^ 2 ≤ 1) : ∑ i, σ i ≤ Real.sqrt D := by
  have hcs := sq_sum_le_card_mul_sum_sq (s := Finset.univ) (f := σ)
  simp only [Finset.card_univ, Fintype.card_fin] at hcs
  have hd : (0 : ℝ) ≤ D := Nat.cast_nonneg _
  have hs := Real.sq_sqrt hd
  have hmul := mul_le_mul_of_nonneg_left htrace hd
  nlinarith [Real.sqrt_nonneg (D : ℝ)]

/-- Proposition 1, sufficient direction. The two damping terms can be chosen independently. -/
theorem cfl_sufficient {D : ℕ} (hD : 0 < D) {ν : ℝ} (hν : 0 ≤ ν)
    (hcfl : ν ≤ 1 / Real.sqrt D) (σ : Fin D → ℝ)
    (htrace : ∑ i, σ i ^ 2 ≤ 1) {damping : ℝ} (hdamp : 0 ≤ damping) :
    ν * (∑ i, σ i - damping) ≤ 1 := by
  have hroot : 0 < Real.sqrt (D : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hD)
  have hn := (le_div_iff₀ hroot).mp hcfl
  have hs := signal_sum_le σ htrace
  have hm := mul_le_mul_of_nonneg_left (show ∑ i, σ i - damping ≤ Real.sqrt D by linarith) hν
  linarith

/-- Equality is attained by equal directional speeds. These are the diagonals of
an M1 free-streaming tensor n ⊗ n with n along the active grid diagonal. -/
theorem equal_speeds {D : ℕ} (hD : 0 < D) :
    (∑ _ : Fin D, (1 / Real.sqrt (D : ℝ)) ^ 2) = 1 ∧
    (∑ _ : Fin D, (1 / Real.sqrt (D : ℝ))) = Real.sqrt D := by
  have hd : (0 : ℝ) < D := by exact_mod_cast hD
  have hr : Real.sqrt (D : ℝ) ≠ 0 := ne_of_gt (Real.sqrt_pos.2 hd)
  have hs := Real.sq_sqrt hd.le
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  constructor <;> field_simp <;> nlinarith

/-- Sharpness in the class of nonnegative, unit-trace directional closures.
Taking damping zero is the optically thin limit used in Proposition 1. -/
theorem cfl_iff {D : ℕ} (hD : 0 < D) {ν : ℝ} (hν : 0 ≤ ν) :
    (∀ σ : Fin D → ℝ, (∀ i, 0 ≤ σ i) → (∑ i, σ i ^ 2 ≤ 1) →
      ∀ damping : ℝ, 0 ≤ damping → ν * (∑ i, σ i - damping) ≤ 1) ↔
    ν ≤ 1 / Real.sqrt D := by
  constructor
  · intro h
    have he := equal_speeds hD
    have hh := h (fun _ => 1 / Real.sqrt (D : ℝ)) (by intro i; positivity) he.1.le 0 le_rfl
    simp only [sub_zero, he.2] at hh
    exact (le_div_iff₀ (Real.sqrt_pos.2 (by exact_mod_cast hD))).mpr hh
  · intro h σ _ ht damping hd
    exact cfl_sufficient hD hν h σ ht hd

/-- The diagonal of the M1 tensor at chi = 1 is exactly n_i squared. -/
theorem m1_streaming_diagonal (n : ℝ) :
    (1 - (1 : ℝ))/2 + (3*(1 : ℝ)-1)/2 * n^2 = n^2 := by ring

/-- A normalized streaming M1 direction attaining the sharp diagonal bound. -/
theorem m1_sharp_witness {D : ℕ} (hD : 0 < D) :
    ∃ n : Fin D → ℝ, (∑ i, n i^2)=1 ∧
      ∀ i, Real.sqrt ((1-(1:ℝ))/2 + (3*(1:ℝ)-1)/2*(n i)^2) = 1/Real.sqrt D := by
  refine ⟨fun _ => 1/Real.sqrt D, (equal_speeds hD).1, ?_⟩
  intro i
  rw [m1_streaming_diagonal, Real.sqrt_sq (by positivity)]

end RSLA
