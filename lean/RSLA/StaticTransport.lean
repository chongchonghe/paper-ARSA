import RSLA.Fourier

/-! The complex Fourier transport symbol of B.3, with one energy variable
and one flux variable per active direction. -/
noncomputable section
namespace RSLA.StaticTransport
open scoped BigOperators
open Complex ComplexConjugate

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

def norm2 (x : Option ι → ℂ) : ℝ := normSq (x none) + ∑ i, normSq (x (some i))
def pairing (x y : Option ι → ℂ) : ℝ :=
  ((x none) * conj (y none)).re + ∑ i, ((x (some i)) * conj (y (some i))).re

def fluxDot (y : ι → ℝ) (x : Option ι → ℂ) : ℂ := ∑ i, (y i : ℂ)*x (some i)

def symbol (s : ℝ) (y : ι → ℝ) (x : Option ι → ℂ) : Option ι → ℂ
  | none => (s : ℂ)*x none + I*fluxDot y x
  | some i => (s : ℂ)*x (some i) + I*((y i : ℂ)*x none)

lemma normSq_real_imag (s : ℝ) (z t : ℂ) :
    normSq ((s : ℂ)*z + I*t) = s^2*normSq z + normSq t +
      2*s*(z.im*t.re-z.re*t.im) := by
  simp [normSq_apply, Complex.mul_re, Complex.mul_im]
  ring

lemma normSq_real_mul (s : ℝ) (z : ℂ) : normSq ((s : ℂ)*z) = s^2*normSq z := by
  rw [normSq_mul, normSq_ofReal]; ring

lemma fluxDot_re (y : ι → ℝ) (x : Option ι → ℂ) :
    (fluxDot y x).re = ∑ i, y i*(x (some i)).re := by simp [fluxDot]
lemma fluxDot_im (y : ι → ℝ) (x : Option ι → ℂ) :
    (fluxDot y x).im = ∑ i, y i*(x (some i)).im := by simp [fluxDot]

/-- Cauchy--Schwarz for a real row acting on a complex vector. -/
theorem fluxDot_bound (y : ι → ℝ) (x : Option ι → ℂ) :
    normSq (fluxDot y x) ≤ (∑ i, (y i)^2) * (∑ i, normSq (x (some i))) := by
  have hr := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ y (fun i => (x (some i)).re)
  have hi := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ y (fun i => (x (some i)).im)
  simp only [normSq_apply, fluxDot_re, fluxDot_im]
  have heq : (∑ i, ((x (some i)).re*(x (some i)).re + (x (some i)).im*(x (some i)).im)) =
      (∑ i, (x (some i)).re^2) + ∑ i, (x (some i)).im^2 := by
    simp [pow_two, Finset.sum_add_distrib]
  rw [heq]
  nlinarith

/-- Self-adjointness of the star-shaped radiation flux coupling cancels all cross terms. -/
theorem symbol_norm_identity (s : ℝ) (y : ι → ℝ) (x : Option ι → ℂ) :
    norm2 (symbol s y x) = s^2*norm2 x + normSq (fluxDot y x) +
      (∑ i, (y i)^2)*normSq (x none) := by
  unfold norm2
  simp only [symbol, normSq_real_imag, normSq_real_mul,
    Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, zero_mul,
    add_zero, sub_zero, fluxDot_re, fluxDot_im]
  have hcross : (∑ i, 2*s*((x (some i)).im*(y i*(x none).re) -
      (x (some i)).re*(y i*(x none).im))) =
      2*s*((x none).re*(∑ i, y i*(x (some i)).im) -
        (x none).im*(∑ i, y i*(x (some i)).re)) := by
    rw [← Finset.mul_sum]
    congr 1
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  simp only [Finset.sum_add_distrib]
  rw [hcross, ← Finset.mul_sum, ← Finset.sum_mul]
  ring

/-- The real quadratic form of W is Σ times the ordinary squared norm. -/
theorem symbol_pairing (s : ℝ) (y : ι → ℝ) (x : Option ι → ℂ) :
    pairing x (symbol s y x) = s*norm2 x := by
  unfold pairing norm2
  simp only [symbol, map_add, map_mul, Complex.conj_ofReal, Complex.conj_I,
    Complex.mul_re, Complex.mul_im, Complex.add_re, Complex.add_im,
    Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im,
    Complex.neg_re, Complex.neg_im, Complex.conj_re, Complex.conj_im,
    normSq_apply, fluxDot_re, fluxDot_im]
  simp only [zero_mul, mul_zero, one_mul, mul_one, neg_zero, neg_neg, add_zero, zero_add,
    sub_zero, zero_sub]
  simp only [mul_add, mul_sub, add_mul, sub_mul, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  ring_nf
  have hr : (∑ i, (x none).im*(x (some i)).re*y i) =
      (x none).im*(∑ i, y i*(x (some i)).re) := by
    rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro i _; ring
  have hi : (∑ i, (x none).re*(x (some i)).im*y i) =
      (x none).re*(∑ i, y i*(x (some i)).im) := by
    rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro i _; ring
  simp only [Finset.sum_neg_distrib, hr, hi, ← Finset.mul_sum]
  ring

/-- The spectral estimate for Y² used in B.8, without invoking eigenvalues. -/
theorem symbol_norm_bound (s : ℝ) (y : ι → ℝ) (x : Option ι → ℂ) :
    norm2 (symbol s y x) ≤ (s^2 + ∑ i, (y i)^2)*norm2 x := by
  rw [symbol_norm_identity]
  have hh := fluxDot_bound y x
  unfold norm2
  nlinarith

noncomputable def W (w θ : ι → ℝ) (x : Option ι → ℂ) : Option ι → ℂ :=
  symbol (∑ i, w i*(1-Real.cos (θ i))) (fun i => w i*Real.sin (θ i)) x

/-- B.3--B.8 combined: the actual multidimensional Fourier symbol satisfies
||Wx||² ≤ a Re<x,Wx>, with a = sum a_d. -/
theorem fourier_accretive (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (x : Option ι → ℂ) :
    norm2 (W w θ x) ≤ (2*∑ i, w i)*pairing x (W w θ x) := by
  unfold W
  rw [symbol_pairing]
  have h₁ := symbol_norm_bound (∑ i, w i*(1-Real.cos (θ i)))
    (fun i => w i*Real.sin (θ i)) x
  have h₂ := Fourier.fourier_trig_bound w θ hw
  have hn : 0 ≤ norm2 x := by
    unfold norm2
    exact add_nonneg (normSq_nonneg _) (Finset.sum_nonneg (fun i _ => normSq_nonneg _))
  simp only [mul_pow] at h₁
  have hh := mul_le_mul_of_nonneg_right h₂ hn
  nlinarith

end RSLA.StaticTransport
