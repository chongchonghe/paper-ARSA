import Mathlib.Analysis.InnerProductSpace.LinearMap
import Mathlib.Analysis.NormedSpace.OperatorNorm.Bilinear
import Mathlib.Analysis.NormedSpace.OperatorNorm.NormedSpace
import Mathlib.Tactic

/-! Explicit Levermore coefficients for the frozen full-M1 Jacobian.
The parameter b is a closure coordinate, not a material velocity. -/
noncomputable section
namespace RSLA.M1Coefficients

def fluxFactor (b : ℝ) : ℝ := 4*b/(3+b^2)
def chi (b : ℝ) : ℝ := (1+3*b^2)/(3+b^2)
def chiPrime (b : ℝ) : ℝ := 4*b/(3-b^2)

theorem parameter_bounds (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1) :
    b ≤ fluxFactor b ∧ fluxFactor b ≤ 1 := by
  have hs : b^2 ≤ 1 := by nlinarith
  have hp : 0 < 3+b^2 := by positivity
  constructor
  · rw [fluxFactor, le_div_iff₀ hp]
    nlinarith [mul_le_mul_of_nonneg_left hs hb]
  · rw [fluxFactor, div_le_iff₀ hp]
    nlinarith [sq_nonneg (b-1)]

theorem chi_difference (b : ℝ) : chi b-1/3 = 8*b^2/(3*(3+b^2)) := by
  have hd : 3+b^2 ≠ 0 := ne_of_gt (by positivity)
  unfold chi
  field_simp
  <;> ring

theorem chi_bounds (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1) :
    1/3 ≤ chi b ∧ chi b ≤ 1 ∧ chi b-1/3 ≤ (fluxFactor b)^2 := by
  have hs : b^2 ≤ 1 := by nlinarith
  have hp : 0 < 3+b^2 := by positivity
  have hf := parameter_bounds b hb hb1
  have hnon : 0 ≤ chi b-1/3 := by rw [chi_difference]; positivity
  have hsmall : chi b-1/3 ≤ b^2 := by
    rw [chi_difference, div_le_iff₀ (by positivity : 0 < 3*(3+b^2))]
    nlinarith [sq_nonneg (b^2)]
  refine ⟨by linarith, ?_, ?_⟩
  · rw [chi, div_le_iff₀ hp]; nlinarith
  · nlinarith [sq_nonneg (fluxFactor b-b)]

theorem chiPrime_bounds (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1) :
    0 ≤ chiPrime b ∧ chiPrime b ≤ 2*fluxFactor b := by
  have hs : b^2 ≤ 1 := by nlinarith
  have hp : 0 < 3-b^2 := by linarith
  have hf := parameter_bounds b hb hb1
  constructor
  · exact div_nonneg (by positivity) hp.le
  · have hh : chiPrime b ≤ 2*b := by
      rw [chiPrime, div_le_iff₀ hp]
      nlinarith [mul_le_mul_of_nonneg_left hs hb]
    linarith

/-- The transverse derivative coefficient simplifies exactly to b. -/
theorem angular_coefficient (b : ℝ) (hb : 0 < b) :
    (3*chi b-1)/(2*fluxFactor b)=b := by
  have hden : 3+b^2 ≠ 0 := ne_of_gt (by positivity)
  unfold chi fluxFactor
  field_simp
  <;> ring

section Tensor
open scoped RealInnerProductSpace
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]

def beamProjector (n : E) : E →L[ℝ] E := (innerSL ℝ n).smulRight n
def tensor (b : ℝ) (n : E) : E →L[ℝ] E :=
  ((1-chi b)/2) • (1 : E →L[ℝ] E) + ((3*chi b-1)/2) • beamProjector n

lemma operator_norm_smul (a : ℝ) (T : E →L[ℝ] E) : ‖a • T‖ = |a| * ‖T‖ := by
  simpa only [Real.norm_eq_abs] using norm_smul a T

theorem projector_norm (n : E) (hn : ‖n‖=1) : ‖beamProjector n‖=1 := by
  simp [beamProjector, ContinuousLinearMap.norm_smulRight_apply, innerSL_apply_norm, hn]

theorem tensor_norm_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n : E) (hn : ‖n‖=1) : ‖tensor b n‖ ≤ 1 := by
  have hχ := chi_bounds b hb hb1
  have ha : 0 ≤ (1-chi b)/2 := by linarith
  have hu : 0 ≤ (3*chi b-1)/2 := by linarith
  calc
    ‖tensor b n‖ ≤ ‖((1-chi b)/2) • (1 : E →L[ℝ] E)‖ +
        ‖((3*chi b-1)/2) • beamProjector n‖ := norm_add_le _ _
    _ = chi b := by
      simp only [operator_norm_smul, abs_of_nonneg ha, abs_of_nonneg hu,
        norm_one, projector_norm n hn]
      ring
    _ ≤ 1 := hχ.2.1

theorem tensor_difference (b : ℝ) (n : E) :
    tensor b n - (1/3 : ℝ) • (1 : E →L[ℝ] E) =
      (chi b-1/3) • ((-1/2 : ℝ) • (1 : E →L[ℝ] E)+(3/2 : ℝ) • beamProjector n) := by
  unfold tensor
  module

/-- A uniform quadratic isotropic-limit bound, independent of beam direction. -/
theorem tensor_difference_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n : E) (hn : ‖n‖=1) :
    ‖tensor b n - (1/3 : ℝ) • (1 : E →L[ℝ] E)‖ ≤ 2*(fluxFactor b)^2 := by
  have hχ := chi_bounds b hb hb1
  have hδ : 0 ≤ chi b-1/3 := by linarith
  have hp : ‖(-1/2 : ℝ) • (1 : E →L[ℝ] E)+(3/2 : ℝ) • beamProjector n‖ ≤ 2 := by
    calc
      _ ≤ ‖(-1/2 : ℝ) • (1 : E →L[ℝ] E)‖+‖(3/2 : ℝ) • beamProjector n‖ := norm_add_le _ _
      _ = 2 := by
        rw [operator_norm_smul, operator_norm_smul, norm_one, projector_norm n hn]
        norm_num
        rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/2),
          abs_of_nonneg (by norm_num : (0:ℝ) ≤ 3/2)]
        norm_num
  rw [tensor_difference, operator_norm_smul, abs_of_nonneg hδ]
  calc
    _ ≤ (chi b-1/3)*2 := mul_le_mul_of_nonneg_left hp hδ
    _ ≤ 2*(fluxFactor b)^2 := by nlinarith [hχ.2.2]

theorem tensor_positive (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1) (n x : E) :
    0 ≤ ⟪x, tensor b n x⟫ := by
  have hχ := chi_bounds b hb hb1
  have ha : 0 ≤ (1-chi b)/2 := by linarith
  have hu : 0 ≤ (3*chi b-1)/2 := by linarith
  have he : ⟪x, tensor b n x⟫ =
      ((1-chi b)/2)*‖x‖^2 + ((3*chi b-1)/2)*⟪n,x⟫^2 := by
    simp [tensor, beamProjector, inner_add_right, inner_smul_right,
      real_inner_self_eq_norm_sq, real_inner_comm x n, pow_two, mul_assoc]
  rw [he]
  positivity

def angularVector (n d : E) : E := (-1/2 : ℝ) • d+(3/2*⟪n,d⟫) • n

theorem angularVector_bound (n d : E) (hn : ‖n‖=1) (hd : ‖d‖=1) :
    ‖angularVector n d‖ ≤ 2 := by
  have hi : |⟪n,d⟫| ≤ 1 := by simpa [hn, hd] using abs_real_inner_le_norm n d
  calc
    _ ≤ ‖(-1/2 : ℝ) • d‖+‖(3/2*⟪n,d⟫) • n‖ := norm_add_le _ _
    _ = 1/2+3/2*|⟪n,d⟫| := by
      rw [norm_smul, norm_smul, hn, hd, Real.norm_eq_abs, Real.norm_eq_abs, abs_mul]
      norm_num
      rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/2),
        abs_of_nonneg (by norm_num : (0:ℝ) ≤ 3/2)]
    _ ≤ 2 := by linarith

def fluxMomentJacobian (b : ℝ) (n d : E) : E →L[ℝ] E :=
  chiPrime b • (innerSL ℝ n).smulRight (angularVector n d) +
  b • ((innerSL ℝ d).comp (1-beamProjector n)).smulRight n +
  (b*⟪n,d⟫) • (1-beamProjector n)

/-- The full M1 flux derivative with respect to J is O(f), uniformly in
beam direction and grid-normal direction. The closure derivatives are included. -/
theorem fluxMomentJacobian_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n d : E) (hn : ‖n‖=1) (hd : ‖d‖=1) :
    ‖fluxMomentJacobian b n d‖ ≤ 8*fluxFactor b := by
  have hc := chiPrime_bounds b hb hb1
  have hf := parameter_bounds b hb hb1
  have hf0 : 0 ≤ fluxFactor b := hb.trans hf.1
  have hbf := hf.1
  have hcp := hc.2
  have hi : |⟪n,d⟫| ≤ 1 := by simpa [hn, hd] using abs_real_inner_le_norm n d
  have ht : ‖(1 : E →L[ℝ] E)-beamProjector n‖ ≤ 2 := by
    exact (norm_sub_le _ _).trans (by norm_num [projector_norm n hn])
  have hr : ‖((innerSL ℝ d).comp (1-beamProjector n)).smulRight n‖ ≤ 2 := by
    rw [ContinuousLinearMap.norm_smulRight_apply, hn, mul_one]
    calc
      _ ≤ ‖innerSL ℝ d‖ * ‖(1 : E →L[ℝ] E)-beamProjector n‖ :=
        ContinuousLinearMap.opNorm_comp_le _ _
      _ ≤ 2 := by simpa [innerSL_apply_norm, hd] using ht
  have hrn : ‖(innerSL ℝ n).smulRight (angularVector n d)‖ ≤ 2 := by
    simpa [ContinuousLinearMap.norm_smulRight_apply, innerSL_apply_norm, hn]
      using angularVector_bound n d hn hd
  calc
    ‖fluxMomentJacobian b n d‖ ≤
      ‖chiPrime b • (innerSL ℝ n).smulRight (angularVector n d)‖ +
      ‖b • ((innerSL ℝ d).comp (1-beamProjector n)).smulRight n‖ +
      ‖(b*⟪n,d⟫) • (1-beamProjector n)‖ := by
        exact (norm_add_le _ _).trans (add_le_add_right (norm_add_le _ _) _)
    _ = chiPrime b*‖(innerSL ℝ n).smulRight (angularVector n d)‖ +
      b*‖((innerSL ℝ d).comp (1-beamProjector n)).smulRight n‖ +
      b*|⟪n,d⟫| * ‖(1 : E →L[ℝ] E)-beamProjector n‖ := by
        rw [operator_norm_smul, operator_norm_smul, operator_norm_smul,
          abs_of_nonneg hc.1, abs_of_nonneg hb, abs_mul, abs_of_nonneg hb]
    _ ≤ (2*fluxFactor b)*2+fluxFactor b*2+fluxFactor b*1*2 := by gcongr
    _ = 8*fluxFactor b := by ring

def fluxEnergyColumn (b : ℝ) (n d : E) : E :=
  tensor b n d - (fluxFactor b*chiPrime b) • angularVector n d

theorem fluxEnergyColumn_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n d : E) (hn : ‖n‖=1) (hd : ‖d‖=1) :
    ‖fluxEnergyColumn b n d-(1/3 : ℝ) • d‖ ≤ 6*(fluxFactor b)^2 := by
  have hc := chiPrime_bounds b hb hb1
  have hf : 0 ≤ fluxFactor b := hb.trans (parameter_bounds b hb hb1).1
  have hT : ‖tensor b n d-(1/3 : ℝ) • d‖ ≤ 2*(fluxFactor b)^2 := by
    have hh := (tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)).le_opNorm d
    simpa [hd] using hh.trans
      (mul_le_mul_of_nonneg_right (tensor_difference_bound b hb hb1 n hn) (norm_nonneg d))
  have hid : fluxEnergyColumn b n d-(1/3 : ℝ) • d =
      (tensor b n d-(1/3 : ℝ) • d)-(fluxFactor b*chiPrime b) • angularVector n d := by
    unfold fluxEnergyColumn
    abel
  rw [hid]
  calc
    _ ≤ ‖tensor b n d-(1/3 : ℝ) • d‖+‖(fluxFactor b*chiPrime b) • angularVector n d‖ := norm_sub_le _ _
    _ = ‖tensor b n d-(1/3 : ℝ) • d‖+(fluxFactor b*chiPrime b)*‖angularVector n d‖ := by
      rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (mul_nonneg hf hc.1)]
    _ ≤ 2*(fluxFactor b)^2+(fluxFactor b*(2*fluxFactor b))*2 := by
      gcongr
      · exact hc.2
      · exact angularVector_bound n d hn hd
    _ = _ := by ring

/-- Full directional radiation Jacobian on (E,J), including closure variation. -/
def fluxJacobian (b : ℝ) (n d : E) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  ((innerSL ℝ d).comp (ContinuousLinearMap.snd ℝ ℝ E)).prod
    ((ContinuousLinearMap.fst ℝ ℝ E).smulRight (fluxEnergyColumn b n d) +
      (fluxMomentJacobian b n d).comp (ContinuousLinearMap.snd ℝ ℝ E))

def isotropicFluxJacobian (d : E) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  ((innerSL ℝ d).comp (ContinuousLinearMap.snd ℝ ℝ E)).prod
    ((ContinuousLinearMap.fst ℝ ℝ E).smulRight ((1/3 : ℝ) • d))

/-- A concrete uniform Lipschitz constant for the entire full-M1 physical
flux Jacobian, in the block maximum norm. Fixed norm changes preserve uniformity. -/
theorem fluxJacobian_difference_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n d : E) (hn : ‖n‖=1) (hd : ‖d‖=1) :
    ‖fluxJacobian b n d-isotropicFluxJacobian d‖ ≤ 14*fluxFactor b := by
  have hfb := parameter_bounds b hb hb1
  have hf : 0 ≤ fluxFactor b := hb.trans hfb.1
  have hf2 : (fluxFactor b)^2 ≤ fluxFactor b := by nlinarith [hfb.2]
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  have hid : (fluxJacobian b n d-isotropicFluxJacobian d) x =
      (0, x.1 • (fluxEnergyColumn b n d-(1/3 : ℝ) • d)+fluxMomentJacobian b n d x.2) := by
    apply Prod.ext
    · change ⟪d,x.2⟫-⟪d,x.2⟫=0
      ring
    · change x.1 • fluxEnergyColumn b n d+fluxMomentJacobian b n d x.2-
        x.1 • ((1/3 : ℝ) • d) = _
      module
  have hcol : ‖fluxEnergyColumn b n d-(1/3 : ℝ) • d‖ ≤ 6*fluxFactor b :=
    (fluxEnergyColumn_bound b hb hb1 n d hn hd).trans (by nlinarith)
  have hmom := fluxMomentJacobian_bound b hb hb1 n d hn hd
  have he := norm_fst_le x
  have hj := norm_snd_le x
  rw [hid]
  simp only [Prod.norm_def, norm_zero, max_eq_right (norm_nonneg _)]
  calc
    _ ≤ ‖x.1 • (fluxEnergyColumn b n d-(1/3 : ℝ) • d)‖+
      ‖fluxMomentJacobian b n d x.2‖ := norm_add_le _ _
    _ ≤ ‖x.1‖*‖fluxEnergyColumn b n d-(1/3 : ℝ) • d‖+
      ‖fluxMomentJacobian b n d‖*‖x.2‖ := by
        rw [norm_smul]
        exact add_le_add_left ((fluxMomentJacobian b n d).le_opNorm x.2) _
    _ ≤ ‖x‖*(6*fluxFactor b)+(8*fluxFactor b)*‖x‖ := by gcongr
    _ = 14*fluxFactor b*‖x‖ := by ring

def coupling (b η : ℝ) (n : E) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  (ContinuousLinearMap.fst ℝ ℝ E +
    fluxFactor b • (innerSL ℝ n).comp (ContinuousLinearMap.snd ℝ ℝ E)).prod
    ((η • (1 : E →L[ℝ] E)+tensor b n).comp (ContinuousLinearMap.snd ℝ ℝ E))

def isotropicCoupling (η : ℝ) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  (ContinuousLinearMap.fst ℝ ℝ E).prod
    (((η+1/3) • (1 : E →L[ℝ] E)).comp (ContinuousLinearMap.snd ℝ ℝ E))

theorem coupling_difference_bound (b η : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n : E) (hn : ‖n‖=1) : ‖coupling b η n-isotropicCoupling (E := E) η‖ ≤ 2*fluxFactor b := by
  have hfb := parameter_bounds b hb hb1
  have hf : 0 ≤ fluxFactor b := hb.trans hfb.1
  have hf2 : (fluxFactor b)^2 ≤ fluxFactor b := by nlinarith [hfb.2]
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  have hid : (coupling b η n-isotropicCoupling (E := E) η) x =
      (fluxFactor b*⟪n,x.2⟫, (tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)) x.2) := by
    apply Prod.ext
    · change x.1+fluxFactor b*⟪n,x.2⟫-x.1=_
      ring
    · change η • x.2+tensor b n x.2-(η+1/3) • x.2 = tensor b n x.2-(1/3 : ℝ) • x.2
      module
  rw [hid, norm_prod_le_iff]
  have hj := norm_snd_le x
  have hi : |⟪n,x.2⟫| ≤ ‖x‖ := (by simpa [hn] using abs_real_inner_le_norm n x.2 :
      |⟪n,x.2⟫| ≤ ‖x.2‖).trans hj
  constructor
  · simp only [Real.norm_eq_abs, abs_mul, abs_of_nonneg hf]
    nlinarith [mul_le_mul_of_nonneg_left hi hf, norm_nonneg x]
  · calc
      _ ≤ ‖tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)‖*‖x.2‖ :=
        ContinuousLinearMap.le_opNorm _ _
      _ ≤ (2*(fluxFactor b)^2)*‖x‖ := by
        gcongr
        exact tensor_difference_bound b hb hb1 n hn
      _ ≤ (2*fluxFactor b)*‖x‖ := by gcongr

end Tensor

/-- The actual HLL speed floor preserves the O(f) isotropic-limit bound.
The hypothesis on t is supplied by the tensor's quadratic difference bound. -/
theorem floored_speed_bound (t f : ℝ) (ht : 0 ≤ t) (hf : 0 ≤ f)
    (hdiff : |t-1/3| ≤ 2*f^2) :
    |max (1/10 : ℝ) (Real.sqrt t) - Real.sqrt (1/3 : ℝ)| ≤ 2*f := by
  have ht2 := Real.sq_sqrt ht
  have hs2 := Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 1/3)
  have htn := Real.sqrt_nonneg t
  have hsn := Real.sqrt_nonneg (1/3 : ℝ)
  have hd := abs_le.mp hdiff
  have hs : |Real.sqrt t-Real.sqrt (1/3 : ℝ)| ≤ 2*f := by
    by_cases horder : Real.sqrt (1/3 : ℝ) ≤ Real.sqrt t
    · rw [abs_of_nonneg (sub_nonneg.mpr horder)]
      have hp := mul_nonneg hsn (sub_nonneg.mpr horder)
      have hsq : (Real.sqrt t-Real.sqrt (1/3 : ℝ))^2 ≤ 2*f^2 := by nlinarith
      nlinarith
    · have horder' := le_of_not_ge horder
      rw [abs_of_nonpos (sub_nonpos.mpr horder')]
      have hp := mul_nonneg htn (sub_nonneg.mpr horder')
      have hsq : (Real.sqrt (1/3 : ℝ)-Real.sqrt t)^2 ≤ 2*f^2 := by nlinarith
      nlinarith
  have hfloor : (1/10 : ℝ) ≤ Real.sqrt (1/3 : ℝ) := by nlinarith
  have hh := abs_max_sub_max_le_max (1/10 : ℝ) (Real.sqrt t)
    (1/10 : ℝ) (Real.sqrt (1/3 : ℝ))
  simpa only [max_eq_right hfloor] using
    hh.trans (max_le (by simp only [sub_self, abs_zero]; positivity) hs)

section Speed
open scoped RealInnerProductSpace
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]

def speed (b : ℝ) (n d : E) : ℝ := max (1/10 : ℝ) (Real.sqrt ⟪d,tensor b n d⟫)

theorem speed_difference_bound (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n d : E) (hn : ‖n‖=1) (hd : ‖d‖=1) :
    |speed b n d-Real.sqrt (1/3 : ℝ)| ≤ 2*fluxFactor b := by
  have hf := parameter_bounds b hb hb1
  apply floored_speed_bound _ _ (tensor_positive b hb hb1 n d) (hb.trans hf.1)
  have heq : ⟪d,tensor b n d⟫-1/3 =
      ⟪d,(tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)) d⟫ := by
    simp [inner_sub_right, inner_smul_right, real_inner_self_eq_norm_sq, hd]
  rw [heq]
  calc
    _ ≤ ‖d‖*‖(tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)) d‖ := abs_real_inner_le_norm _ _
    _ ≤ ‖d‖*(‖tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)‖*‖d‖) := by
      gcongr
      exact ContinuousLinearMap.le_opNorm _ _
    _ ≤ 2*(fluxFactor b)^2 := by
      simpa only [hd, one_mul, mul_one] using tensor_difference_bound b hb hb1 n hn

end Speed
end RSLA.M1Coefficients
