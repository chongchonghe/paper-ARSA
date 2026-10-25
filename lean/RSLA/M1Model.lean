import RSLA.M1Theorem
import RSLA.M1Coordinates

/-! The explicit full-M1 source and HLL symbol in physical coordinates,
with both real and imaginary Fourier amplitudes retained. -/
noncomputable section
namespace RSLA.M1Model
open scoped BigOperators RealInnerProductSpace
open M1Coefficients M1Coordinates
variable {E F : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]
  [NormedAddCommGroup F] [NormedSpace ℝ F]
variable {ι : Type*} [Fintype ι]

abbrev State (E : Type*) := (ℝ × E) × (ℝ × E)
abbrev FourierState (E : Type*) := State E × State E

def source (b η μE μM : ℝ) (n : E) : State E →L[ℝ] State E :=
  M1Source.sourceOperator b η μE μM n

def referenceSource (η μE μM : ℝ) : State E →L[ℝ] State E :=
  M1Source.exchangeOperator (isotropicCoupling η) (M1Source.response μE μM)

def physicalSymbol (b : ℝ) (n : E) (dir : ι → E) (ν θ : ι → ℝ) :
    FourierState E →L[ℝ] FourierState E :=
  M1Fourier.symbol ν (fun i => speed b n (dir i)) θ rotation
    (duplicate (radiationBlock 1))
    (fun i => duplicate (radiationBlock (fluxJacobian b n (dir i))))

def referenceSymbol (dir : ι → E) (ν θ : ι → ℝ) :
    FourierState E →L[ℝ] FourierState E :=
  M1Fourier.symbol ν (fun _ => Real.sqrt (1/3 : ℝ)) θ rotation
    (duplicate (radiationBlock 1))
    (fun i => duplicate (radiationBlock (isotropicFluxJacobian (dir i))))

theorem physical_symbol_difference (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n : E) (hn : ‖n‖=1) (dir : ι → E) (hdir : ∀ i, ‖dir i‖=1)
    (ν θ : ι → ℝ) (hν : ∀ i, 0 ≤ ν i) :
    ‖physicalSymbol b n dir ν θ-referenceSymbol dir ν θ‖ ≤
      18*fluxFactor b*(∑ i, ν i) := by
  have hf := parameter_bounds b hb hb1
  rw [show (18:ℝ)=2*2+14 by norm_num]
  apply M1Fourier.symbol_difference_bound ν _ _ θ _ _ _ _ (fluxFactor b) 2 14
    hν (hb.trans hf.1) (by norm_num) (by norm_num) rotation_bound
  · exact (duplicate_bound _).trans ((radiationBlock_bound _).trans (by simp))
  · intro i
    exact speed_difference_bound b hb hb1 n (dir i) hn (hdir i)
  · intro i
    rw [duplicate_difference, radiationBlock_difference]
    exact (duplicate_bound _).trans ((radiationBlock_bound _).trans
      (fluxJacobian_difference_bound b hb hb1 n (dir i) hn (hdir i)))

/-- Simultaneous real/imaginary lift followed by one fixed coordinate change. -/
def lift (S : FourierState E ≃L[ℝ] F) (A : State E →L[ℝ] State E) : F →L[ℝ] F :=
  conjugate S (duplicate A)

theorem lift_bound (S : FourierState E ≃L[ℝ] F) (A : State E →L[ℝ] State E) :
    ‖lift S A‖ ≤ condition S*‖A‖ :=
  (conjugate_bound S _).trans (mul_le_mul_of_nonneg_left (duplicate_bound A) (condition_nonneg S))

theorem lift_add (S : FourierState E ≃L[ℝ] F) (A B : State E →L[ℝ] State E) :
    lift S (A+B)=lift S A+lift S B := by
  simp only [lift, duplicate_add, conjugate_add]

theorem lift_mul (S : FourierState E ≃L[ℝ] F) (A B : State E →L[ℝ] State E) :
    lift S (A*B)=lift S A*lift S B := by
  simp only [lift, duplicate_mul, conjugate_mul]

theorem lift_one (S : FourierState E ≃L[ℝ] F) : lift S (1 : State E →L[ℝ] State E)=1 := by
  simp only [lift, duplicate_one, conjugate_one]

theorem lift_smul (S : FourierState E ≃L[ℝ] F) (a : ℝ) (A : State E →L[ℝ] State E) :
    lift S (a • A)=a • lift S A := by simp only [lift, duplicate_smul, conjugate_smul]

theorem lift_two (S : FourierState E ≃L[ℝ] F) : lift S (2 : State E →L[ℝ] State E)=2 := by
  have h2 : (2 : State E →L[ℝ] State E)=1+1 := by norm_num
  rw [h2, lift_add, lift_one]
  norm_num

theorem source_difference (S : FourierState E ≃L[ℝ] F) (b η μE μM q : ℝ)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (hq : 0 ≤ q) (n : E) (hn : ‖n‖=1) :
    ‖lift S (q • source b η μE μM n)-lift S (q • referenceSource (E := E) η μE μM)‖ ≤
      (condition S*(2*(‖M1Source.response (E := E) μE μM‖+1)))*q*fluxFactor b := by
  letI := NormedSpace.isBoundedSMul (𝕜 := ℝ) (E := (State E →L[ℝ] State E))
  have hcond := condition_nonneg S
  have hid : lift S (q • source b η μE μM n)-lift S (q • referenceSource (E := E) η μE μM)=
      lift S (q • (source b η μE μM n-referenceSource (E := E) η μE μM)) := by
    simp only [lift, conjugate_difference, duplicate_difference, smul_sub]
  rw [hid]
  calc
    _ ≤ condition S*‖q • (source b η μE μM n-referenceSource (E := E) η μE μM)‖ := lift_bound S _
    _ = condition S*(q*‖source b η μE μM n-referenceSource (E := E) η μE μM‖) := by
      congr 1
      simpa only [Real.norm_eq_abs, abs_of_nonneg hq] using
        norm_smul q (source b η μE μM n-referenceSource (E := E) η μE μM)
    _ ≤ condition S*(q*((2*(‖M1Source.response (E := E) μE μM‖+1))*fluxFactor b)) := by
      gcongr
      exact M1Source.m1_source_difference_bound b η μE μM hb hb1 n hn
    _ = _ := by ring

theorem symbol_difference (S : FourierState E ≃L[ℝ] F)
    (b : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1) (n : E) (hn : ‖n‖=1)
    (dir : ι → E) (hdir : ∀ i, ‖dir i‖=1) (ν θ : ι → ℝ) (hν : ∀ i, 0 ≤ ν i) :
    ‖conjugate S (physicalSymbol b n dir ν θ)-conjugate S (referenceSymbol dir ν θ)‖ ≤
      (condition S*18)*fluxFactor b*(∑ i, ν i) := by
  have hh := conjugate_difference_bound S (physicalSymbol b n dir ν θ) (referenceSymbol dir ν θ)
    18 (fluxFactor b*(∑ i, ν i)) (by
      convert physical_symbol_difference b hb hb1 n hn dir hdir ν θ hν using 1 <;> ring)
  convert hh using 1 <;> ring

end RSLA.M1Model
