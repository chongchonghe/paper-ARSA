import RSLA.StaticTransport
import RSLA.Channels

/-! Assembly of the static multidimensional Fourier stability argument. -/
noncomputable section
namespace RSLA.AllModes
open Complex ComplexConjugate
open scoped BigOperators
open Channels
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The common squared weighted norm, independent of Fourier phase. -/
def weightedNorm (q K : Option ι → ℝ) (x : State (Option ι)) : ℝ :=
  Real.sqrt (totalEnergy q K x)

lemma accretive_sum (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (v : Option ι → ℂ) :
    (∑ k, normSq (StaticTransport.W w θ v k)) ≤
      (2*∑ i, w i)*(∑ k, (v k*conj (StaticTransport.W w θ v k)).re) := by
  simpa [StaticTransport.norm2, StaticTransport.pairing, Fintype.sum_option] using
    StaticTransport.fourier_accretive w θ hw v

/-- Lemma B.2, sufficiency for the explicit complex symbol in any dimension. -/
theorem all_modes_stage_contract (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hb : ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k))
    (θ : ι → ℝ) (x : State (Option ι)) :
    bEnergy q K (stage q K (StaticTransport.W w θ) x) ≤ bEnergy q K x := by
  exact stage_contract q K hq hK _
    (mul_nonneg (by norm_num) (Finset.sum_nonneg (fun i _ => hw i))) hb (accretive_sum w θ hw) x

/-- Proposition B.1, the common-energy bound, with the phase-independent weight explicit. -/
theorem all_modes_energy_contract (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hb : ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k))
    (θ : ι → ℝ) (x : State (Option ι)) :
    totalEnergy q K (fullStep q K (StaticTransport.W w θ) x) ≤ totalEnergy q K x := by
  have hh := full_contract q K hq hK (StaticTransport.W w θ)
    (all_modes_stage_contract q K hq hK w hw hb θ) x
  have hp : 0 ≤ baseEnergy q K (x - fullStep q K (StaticTransport.W w θ) x) :=
    Finset.sum_nonneg (fun k _ => channelBase_nonneg (hq k) (hK k) _ _)
  linarith

/-- The same weighted length is nonincreasing for every Fourier phase. -/
theorem common_weighted_norm (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hb : ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k))
    (θ : ι → ℝ) (x : State (Option ι)) :
    weightedNorm q K (fullStep q K (StaticTransport.W w θ) x) ≤ weightedNorm q K x :=
  Real.sqrt_le_sqrt (all_modes_energy_contract q K hq hK w hw hb θ x)

/-- Proposition B.1, sufficiency for every complex eigenvalue of every mode. -/
theorem all_modes_spectrally_stable (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hb : ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k)) (θ : ι → ℝ) :
    SpectrallyStable (fullStep q K (StaticTransport.W w θ)) :=
  stable_of_energy_contract q K hq hK _ (all_modes_energy_contract q K hq hK w hw hb θ)

lemma W_checkerboard (w : ι → ℝ) (v : Option ι → ℂ) :
    StaticTransport.W w (fun _ => Real.pi) v = fun k => ((2*∑ i, w i : ℝ) : ℂ)*v k := by
  funext k
  cases k <;> simp [StaticTransport.W, StaticTransport.symbol, StaticTransport.fluxDot,
    Finset.sum_mul, mul_comm] <;> norm_num

/-- A vector supported on a single exchange channel. -/
def singleState (k : Option ι) (u v : ℂ) : State (Option ι) :=
  (Pi.single k u, Pi.single k v)

lemma bEnergy_single (q K : Option ι → ℝ) (k : Option ι) (u v : ℂ) :
    bEnergy q K (singleState k u v) = channelB (q k) (K k) u v := by
  classical
  unfold bEnergy singleState
  have hid : ∀ j, channelB (q j) (K j) ((Pi.single k u : Option ι → ℂ) j) ((Pi.single k v : Option ι → ℂ) j) =
      if j=k then channelB (q k) (K k) u v else 0 := by
    intro j
    by_cases hj : j=k
    · subst j; simp
    · simp [Pi.single_apply, hj, Ne.symm hj, channelB, channelBase]
  simp_rw [hid]
  simp

lemma stage_single (q K : Option ι → ℝ) (a : ℝ) (k : Option ι) (u v : ℂ) :
    stage q K (fun v j => (a : ℂ)*v j) (singleState k u v) =
    singleState k (gasSolve (q k) (K k) u (v-(a : ℂ)*v))
      (radSolve (q k) (K k) u (v-(a : ℂ)*v)) := by
  apply Prod.ext <;> funext j <;> by_cases hj : j=k
  all_goals simp_all [stage, singleState, Pi.single_apply, gasSolve, radSolve]

/-- Lemma B.2, both directions, including the checkerboard witness for necessity. -/
theorem all_modes_stage_iff (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i) :
    (∀ θ x, bEnergy q K (stage q K (StaticTransport.W w θ) x) ≤ bEnergy q K x) ↔
      ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k) := by
  constructor
  · intro hs k
    let a : ℝ := 2*∑ i, w i
    let u : ℂ := 1 - (a : ℂ)/((2+q k+K k : ℝ) : ℂ)
    have hh := hs (fun _ => Real.pi) (singleState k u 1)
    have hwc : StaticTransport.W w (fun _ => Real.pi) = fun v j => (a : ℂ)*v j := by
      funext v j
      exact congrFun (W_checkerboard w v) j
    rw [hwc, stage_single, bEnergy_single, bEnergy_single] at hh
    have he := complex_schur (hq k) (hK k) u 1 (a : ℂ)
    have hu : u - 1 + (a : ℂ)/((2+q k+K k : ℝ) : ℂ) = 0 := by dsimp [u]; ring
    rw [hu] at he
    simp only [normSq_zero, mul_zero, add_zero, zero_add, one_mul, Complex.conj_ofReal,
      Complex.ofReal_re, normSq_ofReal, mul_one] at he hh
    have hsign : 0 ≤ 2*a-(2+K k)/(2+q k+K k)*a^2 := by nlinarith
    exact (Fourier.checkerboard_schur_iff (hq k) (hK k)
      (mul_nonneg (by norm_num) (Finset.sum_nonneg (fun i _ => hw i)))).mp hsign
  · intro hb θ x
    exact all_modes_stage_contract q K hq hK w hw hb θ x

end RSLA.AllModes
