import RSLA.M1Model

/-! Full-M1 assembly. The input identification equations specify the fixed
isotropic coordinates and the implicit solves. No source norm, flux norm,
amplification bound, or reference stability estimate is assumed. -/
noncomputable section
namespace RSLA.M1Full
open scoped BigOperators
open M1Coordinates M1Model M1Coefficients M1ReferenceNorm M1Diffusion
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

set_option maxHeartbeats 2000000 in
/-- Opacity-independent finite-time bound for the explicitly defined full M1
source and physical-flux HLL symbol. S is a fixed coordinate equivalence;
its two reference identification equations are algebraic model data. -/
theorem finite_time
    (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)] [Nontrivial (WeightedState ξ)]
    (S : FourierState E ≃L[ℝ] WeightedState ξ)
    (b η μE μM ξmax α B τ dt c : ℝ) (beam : E) (dir : ι → E) (dx θ ν : ι → ℝ)
    (A R : State E →L[ℝ] State E)
    (a₀ r₀ : WeightedState ξ →L[ℝ] WeightedState ξ)
    (hb : 0 ≤ b) (hb1 : b < 1) (hη : 0 ≤ η) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hbeam : ‖beam‖=1) (hdir : ∀ i, ‖dir i‖=1)
    (hmax : ∀ k, ξ k ≤ ξmax)
    (hν : ∀ i, 0 ≤ ν i) (hνcap : ∀ i, ν i ≤ dt*c/dx i)
    (hα : 0 ≤ α) (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hq : 0 < dt*c*α) (hdx : ∀ i, 0 < dx i)
    (hflux : α*fluxFactor b ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (hL : (∑ i, ν i/Real.sqrt 3) ≤ 1+(dt*c*α)/(2+ξmax*(dt*c*α)))
    (hAleft : A*(2+(dt*c*α) • M1Model.source b η μE μM beam)=1)
    (hAright : (2+(dt*c*α) • M1Model.source b η μE μM beam)*A=1)
    (hRleft : R*(1+(dt*c*α) • M1Model.source b η μE μM beam)=1)
    (hRright : (1+(dt*c*α) • M1Model.source b η μE μM beam)*R=1)
    (hQ₀ : ∀ x, lift S ((dt*c*α) • referenceSource (E := E) η μE μM) x=
      M1ReferenceNorm.source ξ (dt*c*α) x)
    (hD₀ : ∀ x, (1-conjugate S (referenceSymbol dir ν θ)) x=
      explicitStep ξ (StaticTransport.W (fun i => ν i/Real.sqrt 3) θ) x)
    (ha₀ : (2+lift S ((dt*c*α) • referenceSource (E := E) η μE μM))*a₀=1)
    (hr₀ : (1+lift S ((dt*c*α) • referenceSource (E := E) η μE μM))*r₀=1)
    (hr₀action : ∀ x, r₀ x=solve ξ (dt*c*α) x)
    (steps : ℕ) :
    let K := condition S*M1Source.sourceBound μE μM (η+1)
    let cQ := condition S*(2*(‖M1Source.response (E := E) μE μM‖+1))
    let cW := condition S*18
    let dbase := 1+2*(1+1/ξmax)
    let dactual := dbase+cW*Real.sqrt 3*(1+1/ξmax)
    let K₀ := (1+dbase^2)/2
    let G := pdars (lift S A) (lift S R)
      (1-conjugate S (physicalSymbol b beam dir ν θ))
    ‖G^steps‖ ≤ K₀*Real.exp
      (K₀*(cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        cW*(K/2*(K*dactual+dbase))*(Fintype.card ι : ℝ)/τ)*c*B*((steps:ℝ)*dt)) := by
  let q := dt*c*α
  let w : ι → ℝ := fun i => ν i/Real.sqrt 3
  let Qf := lift S (q • M1Model.source b η μE μM beam)
  let Q₀ := lift S (q • referenceSource (E := E) η μE μM)
  let df := 1-conjugate S (physicalSymbol b beam dir ν θ)
  let d₀ := 1-conjugate S (referenceSymbol dir ν θ)
  let K := condition S*M1Source.sourceBound μE μM (η+1)
  let cQ := condition S*(2*(‖M1Source.response (E := E) μE μM‖+1))
  let cW := condition S*18
  let Lstar := 1+1/ξmax
  let dbase := 1+2*Lstar
  let dactual := dbase+cW*Real.sqrt 3*Lstar
  have hxp : 0 < ξmax := (hξ.out none).trans_le (hmax none)
  have hLs : 0 ≤ Lstar := by dsimp [Lstar]; positivity
  have hcQ : 0 ≤ cQ := by dsimp [cQ]; exact mul_nonneg (condition_nonneg S) (by positivity)
  have hcW : 0 ≤ cW := mul_nonneg (condition_nonneg S) (by norm_num)
  have hK : 0 ≤ K := mul_nonneg (condition_nonneg S)
    (M1Source.sourceBound_nonneg μE μM (η+1) hμE hμM (by positivity))
  have hd0 : 1 ≤ dbase := by dsimp [dbase]; linarith
  have hdf : 0 ≤ dactual := by dsimp [dactual]; positivity
  have hw : ∀ i, 0 ≤ w i := fun i => div_nonneg (hν i) (Real.sqrt_nonneg _)
  have hcap : (∑ i, w i) ≤ Lstar := M1Reference.arsa_cap q ξmax _ hq.le hxp hL
  have hsqrt : 0 < Real.sqrt (3:ℝ) := Real.sqrt_pos.mpr (by norm_num)
  have hsum : (∑ i, ν i)=Real.sqrt 3*(∑ i, w i) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    dsimp [w]
    field_simp
  have hsumcap : (∑ i, ν i) ≤ Real.sqrt 3*Lstar := by rw [hsum]; gcongr
  have hf := parameter_bounds b hb hb1.le
  have hf0 : 0 ≤ fluxFactor b := hb.trans hf.1
  have hD : ‖df-d₀‖ ≤ cW*fluxFactor b*(∑ i, ν i) := by
    dsimp [df, d₀]
    rw [sub_sub_sub_cancel_left, norm_sub_rev]
    exact M1Model.symbol_difference S b hb hb1.le beam hbeam dir hdir ν θ hν
  have hnd₀ : ‖d₀‖ ≤ dbase := by
    apply ContinuousLinearMap.opNorm_le_bound _ (by linarith)
    intro x
    change ‖(1-conjugate S (referenceSymbol dir ν θ)) x‖ ≤ _
    rw [hD₀]
    exact (explicit_norm_bound ξ w θ hw x).trans (by
      have hn := norm_nonneg x
      dsimp [dbase]
      nlinarith)
  have hsum0 : 0 ≤ ∑ i, ν i := Finset.sum_nonneg (fun i _ => hν i)
  have hf1 := hf.2
  have hnd : ‖df‖ ≤ dactual := by
    calc
      ‖df‖ ≤ ‖df-d₀‖+‖d₀‖ := by
        simpa only [sub_add_cancel] using norm_add_le (df-d₀) d₀
      _ ≤ cW*fluxFactor b*(∑ i, ν i)+dbase := add_le_add hD hnd₀
      _ ≤ cW*1*(Real.sqrt 3*Lstar)+dbase := by gcongr
      _ = dactual := by dsimp [dactual]; ring
  have hresR : ∀ x, R x+q • M1Model.source b η μE μM beam (R x)=x := by
    intro x
    have hh := congrArg (fun T : State E →L[ℝ] State E => T x) hRright
    simpa only [ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, ContinuousLinearMap.smul_apply] using hh
  have hresA : ∀ x, (2:ℝ) • A x+q • M1Model.source b η μE μM beam (A x)=x := by
    intro x
    have hh := congrArg (fun T : State E →L[ℝ] State E => T x) hAright
    have htwo : (2 : State E →L[ℝ] State E)=1+1 := by norm_num
    simpa only [htwo, ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, ContinuousLinearMap.smul_apply, two_smul] using hh
  have hnr : ‖lift S R‖ ≤ K := (lift_bound S R).trans (mul_le_mul_of_nonneg_left
    (M1Source.m1_resolvent_bound b η μE μM q beam hb hb1.le hη hμE hμM hq.le hbeam R hresR)
    (condition_nonneg S))
  have hna : ‖lift S A‖ ≤ K/2 := by
    have hh := (lift_bound S A).trans (mul_le_mul_of_nonneg_left
      (M1Source.m1_half_resolvent_bound b η μE μM q beam hb hb1.le hη hμE hμM hq.le hbeam A hresA)
      (condition_nonneg S))
    convert hh using 1 <;> dsimp [K] <;> ring
  have haf : lift S A*(2+Qf)=1 := by
    have hh := congrArg (lift S) hAleft
    simpa only [lift_mul, lift_add, lift_one, lift_two] using hh
  have hrf : lift S R*(1+Qf)=1 := by
    have hh := congrArg (lift S) hRleft
    simpa only [lift_mul, lift_add, lift_one] using hh
  exact M1Theorem.finite_time_from_arsa ξ Qf Q₀ (lift S A) a₀ (lift S R) r₀ df d₀
    ξmax K dbase dactual cQ cW α (fluxFactor b) B τ dt c dx θ ν steps hmax hν
    hK hd0 hdf hcQ hcW hα hf0 hB hdt hc hτ hq hdx hflux hcell hL
    haf ha₀ hrf hr₀ hna hnr hnd hnd₀
    (M1Model.source_difference S b η μE μM q hb hb1.le hq.le beam hbeam)
    (hD.trans (by
      apply mul_le_mul_of_nonneg_left _ (mul_nonneg hcW hf0)
      exact Finset.sum_le_sum (fun i _ => hνcap i))) hQ₀ hr₀action hD₀

end RSLA.M1Full
