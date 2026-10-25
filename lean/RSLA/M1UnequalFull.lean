import RSLA.M1UnequalSource
import RSLA.M1UnequalReference

/-! Unequal-opacity full-M1 assembly. The input identification equations specify the fixed
isotropic coordinates and the implicit solves. No source norm, flux norm,
amplification bound, or reference stability estimate is assumed. -/
noncomputable section
namespace RSLA.M1UnequalFull
open scoped BigOperators
open M1Coordinates M1Model M1Coefficients M1ReferenceNorm M1Diffusion
open M1UnequalSource
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

def channelSteps (q rateE : ℝ) : Option ι → ℝ
  | none => q*rateE
  | some _ => q

theorem lift_source_difference
    {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (S : FourierState E ≃L[ℝ] F) (p ell b gamma muE muM rateE q : ℝ)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (hq : 0 ≤ q) (n : E) (hn : ‖n‖=1) :
    ‖lift S (q • sourceOperator p ell b gamma muE muM rateE n)-
      lift S (q • referenceOperator (E := E) p gamma muE muM rateE)‖ ≤
      (condition S*(‖M1Source.response (E := E) rateE 1‖*
        (‖M1Source.response (E := E) muE muM‖+1)*(|ell|+2)))*q*fluxFactor b := by
  letI := NormedSpace.isBoundedSMul (𝕜 := ℝ) (E := (State E →L[ℝ] State E))
  have hc := condition_nonneg S
  have hid : lift S (q • sourceOperator p ell b gamma muE muM rateE n)-
      lift S (q • referenceOperator (E := E) p gamma muE muM rateE)=
      lift S (q • (sourceOperator p ell b gamma muE muM rateE n-
        referenceOperator (E := E) p gamma muE muM rateE)) := by
    simp only [lift, conjugate_difference, duplicate_difference, smul_sub]
  rw [hid]
  calc
    _ ≤ condition S*‖q • (sourceOperator p ell b gamma muE muM rateE n-
        referenceOperator (E := E) p gamma muE muM rateE)‖ := lift_bound S _
    _ = condition S*(q*‖sourceOperator p ell b gamma muE muM rateE n-
        referenceOperator (E := E) p gamma muE muM rateE‖) := by
      congr 1
      simpa only [Real.norm_eq_abs, abs_of_nonneg hq] using
        norm_smul q (sourceOperator p ell b gamma muE muM rateE n-
          referenceOperator (E := E) p gamma muE muM rateE)
    _ ≤ condition S*(q*((‖M1Source.response (E := E) rateE 1‖*
        (‖M1Source.response (E := E) muE muM‖+1)*(|ell|+2))*fluxFactor b)) := by
      gcongr
      exact M1UnequalSource.source_difference p ell b gamma muE muM rateE hb hb1 n hn
    _ = _ := by ring

set_option maxHeartbeats 2000000 in
/-- Opacity-independent finite-time bound for the explicitly defined full M1
source and physical-flux HLL symbol. S is a fixed coordinate equivalence;
its reference identification equations are algebraic model data.
For opacity ratios p=kappa_P/kappa_E and r=kappa_F/kappa_E, take
rateE=1/r, ell=2-r, gamma=1+(p*B0/E0-1)/r. The sufficient drag
condition gamma>=0 includes energy balance and all equal-gray cases. -/
theorem finite_time
    (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)] [Nontrivial (WeightedState ξ)]
    (S : FourierState E ≃L[ℝ] WeightedState ξ)
    (p ell b gamma μE μM rateE α B τ dt c : ℝ) (beam : E) (dir : ι → E) (dx θ ν : ι → ℝ)
    (A R : State E →L[ℝ] State E)
    (a₀ r₀ : WeightedState ξ →L[ℝ] WeightedState ξ)
    (hb : 0 ≤ b) (hb1 : b < 1) (hp : 0 ≤ p) (hg : 0 ≤ gamma) (hrate : 0 < rateE) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hbeam : ‖beam‖=1) (hdir : ∀ i, ‖dir i‖=1)
    (hν : ∀ i, 0 ≤ ν i) (hνcap : ∀ i, ν i ≤ dt*c/dx i)
    (hα : 0 ≤ α) (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hq : 0 < dt*c*α) (hdx : ∀ i, 0 < dx i)
    (hflux : α*fluxFactor b ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (hL : ∀ k, (∑ i, ν i/Real.sqrt 3) ≤
      1+channelSteps (dt*c*α) rateE k/(2+channelSteps (dt*c*α) rateE k*ξ k))
    (hAleft : A*(2+(dt*c*α) • sourceOperator p ell b gamma μE μM rateE beam)=1)
    (hAright : (2+(dt*c*α) • sourceOperator p ell b gamma μE μM rateE beam)*A=1)
    (hRleft : R*(1+(dt*c*α) • sourceOperator p ell b gamma μE μM rateE beam)=1)
    (hRright : (1+(dt*c*α) • sourceOperator p ell b gamma μE μM rateE beam)*R=1)
    (hQ₀ : ∀ x, lift S ((dt*c*α) • referenceOperator (E := E) p gamma μE μM rateE) x=
      M1UnequalReference.source ξ (channelSteps (dt*c*α) rateE) x)
    (hD₀ : ∀ x, (1-conjugate S (referenceSymbol dir ν θ)) x=
      explicitStep ξ (StaticTransport.W (fun i => ν i/Real.sqrt 3) θ) x)
    (ha₀ : (2+lift S ((dt*c*α) • referenceOperator (E := E) p gamma μE μM rateE))*a₀=1)
    (hr₀ : (1+lift S ((dt*c*α) • referenceOperator (E := E) p gamma μE μM rateE))*r₀=1)
    (hr₀action : ∀ x, r₀ x=M1UnequalReference.solve ξ (channelSteps (dt*c*α) rateE) x)
    (steps : ℕ) :
    let K := condition S*uniformBound p μE μM ell (gamma+1)
    let cQ := condition S*(‖M1Source.response (E := E) rateE 1‖*(‖M1Source.response (E := E) μE μM‖+1)*(|ell|+2))
    let cW := condition S*18
    let dbase := 1+2*(1+1/ξ none)
    let dactual := dbase+cW*Real.sqrt 3*(1+1/ξ none)
    let K₀ := (1+dbase^2)/2
    let G := pdars (lift S A) (lift S R)
      (1-conjugate S (physicalSymbol b beam dir ν θ))
    ‖G^steps‖ ≤ K₀*Real.exp
      (K₀*(cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        cW*(K/2*(K*dactual+dbase))*(Fintype.card ι : ℝ)/τ)*c*B*((steps:ℝ)*dt)) := by
  let q := dt*c*α
  let w : ι → ℝ := fun i => ν i/Real.sqrt 3
  let Qf := lift S (q • sourceOperator p ell b gamma μE μM rateE beam)
  let Q₀ := lift S (q • referenceOperator (E := E) p gamma μE μM rateE)
  let df := 1-conjugate S (physicalSymbol b beam dir ν θ)
  let d₀ := 1-conjugate S (referenceSymbol dir ν θ)
  let K := condition S*uniformBound p μE μM ell (gamma+1)
  let cQ := condition S*(‖M1Source.response (E := E) rateE 1‖*(‖M1Source.response (E := E) μE μM‖+1)*(|ell|+2))
  let cW := condition S*18
  let Lstar := 1+1/ξ none
  let dbase := 1+2*Lstar
  let dactual := dbase+cW*Real.sqrt 3*Lstar
  have hxp : 0 < ξ none := hξ.out none
  have hchannels : ∀ k : Option ι, 0 < channelSteps q rateE k := by
    intro k; cases k with
    | none => exact mul_pos hq hrate
    | some i => exact hq
  have hLs : 0 ≤ Lstar := by dsimp [Lstar]; positivity
  have hcQ : 0 ≤ cQ := by dsimp [cQ]; exact mul_nonneg (condition_nonneg S) (by positivity)
  have hcW : 0 ≤ cW := mul_nonneg (condition_nonneg S) (by norm_num)
  have hK : 0 ≤ K := mul_nonneg (condition_nonneg S)
    (uniformBound_nonneg p μE μM ell (gamma+1) hp hμE hμM (by positivity))
  have hd0 : 1 ≤ dbase := by dsimp [dbase]; linarith
  have hdf : 0 ≤ dactual := by dsimp [dactual]; positivity
  have hw : ∀ i, 0 ≤ w i := fun i => div_nonneg (hν i) (Real.sqrt_nonneg _)
  have hcap : (∑ i, w i) ≤ Lstar := M1Reference.arsa_cap (q*rateE) (ξ none) _
    (mul_pos hq hrate).le hxp (by simpa only [channelSteps, q, w, mul_comm] using hL none)
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
  have hresR : ∀ x, R x+q • sourceOperator p ell b gamma μE μM rateE beam (R x)=x := by
    intro x
    have hh := congrArg (fun T : State E →L[ℝ] State E => T x) hRright
    simpa only [ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, ContinuousLinearMap.smul_apply] using hh
  have hresA : ∀ x, (2:ℝ) • A x+q • sourceOperator p ell b gamma μE μM rateE beam (A x)=x := by
    intro x
    have hh := congrArg (fun T : State E →L[ℝ] State E => T x) hAright
    have htwo : (2 : State E →L[ℝ] State E)=1+1 := by norm_num
    simpa only [htwo, ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, ContinuousLinearMap.smul_apply, two_smul] using hh
  have hnr : ‖lift S R‖ ≤ K := (lift_bound S R).trans (mul_le_mul_of_nonneg_left
    (resolvent_bound p ell b gamma μE μM rateE q beam hp hb hb1.le hg hμE hμM hrate.le hq.le hbeam R hresR)
    (condition_nonneg S))
  have hna : ‖lift S A‖ ≤ K/2 := by
    have hh := (lift_bound S A).trans (mul_le_mul_of_nonneg_left
      (half_resolvent_bound p ell b gamma μE μM rateE q beam hp hb hb1.le hg hμE hμM hrate.le hq.le hbeam A hresA)
      (condition_nonneg S))
    convert hh using 1 <;> dsimp [K] <;> ring
  have haf : lift S A*(2+Qf)=1 := by
    have hh := congrArg (lift S) hAleft
    simpa only [lift_mul, lift_add, lift_one, lift_two] using hh
  have hrf : lift S R*(1+Qf)=1 := by
    have hh := congrArg (lift S) hRleft
    simpa only [lift_mul, lift_add, lift_one] using hh
  exact M1UnequalReference.finite_time_from_channels ξ (channelSteps q rateE) Qf Q₀ (lift S A) a₀ (lift S R) r₀ df d₀
    K dbase dactual cQ cW α (fluxFactor b) B τ dt c dx θ ν steps hchannels hν
    hK hd0 hdf hcQ hcW hα hf0 hB hdt hc hτ hdx hflux hcell hL
    haf ha₀ hrf hr₀ hna hnr hnd hnd₀
    (lift_source_difference S p ell b gamma μE μM rateE q hb hb1.le hq.le beam hbeam)
    (hD.trans (by
      apply mul_le_mul_of_nonneg_left _ (mul_nonneg hcW hf0)
      exact Finset.sum_le_sum (fun i _ => hνcap i))) hQ₀ hr₀action hD₀

end RSLA.M1UnequalFull
