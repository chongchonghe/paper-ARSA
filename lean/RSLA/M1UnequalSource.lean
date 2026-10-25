import RSLA.M1Model

/-! Unequal opacity source estimates with all zero-background-velocity couplings.
p is kappa_P/kappa_E; ell is 2-kappa_F/kappa_E. -/
noncomputable section
namespace RSLA.M1UnequalSource
set_option maxHeartbeats 1000000
open scoped RealInnerProductSpace
open M1Source
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]

def coupling (p ell b gamma : ℝ) (n : E) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  (p • ContinuousLinearMap.fst ℝ ℝ E+
    (ell*M1Coefficients.fluxFactor b) • (innerSL ℝ n).comp (ContinuousLinearMap.snd ℝ ℝ E)).prod
    ((gamma • (1 : E →L[ℝ] E)+M1Coefficients.tensor b n).comp (ContinuousLinearMap.snd ℝ ℝ E))

def referenceCoupling (p gamma : ℝ) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  (p • ContinuousLinearMap.fst ℝ ℝ E).prod
    (((gamma+1/3) • (1 : E →L[ℝ] E)).comp (ContinuousLinearMap.snd ℝ ℝ E))

theorem coupling_difference (p ell b gamma : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (n : E) (hn : ‖n‖=1) :
    ‖coupling p ell b gamma n-referenceCoupling (E := E) p gamma‖ ≤
      (|ell|+2)*M1Coefficients.fluxFactor b := by
  have hfb := M1Coefficients.parameter_bounds b hb hb1
  have hf : 0 ≤ M1Coefficients.fluxFactor b := hb.trans hfb.1
  have hf2 : (M1Coefficients.fluxFactor b)^2 ≤ M1Coefficients.fluxFactor b := by nlinarith [hfb.2]
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  have hid : (coupling p ell b gamma n-referenceCoupling (E := E) p gamma) x=
      (ell*M1Coefficients.fluxFactor b*⟪n,x.2⟫,
       (M1Coefficients.tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)) x.2) := by
    apply Prod.ext
    · change p*x.1+ell*M1Coefficients.fluxFactor b*⟪n,x.2⟫-p*x.1=_
      ring
    · change gamma • x.2+M1Coefficients.tensor b n x.2-(gamma+1/3) • x.2=
        M1Coefficients.tensor b n x.2-(1/3 : ℝ) • x.2
      module
  rw [hid, norm_prod_le_iff]
  have hi : |⟪n,x.2⟫| ≤ ‖x‖ :=
    (by simpa [hn] using abs_real_inner_le_norm n x.2 : |⟪n,x.2⟫| ≤ ‖x.2‖).trans (norm_snd_le x)
  constructor
  · simp only [Real.norm_eq_abs, abs_mul, abs_of_nonneg hf]
    have hh := mul_le_mul_of_nonneg_left hi (mul_nonneg (abs_nonneg ell) hf)
    nlinarith [mul_nonneg hf (norm_nonneg x)]
  · calc
      _ ≤ ‖M1Coefficients.tensor b n-(1/3 : ℝ) • (1 : E →L[ℝ] E)‖*‖x.2‖ :=
        ContinuousLinearMap.le_opNorm _ _
      _ ≤ (2*(M1Coefficients.fluxFactor b)^2)*‖x‖ := by
        gcongr
        · exact M1Coefficients.tensor_difference_bound b hb hb1 n hn
        · exact norm_snd_le x
      _ ≤ ((|ell|+2)*M1Coefficients.fluxFactor b)*‖x‖ := by
        gcongr
        nlinarith [abs_nonneg ell]

def sourceOperator (p ell b gamma muE muM rateE : ℝ) (n : E) :
    M1Model.State E →L[ℝ] M1Model.State E :=
  M1Coordinates.duplicate (response rateE 1)*exchangeOperator (coupling p ell b gamma n) (response muE muM)

def referenceOperator (p gamma muE muM rateE : ℝ) :
    M1Model.State E →L[ℝ] M1Model.State E :=
  M1Coordinates.duplicate (response rateE 1)*exchangeOperator (referenceCoupling p gamma) (response muE muM)

theorem source_difference_rates {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    (C C₀ M R : V →L[ℝ] V) :
    ‖M1Coordinates.duplicate R*exchangeOperator C M-M1Coordinates.duplicate R*exchangeOperator C₀ M‖ ≤
      ‖R‖*(‖M‖+1)*‖C-C₀‖ := by
  rw [← mul_sub]
  calc
    _ ≤ ‖M1Coordinates.duplicate R‖*‖exchangeOperator C M-exchangeOperator C₀ M‖ := norm_mul_le _ _
    _ ≤ ‖R‖*((‖M‖+1)*‖C-C₀‖) := by
      gcongr
      · exact M1Coordinates.duplicate_bound R
      · exact M1Source.source_difference_bound C C₀ M
    _ = _ := by ring

theorem source_difference (p ell b gamma muE muM rateE : ℝ)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (n : E) (hn : ‖n‖=1) :
    ‖sourceOperator p ell b gamma muE muM rateE n-referenceOperator (E := E) p gamma muE muM rateE‖ ≤
      (‖response (E := E) rateE 1‖*(‖response (E := E) muE muM‖+1)*(|ell|+2))*M1Coefficients.fluxFactor b := by
  have h1 := source_difference_rates (coupling p ell b gamma n)
    (referenceCoupling (E := E) p gamma) (response (E := E) muE muM) (response (E := E) rateE 1)
  have h2 := mul_le_mul_of_nonneg_left (coupling_difference p ell b gamma hb hb1 n hn)
    (show 0 ≤ ‖response (E := E) rateE 1‖*(‖response (E := E) muE muM‖+1) by positivity)
  exact h1.trans (by convert h2 using 1 <;> ring)

def velocityBound (muM : ℝ) := 1+muM
def thermalBound (muE muM ell : ℝ) := 1+muE+muE*|ell| *velocityBound muM
def radiationEnergyBound (p muE muM ell : ℝ) := 1+p*thermalBound muE muM ell+|ell| *velocityBound muM
def radiationFluxBound (muM s : ℝ) := 1+s*velocityBound muM
def uniformBound (p muE muM ell s : ℝ) :=
  velocityBound muM+thermalBound muE muM ell+radiationEnergyBound p muE muM ell+radiationFluxBound muM s

theorem thermal_bound (p ell t mu f : ℝ) (n v : E) (e a b : ℝ)
    (hp : 0 ≤ p) (ht : 0 ≤ t) (ht1 : t ≤ 1) (hmu : 0 ≤ mu)
    (hf : 0 ≤ f) (hf1 : f ≤ 1) (hn : ‖n‖=1)
    (heq : (1+t*mu*p)*e=a+t*mu*b-t*mu*ell*f*⟪n,v⟫) :
    |e| ≤ |a|+mu*|b|+mu*|ell| *‖v‖ := by
  have hi : |⟪n,v⟫| ≤ ‖v‖ := by simpa [hn] using abs_real_inner_le_norm n v
  have htm : t*mu ≤ mu := by nlinarith
  have hs : |e| ≤ |(1+t*mu*p)*e| := by
    rw [abs_mul, abs_of_nonneg (by positivity : 0 ≤ 1+t*mu*p)]
    nlinarith [mul_nonneg (mul_nonneg ht hmu) hp, abs_nonneg e]
  calc
    |e| ≤ |(1+t*mu*p)*e| := hs
    _ = |a+t*mu*b-t*mu*ell*f*⟪n,v⟫| := congrArg abs heq
    _ ≤ |a+t*mu*b|+|t*mu*ell*f*⟪n,v⟫| := abs_sub _ _
    _ ≤ (|a|+|t*mu*b|)+|t*mu*ell*f*⟪n,v⟫| := by gcongr; exact abs_add _ _
    _ = |a|+(t*mu)*|b|+(t*mu)*|ell| *f*|⟪n,v⟫| := by
      simp only [abs_mul, abs_of_nonneg ht, abs_of_nonneg hmu, abs_of_nonneg hf]
    _ ≤ |a|+mu*|b|+mu*|ell| *1*‖v‖ := by gcongr
    _ = _ := by ring

/-- Two independent relaxation fractions tE,tF, both in [0,1]. -/
theorem reduced_bound (S : E →L[ℝ] E) (n : E)
    (p ell tE tF muE muM f s N : ℝ) (me re ae be : ℝ) (mv rv av bv : E)
    (hp : 0 ≤ p) (htE : 0 ≤ tE) (htE1 : tE ≤ 1) (htF : 0 ≤ tF) (htF1 : tF ≤ 1)
    (he : 0 ≤ muE) (hm : 0 ≤ muM) (hf : 0 ≤ f) (hf1 : f ≤ 1)
    (hs : 0 ≤ s) (hN : 0 ≤ N) (hn : ‖n‖=1)
    (hpos : ∀ v, 0 ≤ ⟪v,S v⟫) (hS : ∀ v, ‖S v‖ ≤ s*‖v‖)
    (hae : |ae| ≤ N) (hbe : |be| ≤ N) (hav : ‖av‖ ≤ N) (hbv : ‖bv‖ ≤ N)
    (hmv : mv+(tF*muM) • S mv=av+(tF*muM) • bv)
    (hme : (1+tE*muE*p)*me=ae+tE*muE*be-tE*muE*ell*f*⟪n,mv⟫)
    (hre : re=(1-tE)*be+tE*(p*me+ell*f*⟪n,mv⟫))
    (hrv : rv=(1-tF) • bv+tF • S mv) :
    |me| ≤ uniformBound p muE muM ell s*N ∧ ‖mv‖ ≤ uniformBound p muE muM ell s*N ∧
      |re| ≤ uniformBound p muE muM ell s*N ∧ ‖rv‖ ≤ uniformBound p muE muM ell s*N := by
  have hv : ‖mv‖ ≤ velocityBound muM*N := by
    calc
      _ ≤ ‖av‖+muM*‖bv‖ := momentum_bound S tF muM htF htF1 hm hpos mv av bv hmv
      _ ≤ N+muM*N := by gcongr
      _ = _ := by unfold velocityBound; ring
  have hmeb : |me| ≤ thermalBound muE muM ell*N := by
    calc
      _ ≤ |ae|+muE*|be|+muE*|ell| *‖mv‖ :=
        thermal_bound p ell tE muE f n mv me ae be hp htE htE1 he hf hf1 hn hme
      _ ≤ N+muE*N+muE*|ell| *(velocityBound muM*N) := by gcongr
      _ = _ := by unfold thermalBound; ring
  have hi : |⟪n,mv⟫| ≤ ‖mv‖ := by simpa [hn] using abs_real_inner_le_norm n mv
  have hv0 : 0 ≤ velocityBound muM := by unfold velocityBound; positivity
  have he0 : 0 ≤ thermalBound muE muM ell := by unfold thermalBound; positivity
  have hr0 : 0 ≤ radiationEnergyBound p muE muM ell := by unfold radiationEnergyBound; positivity
  have hj0 : 0 ≤ radiationFluxBound muM s := by unfold radiationFluxBound; positivity
  have hrb : |re| ≤ radiationEnergyBound p muE muM ell*N := by
    rw [hre]
    calc
      _ ≤ |(1-tE)*be|+|tE*(p*me+ell*f*⟪n,mv⟫)| := abs_add _ _
      _ = (1-tE)*|be|+tE*|p*me+ell*f*⟪n,mv⟫| := by
        simp only [abs_mul, abs_of_nonneg (sub_nonneg.mpr htE1), abs_of_nonneg htE]
      _ ≤ (1-tE)*|be|+tE*(|p*me|+|ell*f*⟪n,mv⟫|) := by gcongr; exact abs_add _ _
      _ = (1-tE)*|be|+tE*(p*|me|+|ell| *f*|⟪n,mv⟫|) := by
        simp only [abs_mul, abs_of_nonneg hp, abs_of_nonneg hf]
      _ ≤ 1*N+1*(p*(thermalBound muE muM ell*N)+|ell| *1*(velocityBound muM*N)) := by
        gcongr
        · linarith
        · exact hi.trans hv
      _ = _ := by unfold radiationEnergyBound; ring
  have hjb : ‖rv‖ ≤ radiationFluxBound muM s*N := by
    rw [hrv]
    calc
      _ ≤ ‖(1-tF) • bv‖+‖tF • S mv‖ := norm_add_le _ _
      _ = (1-tF)*‖bv‖+tF*‖S mv‖ := by
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_of_nonneg (sub_nonneg.mpr htF1), abs_of_nonneg htF]
      _ ≤ 1*N+1*(s*(velocityBound muM*N)) := by
        gcongr
        · linarith
        · exact (hS mv).trans (mul_le_mul_of_nonneg_left hv hs)
      _ = _ := by unfold radiationFluxBound; ring
  have hvN := mul_nonneg hv0 hN
  have heN := mul_nonneg he0 hN
  have hrN := mul_nonneg hr0 hN
  have hjN := mul_nonneg hj0 hN
  unfold uniformBound
  exact ⟨by nlinarith, by nlinarith, by nlinarith, by nlinarith⟩

/-- Full unequal-rate source equations reduce to the two bounded fractions. -/
theorem full_bound (S : E →L[ℝ] E) (n : E)
    (p ell qE qF muE muM f s N : ℝ) (me re ae be : ℝ) (mv rv av bv : E)
    (hp : 0 ≤ p) (hqE : 0 ≤ qE) (hqF : 0 ≤ qF) (he : 0 ≤ muE) (hm : 0 ≤ muM)
    (hf : 0 ≤ f) (hf1 : f ≤ 1) (hs : 0 ≤ s) (hN : 0 ≤ N) (hn : ‖n‖=1)
    (hpos : ∀ v, 0 ≤ ⟪v,S v⟫) (hS : ∀ v, ‖S v‖ ≤ s*‖v‖)
    (hae : |ae| ≤ N) (hbe : |be| ≤ N) (hav : ‖av‖ ≤ N) (hbv : ‖bv‖ ≤ N)
    (hmv : mv+(qF*muM) • (S mv-rv)=av)
    (hrv : rv+qF • (rv-S mv)=bv)
    (hme : me+qE*muE*(p*me+ell*f*⟪n,mv⟫-re)=ae)
    (hre : re+qE*(re-(p*me+ell*f*⟪n,mv⟫))=be) :
    |me| ≤ uniformBound p muE muM ell s*N ∧ ‖mv‖ ≤ uniformBound p muE muM ell s*N ∧
      |re| ≤ uniformBound p muE muM ell s*N ∧ ‖rv‖ ≤ uniformBound p muE muM ell s*N := by
  have htE := relaxation_parameter qE hqE
  have htF := relaxation_parameter qF hqF
  have hv := exchange_reduction qF muM hqF mv rv (S mv) av bv hmv hrv
  have ht := exchange_reduction qE muE hqE me re (p*me+ell*f*⟪n,mv⟫) ae be
    (by simpa only [smul_eq_mul] using hme) (by simpa only [smul_eq_mul] using hre)
  simp only [smul_eq_mul] at ht
  have hthermal : (1+qE/(1+qE)*muE*p)*me=
      ae+qE/(1+qE)*muE*be-qE/(1+qE)*muE*ell*f*⟪n,mv⟫ := by nlinarith [ht.1]
  exact reduced_bound S n p ell (qE/(1+qE)) (qF/(1+qF)) muE muM f s N me re ae be mv rv av bv
    hp htE.1 htE.2 htF.1 htF.2 he hm hf hf1 hs hN hn hpos hS hae hbe hav hbv hv.1 hthermal ht.2 hv.2

theorem uniformBound_nonneg (p muE muM ell s : ℝ)
    (hp : 0 ≤ p) (he : 0 ≤ muE) (hm : 0 ≤ muM) (hs : 0 ≤ s) :
    0 ≤ uniformBound p muE muM ell s := by
  unfold uniformBound radiationEnergyBound radiationFluxBound thermalBound velocityBound
  positivity

/-- Actual unequal-source resolvent, uniformly over its implicit step. -/
theorem resolvent_bound (p ell b gamma muE muM rateE q : ℝ) (n : E)
    (hp : 0 ≤ p) (hb : 0 ≤ b) (hb1 : b ≤ 1) (hg : 0 ≤ gamma)
    (he : 0 ≤ muE) (hm : 0 ≤ muM) (hrate : 0 ≤ rateE) (hq : 0 ≤ q) (hn : ‖n‖=1)
    (R : M1Model.State E →L[ℝ] M1Model.State E)
    (hres : ∀ x, R x+q • sourceOperator p ell b gamma muE muM rateE n (R x)=x) :
    ‖R‖ ≤ uniformBound p muE muM ell (gamma+1) := by
  have hs := m1_moment_conditions b gamma hb hb1 hg n hn
  have hf := M1Coefficients.parameter_bounds b hb hb1
  apply ContinuousLinearMap.opNorm_le_bound _ (uniformBound_nonneg p muE muM ell (gamma+1) hp he hm (by positivity))
  intro x
  have hme := congrArg (fun y => y.1.1) (hres x)
  have hmv := congrArg (fun y => y.1.2) (hres x)
  have hre := congrArg (fun y => y.2.1) (hres x)
  have hrv := congrArg (fun y => y.2.2) (hres x)
  change (R x).1.1+q*(rateE*(muE*(p*(R x).1.1+ell*M1Coefficients.fluxFactor b*⟪n,(R x).1.2⟫-(R x).2.1)))=x.1.1 at hme
  change (R x).1.2+q • ((1:ℝ) • (muM • (momentBlock b gamma n (R x).1.2-(R x).2.2)))=x.1.2 at hmv
  change (R x).2.1+q*(rateE*((R x).2.1-(p*(R x).1.1+ell*M1Coefficients.fluxFactor b*⟪n,(R x).1.2⟫)))=x.2.1 at hre
  change (R x).2.2+q • ((1:ℝ) • ((R x).2.2-momentBlock b gamma n (R x).1.2))=x.2.2 at hrv
  rw [one_smul, ← mul_smul] at hmv
  rw [one_smul] at hrv
  have hae : |x.1.1| ≤ ‖x‖ := by simpa only [Real.norm_eq_abs] using (norm_fst_le x.1).trans (norm_fst_le x)
  have hbe : |x.2.1| ≤ ‖x‖ := by simpa only [Real.norm_eq_abs] using (norm_fst_le x.2).trans (norm_snd_le x)
  have hav := (norm_snd_le x.1).trans (norm_fst_le x)
  have hbv := (norm_snd_le x.2).trans (norm_snd_le x)
  have hh := full_bound (momentBlock b gamma n) n p ell (q*rateE) q muE muM
    (M1Coefficients.fluxFactor b) (gamma+1) ‖x‖
    (R x).1.1 (R x).2.1 x.1.1 x.2.1 (R x).1.2 (R x).2.2 x.1.2 x.2.2
    hp (mul_nonneg hq hrate) hq he hm (hb.trans hf.1) hf.2 (by positivity) (norm_nonneg _) hn hs.1 hs.2
    hae hbe hav hbv hmv hrv (by nlinarith [hme]) (by nlinarith [hre])
  rw [norm_prod_le_iff, norm_prod_le_iff, norm_prod_le_iff]
  simpa only [Real.norm_eq_abs] using And.intro (And.intro hh.1 hh.2.1) (And.intro hh.2.2.1 hh.2.2.2)

theorem half_resolvent_bound (p ell b gamma muE muM rateE q : ℝ) (n : E)
    (hp : 0 ≤ p) (hb : 0 ≤ b) (hb1 : b ≤ 1) (hg : 0 ≤ gamma)
    (he : 0 ≤ muE) (hm : 0 ≤ muM) (hrate : 0 ≤ rateE) (hq : 0 ≤ q) (hn : ‖n‖=1)
    (A : M1Model.State E →L[ℝ] M1Model.State E)
    (hres : ∀ x, (2:ℝ) • A x+q • sourceOperator p ell b gamma muE muM rateE n (A x)=x) :
    ‖A‖ ≤ uniformBound p muE muM ell (gamma+1)/2 := by
  have hh := resolvent_bound p ell b gamma muE muM rateE (q/2) n hp hb hb1 hg he hm hrate
    (by positivity) hn ((2:ℝ) • A) (by
      intro x
      change (2:ℝ) • A x+(q/2) • sourceOperator p ell b gamma muE muM rateE n ((2:ℝ) • A x)=x
      rw [map_smul, smul_smul]
      convert hres x using 1 <;> ring)
  have hnorm : ‖(2:ℝ) • A‖=2*‖A‖ := by
    letI := NormedSpace.isBoundedSMul (𝕜 := ℝ) (E := (M1Model.State E →L[ℝ] M1Model.State E))
    simpa using norm_smul (2:ℝ) A
  rw [hnorm] at hh
  linarith

end RSLA.M1UnequalSource
