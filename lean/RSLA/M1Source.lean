import RSLA.M1Energy
import RSLA.M1Coefficients
import Mathlib.Analysis.InnerProductSpace.LinearMap

/-! Uniform bounds for the full gray exchange source, including the
nonzero-background-flux velocity coupling in the energy equation. -/
noncomputable section
namespace RSLA.M1Source
open scoped RealInnerProductSpace
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

def momentBlock (b η : ℝ) (n : E) : E →L[ℝ] E := η • 1+M1Coefficients.tensor b n

theorem m1_moment_conditions [Nontrivial E] (b η : ℝ) (hb : 0 ≤ b) (hb1 : b ≤ 1)
    (hη : 0 ≤ η) (n : E) (hn : ‖n‖=1) :
    (∀ v, 0 ≤ ⟪v,momentBlock b η n v⟫) ∧
      (∀ v, ‖momentBlock b η n v‖ ≤ (η+1)*‖v‖) := by
  constructor
  · intro v
    change 0 ≤ ⟪v,η • v+M1Coefficients.tensor b n v⟫
    rw [inner_add_right, inner_smul_right, real_inner_self_eq_norm_sq]
    exact add_nonneg (mul_nonneg hη (sq_nonneg _))
      (M1Coefficients.tensor_positive b hb hb1 n v)
  · intro v
    change ‖η • v+M1Coefficients.tensor b n v‖ ≤ _
    calc
      _ ≤ ‖η • v‖+‖M1Coefficients.tensor b n v‖ := norm_add_le _ _
      _ ≤ η*‖v‖+‖M1Coefficients.tensor b n‖*‖v‖ := by
        rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg hη]
        exact add_le_add_left ((M1Coefficients.tensor b n).le_opNorm v) _
      _ ≤ η*‖v‖+1*‖v‖ := by
        gcongr
        exact M1Coefficients.tensor_norm_bound b hb hb1 n hn
      _ = _ := by ring

theorem relaxation_parameter (q : ℝ) (hq : 0 ≤ q) :
    0 ≤ q/(1+q) ∧ q/(1+q) ≤ 1 := by
  constructor
  · positivity
  · apply (div_le_iff₀ (by positivity : 0 < 1+q)).mpr
    linarith

/-- Exact elimination of the radiation unknown in a material/radiation channel.
No commutation or diagonalization of the momentum operator is used. -/
theorem exchange_reduction (q μ : ℝ) (hq : 0 ≤ q) (m r w a b : E)
    (hm : m+(q*μ) • (w-r)=a) (hr : r+q • (r-w)=b) :
    m+(q/(1+q)*μ) • w=a+(q/(1+q)*μ) • b ∧
      r=(1-q/(1+q)) • b+(q/(1+q)) • w := by
  have hd : 1+q ≠ 0 := ne_of_gt (by positivity)
  have hri : (1+q) • r=b+q • w := by
    calc
      (1+q) • r=(r+q • (r-w))+q • w := by module
      _ = b+q • w := by rw [hr]
  have hrv := congrArg (fun z => (1+q)⁻¹ • z) hri
  have hi : (1+q)⁻¹=1-q/(1+q) := by field_simp
  have hit : (1+q)⁻¹*q=q/(1+q) := by ring
  simp only [smul_add, smul_smul, inv_mul_cancel₀ hd, one_smul] at hrv
  rw [hit, hi] at hrv
  have hdiff : w-r=(1-q/(1+q)) • (w-b) := by rw [hrv]; module
  have hcoef : (q*μ)*(1-q/(1+q))=q/(1+q)*μ := by field_simp
  rw [hdiff, smul_smul, hcoef] at hm
  constructor
  · calc
      m+(q/(1+q)*μ) • w = (m+(q/(1+q)*μ) • (w-b))+(q/(1+q)*μ) • b := by module
      _ = a+(q/(1+q)*μ) • b := by rw [hm]
  · exact hrv

theorem coercive_solution_bound (S : E →L[ℝ] E) (a : ℝ) (ha : 0 ≤ a)
    (hpos : ∀ v, 0 ≤ ⟪v, S v⟫) (v y : E) (heq : v+a • S v=y) : ‖v‖ ≤ ‖y‖ := by
  have hh := real_inner_le_norm v y
  have hi : ⟪v,y⟫ = ‖v‖^2+a*⟪v,S v⟫ := by
    rw [← heq]
    simp only [inner_add_right, inner_smul_right, real_inner_self_eq_norm_sq]
  rw [hi] at hh
  have hp := mul_nonneg ha (hpos v)
  by_cases hz : ‖v‖=0
  · rw [hz]; exact norm_nonneg _
  · apply (mul_le_mul_left (lt_of_le_of_ne (norm_nonneg v) (Ne.symm hz))).mp
    nlinarith

/-- The momentum part of the source resolvent has a bound independent of q.
t=q/(1+q), so t is between zero and one. -/
theorem momentum_bound (S : E →L[ℝ] E) (t μ : ℝ)
    (ht : 0 ≤ t) (ht1 : t ≤ 1) (hμ : 0 ≤ μ)
    (hpos : ∀ v, 0 ≤ ⟪v,S v⟫) (v a b : E)
    (heq : v+(t*μ) • S v=a+(t*μ) • b) : ‖v‖ ≤ ‖a‖+μ*‖b‖ := by
  have hh := coercive_solution_bound S (t*μ) (mul_nonneg ht hμ) hpos v _ heq
  calc
    ‖v‖ ≤ ‖a+(t*μ) • b‖ := hh
    _ ≤ ‖a‖+‖(t*μ) • b‖ := norm_add_le _ _
    _ = ‖a‖+(t*μ)*‖b‖ := by rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg (mul_nonneg ht hμ)]
    _ ≤ ‖a‖+μ*‖b‖ := by gcongr; nlinarith

/-- The triangular energy solve retains its full f times velocity coupling. -/
theorem thermal_bound (t μ f : ℝ) (n v : E) (e a b : ℝ)
    (ht : 0 ≤ t) (ht1 : t ≤ 1) (hμ : 0 ≤ μ) (hf : 0 ≤ f)
    (hn : ‖n‖=1) (heq : (1+t*μ)*e=a+t*μ*b-t*μ*f*⟪n,v⟫) :
    |e| ≤ |a|+μ*|b|+μ*f*‖v‖ := by
  have hc : 0 ≤ t*μ := mul_nonneg ht hμ
  have hi : |⟪n,v⟫| ≤ ‖v‖ := by simpa [hn] using abs_real_inner_le_norm n v
  have hs : |e| ≤ |(1+t*μ)*e| := by
    rw [abs_mul, abs_of_nonneg (by positivity : 0 ≤ 1+t*μ)]
    nlinarith [abs_nonneg e]
  calc
    |e| ≤ |(1+t*μ)*e| := hs
    _ = |a+t*μ*b-t*μ*f*⟪n,v⟫| := congrArg abs heq
    _ ≤ |a+t*μ*b|+|t*μ*f*⟪n,v⟫| := abs_sub _ _
    _ ≤ (|a|+|t*μ*b|)+|t*μ*f*⟪n,v⟫| := by gcongr; exact abs_add _ _
    _ = |a|+(t*μ)*|b|+(t*μ*f)*|⟪n,v⟫| := by
      simp only [abs_mul, abs_of_nonneg ht, abs_of_nonneg hμ, abs_of_nonneg hf]
    _ ≤ |a|+μ*|b|+μ*f*‖v‖ := by gcongr <;> nlinarith

def sourceBound (μE μM s : ℝ) : ℝ :=
  (3+2*μE+μE*μM+μM)+(1+s*(1+μM))

/-- A q-independent bound on every component of the full gray resolvent.
The four equations are its explicit triangular reduction, with t=q/(1+q).
They include the entire anisotropic moment block S and the thermal f n.v term. -/
theorem reduced_resolvent_bound (S : E →L[ℝ] E) (n : E)
    (t μE μM f s N : ℝ) (me re ae be : ℝ) (mv rv av bv : E)
    (ht : 0 ≤ t) (ht1 : t ≤ 1) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hf : 0 ≤ f) (hf1 : f ≤ 1) (hs : 0 ≤ s) (hN : 0 ≤ N) (hn : ‖n‖=1)
    (hpos : ∀ v, 0 ≤ ⟪v,S v⟫) (hS : ∀ v, ‖S v‖ ≤ s*‖v‖)
    (hae : |ae| ≤ N) (hbe : |be| ≤ N) (hav : ‖av‖ ≤ N) (hbv : ‖bv‖ ≤ N)
    (hmv : mv+(t*μM) • S mv=av+(t*μM) • bv)
    (hme : (1+t*μE)*me=ae+t*μE*be-t*μE*f*⟪n,mv⟫)
    (hre : re=(1-t)*be+t*(me+f*⟪n,mv⟫))
    (hrv : rv=(1-t) • bv+t • S mv) :
    |me| ≤ sourceBound μE μM s*N ∧ ‖mv‖ ≤ sourceBound μE μM s*N ∧
      |re| ≤ sourceBound μE μM s*N ∧ ‖rv‖ ≤ sourceBound μE μM s*N := by
  have hvm : ‖mv‖ ≤ (1+μM)*N := by
    calc
      _ ≤ ‖av‖+μM*‖bv‖ := momentum_bound S t μM ht ht1 hμM hpos mv av bv hmv
      _ ≤ N+μM*N := by gcongr
      _ = _ := by ring
  have hem : |me| ≤ (1+2*μE+μE*μM)*N := by
    calc
      _ ≤ |ae|+μE*|be|+μE*f*‖mv‖ := thermal_bound t μE f n mv me ae be ht ht1 hμE hf hn hme
      _ ≤ N+μE*N+μE*1*((1+μM)*N) := by gcongr
      _ = _ := by ring
  have hin : |⟪n,mv⟫| ≤ ‖mv‖ := by simpa [hn] using abs_real_inner_le_norm n mv
  have her : |re| ≤ (3+2*μE+μE*μM+μM)*N := by
    calc
      |re| = |(1-t)*be+t*(me+f*⟪n,mv⟫)| := congrArg abs hre
      _ ≤ |(1-t)*be|+|t*(me+f*⟪n,mv⟫)| := abs_add _ _
      _ = (1-t)*|be|+t*|me+f*⟪n,mv⟫| := by
        simp only [abs_mul, abs_of_nonneg (sub_nonneg.mpr ht1), abs_of_nonneg ht]
      _ ≤ (1-t)*|be|+t*(|me|+|f*⟪n,mv⟫|) := by gcongr; exact abs_add _ _
      _ = (1-t)*|be|+t*(|me|+f*|⟪n,mv⟫|) := by rw [abs_mul, abs_of_nonneg hf]
      _ ≤ 1*N+1*((1+2*μE+μE*μM)*N+1*((1+μM)*N)) := by
        gcongr
        · linarith
        · exact hin.trans hvm
      _ = _ := by ring
  have hvr : ‖rv‖ ≤ (1+s*(1+μM))*N := by
    calc
      ‖rv‖ = ‖(1-t) • bv+t • S mv‖ := congrArg norm hrv
      _ ≤ ‖(1-t) • bv‖+‖t • S mv‖ := norm_add_le _ _
      _ = (1-t)*‖bv‖+t*‖S mv‖ := by
        rw [norm_smul, norm_smul, Real.norm_eq_abs, Real.norm_eq_abs,
          abs_of_nonneg (sub_nonneg.mpr ht1), abs_of_nonneg ht]
      _ ≤ 1*N+1*(s*((1+μM)*N)) := by
        gcongr
        · linarith
        · exact (hS mv).trans (mul_le_mul_of_nonneg_left hvm hs)
      _ = _ := by ring
  have hE : 0 ≤ (3+2*μE+μE*μM+μM)*N := by positivity
  have hM : 0 ≤ (1+s*(1+μM))*N := by positivity
  have hEN : 0 ≤ (2+μM)*N := by positivity
  have hMN : 0 ≤ (2+2*μE+μE*μM)*N := by positivity
  dsimp [sourceBound]
  constructor
  · nlinarith
  constructor
  · nlinarith
  exact ⟨by nlinarith, by nlinarith⟩

/-- Uniform bound for the actual implicit gray four-force equations. No bound
on q is assumed. The moment block is S=eta_B I+T for M1. -/
theorem full_resolvent_bound (S : E →L[ℝ] E) (n : E)
    (q μE μM f s N : ℝ) (me re ae be : ℝ) (mv rv av bv : E)
    (hq : 0 ≤ q) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hf : 0 ≤ f) (hf1 : f ≤ 1) (hs : 0 ≤ s) (hN : 0 ≤ N) (hn : ‖n‖=1)
    (hpos : ∀ v, 0 ≤ ⟪v,S v⟫) (hS : ∀ v, ‖S v‖ ≤ s*‖v‖)
    (hae : |ae| ≤ N) (hbe : |be| ≤ N) (hav : ‖av‖ ≤ N) (hbv : ‖bv‖ ≤ N)
    (hmv : mv+(q*μM) • (S mv-rv)=av)
    (hrv : rv+q • (rv-S mv)=bv)
    (hme : me+q*μE*(me+f*⟪n,mv⟫-re)=ae)
    (hre : re+q*(re-(me+f*⟪n,mv⟫))=be) :
    |me| ≤ sourceBound μE μM s*N ∧ ‖mv‖ ≤ sourceBound μE μM s*N ∧
      |re| ≤ sourceBound μE μM s*N ∧ ‖rv‖ ≤ sourceBound μE μM s*N := by
  have ht := relaxation_parameter q hq
  have hv := exchange_reduction q μM hq mv rv (S mv) av bv hmv hrv
  have he := exchange_reduction q μE hq me re (me+f*⟪n,mv⟫) ae be
    (by simpa only [smul_eq_mul] using hme) (by simpa only [smul_eq_mul] using hre)
  simp only [smul_eq_mul] at he
  have hthermal : (1+q/(1+q)*μE)*me=ae+q/(1+q)*μE*be-q/(1+q)*μE*f*⟪n,mv⟫ := by
    nlinarith [he.1]
  exact reduced_resolvent_bound S n (q/(1+q)) μE μM f s N me re ae be mv rv av bv
    ht.1 ht.2 hμE hμM hf hf1 hs hN hn hpos hS hae hbe hav hbv hv.1 hthermal he.2 hv.2

def response (μE μM : ℝ) : (ℝ × E) →L[ℝ] (ℝ × E) :=
  (μE • ContinuousLinearMap.fst ℝ ℝ E).prod (μM • ContinuousLinearMap.snd ℝ ℝ E)

def exchangeOperator {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    (C M : V →L[ℝ] V) : (V × V) →L[ℝ] (V × V) :=
  (M.comp (C.comp (ContinuousLinearMap.fst ℝ V V)-ContinuousLinearMap.snd ℝ V V)).prod
    (ContinuousLinearMap.snd ℝ V V-C.comp (ContinuousLinearMap.fst ℝ V V))

/-- Explicit full gray source operator on ((material energy,velocity),(E,J)). -/
def sourceOperator (b η μE μM : ℝ) (n : E) :
    ((ℝ × E) × (ℝ × E)) →L[ℝ] ((ℝ × E) × (ℝ × E)) :=
  exchangeOperator (M1Coefficients.coupling b η n) (response μE μM)

theorem source_difference_bound {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    (C C₀ M : V →L[ℝ] V) :
    ‖exchangeOperator C M-exchangeOperator C₀ M‖ ≤ (‖M‖+1)*‖C-C₀‖ := by
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  have hid : (exchangeOperator C M-exchangeOperator C₀ M) x =
      (M ((C-C₀) x.1), -((C-C₀) x.1)) := by
    apply Prod.ext
    · change M (C x.1-x.2)-M (C₀ x.1-x.2)=M (C x.1-C₀ x.1)
      simp only [map_sub]
      abel
    · change (x.2-C x.1)-(x.2-C₀ x.1)=-(C x.1-C₀ x.1)
      abel
  have hz : ‖(C-C₀) x.1‖ ≤ ‖C-C₀‖*‖x‖ :=
    ((C-C₀).le_opNorm x.1).trans (mul_le_mul_of_nonneg_left (norm_fst_le x) (norm_nonneg _))
  rw [hid, norm_prod_le_iff]
  constructor
  · calc
      _ ≤ ‖M‖*‖(C-C₀) x.1‖ := M.le_opNorm _
      _ ≤ ‖M‖*(‖C-C₀‖*‖x‖) := mul_le_mul_of_nonneg_left hz (norm_nonneg _)
      _ ≤ (‖M‖+1)*‖C-C₀‖*‖x‖ := by nlinarith [mul_nonneg (norm_nonneg (C-C₀)) (norm_nonneg x)]
  · rw [norm_neg]
    nlinarith [mul_nonneg (norm_nonneg M) (mul_nonneg (norm_nonneg (C-C₀)) (norm_nonneg x))]

/-- The source coefficient perturbation is O(f), in an opacity-independent norm. -/
theorem m1_source_difference_bound [Nontrivial E] (b η μE μM : ℝ)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (n : E) (hn : ‖n‖=1) :
    ‖sourceOperator b η μE μM n-
      exchangeOperator (M1Coefficients.isotropicCoupling (E := E) η) (response μE μM)‖ ≤
      (2*(‖response (E := E) μE μM‖+1))*M1Coefficients.fluxFactor b := by
  calc
    _ ≤ (‖response (E := E) μE μM‖+1)*
        ‖M1Coefficients.coupling b η n-M1Coefficients.isotropicCoupling (E := E) η‖ :=
      source_difference_bound _ _ _
    _ ≤ (‖response (E := E) μE μM‖+1)*(2*M1Coefficients.fluxFactor b) := by
      gcongr
      exact M1Coefficients.coupling_difference_bound b η hb hb1 n hn
    _ = _ := by ring

/-- Operator-norm version for the actual M1 source, uniform over q and beam
orientation, with every velocity coupling included. -/
theorem m1_resolvent_bound [Nontrivial E] (b η μE μM q : ℝ) (n : E)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (hη : 0 ≤ η) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hq : 0 ≤ q) (hn : ‖n‖=1)
    (R : ((ℝ × E) × (ℝ × E)) →L[ℝ] ((ℝ × E) × (ℝ × E)))
    (hres : ∀ x, R x+q • sourceOperator b η μE μM n (R x)=x) :
    ‖R‖ ≤ sourceBound μE μM (η+1) := by
  have hs := m1_moment_conditions b η hb hb1 hη n hn
  have hf := M1Coefficients.parameter_bounds b hb hb1
  apply ContinuousLinearMap.opNorm_le_bound _ (by dsimp [sourceBound]; positivity)
  intro x
  have he := congrArg (fun y => y.1.1) (hres x)
  have hv := congrArg (fun y => y.1.2) (hres x)
  have hre := congrArg (fun y => y.2.1) (hres x)
  have hrv := congrArg (fun y => y.2.2) (hres x)
  change (R x).1.1+q*(μE*((R x).1.1+M1Coefficients.fluxFactor b*⟪n,(R x).1.2⟫-(R x).2.1))=x.1.1 at he
  change (R x).1.2+q • (μM • (momentBlock b η n (R x).1.2-(R x).2.2))=x.1.2 at hv
  change (R x).2.1+q*((R x).2.1-((R x).1.1+M1Coefficients.fluxFactor b*⟪n,(R x).1.2⟫))=x.2.1 at hre
  change (R x).2.2+q • ((R x).2.2-momentBlock b η n (R x).1.2)=x.2.2 at hrv
  rw [← mul_smul] at hv
  have hae : |x.1.1| ≤ ‖x‖ := by
    simpa only [Real.norm_eq_abs] using (norm_fst_le x.1).trans (norm_fst_le x)
  have hbe : |x.2.1| ≤ ‖x‖ := by
    simpa only [Real.norm_eq_abs] using (norm_fst_le x.2).trans (norm_snd_le x)
  have hav := (norm_snd_le x.1).trans (norm_fst_le x)
  have hbv := (norm_snd_le x.2).trans (norm_snd_le x)
  have hh := full_resolvent_bound (momentBlock b η n) n q μE μM
    (M1Coefficients.fluxFactor b) (η+1) ‖x‖
    (R x).1.1 (R x).2.1 x.1.1 x.2.1 (R x).1.2 (R x).2.2 x.1.2 x.2.2
    hq hμE hμM (hb.trans hf.1) hf.2 (by positivity) (norm_nonneg _) hn hs.1 hs.2
    hae hbe hav hbv hv hrv (by simpa only [mul_assoc] using he) hre
  rw [norm_prod_le_iff, norm_prod_le_iff, norm_prod_le_iff]
  simpa only [Real.norm_eq_abs] using And.intro (And.intro hh.1 hh.2.1)
    (And.intro hh.2.2.1 hh.2.2.2)

theorem sourceBound_nonneg (μE μM s : ℝ) (he : 0 ≤ μE) (hm : 0 ≤ μM) (hs : 0 ≤ s) :
    0 ≤ sourceBound μE μM s := by unfold sourceBound; positivity

theorem m1_half_resolvent_bound [Nontrivial E] (b η μE μM q : ℝ) (n : E)
    (hb : 0 ≤ b) (hb1 : b ≤ 1) (hη : 0 ≤ η) (hμE : 0 ≤ μE) (hμM : 0 ≤ μM)
    (hq : 0 ≤ q) (hn : ‖n‖=1)
    (A : ((ℝ × E) × (ℝ × E)) →L[ℝ] ((ℝ × E) × (ℝ × E)))
    (hres : ∀ x, (2:ℝ) • A x+q • sourceOperator b η μE μM n (A x)=x) :
    ‖A‖ ≤ sourceBound μE μM (η+1)/2 := by
  have hh := m1_resolvent_bound b η μE μM (q/2) n hb hb1 hη hμE hμM
    (by positivity) hn ((2:ℝ) • A) (by
      intro x
      change (2:ℝ) • A x+(q/2) • sourceOperator b η μE μM n ((2:ℝ) • A x)=x
      rw [map_smul, smul_smul]
      convert hres x using 1 <;> ring)
  have hnorm : ‖(2:ℝ) • A‖=2*‖A‖ := by
    letI := NormedSpace.isBoundedSMul (𝕜 := ℝ)
      (E := (((ℝ × E) × (ℝ × E)) →L[ℝ] ((ℝ × E) × (ℝ × E))))
    simpa using norm_smul (2:ℝ) A
  rw [hnorm] at hh
  linarith

end RSLA.M1Source
