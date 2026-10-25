import RSLA.M1Reference
import Mathlib.Analysis.InnerProductSpace.PiL2

/-! The fixed reference norm has weights 1/xi, rather than weights depending
on opacity q. This connects the concrete P1 channel proof to the Hilbert-space
reference-energy argument. -/
noncomputable section
namespace RSLA.M1ReferenceNorm
open scoped RealInnerProductSpace BigOperators
variable {κ : Type*} [Fintype κ] [DecidableEq κ]

def WeightedState (ξ : κ → ℝ) := Channels.State κ

instance (ξ : κ → ℝ) : AddCommGroup (WeightedState ξ) :=
  inferInstanceAs (AddCommGroup (Channels.State κ))
instance (ξ : κ → ℝ) : Module ℝ (WeightedState ξ) :=
  inferInstanceAs (Module ℝ (Channels.State κ))

def pairing (ξ : κ → ℝ) (x y : WeightedState ξ) : ℝ :=
  ∑ k, ((1/ξ k)*⟪x.1 k,y.1 k⟫+⟪x.2 k,y.2 k⟫)

theorem pairing_self (ξ : κ → ℝ) (x : WeightedState ξ) :
    pairing ξ x x = ∑ k, ((1/ξ k)*‖x.1 k‖^2+‖x.2 k‖^2) := by
  simp only [pairing, real_inner_self_eq_norm_sq]

def innerCore (ξ : κ → ℝ) (hξ : ∀ k, 0 < ξ k) : InnerProductSpace.Core ℝ (WeightedState ξ) where
  inner := pairing ξ
  conj_inner_symm x y := by
    simp only [RCLike.conj_to_real, pairing]
    apply Finset.sum_congr rfl
    intro k _
    rw [real_inner_comm (y.1 k) (x.1 k), real_inner_comm (y.2 k) (x.2 k)]
  re_inner_nonneg x := by
    change 0 ≤ pairing ξ x x
    rw [pairing_self]
    apply Finset.sum_nonneg
    intro k _
    have hk := hξ k
    positivity
  add_left x y z := by
    change pairing ξ (x+y) z=pairing ξ x z+pairing ξ y z
    simp only [pairing, Prod.fst_add, Prod.snd_add, Pi.add_apply, inner_add_left,
      mul_add, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro k _
    change (1/ξ k)*⟪x.1 k+y.1 k,z.1 k⟫+⟪x.2 k+y.2 k,z.2 k⟫ = _
    rw [inner_add_left, inner_add_left]
    ring
  smul_left x y r := by
    change pairing ξ (r • x) y=r*pairing ξ x y
    simp only [pairing, Prod.smul_fst, Prod.smul_snd, Pi.smul_apply,
      real_inner_smul_left, Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro k _
    change (1/ξ k)*⟪r • x.1 k,y.1 k⟫+⟪r • x.2 k,y.2 k⟫ = _
    rw [real_inner_smul_left, real_inner_smul_left]
    ring
  definite x hx := by
    change pairing ξ x x=0 at hx
    rw [pairing_self] at hx
    have hz (k : κ) : (1/ξ k)*‖x.1 k‖^2+‖x.2 k‖^2=0 := by
      have hn : ∀ j ∈ Finset.univ, 0 ≤ (1/ξ j)*‖x.1 j‖^2+‖x.2 j‖^2 := by
        intro j _; have hj := hξ j; positivity
      have hh := Finset.single_le_sum hn (Finset.mem_univ k)
      rw [hx] at hh
      exact le_antisymm hh (hn k (Finset.mem_univ k))
    apply Prod.ext <;> funext k
    · have hh := hz k
      have hk : 0 < 1/ξ k := one_div_pos.mpr (hξ k)
      have hn := sq_nonneg ‖x.2 k‖
      have : ‖x.1 k‖^2=0 := by nlinarith [sq_nonneg ‖x.1 k‖]
      exact norm_eq_zero.mp (by nlinarith [norm_nonneg (x.1 k)])
    · have hh := hz k
      have hk : 0 ≤ (1/ξ k)*‖x.1 k‖^2 := by have hp := hξ k; positivity
      exact norm_eq_zero.mp (by nlinarith [norm_nonneg (x.2 k)])

instance (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)] : NormedAddCommGroup (WeightedState ξ) :=
  @InnerProductSpace.Core.toNormedAddCommGroup ℝ (WeightedState ξ) _ _ _ (innerCore ξ hξ.out)

instance (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)] : InnerProductSpace ℝ (WeightedState ξ) :=
  InnerProductSpace.ofCore (innerCore ξ hξ.out)

theorem norm_sq (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)] (x : WeightedState ξ) :
    ‖x‖^2 = ∑ k, ((1/ξ k)*Complex.normSq (x.1 k)+Complex.normSq (x.2 k)) := by
  rw [← real_inner_self_eq_norm_sq]
  change pairing ξ x x=_
  rw [pairing_self]
  simp only [Complex.normSq_eq_norm_sq]

theorem base_energy_eq (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 < q) (x : WeightedState ξ) :
    Channels.baseEnergy (fun _ => q) (fun k => q*ξ k) x=‖x‖^2 := by
  rw [norm_sq]
  unfold Channels.baseEnergy Channels.channelBase
  apply Finset.sum_congr rfl
  intro k _
  have hq0 : q ≠ 0 := ne_of_gt hq
  have hk : ξ k ≠ 0 := ne_of_gt (hξ.out k)
  field_simp
  ring

def source (ξ : κ → ℝ) (q : ℝ) (x : WeightedState ξ) : WeightedState ξ :=
  (fun k => (q*ξ k) • (x.1 k-x.2 k), fun k => q • (x.2 k-x.1 k))

def solve (ξ : κ → ℝ) (q : ℝ) (x : WeightedState ξ) : WeightedState ξ :=
  (fun k => Channels.gasSolve q (q*ξ k) (x.1 k) (x.2 k),
    fun k => Channels.radSolve q (q*ξ k) (x.1 k) (x.2 k))

def explicitStep (ξ : κ → ℝ) (W : (κ → ℂ) → (κ → ℂ))
    (x : WeightedState ξ) : WeightedState ξ := (x.1, x.2-W x.2)

theorem source_pairing (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (x y : WeightedState ξ) :
    ⟪x,source ξ q y⟫ = ∑ k, q*⟪x.1 k-x.2 k,y.1 k-y.2 k⟫ := by
  change pairing ξ x (source ξ q y)=_
  unfold pairing source
  apply Finset.sum_congr rfl
  intro k _
  simp only [real_inner_smul_right, inner_sub_left, inner_sub_right]
  have hk : ξ k ≠ 0 := ne_of_gt (hξ.out k)
  field_simp
  <;> ring

theorem source_positive (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 ≤ q) (x : WeightedState ξ) : 0 ≤ ⟪x,source ξ q x⟫ := by
  rw [source_pairing]
  apply Finset.sum_nonneg
  intro k _
  exact mul_nonneg hq real_inner_self_nonneg

theorem b_energy_eq (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 < q) (x : WeightedState ξ) :
    Channels.bEnergy (fun _ => q) (fun k => q*ξ k) x = ‖x‖^2+⟪x,source ξ q x⟫ := by
  rw [← base_energy_eq ξ q hq, source_pairing]
  simp only [Channels.bEnergy, Channels.channelB, Channels.baseEnergy,
    real_inner_self_eq_norm_sq, Complex.normSq_eq_norm_sq, Finset.sum_add_distrib]

theorem total_energy_eq (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 < q) (x : WeightedState ξ) :
    Channels.totalEnergy (fun _ => q) (fun k => q*ξ k) x = M1Energy.sourceEnergy (source ξ q) x := by
  unfold Channels.totalEnergy
  rw [base_energy_eq ξ q hq, b_energy_eq ξ q hq]
  unfold M1Energy.sourceEnergy
  ring

theorem solve_identity (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 < q) (x : WeightedState ξ) :
    solve ξ q x+source ξ q (solve ξ q x)=x := by
  apply Prod.ext <;> funext k
  · change Channels.gasSolve q (q*ξ k) (x.1 k) (x.2 k) +
      (q*ξ k) • (Channels.gasSolve q (q*ξ k) (x.1 k) (x.2 k)-
        Channels.radSolve q (q*ξ k) (x.1 k) (x.2 k)) = x.1 k
    have hk := hξ.out k
    have hd : 1+q+q*ξ k ≠ 0 := ne_of_gt (by positivity)
    apply Complex.ext <;> simp only [Complex.real_smul, Complex.add_re, Complex.add_im,
      Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re,
      Complex.ofReal_im, zero_mul, mul_zero, add_zero, sub_zero,
      Channels.gasSolve_re, Channels.gasSolve_im, Channels.radSolve_re, Channels.radSolve_im]
    all_goals field_simp; ring
  · change Channels.radSolve q (q*ξ k) (x.1 k) (x.2 k) +
      q • (Channels.radSolve q (q*ξ k) (x.1 k) (x.2 k)-
        Channels.gasSolve q (q*ξ k) (x.1 k) (x.2 k)) = x.2 k
    have hk := hξ.out k
    have hd : 1+q+q*ξ k ≠ 0 := ne_of_gt (by positivity)
    apply Complex.ext <;> simp only [Complex.real_smul, Complex.add_re, Complex.add_im,
      Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im, Complex.ofReal_re,
      Complex.ofReal_im, zero_mul, mul_zero, add_zero, sub_zero,
      Channels.gasSolve_re, Channels.gasSolve_im, Channels.radSolve_re, Channels.radSolve_im]
    all_goals field_simp; ring
theorem solve_explicit_eq_stage (ξ : κ → ℝ) (q : ℝ) (W : (κ → ℂ) → (κ → ℂ))
    (x : WeightedState ξ) : solve ξ q (explicitStep ξ W x)=
      Channels.stage (fun _ => q) (fun k => q*ξ k) W x := rfl

def full (ξ : κ → ℝ) (q : ℝ) (W : (κ → ℂ) → (κ → ℂ))
    (x : WeightedState ξ) : WeightedState ξ := Channels.fullStep (fun _ => q) (fun k => q*ξ k) W x

def radOnly (ξ : κ → ℝ) (v : κ → ℂ) : WeightedState ξ := (0,v)

theorem full_update (ξ : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q : ℝ) (hq : 0 < q) (W : (κ → ℂ) → (κ → ℂ)) (x : WeightedState ξ) :
    let G := full ξ q W
    (2:ℝ) • G x+source ξ q (G x)=x+explicitStep ξ W (solve ξ q (explicitStep ξ W x)) := by
  let z := solve ξ q (explicitStep ξ W (solve ξ q (explicitStep ξ W x)))
  have hz : z+source ξ q z=explicitStep ξ W (solve ξ q (explicitStep ξ W x)) :=
    solve_identity ξ q hq _
  change (2:ℝ) • full ξ q W x+source ξ q (full ξ q W x) = _
  rw [← hz]
  have hstage : z=Channels.stage (fun _ => q) (fun k => q*ξ k) W
      (Channels.stage (fun _ => q) (fun k => q*ξ k) W x) := rfl
  apply Prod.ext <;> funext k
  · change (2:ℝ) • Channels.finalGas q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k) +
      (q*ξ k) • (Channels.finalGas q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)-
        Channels.finalRad q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)) =
      x.1 k+(z.1 k+(q*ξ k) • (z.1 k-z.2 k))
    have hh := (Channels.final_update hq (mul_pos hq (hξ.out k))
      (x.1 k) (x.2 k) (z.1 k) (z.2 k)).1
    simp only [Complex.real_smul, Complex.ofReal_add, Complex.ofReal_ofNat, Complex.ofReal_one] at *
    linear_combination hh
  · change (2:ℝ) • Channels.finalRad q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k) +
      q • (Channels.finalRad q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)-
        Channels.finalGas q (q*ξ k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)) =
      x.2 k+(z.2 k+q • (z.2 k-z.1 k))
    have hh := (Channels.final_update hq (mul_pos hq (hξ.out k))
      (x.1 k) (x.2 k) (z.1 k) (z.2 k)).2
    simp only [Complex.real_smul, Complex.ofReal_add, Complex.ofReal_ofNat, Complex.ofReal_one] at *
    linear_combination hh

section Fourier
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem radiation_norm_le (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (x : WeightedState ξ) : StaticTransport.norm2 x.2 ≤ ‖x‖^2 := by
  rw [norm_sq]
  have hh : (∑ k, Complex.normSq (x.2 k)) ≤
      ∑ k, ((1/ξ k)*Complex.normSq (x.1 k)+Complex.normSq (x.2 k)) := by
    apply Finset.sum_le_sum
    intro k _
    have hk := hξ.out k
    have hn := Complex.normSq_nonneg (x.1 k)
    exact le_add_of_nonneg_left (mul_nonneg (by positivity) hn)
  simpa only [Fintype.sum_option, StaticTransport.norm2] using hh

theorem radiation_only_norm (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (v : Option ι → ℂ) : ‖radOnly ξ v‖^2=StaticTransport.norm2 v := by
  rw [norm_sq]
  simp [radOnly, Fintype.sum_option, StaticTransport.norm2]

theorem explicit_norm_bound (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i) (x : WeightedState ξ) :
    ‖explicitStep ξ (StaticTransport.W w θ) x‖ ≤ (1+2*(∑ i, w i))*‖x‖ := by
  let L := ∑ i, w i
  let s := ∑ i, w i*(1-Real.cos (θ i))
  have hL : 0 ≤ L := Finset.sum_nonneg (fun i _ => hw i)
  have hs : s ≤ 2*L := by
    dsimp [s, L]
    rw [Finset.mul_sum]
    apply Finset.sum_le_sum
    intro i _
    nlinarith [hw i, Real.neg_one_le_cos (θ i)]
  have hW := StaticTransport.fourier_accretive w θ hw x.2
  have hpair : StaticTransport.pairing x.2 (StaticTransport.W w θ x.2)=
      s*StaticTransport.norm2 x.2 := StaticTransport.symbol_pairing _ _ _
  rw [hpair] at hW
  have hv := radiation_norm_le ξ x
  have hv0 : 0 ≤ StaticTransport.norm2 x.2 := by
    unfold StaticTransport.norm2
    exact add_nonneg (Complex.normSq_nonneg _) (Finset.sum_nonneg (fun i _ => Complex.normSq_nonneg _))
  have hsq : StaticTransport.norm2 (StaticTransport.W w θ x.2) ≤
      (2*L*‖x‖)^2 := by
    calc
      _ ≤ (2*L)*(s*StaticTransport.norm2 x.2) := hW
      _ ≤ (2*L)*((2*L)*‖x‖^2) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        exact (mul_le_mul_of_nonneg_right hs hv0).trans
          (mul_le_mul_of_nonneg_left hv (by positivity))
      _ = _ := by ring
  have hn : ‖radOnly ξ (StaticTransport.W w θ x.2)‖ ≤ 2*L*‖x‖ := by
    rw [← radiation_only_norm ξ] at hsq
    exact (sq_le_sq₀ (norm_nonneg _) (by positivity)).mp hsq
  have heq : explicitStep ξ (StaticTransport.W w θ) x =
      x-radOnly ξ (StaticTransport.W w θ x.2) := by
    apply Prod.ext
    · change x.1=x.1-0
      simp
    · rfl
  rw [heq]
  exact (norm_sub_le _ _).trans (by nlinarith)

/-- Concrete all-phase P1 bound in the fixed, opacity-independent norm. -/
theorem isotropic_power_bound (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q ξmax : ℝ) (hq : 0 < q) (hmax : ∀ k, ξ k ≤ ξmax)
    (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hL : (∑ i, w i) ≤ 1+q/(2+ξmax*q)) (n : ℕ) (x : WeightedState ξ) :
    ‖(full ξ q (StaticTransport.W w θ))^[n] x‖ ≤
      ((1+(1+2*(1+1/ξmax))^2)/2)*‖x‖ := by
  have hxp : 0 < ξmax := (hξ.out none).trans_le (hmax none)
  have hcap := M1Reference.arsa_cap q ξmax (∑ i, w i) hq.le hxp hL
  apply M1Energy.functional_uniform_power (source ξ q) (solve ξ q)
    (explicitStep ξ (StaticTransport.W w θ)) (full ξ q (StaticTransport.W w θ))
    (1+2*(1+1/ξmax)) (by have hp := one_div_pos.mpr hxp; linarith)
    (source_positive ξ q hq.le)
    (M1Energy.resolvent_contract _ _ (source_positive ξ q hq.le) (solve_identity ξ q hq))
    (fun y => (explicit_norm_bound ξ w θ hw y).trans (by
      nlinarith [norm_nonneg y]))
    (full_update ξ q hq (StaticTransport.W w θ))
  intro y
  rw [← total_energy_eq ξ q hq, ← total_energy_eq ξ q hq]
  exact M1Reference.isotropic_full_energy_from_arsa q ξmax ξ hq hξ.out hmax w hw hL θ y

/-- Operator form, requiring only identification with the concrete update. -/
theorem isotropic_operator_power_bound (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (q ξmax : ℝ) (hq : 0 < q) (hmax : ∀ k, ξ k ≤ ξmax)
    (w θ : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hL : (∑ i, w i) ≤ 1+q/(2+ξmax*q))
    (G : WeightedState ξ →L[ℝ] WeightedState ξ)
    (hG : ∀ x, G x=full ξ q (StaticTransport.W w θ) x) (n : ℕ) :
    ‖G^n‖ ≤ (1+(1+2*(1+1/ξmax))^2)/2 := by
  apply ContinuousLinearMap.opNorm_le_bound _ (by positivity)
  intro x
  have heq : ∀ k, (G^k) x=(full ξ q (StaticTransport.W w θ))^[k] x := by
    intro k
    induction k with
    | zero => rfl
    | succ k ih =>
      rw [pow_succ', ContinuousLinearMap.mul_apply, ih, hG, Function.iterate_succ_apply']
  rw [heq]
  exact isotropic_power_bound ξ q ξmax hq hmax w θ hw hL n x

end Fourier
end RSLA.M1ReferenceNorm
