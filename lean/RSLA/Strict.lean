import RSLA.Checkerboard

/-! Strict damping (Lemma B.3), with the necessary positive-step assumption. -/
noncomputable section
namespace RSLA.Strict
open Complex ComplexConjugate Channels AllModes
open scoped BigOperators
variable {κ : Type*} [Fintype κ] [DecidableEq κ]

lemma sum_normSq_zero {v : κ → ℂ} (h : (∑ k, normSq (v k)) ≤ 0) : v = 0 := by
  funext k
  apply normSq_eq_zero.mp
  have hh := Finset.single_le_sum (fun j _ => normSq_nonneg (v j)) (Finset.mem_univ k)
  exact le_antisymm (hh.trans h) (normSq_nonneg _)

lemma baseEnergy_pos (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    {x : State κ} (hx : x ≠ 0) : 0 < baseEnergy q K x := by
  by_contra! h
  have hk : ∀ k, channelBase (q k) (K k) (x.1 k) (x.2 k) ≤ 0 := fun k =>
    (Finset.single_le_sum (fun j _ => channelBase_nonneg (hq j) (hK j) _ _)
      (Finset.mem_univ k)).trans h
  apply hx
  apply Prod.ext <;> funext k <;> apply normSq_eq_zero.mp
  all_goals
    have hh := hk k
    unfold channelBase at hh
    have hp : 0 < q k/K k := div_pos (hq k) (hK k)
    have hn := normSq_nonneg (x.1 k)
    have hm := normSq_nonneg (x.2 k)
    nlinarith

lemma stage_update {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u r : ℂ) :
    ((1+K : ℝ) : ℂ)*gasSolve q K u r - (K : ℂ)*radSolve q K u r = u ∧
    -(q : ℂ)*gasSolve q K u r + ((1+q : ℝ) : ℂ)*radSolve q K u r = r := by
  constructor <;> apply Complex.ext <;>
    (try simp only [gasSolve_re, gasSolve_im, radSolve_re, radSolve_im,
      Complex.mul_re, Complex.mul_im, Complex.sub_re, Complex.sub_im,
      Complex.neg_re, Complex.neg_im, Complex.ofReal_re, Complex.ofReal_im,
      zero_mul, mul_zero, add_zero, sub_zero, zero_sub, neg_zero, zero_add]) <;>
    field_simp <;> ring

lemma exchange_injective {q K : ℝ} (hq : 0 < q) (hK : 0 < K) {u v : ℂ}
    (h₁ : ((1+K : ℝ) : ℂ)*u-(K : ℂ)*v = 0)
    (h₂ : -(q : ℂ)*u+((1+q : ℝ) : ℂ)*v = 0) : u=0 ∧ v=0 := by
  have hn : (1+(q : ℂ)+(K : ℂ)) ≠ 0 := by exact_mod_cast (ne_of_gt (by positivity : 0 < 1+q+K))
  push_cast at h₁ h₂
  constructor
  · have hh : (1+(q : ℂ)+(K : ℂ))*u = 0 := by
      linear_combination (1+(q : ℂ))*h₁ + (K : ℂ)*h₂
    exact (mul_eq_zero.mp hh).resolve_left hn
  · have hh : (1+(q : ℂ)+(K : ℂ))*v = 0 := by
      linear_combination (q : ℂ)*h₁ + (1+(K : ℂ))*h₂
    exact (mul_eq_zero.mp hh).resolve_left hn

lemma full_fixed_twice (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ)) {x : State κ} (hx : fullStep q K W x=x) :
    stage q K W (stage q K W x)=x := by
  let z := stage q K W (stage q K W x)
  have he : ∀ k, z.1 k=x.1 k ∧ z.2 k=x.2 k := by
    intro k
    obtain ⟨h₁,h₂⟩ := final_update (hq k) (hK k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)
    have hg := congrArg (fun y : State κ => y.1 k) hx
    have hr := congrArg (fun y : State κ => y.2 k) hx
    change finalGas (q k) (K k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)=x.1 k at hg
    change finalRad (q k) (K k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)=x.2 k at hr
    rw [hg,hr] at h₁ h₂
    obtain ⟨hu,hv⟩ := exchange_injective (hq k) (hK k)
      (u := z.1 k-x.1 k) (v := z.2 k-x.2 k)
      (by push_cast at h₁ ⊢; linear_combination -h₁)
      (by push_cast at h₂ ⊢; linear_combination -h₂)
    exact ⟨sub_eq_zero.mp hu, sub_eq_zero.mp hv⟩
  exact Prod.ext (funext (fun k => (he k).1)) (funext (fun k => (he k).2))

lemma stage_add (q K : κ → ℝ) (W : (κ → ℂ) → (κ → ℂ))
    (hW : ∀ v r, W (v+r)=W v+W r) (x y : State κ) :
    stage q K W (x+y)=stage q K W x+stage q K W y := by
  apply Prod.ext <;> funext k <;>
    simp [stage, hW, gasSolve, radSolve, mul_add, add_div, sub_div] <;> ring

lemma stage_fixed_zero (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ)) {s : ℝ} (hs : 0 < s)
    (hpair : ∀ v, (∑ k, (v k*conj (W v k)).re)=s*(∑ k, normSq (v k)))
    {x : State κ} (hx : stage q K W x=x) : x=0 := by
  have he : ∀ k, x.1 k=x.2 k ∧ W x.2 k=0 := by
    intro k
    obtain ⟨h₁,h₂⟩ := stage_update (hq k) (hK k) (x.1 k) (x.2 k-W x.2 k)
    have hg := congrArg (fun y : State κ => y.1 k) hx
    have hr := congrArg (fun y : State κ => y.2 k) hx
    change gasSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)=x.1 k at hg
    change radSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)=x.2 k at hr
    rw [hg,hr] at h₁ h₂
    have hkn : (K k : ℂ) ≠ 0 := by exact_mod_cast (ne_of_gt (hK k))
    have heq : x.1 k=x.2 k := by
      have hh : (K k : ℂ)*(x.1 k-x.2 k)=0 := by push_cast at h₁; linear_combination h₁
      exact sub_eq_zero.mp ((mul_eq_zero.mp hh).resolve_left hkn)
    refine ⟨heq, ?_⟩
    push_cast at h₂
    linear_combination h₂ + (q k : ℂ)*heq
  have hh := hpair x.2
  simp only [(he _).2, map_zero, mul_zero, Complex.zero_re, Finset.sum_const_zero] at hh
  have hv : x.2=0 := sum_normSq_zero (le_of_eq ((mul_eq_zero.mp hh.symm).resolve_left hs.ne'))
  exact Prod.ext (by funext k; rw [(he k).1, hv]; rfl) hv

lemma pair_real_mul (r : ℝ) (z : ℂ) : (z*conj ((r : ℂ)*z)).re = r*normSq z := by
  simp [Complex.mul_re, Complex.mul_im, normSq_apply]
  ring

lemma stage_negative_zero (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ)) {s : ℝ}
    (hs : ∀ k, s < EnergyBlock.boundary (q k) (K k))
    (hpair : ∀ v, (∑ k, (v k*conj (W v k)).re)=s*(∑ k, normSq (v k)))
    {x : State κ} (hx : stage q K W x = -x) : x=0 := by
  have he : ∀ k, W x.2 k=(EnergyBlock.boundary (q k) (K k) : ℂ)*x.2 k := by
    intro k
    obtain ⟨h₁,h₂⟩ := stage_update (hq k) (hK k) (x.1 k) (x.2 k-W x.2 k)
    have hg := congrArg (fun y : State κ => y.1 k) hx
    have hr := congrArg (fun y : State κ => y.2 k) hx
    change gasSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)= -x.1 k at hg
    change radSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)= -x.2 k at hr
    rw [hg,hr] at h₁ h₂
    have hn : (2+(K k : ℂ)) ≠ 0 := by exact_mod_cast (ne_of_gt (by linarith [hK k] : 0 < 2+K k))
    rw [EnergyBlock.boundary]
    push_cast at h₁ h₂ ⊢
    field_simp [hn]
    linear_combination (q k : ℂ)*h₁ + (2+(K k : ℂ))*h₂
  have hh := hpair x.2
  simp only [he, pair_real_mul] at hh
  have hz : (∑ k, (EnergyBlock.boundary (q k) (K k)-s)*normSq (x.2 k)) = 0 := by
    simp only [sub_mul, Finset.sum_sub_distrib, ← Finset.mul_sum]
    linarith
  have hv : x.2=0 := by
    funext k
    apply normSq_eq_zero.mp
    have hp : 0 < EnergyBlock.boundary (q k) (K k)-s := sub_pos.mpr (hs k)
    have hle := Finset.single_le_sum (s := Finset.univ)
      (fun j _ => mul_nonneg (sub_pos.mpr (hs j)).le (normSq_nonneg (x.2 j))) (Finset.mem_univ k)
    rw [hz] at hle
    nlinarith [normSq_nonneg (x.2 k)]
  apply Prod.ext _ hv
  funext k
  obtain ⟨h₁,_⟩ := stage_update (hq k) (hK k) (x.1 k) (x.2 k-W x.2 k)
  have hg := congrArg (fun y : State κ => y.1 k) hx
  have hr := congrArg (fun y : State κ => y.2 k) hx
  change gasSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)= -x.1 k at hg
  change radSolve (q k) (K k) (x.1 k) (x.2 k-W x.2 k)= -x.2 k at hr
  rw [hg,hr,hv] at h₁
  simp only [Pi.zero_apply, neg_zero, mul_zero, sub_zero] at h₁
  have hn : (2+(K k : ℂ)) ≠ 0 := by exact_mod_cast (ne_of_gt (by linarith [hK k] : 0 < 2+K k))
  have hz : (2+(K k : ℂ))*x.1 k=0 := by push_cast at h₁; linear_combination -h₁
  exact (mul_eq_zero.mp hz).resolve_left hn

lemma no_full_fixed (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ)) {s : ℝ} (hs : 0 < s)
    (hupper : ∀ k, s < EnergyBlock.boundary (q k) (K k))
    (hpair : ∀ v, (∑ k, (v k*conj (W v k)).re)=s*(∑ k, normSq (v k)))
    (hadd : ∀ v r, W (v+r)=W v+W r) {x : State κ} (hx : fullStep q K W x=x) : x=0 := by
  have htwice := full_fixed_twice q K hq hK W hx
  have hy : stage q K W x+x=0 := stage_fixed_zero q K hq hK W hs hpair (by
    rw [stage_add q K W hadd, htwice, add_comm])
  have hm : stage q K W x = -x := eq_neg_of_add_eq_zero_left hy
  exact stage_negative_zero q K hq hK W hupper hpair hm

/-- General strictness principle: the full update has no unit-circle eigenvalue
when its only fixed point is zero. -/
theorem strict_of_no_fixed (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ))
    (hstage : ∀ x, bEnergy q K (stage q K W x) ≤ bEnergy q K x)
    (hnofix : ∀ x, fullStep q K W x=x → x=0) :
    ∀ g : ℂ, ∀ x : State κ, x ≠ 0 → fullStep q K W x=g • x → ‖g‖ < 1 := by
  classical
  intro g x hx heig
  have hc := full_contract q K hq hK W hstage x
  have hp := totalEnergy_pos q K hq hK hx
  have hb : 0 ≤ baseEnergy q K (x-fullStep q K W x) :=
    Finset.sum_nonneg (fun k _ => channelBase_nonneg (hq k) (hK k) _ _)
  rw [heig] at hb
  rw [heig, totalEnergy_smul, normSq_eq_norm_sq] at hc
  have hsquare : ‖g‖^2 ≤ 1 := by
    apply (mul_le_mul_right hp).mp
    nlinarith only [hc, hb]
  have hle : ‖g‖ ≤ 1 := by nlinarith [norm_nonneg g]
  apply lt_of_le_of_ne hle
  intro hequal
  have hnorm : ‖g‖=1 := hequal
  rw [hnorm] at hc
  have hz : x - fullStep q K W x=0 := by
    apply Classical.byContradiction
    intro hn
    have hpos := baseEnergy_pos q K hq hK hn
    rw [heig] at hpos
    nlinarith
  have hfix : fullStep q K W x=x := (sub_eq_zero.mp hz).symm
  exact hx (hnofix x hfix)

lemma W_add {ι : Type*} [Fintype ι] [DecidableEq ι] (w θ : ι → ℝ)
    (x y : Option ι → ℂ) : StaticTransport.W w θ (x+y)=StaticTransport.W w θ x+StaticTransport.W w θ y := by
  funext k
  cases k <;> simp [StaticTransport.W, StaticTransport.symbol, StaticTransport.fluxDot,
    mul_add, Finset.sum_add_distrib] <;> ring

/-- Corrected Lemma B.3. `a>0` is essential; it was implicit in the paper's proof. -/
theorem lemma_B3 (q K : Option (Fin 1) → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    {a θ : ℝ} (ha : 0 < a) (hθ : 0 < θ) (hθπ : θ < Real.pi)
    (hb : ∀ k, a ≤ EnergyBlock.boundary (q k) (K k)) :
    ∀ g : ℂ, ∀ x : State (Option (Fin 1)), x ≠ 0 →
      fullStep q K (StaticTransport.W (fun _ => a/2) (fun _ => θ)) x=g • x → ‖g‖ < 1 := by
  have hw : ∀ i : Fin 1, 0 ≤ a/2 := fun _ => by positivity
  have hbound : ∀ k, 2*∑ _ : Fin 1, a/2 ≤ EnergyBlock.boundary (q k) (K k) := by
    intro k
    simp only [Fin.sum_univ_one]
    convert hb k using 1 <;> ring
  have hs : 0 < a/2*(1-Real.cos θ) := by
    have hh := Real.cos_lt_cos_of_nonneg_of_le_pi (x := 0) (by norm_num) hθπ.le hθ
    rw [Real.cos_zero] at hh
    exact mul_pos (by positivity) (by linarith)
  have hu : a/2*(1-Real.cos θ) < a := by
    have hh := Real.cos_lt_cos_of_nonneg_of_le_pi hθ.le (y := Real.pi) le_rfl hθπ
    rw [Real.cos_pi] at hh
    nlinarith [mul_pos ha (show 0 < 1+Real.cos θ by linarith)]
  have hp : ∀ v : Option (Fin 1) → ℂ,
      (∑ k, (v k*conj (StaticTransport.W (fun _ => a/2) (fun _ => θ) v k)).re) =
      (a/2*(1-Real.cos θ))*(∑ k, normSq (v k)) := by
    intro v
    simpa [StaticTransport.W, StaticTransport.pairing, StaticTransport.norm2, Fintype.sum_option] using
      StaticTransport.symbol_pairing (a/2*(1-Real.cos θ)) (fun _ : Fin 1 => a/2*Real.sin θ) v
  exact strict_of_no_fixed q K hq hK _
    (all_modes_stage_contract q K hq hK (fun _ => a/2) hw hbound (fun _ => θ))
    (fun x hx => no_full_fixed q K hq hK _ hs (fun k => hu.trans_le (hb k)) hp
      (W_add _ _) hx)

/-- At zero transport, an exchange equilibrium remains neutral at every phase.
This witnesses why the displayed strict-damping statement needs a>0. -/
theorem zero_step_neutral (θ : ℝ) :
    let x : State (Option (Fin 1)) := (fun _ => 1, fun _ => 1)
    x ≠ 0 ∧ fullStep (fun _ => 1) (fun _ => 1)
      (StaticTransport.W (fun _ => 0) (fun _ => θ)) x = (1 : ℂ) • x := by
  dsimp only
  constructor
  · intro h
    have hh := congrArg (fun x : State (Option (Fin 1)) => x.1 none) h
    norm_num at hh
  · apply Prod.ext <;> funext k <;> cases k <;>
      norm_num [fullStep, stage, gasSolve, radSolve, finalGas, finalRad,
        StaticTransport.W, StaticTransport.symbol, StaticTransport.fluxDot]

end RSLA.Strict
