import RSLA.Fourier

/-! Exchange channels in complex arithmetic. The gas coordinate is scaled so
that exchange equilibrium is u=v. The positive weight q/K implements the same
symmetrization as Appendix B, without requiring square roots. -/
noncomputable section
namespace RSLA.Channels
open Complex ComplexConjugate
open scoped BigOperators
variable {κ : Type*} [Fintype κ] [DecidableEq κ]

abbrev State (κ : Type*) := (κ → ℂ) × (κ → ℂ)

def channelBase (q K : ℝ) (u v : ℂ) : ℝ := q/K * normSq u + normSq v
def channelB (q K : ℝ) (u v : ℂ) : ℝ := channelBase q K u v + q*normSq (u-v)
def baseEnergy (q K : κ → ℝ) (x : State κ) : ℝ := ∑ k, channelBase (q k) (K k) (x.1 k) (x.2 k)
def bEnergy (q K : κ → ℝ) (x : State κ) : ℝ := ∑ k, channelB (q k) (K k) (x.1 k) (x.2 k)
def totalEnergy (q K : κ → ℝ) (x : State κ) : ℝ := baseEnergy q K x + bEnergy q K x

def gasSolve (q K : ℝ) (u r : ℂ) : ℂ := (((1+q : ℝ) : ℂ)*u+(K : ℂ)*r)/((1+q+K : ℝ) : ℂ)
def radSolve (q K : ℝ) (u r : ℂ) : ℂ := ((q : ℂ)*u+((1+K : ℝ) : ℂ)*r)/((1+q+K : ℝ) : ℂ)

@[simp] lemma gasSolve_re (q K : ℝ) (u r : ℂ) :
    (gasSolve q K u r).re = ((1+q)*u.re+K*r.re)/(1+q+K) := by
  rw [gasSolve, Complex.div_ofReal_re]
  simp

@[simp] lemma gasSolve_im (q K : ℝ) (u r : ℂ) :
    (gasSolve q K u r).im = ((1+q)*u.im+K*r.im)/(1+q+K) := by
  rw [gasSolve, Complex.div_ofReal_im]
  simp

@[simp] lemma radSolve_re (q K : ℝ) (u r : ℂ) :
    (radSolve q K u r).re = (q*u.re+(1+K)*r.re)/(1+q+K) := by
  rw [radSolve, Complex.div_ofReal_re]
  simp

@[simp] lemma radSolve_im (q K : ℝ) (u r : ℂ) :
    (radSolve q K u r).im = (q*u.im+(1+K)*r.im)/(1+q+K) := by
  rw [radSolve, Complex.div_ofReal_im]
  simp


def stage (q K : κ → ℝ) (W : (κ → ℂ) → (κ → ℂ)) (x : State κ) : State κ :=
  (fun k => gasSolve (q k) (K k) (x.1 k) (x.2 k - W x.2 k),
   fun k => radSolve (q k) (K k) (x.1 k) (x.2 k - W x.2 k))

lemma complex_schur {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v w : ℂ) :
    channelB q K u v - channelB q K (gasSolve q K u (v-w)) (radSolve q K u (v-w)) =
      q*(2+q+K)/(1+q+K)*normSq (u-v+w/((2+q+K : ℝ) : ℂ)) +
      2*(v*conj w).re - (2+K)/(2+q+K)*normSq w := by
  have hr := Fourier.channel_schur_identity hq hK u.re v.re w.re
  have hi := Fourier.channel_schur_identity hq hK u.im v.im w.im
  simp only [channelB, channelBase, normSq_apply, gasSolve_re, gasSolve_im, radSolve_re, radSolve_im,
    Complex.add_re, Complex.add_im, Complex.sub_re, Complex.sub_im, Complex.div_ofReal_re,
    Complex.div_ofReal_im, Complex.mul_re, Complex.conj_re, Complex.conj_im]
  unfold Fourier.channelEnergy Fourier.stageGas Fourier.stageRad at hr hi
  linear_combination hr + hi

lemma channelB_nonneg {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v : ℂ) :
    0 ≤ channelB q K u v := by
  exact add_nonneg (add_nonneg (mul_nonneg (div_pos hq hK).le (normSq_nonneg _))
    (normSq_nonneg _)) (mul_nonneg hq.le (normSq_nonneg _))

lemma channelBase_nonneg {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v : ℂ) :
    0 ≤ channelBase q K u v := by
  exact add_nonneg (mul_nonneg (div_pos hq hK).le (normSq_nonneg _)) (normSq_nonneg _)

/-- Sufficiency in Lemma B.2 for any transport satisfying the Fourier estimate.
The concrete symbol satisfies this hypothesis in `StaticTransport.fourier_accretive`. -/
theorem stage_contract (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ)) {a : ℝ} (ha : 0 ≤ a)
    (hbound : ∀ k, a ≤ EnergyBlock.boundary (q k) (K k))
    (hW : ∀ v, (∑ k, normSq (W v k)) ≤ a*(∑ k, (v k*conj (W v k)).re))
    (x : State κ) : bEnergy q K (stage q K W x) ≤ bEnergy q K x := by
  let V : κ → ℝ := fun k => (2+K k)/(2+q k+K k)
  have hV : ∀ k, a*V k ≤ 2 := by
    intro k
    have hqk := hq k
    have hKk := hK k
    have hh := hbound k
    rw [EnergyBlock.boundary, le_div_iff₀ (by positivity : 0 < 2+K k)] at hh
    dsimp [V]
    rw [← mul_div_assoc, div_le_iff₀ (by positivity : 0 < 2+q k+K k)]
    nlinarith
  have hc : (∑ k, V k*normSq (W x.2 k)) ≤ 2*(∑ k, (x.2 k*conj (W x.2 k)).re) := by
    by_cases haz : a=0
    · have hh := hW x.2
      rw [haz, zero_mul] at hh
      have hzero : ∀ k, W x.2 k = 0 := by
        intro k
        have hs := Finset.single_le_sum (fun j _ => normSq_nonneg (W x.2 j)) (Finset.mem_univ k)
        have : normSq (W x.2 k) = 0 := le_antisymm (hs.trans hh) (normSq_nonneg _)
        exact normSq_eq_zero.mp this
      simp [hzero]
    · have hap : 0 < a := lt_of_le_of_ne ha (Ne.symm haz)
      have hs := Finset.sum_le_sum (s := Finset.univ) (fun k _ =>
        mul_le_mul_of_nonneg_right (hV k) (normSq_nonneg (W x.2 k)))
      have hid : (∑ k, a*V k*normSq (W x.2 k)) = a*(∑ k, V k*normSq (W x.2 k)) := by
        rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro k _; ring
      rw [hid, ← Finset.mul_sum] at hs
      have hww := hW x.2
      have hmul : a*(∑ k, V k*normSq (W x.2 k)) ≤ a*(2*(∑ k, (x.2 k*conj (W x.2 k)).re)) := by
        nlinarith
      exact (mul_le_mul_left hap).mp hmul
  have hid : bEnergy q K x - bEnergy q K (stage q K W x) =
      (∑ k, (q k)*(2+q k+K k)/(1+q k+K k)*
        normSq (x.1 k-x.2 k+W x.2 k/((2+q k+K k : ℝ) : ℂ))) +
      2*(∑ k, (x.2 k*conj (W x.2 k)).re) - ∑ k, V k*normSq (W x.2 k) := by
    unfold bEnergy
    rw [← Finset.sum_sub_distrib]
    simp only [stage, complex_schur (hq _) (hK _), Finset.sum_sub_distrib,
      Finset.sum_add_distrib, ← Finset.mul_sum, V]
  have hp : 0 ≤ (∑ k, (q k)*(2+q k+K k)/(1+q k+K k)*
        normSq (x.1 k-x.2 k+W x.2 k/((2+q k+K k : ℝ) : ℂ))) := by
    apply Finset.sum_nonneg
    intro k _
    have hqk := hq k
    have hKk := hK k
    exact mul_nonneg (by positivity) (normSq_nonneg _)
  linarith

/-- Solve (I+B)y = x+Bz channel by channel, where z is the result of two stages. -/
def finalGas (q K : ℝ) (u v zu zv : ℂ) : ℂ :=
  (((2+q : ℝ) : ℂ)*u + (K : ℂ)*v + ((2+2*K+q : ℝ) : ℂ)*zu - (K : ℂ)*zv) /
    ((2*(2+q+K) : ℝ) : ℂ)
def finalRad (q K : ℝ) (u v zu zv : ℂ) : ℂ :=
  ((q : ℂ)*u + ((2+K : ℝ) : ℂ)*v - (q : ℂ)*zu + ((2+K+2*q : ℝ) : ℂ)*zv) /
    ((2*(2+q+K) : ℝ) : ℂ)

@[simp] lemma finalGas_re (q K : ℝ) (u v zu zv : ℂ) :
    (finalGas q K u v zu zv).re = ((2+q)*u.re+K*v.re+(2+2*K+q)*zu.re-K*zv.re)/(2*(2+q+K)) := by
  rw [finalGas, Complex.div_ofReal_re]
  simp

@[simp] lemma finalGas_im (q K : ℝ) (u v zu zv : ℂ) :
    (finalGas q K u v zu zv).im = ((2+q)*u.im+K*v.im+(2+2*K+q)*zu.im-K*zv.im)/(2*(2+q+K)) := by
  rw [finalGas, Complex.div_ofReal_im]
  simp

@[simp] lemma finalRad_re (q K : ℝ) (u v zu zv : ℂ) :
    (finalRad q K u v zu zv).re = (q*u.re+(2+K)*v.re-q*zu.re+(2+K+2*q)*zv.re)/(2*(2+q+K)) := by
  rw [finalRad, Complex.div_ofReal_re]
  simp

@[simp] lemma finalRad_im (q K : ℝ) (u v zu zv : ℂ) :
    (finalRad q K u v zu zv).im = (q*u.im+(2+K)*v.im-q*zu.im+(2+K+2*q)*zv.im)/(2*(2+q+K)) := by
  rw [finalRad, Complex.div_ofReal_im]
  simp


def fullStep (q K : κ → ℝ) (W : (κ → ℂ) → (κ → ℂ)) (x : State κ) : State κ :=
  let z := stage q K W (stage q K W x)
  (fun k => finalGas (q k) (K k) (x.1 k) (x.2 k) (z.1 k) (z.2 k),
   fun k => finalRad (q k) (K k) (x.1 k) (x.2 k) (z.1 k) (z.2 k))

/-- The final solve used in the definition is exactly B.4. -/
theorem final_update {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v zu zv : ℂ) :
    ((2+K : ℝ) : ℂ)*finalGas q K u v zu zv - (K : ℂ)*finalRad q K u v zu zv =
      u + ((1+K : ℝ) : ℂ)*zu - (K : ℂ)*zv ∧
    -(q : ℂ)*finalGas q K u v zu zv + ((2+q : ℝ) : ℂ)*finalRad q K u v zu zv =
      v - (q : ℂ)*zu + ((1+q : ℝ) : ℂ)*zv := by
  constructor <;> apply Complex.ext <;>
    simp only [finalGas_re, finalGas_im, finalRad_re, finalRad_im, Complex.mul_re, Complex.mul_im,
      Complex.add_re, Complex.add_im, Complex.sub_re, Complex.sub_im,
      Complex.neg_re, Complex.neg_im, Complex.ofReal_re, Complex.ofReal_im,
      zero_mul, mul_zero, add_zero, sub_zero, zero_sub, neg_zero, zero_add] <;>
    field_simp <;> ring

lemma real_final_identity {q K : ℝ} (hK : 0 < K) (u v zu zv yu yv : ℝ)
    (h₁ : (2+K)*yu-K*yv = u+(1+K)*zu-K*zv)
    (h₂ : -q*yu+(2+q)*yv = v-q*zu+(1+q)*zv) :
    (q/K*u^2+v^2) + Fourier.channelEnergy q K u v -
      ((q/K*yu^2+yv^2) + Fourier.channelEnergy q K yu yv) =
    (q/K*(u-yu)^2+(v-yv)^2) + Fourier.channelEnergy q K (zu-yu) (zv-yv) +
      Fourier.channelEnergy q K u v - Fourier.channelEnergy q K zu zv := by
  unfold Fourier.channelEnergy
  field_simp
  linear_combination -2*q*yu*h₁ - 2*K*yv*h₂

/-- Energy identity behind the two-stage contraction, retaining its strictness term. -/
theorem final_energy_identity {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (u v zu zv : ℂ) :
    let yu := finalGas q K u v zu zv
    let yv := finalRad q K u v zu zv
    channelBase q K u v + channelB q K u v -
      (channelBase q K yu yv + channelB q K yu yv) =
    channelBase q K (u-yu) (v-yv) + channelB q K (zu-yu) (zv-yv) +
      channelB q K u v - channelB q K zu zv := by
  dsimp only
  generalize hyu : finalGas q K u v zu zv = yu
  generalize hyv : finalRad q K u v zu zv = yv
  obtain ⟨h₁, h₂⟩ := final_update hq hK u v zu zv
  rw [hyu, hyv] at h₁ h₂
  have hr₁ := congrArg Complex.re h₁
  have hi₁ := congrArg Complex.im h₁
  have hr₂ := congrArg Complex.re h₂
  have hi₂ := congrArg Complex.im h₂
  simp only [Complex.mul_re, Complex.mul_im, Complex.add_re, Complex.add_im,
    Complex.sub_re, Complex.sub_im, Complex.neg_re, Complex.neg_im,
    Complex.ofReal_re, Complex.ofReal_im, zero_mul, mul_zero, add_zero,
    sub_zero, zero_sub, neg_zero, zero_add] at hr₁ hi₁ hr₂ hi₂
  have hr := real_final_identity hK u.re v.re zu.re zv.re yu.re yv.re hr₁ hr₂
  have hi := real_final_identity hK u.im v.im zu.im zv.im yu.im yv.im hi₁ hi₂
  simp only [channelB, channelBase, normSq_apply, Complex.sub_re, Complex.sub_im]
  unfold Fourier.channelEnergy at hr hi
  linear_combination hr + hi

/-- Concrete-channel version of Lemma B.1, including an energy-loss term. -/
theorem full_contract (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (W : (κ → ℂ) → (κ → ℂ))
    (hs : ∀ x, bEnergy q K (stage q K W x) ≤ bEnergy q K x) (x : State κ) :
    totalEnergy q K (fullStep q K W x) + baseEnergy q K (x-fullStep q K W x) ≤
      totalEnergy q K x := by
  let z := stage q K W (stage q K W x)
  let y := fullStep q K W x
  have hz : bEnergy q K z ≤ bEnergy q K x := (hs (stage q K W x)).trans (hs x)
  have hid : totalEnergy q K x - totalEnergy q K y =
      baseEnergy q K (x-y) + bEnergy q K (z-y) + bEnergy q K x - bEnergy q K z := by
    unfold totalEnergy baseEnergy bEnergy
    simp only [← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro k _
    exact final_energy_identity (hq k) (hK k) (x.1 k) (x.2 k) (z.1 k) (z.2 k)
  have hp : 0 ≤ bEnergy q K (z-y) :=
    Finset.sum_nonneg (fun k _ => channelB_nonneg (hq k) (hK k) _ _)
  change totalEnergy q K y + baseEnergy q K (x-y) ≤ totalEnergy q K x
  linarith

lemma totalEnergy_pos (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    {x : State κ} (hx : x ≠ 0) : 0 < totalEnergy q K x := by
  by_contra! h
  have hbase : baseEnergy q K x ≤ 0 := by
    have hp : 0 ≤ bEnergy q K x := Finset.sum_nonneg (fun k _ => channelB_nonneg (hq k) (hK k) _ _)
    unfold totalEnergy at h
    linarith
  have hk : ∀ k, channelBase (q k) (K k) (x.1 k) (x.2 k) ≤ 0 := by
    intro k
    exact (Finset.single_le_sum (fun j _ => channelBase_nonneg (hq j) (hK j) _ _)
      (Finset.mem_univ k)).trans hbase
  apply hx
  apply Prod.ext
  · funext k
    apply normSq_eq_zero.mp
    have ht := hk k
    unfold channelBase at ht
    have hp : 0 < q k / K k := div_pos (hq k) (hK k)
    have hn := normSq_nonneg (x.2 k)
    have hn' := normSq_nonneg (x.1 k)
    nlinarith
  · funext k
    apply normSq_eq_zero.mp
    have ht := hk k
    unfold channelBase at ht
    have hp : 0 ≤ q k/K k*normSq (x.1 k) := mul_nonneg (div_pos (hq k) (hK k)).le (normSq_nonneg _)
    nlinarith [normSq_nonneg (x.2 k)]

lemma totalEnergy_smul (q K : κ → ℝ) (g : ℂ) (x : State κ) :
    totalEnergy q K (g • x) = normSq g * totalEnergy q K x := by
  simp only [totalEnergy, baseEnergy, bEnergy, channelB, channelBase,
    Prod.smul_fst, Prod.smul_snd, Pi.smul_apply, smul_eq_mul, ← mul_sub, normSq_mul]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k _
  ring

/-- Complex eigenvalue definition for the concrete update. -/
def SpectrallyStable (F : State κ → State κ) : Prop :=
  ∀ g : ℂ, ∀ x : State κ, x ≠ 0 → F x = g • x → ‖g‖ ≤ 1

lemma stable_of_energy_contract (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (F : State κ → State κ) (hF : ∀ x, totalEnergy q K (F x) ≤ totalEnergy q K x) :
    SpectrallyStable F := by
  intro g x hx heig
  have hp := totalEnergy_pos q K hq hK hx
  have hh := hF x
  rw [heig, totalEnergy_smul] at hh
  have hg : normSq g ≤ 1 := by nlinarith
  rw [normSq_eq_norm_sq] at hg
  nlinarith [norm_nonneg g]

end RSLA.Channels
