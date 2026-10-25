import RSLA.M1Theorem

/-! Distinct energy and momentum relaxation rates in the fixed P1 norm. -/
noncomputable section
namespace RSLA.M1UnequalReference
open scoped RealInnerProductSpace BigOperators
open M1ReferenceNorm
variable {κ : Type*} [Fintype κ] [DecidableEq κ]

def source (ξ q : κ → ℝ) (x : WeightedState ξ) : WeightedState ξ :=
  (fun k => (q k*ξ k) • (x.1 k-x.2 k), fun k => q k • (x.2 k-x.1 k))

def solve (ξ q : κ → ℝ) (x : WeightedState ξ) : WeightedState ξ :=
  (fun k => Channels.gasSolve (q k) (q k*ξ k) (x.1 k) (x.2 k),
   fun k => Channels.radSolve (q k) (q k*ξ k) (x.1 k) (x.2 k))

theorem base_energy_eq (ξ q : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (hq : ∀ k, 0 < q k) (x : WeightedState ξ) :
    Channels.baseEnergy q (fun k => q k*ξ k) x=‖x‖^2 := by
  rw [M1ReferenceNorm.norm_sq]
  unfold Channels.baseEnergy Channels.channelBase
  apply Finset.sum_congr rfl
  intro k _
  have hq0 : q k ≠ 0 := ne_of_gt (hq k)
  have hk : ξ k ≠ 0 := ne_of_gt (hξ.out k)
  field_simp
  ring

theorem source_pairing (ξ q : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (x y : WeightedState ξ) :
    ⟪x,source ξ q y⟫=∑ k, q k*⟪x.1 k-x.2 k,y.1 k-y.2 k⟫ := by
  change pairing ξ x (source ξ q y)=_
  unfold pairing source
  apply Finset.sum_congr rfl
  intro k _
  simp only [real_inner_smul_right, inner_sub_left, inner_sub_right]
  have hk : ξ k ≠ 0 := ne_of_gt (hξ.out k)
  field_simp
  ring

theorem b_energy_eq (ξ q : κ → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    (hq : ∀ k, 0 < q k) (x : WeightedState ξ) :
    Channels.bEnergy q (fun k => q k*ξ k) x=‖x‖^2+⟪x,source ξ q x⟫ := by
  rw [← base_energy_eq ξ q hq, source_pairing]
  simp only [Channels.bEnergy, Channels.channelB, Channels.baseEnergy,
    real_inner_self_eq_norm_sq, Complex.normSq_eq_norm_sq, Finset.sum_add_distrib]

theorem channel_condition (q ξ L : ℝ) (hq : 0 < q) (hξ : 0 < ξ)
    (hL : L ≤ 1+q/(2+q*ξ)) : 2*L ≤ EnergyBlock.boundary q (q*ξ) := by
  have hh : 2*(1+q/(2+q*ξ))=EnergyBlock.boundary q (q*ξ) := by
    unfold EnergyBlock.boundary
    field_simp
    ring
  rw [← hh]
  linarith

section Assembly
variable {ι : Type*} [Fintype ι] [DecidableEq ι]
open M1Diffusion

/-- Unequal-channel diffusion theorem. All source and transport estimates are
primitive certificates; the all-mode reference stability conclusion is proved
from the individual ARSA channel conditions, not assumed. -/
theorem finite_time_from_channels
    (ξ q : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    [Nontrivial (WeightedState ξ)]
    (Qf Q₀ af a₀ rf r₀ df d₀ : WeightedState ξ →L[ℝ] WeightedState ξ)
    (K dbase dactual cQ cW α f B τ dt c : ℝ) (dx θ ν : ι → ℝ) (n : ℕ)
    (hq : ∀ k, 0 < q k) (hν : ∀ i, 0 ≤ ν i)
    (hK : 0 ≤ K) (hd0 : 1 ≤ dbase) (hdf : 0 ≤ dactual)
    (hcQ : 0 ≤ cQ) (hcW : 0 ≤ cW) (hα : 0 ≤ α) (hf : 0 ≤ f)
    (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hdx : ∀ i, 0 < dx i) (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (hL : ∀ k, (∑ i, ν i/Real.sqrt 3) ≤ 1+q k/(2+q k*ξ k))
    (haf : af*(2+Qf)=1) (ha₀ : (2+Q₀)*a₀=1)
    (hrf : rf*(1+Qf)=1) (hr₀ : (1+Q₀)*r₀=1)
    (hna : ‖af‖ ≤ K/2) (hnr : ‖rf‖ ≤ K)
    (hnd : ‖df‖ ≤ dactual) (hnd₀ : ‖d₀‖ ≤ dbase)
    (hQ : ‖Qf-Q₀‖ ≤ cQ*(dt*c*α)*f)
    (hD : ‖df-d₀‖ ≤ cW*f*(∑ i, dt*c/dx i))
    (hQ₀ : ∀ x, Q₀ x=source ξ q x)
    (hr₀action : ∀ x, r₀ x=solve ξ q x)
    (hd₀action : ∀ x, d₀ x=explicitStep ξ
      (StaticTransport.W (fun i => ν i/Real.sqrt 3) θ) x) :
    ‖(pdars af rf df)^n‖ ≤ ((1+dbase^2)/2)*Real.exp
      (((1+dbase^2)/2)*(cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        cW*(K/2*(K*dactual+dbase))*(Fintype.card ι : ℝ)/τ)*c*B*((n:ℝ)*dt)) := by
  apply M1Theorem.finite_time_from_reference_stage Qf Q₀ af a₀ rf r₀ df d₀ K dbase dactual
    cQ cW α f B τ dt c dx n hK hd0 hdf hcQ hcW hα hf hB hdt hc hτ hdx
    hflux hcell haf ha₀ hrf hr₀ hna hnr hnd hnd₀ hQ hD
  · intro x y
    rw [hQ₀, hQ₀, source_pairing, source_pairing]
    apply Finset.sum_congr rfl
    intro k _
    rw [real_inner_comm]
  · intro x
    rw [hQ₀, source_pairing]
    exact Finset.sum_nonneg (fun k _ => mul_nonneg (hq k).le real_inner_self_nonneg)
  · intro x
    simp only [hQ₀, hr₀action, hd₀action]
    rw [← b_energy_eq ξ q hq, ← b_energy_eq ξ q hq]
    exact AllModes.all_modes_stage_contract q (fun k => q k*ξ k) hq
      (fun k => mul_pos (hq k) (hξ.out k)) (fun i => ν i/Real.sqrt 3)
      (fun i => div_nonneg (hν i) (Real.sqrt_nonneg _))
      (fun k => channel_condition (q k) (ξ k) _ (hq k) (hξ.out k) (hL k)) θ x

end Assembly
end RSLA.M1UnequalReference
