import RSLA.AllModes
import RSLA.M1Energy

/-! The enlarged reference CFL is inherited from the already formalized
P1 channel proof. It is not a postulated power bound for full M1. -/
noncomputable section
namespace RSLA.M1Reference
open scoped BigOperators
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

theorem arsa_cap (q ξ L : ℝ) (hq : 0 ≤ q) (hξ : 0 < ξ)
    (hL : L ≤ 1+q/(2+ξ*q)) : L ≤ 1+1/ξ := by
  apply hL.trans
  apply add_le_add_left
  apply (div_le_div_iff₀ (by positivity : 0 < 2+ξ*q) hξ).mpr
  nlinarith

/-- The exact isotropic ARSA bound implies every source-channel stage bound. -/
theorem channel_bounds (q ξmax L : ℝ) (ξ : Option ι → ℝ)
    (hq : 0 < q) (hξ : ∀ k, 0 < ξ k) (hmax : ∀ k, ξ k ≤ ξmax)
    (hL : L ≤ 1+q/(2+ξmax*q)) :
    ∀ k, 2*L ≤ EnergyBlock.boundary q (q*ξ k) := by
  intro k
  have hk := hξ k
  have hm := hmax k
  have hxp : 0 < ξmax := hk.trans_le hm
  have hden : 2+q*ξ k ≤ 2+ξmax*q := by nlinarith
  have hfrac : q/(2+ξmax*q) ≤ q/(2+q*ξ k) :=
    div_le_div_of_nonneg_left hq.le (by positivity) hden
  calc
    2*L ≤ 2*(1+q/(2+q*ξ k)) := by linarith
    _ = EnergyBlock.boundary q (q*ξ k) := by
      unfold EnergyBlock.boundary
      field_simp
      <;> ring

/-- All Fourier phases of the actual P1 reference symbol, at the enlarged CFL. -/
theorem isotropic_stage_from_arsa (q ξmax : ℝ) (ξ : Option ι → ℝ)
    (hq : 0 < q) (hξ : ∀ k, 0 < ξ k) (hmax : ∀ k, ξ k ≤ ξmax)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hL : (∑ i, w i) ≤ 1+q/(2+ξmax*q))
    (θ : ι → ℝ) (x : Channels.State (Option ι)) :
    Channels.bEnergy (fun _ => q) (fun k => q*ξ k)
      (Channels.stage (fun _ => q) (fun k => q*ξ k) (StaticTransport.W w θ) x) ≤
      Channels.bEnergy (fun _ => q) (fun k => q*ξ k) x := by
  exact AllModes.all_modes_stage_contract (fun _ => q) (fun k => q*ξ k)
    (fun _ => hq) (fun k => mul_pos hq (hξ k)) w hw
    (channel_bounds q ξmax (∑ i, w i) ξ hq hξ hmax hL) θ x

theorem isotropic_full_energy_from_arsa (q ξmax : ℝ) (ξ : Option ι → ℝ)
    (hq : 0 < q) (hξ : ∀ k, 0 < ξ k) (hmax : ∀ k, ξ k ≤ ξmax)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i)
    (hL : (∑ i, w i) ≤ 1+q/(2+ξmax*q))
    (θ : ι → ℝ) (x : Channels.State (Option ι)) :
    Channels.totalEnergy (fun _ => q) (fun k => q*ξ k)
      (Channels.fullStep (fun _ => q) (fun k => q*ξ k) (StaticTransport.W w θ) x) ≤
      Channels.totalEnergy (fun _ => q) (fun k => q*ξ k) x := by
  exact AllModes.all_modes_energy_contract (fun _ => q) (fun k => q*ξ k)
    (fun _ => hq) (fun k => mul_pos hq (hξ k)) w hw
    (channel_bounds q ξmax (∑ i, w i) ξ hq hξ hmax hL) θ x

end RSLA.M1Reference
