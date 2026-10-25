import RSLA.M1Diffusion
import RSLA.M1Energy
import RSLA.M1ReferenceNorm

/-! Final finite-time theorem with the reference power bound derived, rather
than assumed. The input stage certificate is supplied by M1Reference's
checked P1 ARSA inequality; M1Coefficients, M1Source and M1Fourier supply the
primitive bounds for the actual full-M1 operators. -/
noncomputable section
namespace RSLA.M1Theorem
open scoped RealInnerProductSpace BigOperators
open M1Diffusion
variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E] [Nontrivial E]
variable {ι : Type*} [Fintype ι]

/-- The full PD-ARS diffusion-scaling stability theorem. Its hypotheses are
implicit-solve identities, the isotropic reference stage estimate, and
primitive coefficient bounds. No full-M1 amplification or power estimate is
assumed. The exponent has no opacity or mesh-spacing dependence. -/
theorem finite_time_from_reference_stage
    (Qf Q₀ af a₀ rf r₀ df d₀ : E →L[ℝ] E)
    (K dbase dactual cQ cW α f B τ dt c : ℝ) (dx : ι → ℝ) (n : ℕ)
    (hK : 0 ≤ K) (hd0 : 1 ≤ dbase) (hdf : 0 ≤ dactual)
    (hcQ : 0 ≤ cQ) (hcW : 0 ≤ cW) (hα : 0 ≤ α) (hf : 0 ≤ f)
    (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hdx : ∀ i, 0 < dx i) (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (haf : af*(2+Qf)=1) (ha₀ : (2+Q₀)*a₀=1)
    (hrf : rf*(1+Qf)=1) (hr₀ : (1+Q₀)*r₀=1)
    (hna : ‖af‖ ≤ K/2) (hnr : ‖rf‖ ≤ K)
    (hnd : ‖df‖ ≤ dactual) (hnd₀ : ‖d₀‖ ≤ dbase)
    (hQ : ‖Qf-Q₀‖ ≤ cQ*(dt*c*α)*f)
    (hD : ‖df-d₀‖ ≤ cW*f*(∑ i, dt*c/dx i))
    (hsym : ∀ x y, ⟪x,Q₀ y⟫=⟪y,Q₀ x⟫)
    (hpos : ∀ x, 0 ≤ ⟪x,Q₀ x⟫)
    (hstage : ∀ x, ‖r₀ (d₀ x)‖^2+⟪r₀ (d₀ x),Q₀ (r₀ (d₀ x))⟫ ≤
      ‖x‖^2+⟪x,Q₀ x⟫) :
    ‖(pdars af rf df)^n‖ ≤ ((1+dbase^2)/2) * Real.exp
      (((1+dbase^2)/2)*
        (cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        (cW*(K/2*(K*dactual+dbase)))*(Fintype.card ι : ℝ)/τ)*
        c*B*((n : ℝ)*dt)) := by
  have hres₀ (x : E) : r₀ x+Q₀ (r₀ x)=x := by
    have hh := congrArg (fun T : E →L[ℝ] E => T x) hr₀
    simpa only [ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply] using hh
  have htwo : (2 : E →L[ℝ] E)=1+1 := by norm_num
  have hares (x : E) : (2:ℝ) • a₀ x+Q₀ (a₀ x)=x := by
    have hh := congrArg (fun T : E →L[ℝ] E => T x) ha₀
    simpa only [htwo, ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, two_smul] using hh
  have hupdate (x : E) : (2:ℝ) • pdars a₀ r₀ d₀ x+Q₀ (pdars a₀ r₀ d₀ x)=
      x+d₀ (r₀ (d₀ x)) := by
    have hh : (2+Q₀)*pdars a₀ r₀ d₀=1+d₀*r₀*d₀ := by
      unfold pdars
      rw [← mul_assoc, ha₀, one_mul]
    have hx := congrArg (fun T : E →L[ℝ] E => T x) hh
    simpa only [htwo, ContinuousLinearMap.mul_apply, ContinuousLinearMap.add_apply,
      ContinuousLinearMap.one_apply, two_smul] using hx
  have hnr₀ : ‖r₀‖ ≤ 1 := by
    apply ContinuousLinearMap.opNorm_le_bound _ zero_le_one
    intro x
    simpa using M1Energy.resolvent_contract Q₀ r₀ hpos hres₀ x
  have hna₀ := M1Energy.half_resolvent_bound Q₀ a₀ hpos hares
  have href (j : ℕ) : ‖(pdars a₀ r₀ d₀)^j‖ ≤ (1+dbase^2)/2 :=
    M1Energy.reference_power_bound Q₀ r₀ d₀ (pdars a₀ r₀ d₀) dbase hd0 hsym hpos hres₀
      (fun x => d₀.le_of_opNorm_le hnd₀ x) hupdate hstage j
  exact pdars_diffusion_finite_time Qf Q₀ af a₀ rf r₀ df d₀ K ((1+dbase^2)/2)
    dbase dactual cQ cW α f B τ dt c dx n hK (by positivity) (by linarith) hdf
    hcQ hcW hα hf hB hdt hc hτ hdx hflux hcell haf ha₀ hrf hr₀
    hna hna₀ hnr hnr₀ hnd hnd₀ hQ hD href

section ConcreteReference
open M1ReferenceNorm
variable [DecidableEq ι]

/-- Explicit-CFL version: the reference stage inequality is proved for the
actual P1 symbol. The remaining norm hypotheses are local coefficient and
implicit-source estimates, not amplification or power-boundedness assumptions. -/
theorem finite_time_from_arsa
    (ξ : Option ι → ℝ) [hξ : Fact (∀ k, 0 < ξ k)]
    [Nontrivial (WeightedState ξ)]
    (Qf Q₀ af a₀ rf r₀ df d₀ : WeightedState ξ →L[ℝ] WeightedState ξ)
    (ξmax K dbase dactual cQ cW α f B τ dt c : ℝ) (dx θ ν : ι → ℝ) (n : ℕ)
    (hmax : ∀ k, ξ k ≤ ξmax) (hν : ∀ i, 0 ≤ ν i)
    (hK : 0 ≤ K) (hd0 : 1 ≤ dbase) (hdf : 0 ≤ dactual)
    (hcQ : 0 ≤ cQ) (hcW : 0 ≤ cW) (hα : 0 ≤ α) (hf : 0 ≤ f)
    (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hq : 0 < dt*c*α) (hdx : ∀ i, 0 < dx i)
    (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (hL : (∑ i, ν i/Real.sqrt 3) ≤ 1+(dt*c*α)/(2+ξmax*(dt*c*α)))
    (haf : af*(2+Qf)=1) (ha₀ : (2+Q₀)*a₀=1)
    (hrf : rf*(1+Qf)=1) (hr₀ : (1+Q₀)*r₀=1)
    (hna : ‖af‖ ≤ K/2) (hnr : ‖rf‖ ≤ K)
    (hnd : ‖df‖ ≤ dactual) (hnd₀ : ‖d₀‖ ≤ dbase)
    (hQ : ‖Qf-Q₀‖ ≤ cQ*(dt*c*α)*f)
    (hD : ‖df-d₀‖ ≤ cW*f*(∑ i, dt*c/dx i))
    (hQ₀ : ∀ x, Q₀ x=source ξ (dt*c*α) x)
    (hr₀action : ∀ x, r₀ x=solve ξ (dt*c*α) x)
    (hd₀action : ∀ x, d₀ x=explicitStep ξ
      (StaticTransport.W (fun i => ν i/Real.sqrt 3) θ) x) :
    ‖(pdars af rf df)^n‖ ≤ ((1+dbase^2)/2) * Real.exp
      (((1+dbase^2)/2)*
        (cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        (cW*(K/2*(K*dactual+dbase)))*(Fintype.card ι : ℝ)/τ)*
        c*B*((n : ℝ)*dt)) := by
  apply finite_time_from_reference_stage Qf Q₀ af a₀ rf r₀ df d₀ K dbase dactual
    cQ cW α f B τ dt c dx n hK hd0 hdf hcQ hcW hα hf hB hdt hc hτ hdx
    hflux hcell haf ha₀ hrf hr₀ hna hnr hnd hnd₀ hQ hD
  · intro x y
    rw [hQ₀, hQ₀, source_pairing, source_pairing]
    apply Finset.sum_congr rfl
    intro k _
    rw [real_inner_comm]
  · intro x
    rw [hQ₀]
    exact source_positive ξ (dt*c*α) hq.le x
  · intro x
    simp only [hQ₀, hr₀action, hd₀action, solve_explicit_eq_stage]
    rw [← b_energy_eq ξ (dt*c*α) hq, ← b_energy_eq ξ (dt*c*α) hq]
    apply M1Reference.isotropic_stage_from_arsa (dt*c*α) ξmax ξ hq hξ.out hmax
      (fun i => ν i/Real.sqrt 3) _ hL θ x
    intro i
    exact div_nonneg (hν i) (Real.sqrt_nonneg _)

end ConcreteReference
end RSLA.M1Theorem
