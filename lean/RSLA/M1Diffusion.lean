import Mathlib.Analysis.Normed.Ring.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

/-! Opacity-uniform finite-time perturbation bounds for the M1 extension.
The noncommutative argument does not assume that the two amplification
operators commute, are normal, or are individually contractive. -/
noncomputable section
namespace RSLA.M1Diffusion

variable {A : Type*} [NormedRing A]

/-- A power-bounded reference operator controls a nearby operator, including
non-normal perturbations. The auxiliary left factor avoids any diagonalization. -/
theorem power_perturbation_aux (G G₀ : A) (K ε : ℝ)
    (hK : 0 ≤ K) (hε : 0 ≤ ε)
    (hbase : ∀ j : ℕ, ‖G₀ ^ j‖ ≤ K) (hdiff : ‖G - G₀‖ ≤ ε) :
    ∀ n j : ℕ, ‖G₀ ^ j * G ^ n‖ ≤ K * (1 + K * ε) ^ n := by
  intro n
  induction n with
  | zero => intro j; simpa using hbase j
  | succ n ih =>
    intro j
    have hn : ‖G ^ n‖ ≤ K * (1 + K * ε) ^ n := by simpa using ih 0
    have hid : G₀ ^ j * G ^ (n+1) =
        G₀ ^ (j+1) * G ^ n + G₀ ^ j * (G-G₀) * G ^ n := by
      rw [pow_succ' G, pow_succ G₀]
      noncomm_ring
    rw [hid]
    calc
      ‖G₀ ^ (j+1) * G ^ n + G₀ ^ j * (G-G₀) * G ^ n‖
          ≤ ‖G₀ ^ (j+1) * G ^ n‖ + ‖G₀ ^ j * (G-G₀) * G ^ n‖ := norm_add_le _ _
      _ ≤ K * (1+K*ε)^n + K*ε*(K*(1+K*ε)^n) := by
        apply add_le_add (ih (j+1))
        calc
          ‖G₀ ^ j * (G-G₀) * G ^ n‖ ≤ ‖G₀ ^ j‖ * ‖G-G₀‖ * ‖G ^ n‖ := by
            exact (norm_mul_le _ _).trans
              (mul_le_mul_of_nonneg_right (norm_mul_le _ _) (norm_nonneg _))
          _ ≤ K*ε*(K*(1+K*ε)^n) := by gcongr; exact hbase j
      _ = K * (1+K*ε)^(n+1) := by rw [pow_succ]; ring

/-- The finite-time exponential bound, without a one-step contraction assumption. -/
theorem power_perturbation (G G₀ : A) (K ε : ℝ)
    (hK : 0 ≤ K) (hε : 0 ≤ ε)
    (hbase : ∀ j : ℕ, ‖G₀ ^ j‖ ≤ K) (hdiff : ‖G - G₀‖ ≤ ε) (n : ℕ) :
    ‖G ^ n‖ ≤ K * Real.exp (K*ε*(n : ℝ)) := by
  have h := power_perturbation_aux G G₀ K ε hK hε hbase hdiff n 0
  simp only [pow_zero, one_mul] at h
  refine h.trans ?_
  apply mul_le_mul_of_nonneg_left _ hK
  calc
    (1+K*ε)^n ≤ (Real.exp (K*ε))^n := by
      apply pow_le_pow_left₀ (by positivity)
      simpa [add_comm] using Real.add_one_le_exp (K*ε)
    _ = Real.exp (K*ε*(n : ℝ)) := by rw [← Real.exp_nat_mul]; congr 1; ring

section PDARS

/-- The actual PD-ARS algebra, with A=(2I+Q)⁻¹ and R=(I+Q)⁻¹. -/
def pdars (a r d : A) : A := a * (1 + d * r * d)

theorem pdars_difference (af a₀ rf r₀ df d₀ : A) :
    pdars af rf df - pdars a₀ r₀ d₀ =
      (af-a₀)*(1+d₀*r₀*d₀) +
      af*((df-d₀)*rf*df + d₀*(rf-r₀)*df + d₀*r₀*(df-d₀)) := by
  unfold pdars
  noncomm_ring

theorem inverse_difference (a b ia ib : A)
    (ha : ia*a=1) (hb : b*ib=1) : ia-ib = ia*(b-a)*ib := by
  calc
    ia-ib = ia*(b*ib)-(ia*a)*ib := by rw [ha, hb]; simp
    _ = ia*(b-a)*ib := by noncomm_ring

lemma norm_triple (a b c : A) : ‖a*b*c‖ ≤ ‖a‖*‖b‖*‖c‖ :=
  (norm_mul_le _ _).trans
    (mul_le_mul_of_nonneg_right (norm_mul_le _ _) (norm_nonneg _))

variable [NormOneClass A]

theorem inverse_difference_bound (a b ia ib : A) (ka kb δ : ℝ)
    (hka : 0 ≤ ka) (hkb : 0 ≤ kb) (hδ : 0 ≤ δ)
    (ha : ia*a=1) (hb : b*ib=1)
    (hia : ‖ia‖ ≤ ka) (hib : ‖ib‖ ≤ kb) (hdiff : ‖b-a‖ ≤ δ) :
    ‖ia-ib‖ ≤ ka*δ*kb := by
  rw [inverse_difference a b ia ib ha hb]
  exact (norm_triple _ _ _).trans (by gcongr)

theorem pdars_difference_norm (af a₀ rf r₀ df d₀ : A) :
    ‖pdars af rf df - pdars a₀ r₀ d₀‖ ≤
      ‖af-a₀‖*(1+‖d₀‖*‖r₀‖*‖d₀‖) +
      ‖af‖*(‖df-d₀‖*‖rf‖*‖df‖ + ‖d₀‖*‖rf-r₀‖*‖df‖ +
        ‖d₀‖*‖r₀‖*‖df-d₀‖) := by
  rw [pdars_difference]
  have hsum : ‖1+d₀*r₀*d₀‖ ≤ 1+‖d₀‖*‖r₀‖*‖d₀‖ := by
    exact (norm_add_le _ _).trans (by simpa using add_le_add_left (norm_triple d₀ r₀ d₀) ‖(1:A)‖)
  have hthree : ‖(df-d₀)*rf*df + d₀*(rf-r₀)*df + d₀*r₀*(df-d₀)‖ ≤
      ‖df-d₀‖*‖rf‖*‖df‖ + ‖d₀‖*‖rf-r₀‖*‖df‖ + ‖d₀‖*‖r₀‖*‖df-d₀‖ := by
    calc
      _ ≤ ‖(df-d₀)*rf*df + d₀*(rf-r₀)*df‖ + ‖d₀*r₀*(df-d₀)‖ := norm_add_le _ _
      _ ≤ (‖(df-d₀)*rf*df‖ + ‖d₀*(rf-r₀)*df‖) + ‖d₀*r₀*(df-d₀)‖ := by
        gcongr; exact norm_add_le _ _
      _ ≤ _ := add_le_add (add_le_add (norm_triple _ _ _) (norm_triple _ _ _))
        (norm_triple _ _ _)
  calc
    _ ≤ ‖(af-a₀)*(1+d₀*r₀*d₀)‖ +
        ‖af*((df-d₀)*rf*df + d₀*(rf-r₀)*df + d₀*r₀*(df-d₀))‖ := norm_add_le _ _
    _ ≤ ‖af-a₀‖*‖1+d₀*r₀*d₀‖ +
        ‖af‖*‖(df-d₀)*rf*df + d₀*(rf-r₀)*df + d₀*r₀*(df-d₀)‖ :=
      add_le_add (norm_mul_le _ _) (norm_mul_le _ _)
    _ ≤ _ := by gcongr

/-- Quantitative PD-ARS perturbation bound. The hypotheses are separate
resolvent and transport estimates, not an amplification or stability bound. -/
theorem pdars_perturbation_bound (af a₀ rf r₀ df d₀ : A)
    (K dbase dactual δQ δD : ℝ)
    (hK : 0 ≤ K) (hd0 : 0 ≤ dbase)
    (hQ : 0 ≤ δQ) (hD : 0 ≤ δD)
    (ha : ‖af‖ ≤ K/2) (hr : ‖rf‖ ≤ K) (hr0 : ‖r₀‖ ≤ 1)
    (ha_diff : ‖af-a₀‖ ≤ K*δQ/4) (hr_diff : ‖rf-r₀‖ ≤ K*δQ)
    (hd_base : ‖d₀‖ ≤ dbase) (hd_actual : ‖df‖ ≤ dactual)
    (hd_diff : ‖df-d₀‖ ≤ δD) :
    ‖pdars af rf df - pdars a₀ r₀ d₀‖ ≤
      δQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual) +
      δD*(K/2*(K*dactual+dbase)) := by
  calc
    _ ≤ ‖af-a₀‖*(1+‖d₀‖*‖r₀‖*‖d₀‖) +
      ‖af‖*(‖df-d₀‖*‖rf‖*‖df‖ + ‖d₀‖*‖rf-r₀‖*‖df‖ +
        ‖d₀‖*‖r₀‖*‖df-d₀‖) := pdars_difference_norm _ _ _ _ _ _
    _ ≤ (K*δQ/4)*(1+dbase*1*dbase) +
      (K/2)*(δD*K*dactual+dbase*(K*δQ)*dactual+dbase*1*δD) := by gcongr
    _ = _ := by ring

end PDARS

section Scaling
open scoped BigOperators
variable {ι : Type*} [Fintype ι]

/-- Optical thickness and a macroscopic flux bound jointly control f/dx.
Neither optical thickness nor f<1 alone would give this estimate. -/
theorem flux_over_cell_bound (α f B dx τ : ℝ)
    (hf : 0 ≤ f) (hdx : 0 < dx) (hτ : 0 < τ)
    (hflux : α*f ≤ B) (hcell : τ ≤ α*dx) : f/dx ≤ B/τ := by
  apply (div_le_div_iff₀ hdx hτ).mpr
  calc
    f*τ ≤ f*(α*dx) := mul_le_mul_of_nonneg_left hcell hf
    _ = (α*f)*dx := by ring
    _ ≤ B*dx := mul_le_mul_of_nonneg_right hflux hdx.le

theorem diffusion_transport_bound (α f B τ dt c : ℝ) (dx : ι → ℝ)
    (hf : 0 ≤ f) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hdx : ∀ i, 0 < dx i) (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i) :
    f * (∑ i, dt*c/dx i) ≤ dt*c*B*((Fintype.card ι : ℝ)/τ) := by
  calc
    f * (∑ i, dt*c/dx i) = ∑ i, dt*c*(f/dx i) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ ≤ ∑ i : ι, dt*c*(B/τ) := by
      apply Finset.sum_le_sum
      intro i _
      exact mul_le_mul_of_nonneg_left
        (flux_over_cell_bound α f B (dx i) τ hf (hdx i) hτ hflux (hcell i))
        (mul_nonneg hdt hc)
    _ = dt*c*B*((Fintype.card ι : ℝ)/τ) := by simp; ring

/-- The opacity-independent finite-time conclusion from the two primitive
PD-ARS perturbation estimates. All dependence on opacity and cell widths is
explicit and eliminated in the conclusion. -/
theorem diffusion_finite_time (G G₀ : A) (K c₁ c₂ α f B τ dt c : ℝ)
    (dx : ι → ℝ) (n : ℕ)
    (hK : 0 ≤ K) (hc₁ : 0 ≤ c₁) (hc₂ : 0 ≤ c₂)
    (hf : 0 ≤ f) (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hdx : ∀ i, 0 < dx i) (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (hbase : ∀ j : ℕ, ‖G₀ ^ j‖ ≤ K)
    (hprimitive : ‖G-G₀‖ ≤ c₁*(dt*c*α)*f + c₂*f*(∑ i, dt*c/dx i)) :
    ‖G ^ n‖ ≤ K * Real.exp
      (K*(c₁+c₂*(Fintype.card ι : ℝ)/τ)*c*B*((n : ℝ)*dt)) := by
  have hq : (dt*c*α)*f ≤ dt*c*B := by
    calc
      (dt*c*α)*f = (dt*c)*(α*f) := by ring
      _ ≤ dt*c*B := mul_le_mul_of_nonneg_left hflux (mul_nonneg hdt hc)
  have hw := diffusion_transport_bound α f B τ dt c dx hf hdt hc hτ hdx hflux hcell
  let ε := dt*c*B*(c₁+c₂*(Fintype.card ι : ℝ)/τ)
  have hε : 0 ≤ ε := by dsimp [ε]; positivity
  have hd : ‖G-G₀‖ ≤ ε := by
    calc
      _ ≤ c₁*((dt*c*α)*f)+c₂*(f*(∑ i, dt*c/dx i)) := by
        convert hprimitive using 1 <;> ring
      _ ≤ c₁*(dt*c*B)+c₂*(dt*c*B*((Fintype.card ι : ℝ)/τ)) := by gcongr
      _ = ε := by dsimp [ε]; ring
  have hh := power_perturbation G G₀ K ε hK hε hbase hd n
  convert hh using 1 <;> congr 1 <;> dsimp [ε] <;> ring

end Scaling

section Assembly
variable [NormOneClass A]
open scoped BigOperators
variable {ι : Type*} [Fintype ι]

/-- Complete PD-ARS diffusion-scaling assembly from source/transport data.
The actual implicit-solve identities and primitive coefficient estimates are
explicit; no spectral-radius or full-M1 stability hypothesis is assumed. -/
theorem pdars_diffusion_finite_time
    (Qf Q₀ af a₀ rf r₀ df d₀ : A)
    (K K₀ dbase dactual cQ cW α f B τ dt c : ℝ) (dx : ι → ℝ) (n : ℕ)
    (hK : 0 ≤ K) (hK₀ : 0 ≤ K₀) (hd0 : 0 ≤ dbase) (hdf : 0 ≤ dactual)
    (hcQ : 0 ≤ cQ) (hcW : 0 ≤ cW) (hα : 0 ≤ α) (hf : 0 ≤ f)
    (hB : 0 ≤ B) (hdt : 0 ≤ dt) (hc : 0 ≤ c) (hτ : 0 < τ)
    (hdx : ∀ i, 0 < dx i) (hflux : α*f ≤ B) (hcell : ∀ i, τ ≤ α*dx i)
    (haf : af*(2+Qf)=1) (ha₀ : (2+Q₀)*a₀=1)
    (hrf : rf*(1+Qf)=1) (hr₀ : (1+Q₀)*r₀=1)
    (hna : ‖af‖ ≤ K/2) (hna₀ : ‖a₀‖ ≤ 1/2)
    (hnr : ‖rf‖ ≤ K) (hnr₀ : ‖r₀‖ ≤ 1)
    (hnd : ‖df‖ ≤ dactual) (hnd₀ : ‖d₀‖ ≤ dbase)
    (hQ : ‖Qf-Q₀‖ ≤ cQ*(dt*c*α)*f)
    (hD : ‖df-d₀‖ ≤ cW*f*(∑ i, dt*c/dx i))
    (href : ∀ j : ℕ, ‖(pdars a₀ r₀ d₀)^j‖ ≤ K₀) :
    ‖(pdars af rf df)^n‖ ≤ K₀ * Real.exp
      (K₀*(cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual)+
        (cW*(K/2*(K*dactual+dbase)))*(Fintype.card ι : ℝ)/τ)*
        c*B*((n : ℝ)*dt)) := by
  let δQ := cQ*(dt*c*α)*f
  let δD := cW*f*(∑ i, dt*c/dx i)
  have hδQ : 0 ≤ δQ := by dsimp [δQ]; positivity
  have hδD : 0 ≤ δD := by
    dsimp [δD]
    apply mul_nonneg (mul_nonneg hcW hf)
    apply Finset.sum_nonneg
    intro i _
    exact div_nonneg (mul_nonneg hdt hc) (hdx i).le
  have hQa : ‖(2+Q₀)-(2+Qf)‖ ≤ δQ := by
    simpa only [add_sub_add_left_eq_sub, norm_sub_rev] using hQ
  have hQr : ‖(1+Q₀)-(1+Qf)‖ ≤ δQ := by
    simpa only [add_sub_add_left_eq_sub, norm_sub_rev] using hQ
  have haDiff : ‖af-a₀‖ ≤ K*δQ/4 := by
    convert inverse_difference_bound (2+Qf) (2+Q₀) af a₀ (K/2) (1/2) δQ
      (by positivity) (by norm_num) hδQ haf ha₀ hna hna₀ hQa using 1 <;> ring
  have hrDiff : ‖rf-r₀‖ ≤ K*δQ := by
    simpa using inverse_difference_bound (1+Qf) (1+Q₀) rf r₀ K 1 δQ
      hK (by norm_num) hδQ hrf hr₀ hnr hnr₀ hQr
  have hp := pdars_perturbation_bound af a₀ rf r₀ df d₀ K dbase dactual δQ δD
    hK hd0 hδQ hδD hna hnr hnr₀ haDiff hrDiff hnd₀ hnd hD
  apply diffusion_finite_time (pdars af rf df) (pdars a₀ r₀ d₀) K₀
    (cQ*(K/4*(1+dbase^2)+K^2/2*dbase*dactual))
    (cW*(K/2*(K*dactual+dbase))) α f B τ dt c dx n
    hK₀ (by positivity) (by positivity) hf hB hdt hc hτ hdx hflux hcell href
  convert hp using 1 <;> dsimp [δQ, δD] <;> ring

end Assembly

end RSLA.M1Diffusion
