import RSLA.AllModes
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv

/-! The converse in Proposition B.1: embedding an unstable checkerboard channel. -/
noncomputable section
namespace RSLA.Checkerboard
open Complex ComplexConjugate Matrix Channels AllModes
open scoped BigOperators

/-- The checkerboard block in gas coordinates scaled by four relative to (61). -/
def block (q K a : ℝ) : Mat2 :=
  (2*(1+q+K)*(2+q+K))⁻¹ •
    !![-K*a*q+2*K*q+2*K+2*q^2+6*q+4,
       K*(K*a^2-2*K*a+2*K+a^2-a*q-4*a+2*q+4);
       q*(-K*a+2*K-2*a+2*q+4),
       K^2*a^2-2*K^2*a+2*K^2+3*K*a^2-K*a*q-6*K*a+2*K*q+6*K+2*a^2-4*a+2*q+4]

lemma block_original {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) :
    EnergyBlock.G q K a =
      !![block q K a 0 0, block q K a 0 1 / 4;
         4*block q K a 1 0, block q K a 1 1] := by
  rw [EnergyBlock.G_entries]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [block, EnergyBlock.D, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, Matrix.smul_apply, Matrix.add_apply, Matrix.one_apply] <;>
    field_simp <;> ring

lemma block_stable_iff {q K a : ℝ} (hq : 0 < q) (hK : 0 < K) (ha : 0 ≤ a) :
    MatrixStable (block q K a) ↔ a ≤ EnergyBlock.boundary q K := by
  have hm : MatrixStable (EnergyBlock.G q K a) ↔ MatrixStable (block q K a) := by
    rw [block_original hq hK]
    unfold MatrixStable
    simp only [characteristic_two]
    have ht : trace2 (!![block q K a 0 0, block q K a 0 1 / 4;
        4*block q K a 1 0, block q K a 1 1] : Mat2) = trace2 (block q K a) := by rfl
    have hd : (!![block q K a 0 0, block q K a 0 1 / 4;
        4*block q K a 1 0, block q K a 1 1] : Mat2).det = (block q K a).det := by
      simp [Matrix.det_fin_two]; ring
    rw [ht, hd]
  exact hm.symm.trans (EnergyBlock.stability_iff hq hK ha)

set_option maxHeartbeats 1000000 in
/-- Two explicit stages and the final solve agree with this matrix, on complex inputs. -/
theorem local_matrix_action {q K : ℝ} (hq : 0 < q) (hK : 0 < K) (a : ℝ) (u v : ℂ) :
    let u₁ := gasSolve q K u (v-(a : ℂ)*v)
    let v₁ := radSolve q K u (v-(a : ℂ)*v)
    let u₂ := gasSolve q K u₁ (v₁-(a : ℂ)*v₁)
    let v₂ := radSolve q K u₁ (v₁-(a : ℂ)*v₁)
    finalGas q K u v u₂ v₂ = (block q K a 0 0 : ℂ)*u + (block q K a 0 1 : ℂ)*v ∧
    finalRad q K u v u₂ v₂ = (block q K a 1 0 : ℂ)*u + (block q K a 1 1 : ℂ)*v := by
  dsimp only
  constructor <;> apply Complex.ext <;>
    simp only [finalGas_re, finalGas_im, finalRad_re, finalRad_im,
      gasSolve_re, gasSolve_im, radSolve_re, radSolve_im,
      Complex.add_re, Complex.add_im, Complex.sub_re, Complex.sub_im,
      Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
      mul_zero, zero_mul, add_zero, sub_zero] <;>
    simp [block, Matrix.smul_apply] <;> field_simp <;> ring

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

lemma full_single (q K : Option ι → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (a : ℝ) (k : Option ι) (u v : ℂ) :
    fullStep q K (fun v j => (a : ℂ)*v j) (singleState k u v) =
      singleState k ((block (q k) (K k) a 0 0 : ℂ)*u + (block (q k) (K k) a 0 1 : ℂ)*v)
        ((block (q k) (K k) a 1 0 : ℂ)*u + (block (q k) (K k) a 1 1 : ℂ)*v) := by
  have hh := local_matrix_action (hq k) (hK k) a u v
  unfold fullStep
  rw [stage_single, stage_single]
  apply Prod.ext <;> funext j <;> by_cases hj : j=k
  · subst j; simpa [singleState] using hh.1
  · simp [singleState, Pi.single_apply, hj, Ne.symm hj, finalGas]
  · subst j; simpa [singleState] using hh.2
  · simp [singleState, Pi.single_apply, hj, Ne.symm hj, finalRad]

/-- Proposition B.1: exact all-mode complex spectral stability, both directions. -/
theorem all_modes_stability_iff (q K : Option ι → ℝ)
    (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (w : ι → ℝ) (hw : ∀ i, 0 ≤ w i) :
    (∀ θ, SpectrallyStable (fullStep q K (StaticTransport.W w θ))) ↔
      ∀ k, 2*∑ i, w i ≤ EnergyBlock.boundary (q k) (K k) := by
  constructor
  · intro hs k
    let a : ℝ := 2*∑ i, w i
    have ha : 0 ≤ a := mul_nonneg (by norm_num) (Finset.sum_nonneg (fun i _ => hw i))
    apply (block_stable_iff (hq k) (hK k) ha).mp
    intro g hg
    obtain ⟨v, hv, hz⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hg
    have he : (block (q k) (K k) a).map (algebraMap ℝ ℂ) *ᵥ v = g • v := by
      have hh : g • v - (block (q k) (K k) a).map (algebraMap ℝ ℂ) *ᵥ v = 0 := by
        simpa only [Matrix.sub_mulVec, Matrix.smul_mulVec_assoc, Matrix.one_mulVec] using hz
      exact (sub_eq_zero.mp hh).symm
    have he₀ := congrFun he 0
    have he₁ := congrFun he 1
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.map_apply] at he₀ he₁
    let x := singleState k (v 0) (v 1)
    have hx : x ≠ 0 := by
      intro hx
      apply hv
      have hx₀ := congrArg (fun z : State (Option ι) => z.1 k) hx
      have hx₁ := congrArg (fun z : State (Option ι) => z.2 k) hx
      simp [x, singleState] at hx₀ hx₁
      ext j
      fin_cases j <;> assumption
    have hwc : StaticTransport.W w (fun _ => Real.pi) = fun v j => (a : ℂ)*v j := by
      funext v j
      exact congrFun (W_checkerboard w v) j
    apply hs (fun _ => Real.pi) g x hx
    rw [hwc]
    change fullStep q K (fun v j => (a : ℂ)*v j) (singleState k (v 0) (v 1)) = _
    rw [full_single q K hq hK, he₀, he₁]
    apply Prod.ext <;> funext j <;> by_cases hj : j=k <;>
      simp_all [x, singleState, Pi.single_apply]
  · intro hb θ
    exact all_modes_spectrally_stable q K hq hK w hw hb θ

end RSLA.Checkerboard
