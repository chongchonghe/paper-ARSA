import RSLA.AllModes
import Mathlib.Analysis.InnerProductSpace.PiL2

/-! The phase-independent weighted length in Proposition B.1 is a genuine norm:
it is the pullback of the Euclidean norm by an injective linear embedding. -/
noncomputable section
namespace RSLA.WeightedNorm
open Complex Channels
open scoped BigOperators
variable {κ : Type*} [Fintype κ] [DecidableEq κ]

def embedding (q K : κ → ℝ) : State κ →ₗ[ℂ] EuclideanSpace ℂ (κ × Fin 3) where
  toFun x := (WithLp.equiv 2 _).symm (fun p =>
    ![(Real.sqrt (2*q p.1/K p.1) : ℂ)*x.1 p.1,
      (Real.sqrt 2 : ℂ)*x.2 p.1,
      (Real.sqrt (q p.1) : ℂ)*(x.1 p.1-x.2 p.1)] p.2)
  map_add' x y := by
    ext ⟨k,j⟩
    fin_cases j <;> simp <;> ring
  map_smul' g x := by
    ext ⟨k,j⟩
    fin_cases j <;> simp [smul_eq_mul] <;> ring

lemma embedding_norm_sq (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (x : State κ) : ‖embedding q K x‖^2 = totalEnergy q K x := by
  rw [PiLp.norm_sq_eq_of_L2]
  simp_rw [← Complex.normSq_eq_norm_sq]
  change (∑ p : κ × Fin 3,
    normSq (![ (Real.sqrt (2*q p.1/K p.1) : ℂ)*x.1 p.1,
      (Real.sqrt 2 : ℂ)*x.2 p.1,
      (Real.sqrt (q p.1) : ℂ)*(x.1 p.1-x.2 p.1)] p.2)) = _
  rw [Fintype.sum_prod_type]
  simp only [Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.vecHead, Matrix.vecTail]
  unfold totalEnergy baseEnergy bEnergy
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro k _
  change normSq ((Real.sqrt (2*q k/K k) : ℂ)*x.1 k) +
    normSq ((Real.sqrt 2 : ℂ)*x.2 k) +
    normSq ((Real.sqrt (q k) : ℂ)*(x.1 k-x.2 k)) =
    channelBase (q k) (K k) (x.1 k) (x.2 k) + channelB (q k) (K k) (x.1 k) (x.2 k)
  have hqk := hq k
  have hKk := hK k
  rw [StaticTransport.normSq_real_mul, StaticTransport.normSq_real_mul,
    StaticTransport.normSq_real_mul, Real.sq_sqrt (by positivity : 0 ≤ 2*q k/K k),
    Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 2), Real.sq_sqrt (hq k).le]
  unfold channelB channelBase
  ring

lemma embedding_norm (q K : κ → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (x : State κ) : ‖embedding q K x‖ = Real.sqrt (totalEnergy q K x) := by
  rw [← embedding_norm_sq q K hq hK, Real.sqrt_sq (norm_nonneg _)]

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Positive definiteness of the common weighted norm. -/
theorem eq_zero_iff (q K : Option ι → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (x : State (Option ι)) : AllModes.weightedNorm q K x=0 ↔ x=0 := by
  classical
  constructor
  · intro h
    by_contra hx
    have hp := Real.sqrt_pos.2 (totalEnergy_pos q K hq hK hx)
    change Real.sqrt (totalEnergy q K x)=0 at h
    linarith
  · intro h
    subst x
    simp [AllModes.weightedNorm, totalEnergy, baseEnergy, bEnergy, channelB, channelBase]

/-- Absolute homogeneity over complex scalars. -/
theorem smul (q K : Option ι → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (g : ℂ) (x : State (Option ι)) :
    AllModes.weightedNorm q K (g • x)=‖g‖*AllModes.weightedNorm q K x := by
  unfold AllModes.weightedNorm
  rw [← embedding_norm q K hq hK, ← embedding_norm q K hq hK, map_smul, norm_smul]

/-- Triangle inequality for the norm used in Proposition B.1. -/
theorem triangle (q K : Option ι → ℝ) (hq : ∀ k, 0 < q k) (hK : ∀ k, 0 < K k)
    (x y : State (Option ι)) :
    AllModes.weightedNorm q K (x+y) ≤ AllModes.weightedNorm q K x+AllModes.weightedNorm q K y := by
  unfold AllModes.weightedNorm
  rw [← embedding_norm q K hq hK, ← embedding_norm q K hq hK, ← embedding_norm q K hq hK, map_add]
  exact norm_add_le _ _

end RSLA.WeightedNorm
