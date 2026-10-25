import RSLA.M1Source
import RSLA.M1Fourier

/-! Fixed coordinate changes and real/imaginary duplication. These estimates
make explicit that changing from physical block norms to the reference norm
only changes fixed constants, not their opacity dependence. -/
noncomputable section
namespace RSLA.M1Coordinates
variable {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
  [NormedAddCommGroup F] [NormedSpace ℝ F]

/-- The same real operator acts on real and imaginary Fourier amplitudes. -/
def duplicate (A : E →L[ℝ] E) : (E × E) →L[ℝ] (E × E) := A.prodMap A

theorem duplicate_bound (A : E →L[ℝ] E) : ‖duplicate A‖ ≤ ‖A‖ := by
  apply ContinuousLinearMap.opNorm_le_bound _ (norm_nonneg _)
  intro x
  change ‖(A x.1,A x.2)‖ ≤ _
  rw [norm_prod_le_iff]
  exact ⟨(A.le_opNorm x.1).trans (mul_le_mul_of_nonneg_left (norm_fst_le x) (norm_nonneg _)),
    (A.le_opNorm x.2).trans (mul_le_mul_of_nonneg_left (norm_snd_le x) (norm_nonneg _))⟩

theorem duplicate_difference (A B : E →L[ℝ] E) :
    duplicate A-duplicate B=duplicate (A-B) := by ext x <;> rfl

def conjugate (S : E ≃L[ℝ] F) (A : E →L[ℝ] E) : F →L[ℝ] F :=
  S.toContinuousLinearMap.comp (A.comp S.symm.toContinuousLinearMap)

def condition (S : E ≃L[ℝ] F) : ℝ := ‖S.toContinuousLinearMap‖*‖S.symm.toContinuousLinearMap‖

theorem condition_nonneg (S : E ≃L[ℝ] F) : 0 ≤ condition S := by unfold condition; positivity

theorem conjugate_bound (S : E ≃L[ℝ] F) (A : E →L[ℝ] E) :
    ‖conjugate S A‖ ≤ condition S*‖A‖ := by
  calc
    _ ≤ ‖S.toContinuousLinearMap‖*‖A.comp S.symm.toContinuousLinearMap‖ :=
      ContinuousLinearMap.opNorm_comp_le _ _
    _ ≤ ‖S.toContinuousLinearMap‖*(‖A‖*‖S.symm.toContinuousLinearMap‖) := by
      gcongr
      exact ContinuousLinearMap.opNorm_comp_le _ _
    _ = _ := by unfold condition; ring

theorem conjugate_difference (S : E ≃L[ℝ] F) (A B : E →L[ℝ] E) :
    conjugate S A-conjugate S B=conjugate S (A-B) := by
  ext x
  simp [conjugate, map_sub]

theorem conjugate_difference_bound (S : E ≃L[ℝ] F) (A B : E →L[ℝ] E)
    (c f : ℝ) (h : ‖A-B‖ ≤ c*f) :
    ‖conjugate S A-conjugate S B‖ ≤ (condition S*c)*f := by
  rw [conjugate_difference]
  exact (conjugate_bound S (A-B)).trans (by
    have hh := mul_le_mul_of_nonneg_left h (condition_nonneg S)
    nlinarith)

theorem conjugate_resolvent (S : E ≃L[ℝ] F) (Q R : E →L[ℝ] E) (q : ℝ)
    (h : ∀ x, R x+q • Q (R x)=x) :
    ∀ y, conjugate S R y+q • conjugate S Q (conjugate S R y)=y := by
  intro y
  change S (R (S.symm y))+q • S (Q (S.symm (S (R (S.symm y)))))=y
  rw [S.symm_apply_apply, ← map_smul, ← map_add, h, S.apply_symm_apply]

theorem duplicate_add (A B : E →L[ℝ] E) : duplicate (A+B)=duplicate A+duplicate B := by
  ext x <;> rfl

theorem duplicate_mul (A B : E →L[ℝ] E) : duplicate (A*B)=duplicate A*duplicate B := by
  ext x <;> rfl

theorem duplicate_one : duplicate (1 : E →L[ℝ] E)=1 := by ext x <;> rfl

theorem duplicate_smul (a : ℝ) (A : E →L[ℝ] E) : duplicate (a • A)=a • duplicate A := by
  ext x <;> rfl

theorem conjugate_add (S : E ≃L[ℝ] F) (A B : E →L[ℝ] E) :
    conjugate S (A+B)=conjugate S A+conjugate S B := by ext x; simp [conjugate]

theorem conjugate_mul (S : E ≃L[ℝ] F) (A B : E →L[ℝ] E) :
    conjugate S (A*B)=conjugate S A*conjugate S B := by ext x; simp [conjugate]

theorem conjugate_one (S : E ≃L[ℝ] F) : conjugate S (1 : E →L[ℝ] E)=1 := by
  ext x; simp [conjugate]

theorem conjugate_smul (S : E ≃L[ℝ] F) (a : ℝ) (A : E →L[ℝ] E) :
    conjugate S (a • A)=a • conjugate S A := by ext x; simp [conjugate]

/-- Embed a radiation operator with a zero material block. -/
def radiationBlock (A : E →L[ℝ] E) : (E × E) →L[ℝ] (E × E) :=
  (0 : E →L[ℝ] E).prodMap A

theorem radiationBlock_bound (A : E →L[ℝ] E) : ‖radiationBlock A‖ ≤ ‖A‖ := by
  apply ContinuousLinearMap.opNorm_le_bound _ (norm_nonneg _)
  intro x
  change ‖(0,A x.2)‖ ≤ _
  rw [norm_prod_le_iff]
  exact ⟨by simp; positivity,
    (A.le_opNorm x.2).trans (mul_le_mul_of_nonneg_left (norm_snd_le x) (norm_nonneg _))⟩

theorem radiationBlock_difference (A B : E →L[ℝ] E) :
    radiationBlock A-radiationBlock B=radiationBlock (A-B) := by
  apply ContinuousLinearMap.ext
  intro x
  change ((0 : E),A x.2)-((0 : E),B x.2)=((0 : E),(A-B) x.2)
  simp

/-- Multiplication by i in real/imaginary coordinates. -/
def rotation : (E × E) →L[ℝ] (E × E) :=
  (-ContinuousLinearMap.snd ℝ E E).prod (ContinuousLinearMap.fst ℝ E E)

theorem rotation_bound : ‖rotation (E := E)‖ ≤ 1 := by
  apply ContinuousLinearMap.opNorm_le_bound _ zero_le_one
  intro x
  change ‖(-x.2,x.1)‖ ≤ _
  rw [norm_prod_le_iff]
  simpa using And.intro (norm_snd_le x) (norm_fst_le x)

end RSLA.M1Coordinates
