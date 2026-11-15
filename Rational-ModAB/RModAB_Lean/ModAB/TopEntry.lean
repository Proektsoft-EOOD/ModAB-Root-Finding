import ModAB.ExamplesBase

/-! Exact-arithmetic bridge for top switching. The switching point is still a
midpoint. Recording its moved side changes the next correction, not the bracket
or the phase width policy. This is not a verification of the C# loop. -/
namespace ModAB.Convergence
noncomputable section

inductive EndpointSide where
  | left | right
  deriving DecidableEq

/-- `GetABFactor`, with its positive fallback. -/
def abFactor (new old : ℝ) : ℝ :=
  if 0 < 1-new/old then 1-new/old else 1/2

theorem abFactor_positive (new old : ℝ) : 0 < abFactor new old := by
  unfold abFactor
  split_ifs with h
  · exact h
  · norm_num

structure OrdinateState where
  leftTrue : ℝ
  rightTrue : ℝ
  leftWork : ℝ
  rightWork : ℝ
  previous : Option EndpointSide

/-- The side is chosen using true function signs outside this ordinate update. -/
def recordUpdate (s : OrdinateState) (moved : EndpointSide) (y : ℝ) : OrdinateState :=
  match moved with
  | .left => ⟨y,s.rightTrue,y,
      (if s.previous = some .left then abFactor y s.leftWork else 1)*s.rightWork,
      some .left⟩
  | .right => ⟨s.leftTrue,y,
      (if s.previous = some .right then abFactor y s.rightWork else 1)*s.leftWork,y,
      some .right⟩

theorem recorded_side (s : OrdinateState) (moved : EndpointSide) (y : ℝ) :
    (recordUpdate s moved y).previous = some moved := by cases moved <;> rfl

/-- Opposite-side first secant: no scaling; the actual new side is recorded. -/
theorem opposite_side_no_correction (a b u v y : ℝ) :
    recordUpdate ⟨a,b,u,v,some .left⟩ .right y = ⟨a,y,u,y,some .right⟩ ∧
    recordUpdate ⟨a,b,u,v,some .right⟩ .left y = ⟨y,b,y,v,some .left⟩ := by
  simp [recordUpdate]

/-- Same-side first secant: exactly the usual positive AB factor is applied. -/
theorem same_side_correction (a b u v y : ℝ) :
    (recordUpdate ⟨a,b,u,v,some .left⟩ .left y).rightWork = abFactor y u*v ∧
    (recordUpdate ⟨a,b,u,v,some .right⟩ .right y).leftWork = abFactor y v*u := by
  simp [recordUpdate]

theorem recordUpdate_working_signs (s : OrdinateState) (moved : EndpointSide) (y : ℝ)
    (hu : SameSign s.leftWork s.leftTrue) (hv : SameSign s.rightWork s.rightTrue)
    (hy : y ≠ 0) :
    SameSign (recordUpdate s moved y).leftWork (recordUpdate s moved y).leftTrue ∧
    SameSign (recordUpdate s moved y).rightWork (recordUpdate s moved y).rightTrue := by
  have hp (old : ℝ) (side : EndpointSide) :
      0 < (if s.previous = some side then abFactor y old else 1) := by
    split_ifs
    · exact abFactor_positive y old
    · norm_num
  cases moved with
  | left => exact ⟨inserted_sign hy,positive_correction hv (hp _ _)⟩
  | right => exact ⟨positive_correction hu (hp _ _),inserted_sign hy⟩

/-- Before the switching bisection, the previous side is cleared and the
working values are reset to the true endpoint residuals. -/
theorem top_entry_ordinates (a b y : ℝ) :
    recordUpdate ⟨a,b,a,b,none⟩ .left y = ⟨y,b,y,b,some .left⟩ ∧
    recordUpdate ⟨a,b,a,b,none⟩ .right y = ⟨a,y,a,y,some .right⟩ := by
  simp [recordUpdate]

/-- At entry the width test succeeds, independently of residuals. Halving
its stored threshold gives exactly the threshold of the first true secant. -/
theorem top_entry_width_policy {B W : ℝ} (hB : (1/2:ℝ) ≤ B) (hW : 0 ≤ W) (M : ℕ) :
    W/2 ≤ B*W ∧ (B*W)/2 = threshold B (W/2) 1 ∧
    Step B M W .bisection (W/2) (.ab (W/2) 0 0) := by
  refine ⟨?_,?_,Step.bisectStart W⟩
  · nlinarith [mul_nonneg (sub_nonneg.mpr hB) hW]
  · simp [threshold]; ring

/-- Concrete common entry bracket for the examples starting from [-3,1]. -/
theorem top_entry_midpoint :
    ((-3:ℝ)+1)/2 = -1 ∧ (1:ℝ)-(-1)=2 := by norm_num

/-- Concrete left-side top initialization used by the example orbit models. -/
theorem top_left_seed (a v : ℝ) :
    recordUpdate ⟨a,v,a,v,none⟩ .left (-1) = ⟨-1,v,-1,v,some .left⟩ :=
  (top_entry_ordinates a v (-1)).1

end
end ModAB.Convergence

namespace ModAB.Examples
noncomputable section
/-- Projection of the recorded top-entry state onto the one-sided orbit model. -/
def topLeftSeed (a v : ℝ) : LeftState :=
  let s := ModAB.Convergence.recordUpdate ⟨a,v,a,v,none⟩ .left (-1)
  ⟨1,s.rightWork,decide (s.previous = some .left)⟩

theorem topLeftSeed_eq (a v : ℝ) : topLeftSeed a v = ⟨1,v,true⟩ := by
  simp [topLeftSeed,ModAB.Convergence.recordUpdate]
end
end ModAB.Examples
