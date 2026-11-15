import ModAB.RationalLocalState

/-!
The intermediate same-side update in a rational pair. All formulas follow
from the literal source ordinate/history update; the AB multiplier and the
selected discarded point are not assumed as a recurrence interface.
-/
namespace ModAB.Rational.CycleIntermediate
open Set Filter
open scoped Topology
open ModAB.Convergence LocalState
noncomputable section

def latestWork (s : State) : ℝ :=
  match s.recent with | .left => s.ordinates.leftWork | .right => s.ordinates.rightWork

def oppositeWork (s : State) : ℝ :=
  match s.recent with | .left => s.ordinates.rightWork | .right => s.ordinates.leftWork

theorem same_side_update_left {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .left) (hx : x ≠ s.a)
    (hbranch : 0 < f s.a * f x) :
    update f s x =
      ⟨x,s.b,s.a,s.rightHistory,
        ⟨f x,f s.b,f x,abFactor (f x) (f s.a)*f s.b,some .left⟩,
        .left,s.scale,s.clock+1⟩ := by
  simp [update,hbranch,hx,h.fresh,hs,freshOrdinates,recordUpdate]

theorem same_side_update_right {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .right) (hx : x ≠ s.b)
    (hbranch : ¬ 0 < f s.a * f x) :
    update f s x =
      ⟨s.a,x,s.leftHistory,s.b,
        ⟨f s.a,f x,abFactor (f x) (f s.b)*f s.a,f x,some .right⟩,
        .right,s.scale,s.clock+1⟩ := by
  simp [update,hbranch,hx,h.fresh,hs,freshOrdinates,recordUpdate]

/-- Orientation-independent source update. In particular, the model's next
third node is exactly the previous latest endpoint with its true ordinate. -/
theorem same_side_update {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    (update f s x).recent = s.recent ∧
    latest (update f s x) = x ∧
    opposite (update f s x) = opposite s ∧
    selectedHistory (update f s x) = latest s ∧
    latestWork (update f s x) = f x ∧
    oppositeWork (update f s x) = abFactor (f x) (f (latest s))*f (opposite s) ∧
    (update f s x).ordinates.previous = some s.recent ∧
    (update f s x).scale = s.scale ∧
    (update f s x).clock = s.clock+1 := by
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame
    have hbranch := (hsgn _ _ h.near_left hxn).mpr hsame
    rw [same_side_update_left h hs (ne_of_gt hxi.1) hbranch]
    simp [latest,opposite,selectedHistory,latestWork,oppositeWork,hs]
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame
    have hrx : r < x := by nlinarith [h.right]
    have hbranch : ¬ 0 < f s.a*f x := by
      rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg
        (by linarith [h.left]) (by linarith))
    rw [same_side_update_right h hs (ne_of_lt hxi.2) hbranch]
    simp [latest,opposite,selectedHistory,latestWork,oppositeWork,hs]

theorem abFactor_eq_one_sub {new old : ℝ} (h : new/old < 1) :
    abFactor new old = 1-new/old := by
  simp [abFactor,show 0 < 1-new/old by linarith]

/-- For the local pair the executed AB branch is the untruncated factor
1-f(new)/f(old), rather than its fallback value 1/2. -/
theorem same_side_working_ordinates {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r)))
    (hratio : f x/f (latest s) < 1) :
    latestWork (update f s x) = f x ∧
    oppositeWork (update f s x) = (1-f x/f (latest s))*f (opposite s) ∧
    selectedHistory (update f s x) = latest s := by
  have hu := same_side_update h hxi hsame hsgn
  refine ⟨hu.2.2.2.2.1,?_,hu.2.2.2.1⟩
  rw [hu.2.2.2.2.2.1,abFactor_eq_one_sub hratio]


/-- The same-side intermediate state keeps the strict bracket and both
outside histories in the local neighbourhood. Freshness is deliberately
not included: its retained working ordinate has just been rescaled. -/
theorem same_side_geometry {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    let t := update f s x
    t.a < r ∧ r < t.b ∧ t.leftHistory < t.a ∧ t.b < t.rightHistory ∧
    |t.a-r| ≤ ε ∧ |t.b-r| ≤ ε ∧
    |t.leftHistory-r| ≤ ε ∧ |t.rightHistory-r| ≤ ε := by
  dsimp only
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame
    have hxr : x < r := by nlinarith [h.left]
    have hbranch := (hsgn _ _ h.near_left hxn).mpr hsame
    rw [same_side_update_left h hs (ne_of_gt hxi.1) hbranch]
    exact ⟨hxr,h.right,hxi.1,h.right_history,hxn,h.near_right,
      h.near_left,h.near_right_history⟩
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame
    have hrx : r < x := by nlinarith [h.right]
    have hbranch : ¬ 0 < f s.a*f x := by
      rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg
        (by linarith [h.left]) (by linarith))
    rw [same_side_update_right h hs (ne_of_lt hxi.2) hbranch]
    exact ⟨h.left,hrx,h.left_history,hxi.2,h.near_left,hxn,
      h.near_left_history,h.near_right⟩

end
end ModAB.Rational.CycleIntermediate
