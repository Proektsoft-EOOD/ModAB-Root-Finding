import ModAB.RationalLocalState
import ModAB.RationalCycleState

/-!
Two successive genuine source updates in the negative-coefficient cycle.
The first update applies the actual Anderson--Björck factor to the retained
ordinate. The crossing update restores true ordinates and the correct two
discarded-side histories; no state-invariance premise is used.
-/
namespace ModAB.Rational.CycleStateUpdate
open Set Filter
open scoped Topology
open ModAB.Convergence LocalState
noncomputable section

theorem pair_update_left {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .left)
    (hax : s.a < x) (hxr : x < r) (hry : r < y) (hyb : y < s.b)
    (hfirst : 0 < f s.a * f x) (hsecond : ¬ 0 < f x * f y) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxe : x ≠ s.a := ne_of_gt hax
  have hye : y ≠ s.b := ne_of_lt hyb
  have hxnear := interior_near h.left h.right h.near_left h.near_right
    (show x ∈ Ioo s.a s.b from ⟨hax,by linarith [h.right]⟩)
  have hynear := interior_near h.left h.right h.near_left h.near_right
    (show y ∈ Ioo s.a s.b from ⟨by linarith [h.left],hyb⟩)
  have hu : update f (update f s x) y =
      ⟨x,y,s.a,s.b,freshOrdinates f x y .right,.right,s.scale,s.clock+2⟩ := by
    simp [update,hfirst,hsecond,hxe,hye,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨hxr,hry,hax,hyb,hxnear,hynear,h.near_left,h.near_right,rfl⟩,
    rfl,rfl,?_⟩
  simp [opposite,hs,selectedHistory]

theorem pair_update_right {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .right)
    (hrx : r < x) (hxb : x < s.b) (hay : s.a < y) (hyr : y < r)
    (hfirst : ¬ 0 < f s.a * f x) (hsecond : 0 < f s.a * f y) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxe : x ≠ s.b := ne_of_lt hxb
  have hye : y ≠ s.a := ne_of_gt hay
  have hxnear := interior_near h.left h.right h.near_left h.near_right
    (show x ∈ Ioo s.a s.b from ⟨by linarith [h.left],hxb⟩)
  have hynear := interior_near h.left h.right h.near_left h.near_right
    (show y ∈ Ioo s.a s.b from ⟨hay,by linarith [h.right]⟩)
  have hu : update f (update f s x) y =
      ⟨y,x,s.a,s.b,freshOrdinates f y x .left,.left,s.scale,s.clock+2⟩ := by
    simp [update,hfirst,hsecond,hxe,hye,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨hyr,hrx,hay,hxb,hynear,hxnear,h.near_left,h.near_right,rfl⟩,
    rfl,rfl,?_⟩
  simp [opposite,hs,selectedHistory]

/-- The pair update needs only strict interiority and the two proved step
signs. True signs select the source branches. Both history records are
updated by the actual guarded assignments. -/
theorem pair_update_preserves {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s)
    (hxi : x ∈ Ioo s.a s.b) (hyi : y ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hcross : (latest s-r)*(y-r) < 0)
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  have hyn := interior_near h.left h.right h.near_left h.near_right hyi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame hcross
    have hxr : x < r := by nlinarith [h.left]
    have hry : r < y := by nlinarith [h.left]
    apply pair_update_left h hs hxi.1 hxr hry hyi.2
    · exact (hsgn _ _ h.near_left hxn).mpr hsame
    · rw [hsgn _ _ hxn hyn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg (by linarith) (by linarith))
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame hcross
    have hrx : r < x := by nlinarith [h.right]
    have hyr : y < r := by nlinarith [h.right]
    apply pair_update_right h hs hrx hxi.2 hyi.1 hyr
    · rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg (by linarith [h.left]) (by linarith))
    · exact (hsgn _ _ h.near_left hyn).mpr
        (mul_pos_of_neg_of_neg (sub_neg.mpr h.left) (sub_neg.mpr hyr))

end
end ModAB.Rational.CycleStateUpdate
