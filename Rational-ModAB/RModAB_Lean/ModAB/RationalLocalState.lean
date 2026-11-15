import ModAB.RationalLocalNeighborhood
import ModAB.RationalLocalSafeguards
import ModAB.RationalOrder

/-!
The two-sided history and true-sign update in the local interpolation phase.
This module does not assume that an orbit remains local or keeps interpolating.
Those conclusions are obtained by induction from pointwise neighbourhood bounds.
The final C³ theorem instantiates the bounds independently of the orbit.
-/
namespace ModAB.Rational.LocalState
open Set Filter
open scoped Topology
open ModAB.Convergence ModAB.RationalLocalSafeguards
open LocalNeighborhood
noncomputable section

structure State where
  a : ℝ
  b : ℝ
  leftHistory : ℝ
  rightHistory : ℝ
  ordinates : OrdinateState
  recent : EndpointSide
  scale : ℝ
  clock : ℕ

def freshOrdinates (f : ℝ → ℝ) (a b : ℝ) (side : EndpointSide) : OrdinateState :=
  ⟨f a, f b, f a, f b, some side⟩

def latest (s : State) : ℝ :=
  match s.recent with | .left => s.a | .right => s.b

def opposite (s : State) : ℝ :=
  match s.recent with | .left => s.b | .right => s.a

def selectedHistory (s : State) : ℝ :=
  match s.recent with | .left => s.leftHistory | .right => s.rightHistory

structure Valid (f : ℝ → ℝ) (r ε : ℝ) (s : State) : Prop where
  left : s.a < r
  right : r < s.b
  left_history : s.leftHistory < s.a
  right_history : s.b < s.rightHistory
  near_left : |s.a-r| ≤ ε
  near_right : |s.b-r| ≤ ε
  near_left_history : |s.leftHistory-r| ≤ ε
  near_right_history : |s.rightHistory-r| ≤ ε
  fresh : s.ordinates = freshOrdinates f s.a s.b s.recent

def proposal (f : ℝ → ℝ) (allow : State → Bool) (s : State) : ℝ :=
  candidate s.a s.b s.ordinates.leftWork s.ordinates.rightWork
    (selectedHistory s) (f (selectedHistory s)) (allow s)

def trial (f : ℝ → ℝ) (allow : State → Bool) (s : State) : ℝ :=
  closeCandidate s.a s.b 0 (proposal f allow s) s.recent

/-- Genuine interpolation update, including both guarded discarded-side
records and the established AB ordinate update. It has no root-sign oracle. -/
def update (f : ℝ → ℝ) (s : State) (x : ℝ) : State :=
  if 0 < f s.a * f x then
    ⟨x, s.b, if x = s.a then s.leftHistory else s.a, s.rightHistory,
      recordUpdate s.ordinates .left (f x), .left, s.scale, s.clock+1⟩
  else
    ⟨s.a, x, s.leftHistory, if x = s.b then s.rightHistory else s.b,
      recordUpdate s.ordinates .right (f x), .right, s.scale, s.clock+1⟩

def step (f : ℝ → ℝ) (allow : State → Bool) (s : State) : State :=
  update f s (trial f allow s)

/-- The literal uncapped continuation decision, with the new bracket width. -/
def nextMode (B : ℝ) (f : ℝ → ℝ) (allow : State → Bool) (s : State) :
    ModAB.RationalUncapped.Mode :=
  let s' := step f allow s
  ModAB.RationalUncapped.nextMode B (1/42) s.scale
    (min |f s.a| |f s.b|) (s'.b-s'.a) |f (trial f allow s)| s.clock

/-- Direct local analytic estimates; no state invariance, recurrence, orbit
convergence or future guard decision is an assumption of this structure. -/
structure Neighborhood (f : ℝ → ℝ) (r ε : ℝ) (allow : State → Bool) where
  bound : ℝ
  secantFactor : ℝ
  lower : ℝ
  upper : ℝ
  radius_pos : 0 < ε
  bound_nonneg : 0 ≤ bound
  upper_nonneg : 0 ≤ upper
  contraction : bound*ε^2 < 1
  clamp : bound*ε^2 + secantFactor ≤ (7/10:ℝ)*(1-secantFactor)
  residual : 42*upper*(bound*ε^2) < lower
  raw : ∀ a b d, Admissible r a b d → tripleRadius r a b d ≤ ε →
    denominator a b (f a) (f b) d (f d) ≠ 0 ∧
    0 < rawCoefficient f r a b d ∧ |rawCoefficient f r a b d| ≤ bound
  secant : ∀ a b, a < r → r < b → |a-r| ≤ ε → |b-r| ≤ ε →
    |ModAB.Rational.secant a b (f a) (f b)-r| ≤
      secantFactor * min |a-r| |b-r|
  compare : ∀ x, |x-r| ≤ ε →
    lower*|x-r| ≤ |f x| ∧ |f x| ≤ upper*|x-r|
  same_sign : ∀ x y, |x-r| ≤ ε → |y-r| ≤ ε →
    (0 < f x*f y ↔ 0 < (x-r)*(y-r))
  curvature : ∀ s, Valid f r ε s → allow s = true

theorem Valid.admissible {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : Admissible r s.a s.b (selectedHistory s) := by
  cases hs : s.recent <;> simp only [selectedHistory, hs, Admissible]
  · exact ⟨by linarith [h.left,h.right], by linarith [h.left_history,h.left,h.right],
      by linarith [h.left_history], by linarith [h.left], by linarith [h.right],
      by linarith [h.left_history,h.left]⟩
  · exact ⟨by linarith [h.left,h.right], by linarith [h.right_history],
      by linarith [h.left,h.right,h.right_history], by linarith [h.left],
      by linarith [h.right], by linarith [h.right_history,h.right]⟩

theorem Valid.radius {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : tripleRadius r s.a s.b (selectedHistory s) ≤ ε := by
  simp only [tripleRadius, max_le_iff]
  refine ⟨h.near_left,h.near_right,?_⟩
  cases hs : s.recent <;> simp only [selectedHistory, hs]
  · exact h.near_left_history
  · exact h.near_right_history

theorem latest_update (f : ℝ → ℝ) (s : State) (x : ℝ) :
    latest (update f s x) = x := by
  unfold update
  split_ifs <;> rfl

theorem interior_near {a b x r ε : ℝ} (_ha : a < r) (_hb : r < b)
    (hna : |a-r| ≤ ε) (hnb : |b-r| ≤ ε) (hx : x ∈ Ioo a b) : |x-r| ≤ ε := by
  rcases abs_le.mp hna with ⟨hal,har⟩
  rcases abs_le.mp hnb with ⟨hbl,hbr⟩
  exact abs_le.mpr ⟨by linarith [hx.1], by linarith [hx.2]⟩

theorem fresh_update_right {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) (hs : s.recent = .left) {x : ℝ}
    (hx : r < x) (hxb : x < s.b) (hnear : |x-r| ≤ ε)
    (hsgn : ¬ 0 < f s.a*f x) :
    Valid f r ε (update f s x) ∧
    opposite (update f s x) = latest s ∧
    selectedHistory (update f s x) = opposite s := by
  have hxe : x ≠ s.b := ne_of_lt hxb
  have hu : update f s x =
      ⟨s.a,x,s.leftHistory,s.b,freshOrdinates f s.a x .right,
        .right,s.scale,s.clock+1⟩ := by
    simp [update,hsgn,hxe,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨h.left,hx,h.left_history,hxb,h.near_left,hnear,
    h.near_left_history,h.near_right,rfl⟩,?_,?_⟩ <;>
      simp [opposite,latest,selectedHistory,hs]

theorem fresh_update_left {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) (hs : s.recent = .right) {x : ℝ}
    (hx : x < r) (hax : s.a < x) (hnear : |x-r| ≤ ε)
    (hsgn : 0 < f s.a*f x) :
    Valid f r ε (update f s x) ∧
    opposite (update f s x) = latest s ∧
    selectedHistory (update f s x) = opposite s := by
  have hxe : x ≠ s.a := ne_of_gt hax
  have hu : update f s x =
      ⟨x,s.b,s.a,s.rightHistory,freshOrdinates f x s.b .left,
        .left,s.scale,s.clock+1⟩ := by
    simp [update,hsgn,hxe,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨hx,h.right,hax,h.right_history,hnear,h.near_right,
    h.near_left,h.near_right_history,rfl⟩,?_,?_⟩ <;>
      simp [opposite,latest,selectedHistory,hs]

/-- A single admissible local state forces the complete accepted crossing:
both history records remain correct, the AB ordinates remain true, the actual
clamp and closing are inactive, and the actual continuation test passes. -/
theorem step_preserves {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (B : ℝ) :
    Valid f r ε (step f allow s) ∧
    proposal f allow s = rawPoint s.a s.b (f s.a) (f s.b)
      (selectedHistory s) (f (selectedHistory s)) ∧
    trial f allow s = proposal f allow s ∧
    nextMode B f allow s = .interpolation s.scale (s.clock+1) ∧
    opposite (step f allow s) = latest s ∧
    selectedHistory (step f allow s) = opposite s ∧
    |latest (step f allow s)-r| ≤ (N.bound*ε^2)*|latest s-r| := by
  let d := selectedHistory s
  let x := rawPoint s.a s.b (f s.a) (f s.b) d (f d)
  let H := rawCoefficient f r s.a s.b d
  have ha := h.admissible
  have hn := N.raw s.a s.b d ha h.radius
  have he : x-r = H*(s.a-r)*(s.b-r)*(d-r) :=
    raw_error_eq_coefficient ha.2.2.2.1 ha.2.2.2.2.1 ha.2.2.2.2.2
  have hb : |x-r| ≤ N.bound*|s.a-r| *|s.b-r| *|d-r| := by
    rw [he,abs_mul,abs_mul,abs_mul]
    gcongr
    exact hn.2.2
  have hd : |d-r| ≤ ε := (le_max_right _ _).trans
    ((le_max_right _ _).trans h.radius)
  have hemin := triple_error_bound N.bound_nonneg N.radius_pos.le
    h.near_left h.near_right hd hb
  have hxi := interior_of_local_error h.left h.right N.contraction hemin
  have hxnear := interior_near h.left h.right h.near_left h.near_right hxi
  have hflower : N.lower*min |s.a-r| |s.b-r| ≤ min |f s.a| |f s.b| := by
    have hl := (N.compare s.a h.near_left).1
    have hr := (N.compare s.b h.near_right).1
    have hlo : 0 ≤ N.lower := by
      have hq : 0 ≤ N.bound*ε^2 := mul_nonneg N.bound_nonneg (sq_nonneg _)
      nlinarith [N.residual,mul_nonneg N.upper_nonneg hq]
    rw [mul_min_of_nonneg _ _ hlo]
    exact min_le_min hl hr
  have hop := local_operational_bridge (r:=r) (q:=N.bound*ε^2)
    (t:=N.secantFactor) (lower:=N.lower) (upper:=N.upper)
    (B:=B) (S:=s.scale) (w:=(step f allow s).b-(step f allow s).a)
    s.clock s.recent h.left h.right N.contraction hn.1 hemin
    (N.secant s.a s.b h.left h.right h.near_left h.near_right)
    N.clamp N.upper_nonneg (N.compare x hxnear).2 hflower N.residual
  have hp : proposal f allow s = x := by
    simpa only [proposal,h.fresh,freshOrdinates,N.curvature s h] using hop.1
  have ht : trial f allow s = x := by
    simpa only [trial,proposal,h.fresh,freshOrdinates,N.curvature s h] using hop.2.1
  have hvalid : Valid f r ε (step f allow s) ∧
      opposite (step f allow s)=latest s ∧
      selectedHistory (step f allow s)=opposite s := by
    rw [step,ht]
    cases hs : s.recent with
    | left =>
      have hdr : d-r < 0 := by
        dsimp [d,selectedHistory]; rw [hs]; linarith [h.left_history,h.left]
      have hp0 : 0 < H*(s.a-r)*(s.b-r)*(d-r) := by
        exact mul_pos_of_neg_of_neg
          (mul_neg_of_neg_of_pos (mul_neg_of_pos_of_neg hn.2.1 (sub_neg.mpr h.left))
            (sub_pos.mpr h.right)) hdr
      have hxr : r < x := by linarith [he]
      have hsgn : ¬ 0 < f s.a*f x := by
        rw [N.same_sign s.a x h.near_left hxnear]
        exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg (by linarith [h.left])
          (by linarith))
      exact fresh_update_right h hs hxr hxi.2 hxnear hsgn
    | right =>
      have hdr : 0 < d-r := by
        dsimp [d,selectedHistory]; rw [hs]; linarith [h.right_history,h.right]
      have hp0 : H*(s.a-r)*(s.b-r)*(d-r) < 0 := by
        exact mul_neg_of_neg_of_pos
          (mul_neg_of_neg_of_pos (mul_neg_of_pos_of_neg hn.2.1 (sub_neg.mpr h.left))
            (sub_pos.mpr h.right)) hdr
      have hxr : x < r := by linarith [he]
      have hsgn : 0 < f s.a*f x := by
        rw [N.same_sign s.a x h.near_left hxnear]
        exact mul_pos_of_neg_of_neg (sub_neg.mpr h.left) (sub_neg.mpr hxr)
      exact fresh_update_left h hs hxr hxi.1 hxnear hsgn
  refine ⟨hvalid.1,hp,ht.trans hp.symm,?_,hvalid.2.1,hvalid.2.2,?_⟩
  · dsimp only [nextMode]
    rw [ht]
    simpa only [hop.1] using hop.2.2
  · rw [step,latest_update,ht]
    apply hemin.trans
    apply mul_le_mul_of_nonneg_left _ (mul_nonneg N.bound_nonneg (sq_nonneg _))
    cases hs : s.recent <;> simp only [latest,hs]
    · exact min_le_left _ _
    · exact min_le_right _ _

end
end ModAB.Rational.LocalState
