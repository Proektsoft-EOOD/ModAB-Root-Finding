import ModAB.RationalCycleIntermediate
import ModAB.RationalCycleSafeguards
import ModAB.RationalCycleUniformAlgebra
import ModAB.RationalLocalComplete

/-! Ordered source-state guards for a possibly scaled local pair state.
The raw local formulas name the latest endpoint first; the implementation
always orders the bracket. Both the denominator and raw value are transported
under endpoint exchange before applying the actual clamp and continuation.
-/
namespace ModAB.Rational.CycleStateGuards
open Set Filter
open scoped Topology
open ModAB.Convergence LocalState CycleIntermediate CycleSafeguards
noncomputable section

structure Geometry (r ε : ℝ) (s : State) : Prop where
  left : s.a < r
  right : r < s.b
  left_history : s.leftHistory < s.a
  right_history : s.b < s.rightHistory
  near_left : |s.a-r| ≤ ε
  near_right : |s.b-r| ≤ ε
  near_left_history : |s.leftHistory-r| ≤ ε
  near_right_history : |s.rightHistory-r| ≤ ε

theorem Geometry.of_valid {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : Geometry r ε s :=
  ⟨h.left,h.right,h.left_history,h.right_history,h.near_left,h.near_right,
    h.near_left_history,h.near_right_history⟩

theorem Geometry.of_same_side {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) : Geometry r ε (update f s x) := by
  obtain ⟨ha,hb,hl,hr,hna,hnb,hnl,hnr⟩ := same_side_geometry h hxi hsame hsgn
  exact ⟨ha,hb,hl,hr,hna,hnb,hnl,hnr⟩

theorem function_ne_zero {f : ℝ → ℝ} {r C ε δ x : ℝ}
    (N : ScaledNeighborhood f r C ε δ) (hx : x ≠ r) (hn : |x-r| ≤ ε) :
    f x ≠ 0 := by
  have hp : 0 < |deriv f r|/2*|x-r| :=
    mul_pos (by positivity [abs_pos.mpr N.deriv_nonzero])
      (abs_pos.mpr (sub_ne_zero.mpr hx))
  have hh := hp.trans_le (N.compare x hn).1
  exact abs_pos.mp hh

/-- Once a true residual passes factor 42, the source continuation decision
cannot fall back, irrespective of its current width threshold. -/
theorem nextMode_of_trial_residual {f : ℝ → ℝ} {s : State} {raw : ℝ}
    (ht : trial f (LocalComplete.allowRational f) s = raw)
    (hres : 42*|f raw| < min |f s.a| |f s.b|) :
    ∀ B, nextMode B f (LocalComplete.allowRational f) s =
      .interpolation s.scale (s.clock+1) := by
  intro B
  unfold nextMode
  rw [ht]
  unfold ModAB.RationalUncapped.nextMode
  have hr : |f raw| < (1/42:ℝ)*min |f s.a| |f s.b| := by nlinarith
  split_ifs <;> simp_all

/-- The actual ordered candidate and continuation agree with the latest-first
local rational formula. All function-value guards needed for exchanging the
endpoints follow from local simple-root residual estimates. -/
theorem guards_latest {f : ℝ → ℝ} {r C ε δ μ : ℝ} {s : State}
    (N : ScaledNeighborhood f r C ε δ) (h : Geometry r ε s)
    (hl : latestWork s = f (latest s))
    (ho : oppositeWork s = μ*f (opposite s))
    (hμ : μ ≠ 0) (hμnear : |μ-1| ≤ δ)
    (hden : denominator (latest s) (opposite s) (f (latest s))
      (μ*f (opposite s)) (selectedHistory s) (f (selectedHistory s)) ≠ 0)
    (hraw : |rawPoint (latest s) (opposite s) (f (latest s))
      (μ*f (opposite s)) (selectedHistory s) (f (selectedHistory s))-r| ≤
      (C*ε^2)*min |latest s-r| |opposite s-r|) :
    let raw := rawPoint (latest s) (opposite s) (f (latest s))
      (μ*f (opposite s)) (selectedHistory s) (f (selectedHistory s))
    trial f (LocalComplete.allowRational f) s = raw ∧
      ∀ B, nextMode B f (LocalComplete.allowRational f) s =
        .interpolation s.scale (s.clock+1) := by
  dsimp only
  have hfa := function_ne_zero N (ne_of_lt h.left) h.near_left
  have hfb := function_ne_zero N (ne_of_gt h.right) h.near_right
  cases hs : s.recent with
  | left =>
    have hden' : denominator s.a s.b (1*f s.a) (μ*f s.b)
        s.leftHistory (f s.leftHistory) ≠ 0 := by
      simpa only [latest,opposite,selectedHistory,hs,one_mul] using hden
    have hraw' : |rawPoint s.a s.b (1*f s.a) (μ*f s.b)
        s.leftHistory (f s.leftHistory)-r| ≤ (C*ε^2)*min |s.a-r| |s.b-r| := by
      simpa only [latest,opposite,selectedHistory,hs,one_mul] using hraw
    have hg := N.guards s.recent h.left h.right h.left_history h.right_history
      h.near_left_history h.near_right_history h.near_left h.near_right
      (by simpa using N.factor_pos.le) hμnear hden' hraw'
    have hl' : s.ordinates.leftWork = f s.a := by simpa [latestWork,latest,hs] using hl
    have ho' : s.ordinates.rightWork = μ*f s.b := by simpa [oppositeWork,opposite,hs] using ho
    have ht : trial f (LocalComplete.allowRational f) s =
        rawPoint s.a s.b (f s.a) (μ*f s.b) s.leftHistory (f s.leftHistory) := by
      simpa only [trial,proposal,LocalComplete.allowRational,selectedHistory,hs,
        hl',ho',one_mul] using hg.2.1
    have hres : 42*|f (rawPoint s.a s.b (f s.a) (μ*f s.b)
        s.leftHistory (f s.leftHistory))| < min |f s.a| |f s.b| := by
      simpa only [one_mul] using hg.2.2
    simpa only [latest,opposite,selectedHistory,hs] using
      And.intro ht (nextMode_of_trial_residual ht hres)
  | right =>
    have hhist := function_ne_zero N
      (ne_of_gt (h.right.trans h.right_history)) h.near_right_history
    have hden' : denominator s.b s.a (f s.b) (μ*f s.a)
        s.rightHistory (f s.rightHistory) ≠ 0 := by
      simpa only [latest,opposite,selectedHistory,hs] using hden
    have he : rawPoint s.b s.a (f s.b) (μ*f s.a) s.rightHistory (f s.rightHistory) =
        rawPoint s.a s.b (μ*f s.a) (f s.b) s.rightHistory (f s.rightHistory) :=
      CycleUniformAlgebra.rawPoint_swap
        (by linarith [h.left,h.right]) (by linarith [h.left,h.right,h.right_history])
        (by linarith [h.right_history]) hfb (mul_ne_zero hμ hfa) hhist hden'
    have hswapden := CycleUniformAlgebra.denominator_swap_ne_zero
      (by linarith [h.left,h.right] : s.a-s.b ≠ 0)
      (by linarith [h.left,h.right,h.right_history] : s.rightHistory-s.a ≠ 0)
      (by linarith [h.right_history] : s.rightHistory-s.b ≠ 0)
      hfb (mul_ne_zero hμ hfa) hhist hden'
    have horderedDen : denominator s.a s.b (μ*f s.a) (1*f s.b)
        s.rightHistory (f s.rightHistory) ≠ 0 := by simpa only [one_mul] using hswapden
    have hraw' : |rawPoint s.a s.b (μ*f s.a) (1*f s.b)
        s.rightHistory (f s.rightHistory)-r| ≤ (C*ε^2)*min |s.a-r| |s.b-r| := by
      simpa only [latest,opposite,selectedHistory,hs,he,one_mul,min_comm] using hraw
    have hg := N.guards s.recent h.left h.right h.left_history h.right_history
      h.near_left_history h.near_right_history h.near_left h.near_right
      hμnear (by simpa using N.factor_pos.le) horderedDen hraw'
    have hl' : s.ordinates.rightWork = f s.b := by simpa [latestWork,latest,hs] using hl
    have ho' : s.ordinates.leftWork = μ*f s.a := by simpa [oppositeWork,opposite,hs] using ho
    have ht : trial f (LocalComplete.allowRational f) s =
        rawPoint s.a s.b (μ*f s.a) (f s.b) s.rightHistory (f s.rightHistory) := by
      simpa only [trial,proposal,LocalComplete.allowRational,selectedHistory,hs,
        hl',ho',one_mul] using hg.2.1
    have hres : 42*|f (rawPoint s.a s.b (μ*f s.a) (f s.b)
        s.rightHistory (f s.rightHistory))| < min |f s.a| |f s.b| := by
      simpa only [one_mul] using hg.2.2
    simpa only [latest,opposite,selectedHistory,hs,he] using
      And.intro ht (nextMode_of_trial_residual ht hres)

end
end ModAB.Rational.CycleStateGuards
