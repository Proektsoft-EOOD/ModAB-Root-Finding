import ModAB.RationalAlternating
import ModAB.RationalSecantLocal
import ModAB.RationalCurvature

/-!
The complete C³-to-operational-trajectory theorem for the positive-K regime.
All uniform estimates and the actual curvature decision are derived here.
The initial hypotheses concern only the given local state, never its future
trajectory, a recurrence, a coefficient limit or an asymptotic order.
-/
namespace ModAB.Rational.LocalComplete
open Set Filter
open scoped Topology
open LocalState LocalNeighborhood SecantLocal
noncomputable section

/-- The rational-branch selector uses the true residuals of both endpoint and
discarded-side records, exactly as in the inspected implementation. -/
def allowRational (f : ℝ → ℝ) (s : State) : Bool :=
  !(Curvature.curvatureChangesSign f s.leftHistory s.a s.b s.rightHistory)

/-- Finite C³ regularity supplies every analytic premise of the operational
neighbourhood, including the actual curvature test and factor-42 decision. -/
theorem exists_operational_neighborhood {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : 0 < AnalyticBridge.errorConstant f r) :
    ∃ ε > 0, Nonempty (Neighborhood f r ε (allowRational f)) := by
  have hdiff := hf.differentiableAt (by norm_num)
  obtain ⟨ε₀,hε₀,B,hB,hraw⟩ := exists_positive_raw_bounds hf hz hd hK
  obtain ⟨ε₁,hε₁,hsec⟩ := exists_radius_secant_error hdiff hz hd (t:=(1/10:ℝ)) (by norm_num)
  obtain ⟨ε₂,hε₂,hcompare⟩ := exists_radius_simple_root_bounds hdiff hz hd
  obtain ⟨ε₃,hε₃,hcurv⟩ := Curvature.exists_radius_no_curvature_change
    (hf.of_le (by norm_num)) h₂
  let E := min ε₀ (min ε₁ (min ε₂ (min ε₃ (min (1:ℝ) (1/(1000*(B+1)))))))
  have hE : 0 < E := by dsimp [E]; positivity
  let ε := E/2
  have hε : 0 < ε := by dsimp [ε]; positivity
  have hεE : ε < E := by dsimp [ε]; linarith
  have e₀ : ε < ε₀ := hεE.trans_le (by simp [E,min_le_iff])
  have e₁ : ε < ε₁ := hεE.trans_le (by simp [E,min_le_iff])
  have e₂ : ε < ε₂ := hεE.trans_le (by simp [E,min_le_iff])
  have e₃ : ε < ε₃ := hεE.trans_le (by simp [E,min_le_iff])
  have eone : ε ≤ 1 := (hεE.trans_le (show E ≤ 1 by simp [E,min_le_iff])).le
  have ebound : ε ≤ 1/(1000*(B+1)) :=
    (hεE.trans_le (show E ≤ 1/(1000*(B+1)) by simp [E,min_le_iff])).le
  have heB : ε*(1000*(B+1)) ≤ 1 := (le_div_iff₀ (by positivity)).mp ebound
  have hq : B*ε^2 < (1/1000:ℝ) := by
    have hs := mul_nonneg (mul_nonneg hB.le hε.le) (sub_nonneg.mpr eone)
    nlinarith
  have hder : 0 < |deriv f r| := abs_pos.mpr hd
  refine ⟨ε,hε,⟨{
    bound := B
    secantFactor := 1/10
    lower := |deriv f r|/2
    upper := 2*|deriv f r|
    radius_pos := hε
    bound_nonneg := hB.le
    upper_nonneg := by positivity
    contraction := by linarith
    clamp := by linarith
    residual := ?_
    raw := ?_
    secant := ?_
    compare := ?_
    same_sign := ?_
    curvature := ?_
  }⟩⟩
  · have hs := mul_lt_mul_of_pos_left hq hder
    nlinarith
  · intro a b d hp he
    have hh := hraw a b d hp (he.trans e₀.le)
    exact hh.2.2.2
  · intro a b ha hb hea heb
    exact hsec a b ha hb (hea.trans_lt e₁) (heb.trans_lt e₁)
  · intro x hx
    exact (hcompare x (hx.trans_lt e₂)).2
  · intro x y hx hy
    exact (sign_product_of_normalizedSlope hz hd
      (hcompare x (hx.trans_lt e₂)).1 (hcompare y (hy.trans_lt e₂)).1).1
  · intro s hs
    have hl := (abs_le.mp hs.near_left_history).1
    have hh := (abs_le.mp hs.near_right_history).2
    have hc := hcurv s.leftHistory s.a s.b s.rightHistory
      (by linarith) hs.left_history (hs.left.trans hs.right) hs.right_history (by linarith)
    simp only [allowRational,hc,Bool.not_false]

/-- Theorem 9: a single radius chosen from f makes every valid fresh local
state preserve the actual rational branch, clamp, closing and continuation
tests, and gives the tribonacci Q-order to the genuine state iteration. -/
theorem local_alternating_complete {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : 0 < AnalyticBridge.errorConstant f r) :
    ∃ ε > 0, ∀ s : State, Valid f r ε s →
      (∀ n, Valid f r ε (Alternating.orbit f (allowRational f) s n)) ∧
      (∀ B n, nextMode B f (allowRational f) (Alternating.orbit f (allowRational f) s n) =
        .interpolation (Alternating.orbit f (allowRational f) s n).scale
          ((Alternating.orbit f (allowRational f) s n).clock+1)) ∧
      ModAB.RationalOrder.HasQOrder
        (fun n => latest (Alternating.orbit f (allowRational f) s n)-r)
        ModAB.RationalOrder.tau
        (|AnalyticBridge.errorConstant f r| ^ ((ModAB.RationalOrder.tau-1)/2)) := by
  obtain ⟨ε,hε,⟨N⟩⟩ := exists_operational_neighborhood hf hz hd h₂ hK
  refine ⟨ε,hε,fun s hs => ?_⟩
  exact ⟨Alternating.valid_orbit N hs,
    fun B n => Alternating.continue_orbit N hs B n,
    Alternating.hasQOrder_of_neighborhood N hs hf hz hd hK.ne'⟩

end
end ModAB.Rational.LocalComplete
