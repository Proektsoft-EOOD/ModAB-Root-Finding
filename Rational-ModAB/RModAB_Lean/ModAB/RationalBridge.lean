import ModAB.Convergence

/-! Reuse of the original convergence proof for closed-interval proposals.
A degenerate rational proposal may be an existing endpoint. Such a proposal
does not make progress, but it preserves the bracket and consumes a scheduled
step. No convergence or width estimate is assumed as a field below.
This module does not model floating-point evaluation of the new rational formula.
-/
namespace ModAB.RationalBridge
open Set Filter Topology ModAB.Convergence
noncomputable section

theorem updated_invariants_closed {f : ℝ → ℝ} {a b x a' b' : ℝ}
    (hab : a < b) (hx : x ∈ Icc a b) (hs : f a * f b < 0)
    (hz : f x ≠ 0) (hu : Updated f a b x a' b') :
    a' < b' ∧ f a' * f b' < 0 ∧ Icc a' b' ⊆ Icc a b := by
  unfold Updated at hu
  have hi := interval_update hx
  by_cases hax : f a * f x < 0
  · have hlt : a < x := by
      by_contra hn
      have he : x = a := le_antisymm (le_of_not_gt hn) hx.1
      rw [he] at hax
      nlinarith [sq_nonneg (f a)]
    rw [if_pos hax, Prod.mk.injEq] at hu
    rcases hu with ⟨rfl, rfl⟩
    exact ⟨hlt, hax, hi.1⟩
  · have hxb : f x * f b < 0 := by
      rcases sign_update hs hz with h | h
      · exact False.elim (hax h.1)
      · exact h.1
    have hlt : x < b := by
      by_contra hn
      have he : x = b := le_antisymm hx.2 (le_of_not_gt hn)
      rw [he] at hxb
      nlinarith [sq_nonneg (f b)]
    rw [if_neg hax, Prod.mk.injEq] at hu
    rcases hu with ⟨rfl, rfl⟩
    exact ⟨hlt, hxb, hi.2.1⟩

structure ClosedExactRun (B : ℝ) (M : ℕ) (f : ℝ → ℝ) where
  a : ℕ → ℝ
  b : ℕ → ℝ
  point : ℕ → ℝ
  mode : ℕ → Mode
  initial_order : a 0 < b 0
  initial_sign : f (a 0) * f (b 0) < 0
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  initial_mode : mode 0 = .bisection
  candidate_mem : ∀ n, point n ∈ Icc (a n) (b n)
  nonzero : ∀ n, f (point n) ≠ 0
  update : ∀ n, Updated f (a n) (b n) (point n) (a (n+1)) (b (n+1))
  width_step : ∀ n, Step B M (b n-a n) (mode n)
    (b (n+1)-a (n+1)) (mode (n+1))

theorem closed_invariants {B : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (r : ClosedExactRun B M f) :
    ∀ n, r.a n < r.b n ∧ f (r.a n) * f (r.b n) < 0 ∧
      Icc (r.a (n+1)) (r.b (n+1)) ⊆ Icc (r.a n) (r.b n) := by
  have hs : ∀ n, r.a n < r.b n ∧ f (r.a n) * f (r.b n) < 0 := by
    intro n
    induction n with
    | zero => exact ⟨r.initial_order, r.initial_sign⟩
    | succ n ih =>
      have h := updated_invariants_closed ih.1 (r.candidate_mem n) ih.2
        (r.nonzero n) (r.update n)
      exact ⟨h.1, h.2.1⟩
  intro n
  exact ⟨(hs n).1, (hs n).2,
    (updated_invariants_closed (hs n).1 (r.candidate_mem n) (hs n).2
      (r.nonzero n) (r.update n)).2.2⟩

def ClosedExactRun.widthRun {B : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (r : ClosedExactRun B M f) : Run B M where
  width n := r.b n-r.a n
  mode := r.mode
  nonneg_initial := by linarith [r.initial_order]
  initial := r.initial_mode
  next := r.width_step

theorem closed_zeros {B : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (r : ClosedExactRun B M f) :
    ∀ n, ∃ z ∈ Icc (r.a n) (r.b n), f z = 0 := by
  have hi := closed_invariants r
  have hm : Antitone (fun n => Icc (r.a n) (r.b n)) :=
    antitone_nat_of_succ_le (fun n => (hi n).2.2)
  intro n
  exact bracket_contains_zero (hi n).1.le
    (r.continuous.mono (hm (Nat.zero_le n))) (hi n).2.1

theorem capped_global_convergence {B : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (hB : (1/2 : ℝ) ≤ B) (r : ClosedExactRun B M f) :
    Tendsto (fun n => r.b n-r.a n) atTop (𝓝 0) ∧
    ∃ x, (⋂ n, Icc (r.a n) (r.b n)) = {x} ∧ f x = 0 := by
  have hw := width_tendsto_zero hB r.widthRun
  exact ⟨hw, nested_intervals_root (fun n => (closed_invariants r n).2.2)
    (closed_zeros r) r.continuous hw⟩

theorem capped_global_width {B : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (hB : (1/2 : ℝ) ≤ B) (r : ClosedExactRun B M f) (n : ℕ) :
    r.b n-r.a n ≤ (r.b 0-r.a 0) /
      (2 : ℝ)^(n/(dyadicExponent (envelopeConstant B M)+2)) :=
  global_width_bound_log hB r.widthRun.nonneg_initial (history_of_run r.widthRun n)

theorem capped_width_seven {f : ℝ → ℝ}
    (r : ClosedExactRun 2 3 f) (n : ℕ) :
    r.b n-r.a n ≤ (r.b 0-r.a 0) / (2 : ℝ)^(n/7) := by
  exact C32_global_width r.widthRun.nonneg_initial (history_of_run r.widthRun n)

theorem capped_finite_tolerance {B τ : ℝ} {M : ℕ} {f : ℝ → ℝ}
    (hB : (1/2 : ℝ) ≤ B) (hτ : 0 < τ) (r : ClosedExactRun B M f) :
    let q := dyadicExponent (envelopeConstant B M)+2
    let K := toleranceExponent (r.b 0-r.a 0) τ
    r.b (q*K)-r.a (q*K) ≤ τ := by
  dsimp
  exact finite_termination hB r.widthRun.nonneg_initial
    (ceil_log_bound (lt_of_lt_of_le (by norm_num) (C_at_least_one hB)))
    (tolerance_exponent_sufficient r.widthRun.nonneg_initial hτ)
    (history_of_run r.widthRun _)

end
end ModAB.RationalBridge
