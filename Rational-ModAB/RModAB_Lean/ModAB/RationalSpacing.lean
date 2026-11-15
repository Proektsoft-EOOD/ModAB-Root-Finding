import Mathlib.Analysis.SpecificLimits.Basic

/-!
The closing-step obligation is separate from interpolation accuracy and from
residual convergence.  These results isolate exactly what a positive lower
bound on successive trial-point distances contributes.  They do not assert
that an implementation satisfies that lower bound.
-/

namespace ModAB.RationalUncapped
open Filter Topology

/-- Any convergent real sequence eventually has arbitrarily small successive
differences.  This applies to the actual evaluated points, not an interpolating
model's unaccepted proposals. -/
theorem eventually_short_steps {x : ℕ → ℝ} {z δ : ℝ}
    (hx : Tendsto x atTop (𝓝 z)) (hδ : 0 < δ) :
    ∀ᶠ n in atTop, |x (n + 1) - x n| < δ := by
  have hshift : Tendsto (fun n => x (n + 1)) atTop (𝓝 z) :=
    (tendsto_add_atTop_iff_nat 1).mpr hx
  have hzero : Tendsto (fun n => |x (n + 1) - x n|) atTop (𝓝 (0 : ℝ)) := by
    simpa using (hshift.sub hx).abs
  exact (tendsto_order.mp hzero).2 δ hδ

/-- A positive closing distance on an entire infinite suffix is incompatible
with point convergence.  A finite number of earlier bisections is immaterial. -/
theorem no_convergent_separated_tail {x : ℕ → ℝ} {z δ : ℝ}
    (hx : Tendsto x atTop (𝓝 z)) (hδ : 0 < δ)
    (hsep : ∀ᶠ n in atTop, δ ≤ |x (n + 1) - x n|) : False := by
  obtain ⟨n, hsmall, hlarge⟩ := ((eventually_short_steps hx hδ).and hsep).exists
  exact (not_lt_of_ge hlarge) hsmall

/-- Conditional width termination for a convergent sequence whose accepted
closing steps enforce a fixed positive spacing whenever the width test has
not yet passed.  Proving the spacing premise is an implementation obligation. -/
theorem exists_width_stop_of_converges {x width : ℕ → ℝ} {z δ τ : ℝ}
    (hx : Tendsto x atTop (𝓝 z)) (hδ : 0 < δ)
    (hclose : ∀ᶠ n in atTop, τ < width n → δ ≤ |x (n + 1) - x n|) :
    ∃ n, width n ≤ τ := by
  by_contra h
  push Not at h
  apply no_convergent_separated_tail hx hδ
  filter_upwards [hclose] with n hn
  exact hn (h n)

end ModAB.RationalUncapped
