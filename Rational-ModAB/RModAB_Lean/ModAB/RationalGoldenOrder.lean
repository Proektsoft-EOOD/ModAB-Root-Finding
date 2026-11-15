import ModAB.RationalOrder
import Mathlib.NumberTheory.Real.GoldenRatio

/-!
The golden-ratio Q-order of a two-error product recurrence. This is conditional
on the recurrence; it does not assume or prove entry into a one-sided solver
regime, nor convergence of its accumulated Anderson--Björck multiplier.
-/
namespace ModAB.RationalGoldenOrder
open Filter Topology
noncomputable section

/-- Stability of a scalar recurrence with a vanishing additive perturbation. -/
theorem tendsto_zero_of_first_order {q : ℝ} {v ε : ℕ → ℝ}
    (hq0 : 0 ≤ q) (hq1 : q < 1)
    (hrec : ∀ n, v (n + 1) + q * v n = ε n)
    (hε : Tendsto ε atTop (𝓝 0)) : Tendsto v atTop (𝓝 0) := by
  have hstep : ∀ n, |v (n + 1)| ≤ q * |v n| + |ε n| := by
    intro n
    have hv : v (n + 1) = ε n - q * v n := by linarith [hrec n]
    rw [hv]
    calc
      |ε n - q * v n| ≤ |ε n| + |q * v n| := by
        simpa using abs_sub_le (ε n) 0 (q * v n)
      _ = q * |v n| + |ε n| := by rw [abs_mul, abs_of_nonneg hq0]; ring
  have ha : Tendsto (fun n => |v n|) atTop (𝓝 0) :=
    ModAB.RationalOrderStability.tendsto_zero_of_contraction hq0 hq1
      (fun n => abs_nonneg _) hstep (by simpa using hε.abs)
  exact (tendsto_zero_iff_abs_tendsto_zero v).mpr ha

/-- Logarithmic form of a nonvanishing two-error product recurrence. -/
theorem log_abs_recurrence {e C : ℕ → ℝ} (he : ∀ n, e n ≠ 0)
    (hrec : ∀ n, e (n + 2) = C n * e (n + 1) * e n) (n : ℕ) :
    Real.log |e (n + 2)| = Real.log |C n| + Real.log |e (n + 1)| + Real.log |e n| := by
  have hC : C n ≠ 0 := by
    intro hz
    have hz' : e (n + 2) = 0 := by rw [hrec n, hz]; ring
    exact he (n + 2) hz'
  simp only [Real.log_abs]
  rw [hrec n]
  simp [Real.log_mul, hC, he]

/-- The exact logarithmic difference limit for the golden-ratio order. -/
theorem two_error_recurrence_log_limit {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0) (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 2) = C n * e (n + 1) * e n) :
    Tendsto (fun n => Real.log |e (n + 1)| - Real.goldenRatio * Real.log |e n|)
      atTop (𝓝 ((Real.goldenRatio - 1) * Real.log |K|)) := by
  let u : ℕ → ℝ := fun n => Real.log |e n| + Real.log |K|
  let ε : ℕ → ℝ := fun n => Real.log |C n| - Real.log |K|
  let v : ℕ → ℝ := fun n => u (n + 1) - Real.goldenRatio * u n
  have hε : Tendsto ε atTop (𝓝 0) := by
    have hh := (hC.abs.log (abs_ne_zero.mpr hK)).sub_const (Real.log |K|)
    simpa only [sub_self] using hh
  have hu : ∀ n, u (n + 2) = u (n + 1) + u n + ε n := by
    intro n
    dsimp [u, ε]
    linear_combination log_abs_recurrence he hrec n
  have hp : 0 < Real.goldenRatio := Real.goldenRatio_pos
  have hp0 : Real.goldenRatio ≠ 0 := ne_of_gt hp
  have hcoeff : Real.goldenRatio - 1 / Real.goldenRatio = 1 := by
    have hi : 1 / Real.goldenRatio = Real.goldenRatio - 1 :=
      (div_eq_iff hp0).mpr (by nlinarith [Real.goldenRatio_sq])
    linarith
  have hinv : (1 / Real.goldenRatio) * Real.goldenRatio = 1 := by field_simp
  have hvrec : ∀ n, v (n + 1) + (1 / Real.goldenRatio) * v n = ε n := by
    intro n
    dsimp [v]
    rw [show n + 1 + 1 = n + 2 by omega, hu n]
    linear_combination -(u (n + 1)) * hcoeff - u n * hinv
  have hv : Tendsto v atTop (𝓝 0) := tendsto_zero_of_first_order
    (by positivity) ((div_lt_one hp).mpr Real.one_lt_goldenRatio) hvrec hε
  have heq : (fun n => Real.log |e (n + 1)| - Real.goldenRatio * Real.log |e n|) =
      fun n => v n + (Real.goldenRatio - 1) * Real.log |K| := by
    funext n
    dsimp [v, u]
    ring
  rw [heq]
  simpa only [zero_add] using hv.add_const ((Real.goldenRatio - 1) * Real.log |K|)

/-- The genuine power-quotient limit, including the asymptotic constant. -/
theorem two_error_recurrence_q_ratio {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0) (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 2) = C n * e (n + 1) * e n) :
    Tendsto (fun n => |e (n + 1)| / |e n| ^ Real.goldenRatio)
      atTop (𝓝 (|K| ^ (Real.goldenRatio - 1))) := by
  have hlog := two_error_recurrence_log_limit he hK hC hrec
  have hh := Real.continuous_exp.continuousAt.tendsto.comp hlog
  have heq : ∀ n, |e (n + 1)| / |e n| ^ Real.goldenRatio =
      Real.exp (Real.log |e (n + 1)| - Real.goldenRatio * Real.log |e n|) := by
    intro n
    rw [Real.exp_sub, Real.exp_log (abs_pos.mpr (he (n + 1))),
      Real.rpow_def_of_pos (abs_pos.mpr (he n))]
    congr 2
    ring
  change Tendsto (fun n => Real.exp (Real.log |e (n + 1)| -
      Real.goldenRatio * Real.log |e n|)) atTop
    (𝓝 (Real.exp ((Real.goldenRatio - 1) * Real.log |K|))) at hh
  simp_rw [← heq] at hh
  simpa only [Real.rpow_def_of_pos (abs_pos.mpr hK), mul_comm] using hh

/-- Q-order of a convergent two-error recurrence, with no assumed quotient limit. -/
theorem two_error_recurrence_hasQOrder {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0) (hlim : Tendsto e atTop (𝓝 0))
    (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 2) = C n * e (n + 1) * e n) :
    ModAB.RationalOrder.HasQOrder e Real.goldenRatio (|K| ^ (Real.goldenRatio - 1)) := by
  exact ⟨Real.one_lt_goldenRatio, Real.rpow_pos_of_pos (abs_pos.mpr hK) _, hlim,
    two_error_recurrence_q_ratio he hK hC hrec⟩

/-- An arbitrary finite prefix is permitted in the one-sided recurrence theorem. -/
theorem eventual_two_error_recurrence_hasQOrder {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0) (hlim : Tendsto e atTop (𝓝 0))
    (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ᶠ n in atTop, e (n + 2) = C n * e (n + 1) * e n) :
    ModAB.RationalOrder.HasQOrder e Real.goldenRatio (|K| ^ (Real.goldenRatio - 1)) := by
  obtain ⟨N, hN⟩ := eventually_atTop.mp hrec
  have hshift : Tendsto (fun n => C (n + N)) atTop (𝓝 K) :=
    (tendsto_add_atTop_iff_nat N).mpr hC
  have hr : ∀ n, e (n + 2 + N) = C (n + N) * e (n + 1 + N) * e (n + N) := by
    intro n
    simpa only [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using
      hN (n + N) (Nat.le_add_left N n)
  have hq := two_error_recurrence_q_ratio (e := fun n => e (n + N))
    (C := fun n => C (n + N)) (fun n => he (n + N)) hK hshift hr
  refine ⟨Real.one_lt_goldenRatio, Real.rpow_pos_of_pos (abs_pos.mpr hK) _, hlim, ?_⟩
  apply (tendsto_add_atTop_iff_nat N).mp
  simpa only [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hq

end
end ModAB.RationalGoldenOrder
