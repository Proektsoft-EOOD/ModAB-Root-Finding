import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Stability of the logarithmic two-step recurrence

This module supplies the analytic step in the tribonacci Q-order proof.  It uses
an explicit positive quadratic energy and a perturbed scalar contraction; no
asymptotic conclusion is assumed as a hypothesis.
-/

namespace ModAB.RationalOrderStability

open Filter
open scoped Topology

/-- A nonnegative sequence subject to a strict scalar contraction and a forcing
term tending to zero itself tends to zero. -/
theorem tendsto_zero_of_contraction {u d : ℕ → ℝ} {q : ℝ}
    (hq0 : 0 ≤ q) (hq1 : q < 1) (hu : ∀ n, 0 ≤ u n)
    (hstep : ∀ n, u (n + 1) ≤ q * u n + d n)
    (hd : Tendsto d atTop (𝓝 0)) : Tendsto u atTop (𝓝 0) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  have hqpos : 0 < 1 - q := sub_pos.mpr hq1
  have hδ : 0 < ε * (1 - q) / 2 := by positivity
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp hd (ε * (1 - q) / 2) hδ
  have hforce : ∀ n, N ≤ n → d n ≤ ε * (1 - q) / 2 := by
    intro n hn
    have h := hN n hn
    rw [Real.dist_eq, sub_zero] at h
    exact (le_abs_self _).trans h.le
  have hbound : ∀ k, u (N + k) ≤ q ^ k * u N + ε / 2 := by
    intro k
    induction k with
    | zero => simp only [Nat.add_zero, pow_zero, one_mul]; linarith
    | succ k ih =>
      have hs := hstep (N + k)
      have hf := hforce (N + k) (by omega)
      have hi := mul_le_mul_of_nonneg_left ih hq0
      rw [← Nat.add_assoc, pow_succ]
      nlinarith
  have hpow : Tendsto (fun k : ℕ => q ^ k * u N) atTop (𝓝 0) := by
    simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one hq0 hq1).mul_const (u N)
  obtain ⟨L, hL⟩ := Metric.tendsto_atTop.mp hpow (ε / 2) (by positivity)
  refine ⟨N + L, fun n hn => ?_⟩
  have hnl : L ≤ n - N := by omega
  have hg := hL (n - N) hnl
  rw [Real.dist_eq, sub_zero] at hg
  have hb := hbound (n - N)
  have hNn : N + (n - N) = n := by omega
  rw [hNn] at hb
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (hu n)]
  have := (le_abs_self (q ^ (n - N) * u N)).trans_lt hg
  linarith

/-- The positive quadratic energy of the companion recurrence. -/
def energy (a b x y : ℝ) : ℝ := x ^ 2 + a * x * y + b * y ^ 2

/-- Uniform coercivity on the coefficient rectangle needed below. -/
theorem energy_coercive {a b x y : ℝ} (ha0 : 0 ≤ a)
    (ha1 : a ≤ 84 / 100) (hb0 : 54 / 100 ≤ b) :
    (x ^ 2 + y ^ 2) / 10 ≤ energy a b x y := by
  have h1 : 0 ≤ a / 2 * (x + y) ^ 2 := mul_nonneg (by positivity) (sq_nonneg _)
  have h2 : 0 ≤ (9 / 10 - a / 2) * x ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  have h3 : 0 ≤ (b - 1 / 10 - a / 2) * y ^ 2 :=
    mul_nonneg (by linarith) (sq_nonneg _)
  unfold energy
  nlinarith

/-- The forcing term is controlled by the same energy. -/
theorem forcing_sq_le_energy {a b x y : ℝ} (ha0 : 0 ≤ a)
    (ha1 : a ≤ 84 / 100) (hb0 : 54 / 100 ≤ b) (hb1 : b ≤ 55 / 100) :
    (-a * x - 2 * b * y) ^ 2 ≤ 80 * energy a b x y := by
  have he := energy_coercive (x := x) (y := y) ha0 ha1 hb0
  have ha : a ^ 2 ≤ 1 := by nlinarith
  have hb : b ^ 2 ≤ 1 := by nlinarith
  have hx := mul_le_mul_of_nonneg_right ha (sq_nonneg x)
  have hy := mul_le_mul_of_nonneg_right hb (sq_nonneg y)
  have hs := sq_nonneg (a * x - 2 * b * y)
  nlinarith [sq_nonneg x, sq_nonneg y]

/-- Exact energy balance for one forced companion update. -/
theorem energy_update (a b x y e : ℝ) :
    energy a b (-a * x - b * y + e) x =
      b * energy a b x y + e * (-a * x - 2 * b * y) + e ^ 2 := by
  unfold energy
  ring

/-- A uniform contraction estimate, with an explicit forcing bound. -/
theorem energy_contracts {a b x y e : ℝ} (ha0 : 0 ≤ a)
    (ha1 : a ≤ 84 / 100) (hb0 : 54 / 100 ≤ b) (hb1 : b ≤ 55 / 100) :
    energy a b (-a * x - b * y + e) x ≤
      (4 / 5) * energy a b x y + 81 * e ^ 2 := by
  have he := energy_coercive (x := x) (y := y) ha0 ha1 hb0
  have he0 : 0 ≤ energy a b x y := by nlinarith [sq_nonneg x, sq_nonneg y]
  have hf := forcing_sq_le_energy (x := x) (y := y) ha0 ha1 hb0 hb1
  have hy := sq_nonneg ((-a * x - 2 * b * y) - 160 * e)
  have hb := mul_le_mul_of_nonneg_right hb1 he0
  rw [energy_update]
  nlinarith

/-- Stability of the second-order recurrence on a rectangle containing the
tribonacci quadratic factor coefficients. -/
theorem tendsto_zero_of_second_order {a b : ℝ} {v eps : ℕ → ℝ}
    (ha0 : 0 ≤ a) (ha1 : a ≤ 84 / 100)
    (hb0 : 54 / 100 ≤ b) (hb1 : b ≤ 55 / 100)
    (hrec : ∀ n, v (n + 2) + a * v (n + 1) + b * v n = eps n)
    (heps : Tendsto eps atTop (𝓝 0)) : Tendsto v atTop (𝓝 0) := by
  let E : ℕ → ℝ := fun n => energy a b (v (n + 1)) (v n)
  have hcoerce : ∀ n, ((v (n + 1)) ^ 2 + (v n) ^ 2) / 10 ≤ E n := by
    intro n
    exact energy_coercive ha0 ha1 hb0
  have hnonneg : ∀ n, 0 ≤ E n := by
    intro n
    have := hcoerce n
    nlinarith [sq_nonneg (v (n + 1)), sq_nonneg (v n)]
  have hstep : ∀ n, E (n + 1) ≤ (4 / 5) * E n + 81 * (eps n) ^ 2 := by
    intro n
    have hr : v (n + 2) = -a * v (n + 1) - b * v n + eps n := by
      linarith [hrec n]
    dsimp [E]
    rw [show n + 1 + 1 = n + 2 by omega, hr]
    exact energy_contracts ha0 ha1 hb0 hb1
  have hforcing : Tendsto (fun n => 81 * (eps n) ^ 2) atTop (𝓝 0) := by
    simpa using (heps.pow 2).const_mul 81
  have hE : Tendsto E atTop (𝓝 0) :=
    tendsto_zero_of_contraction (by norm_num) (by norm_num) hnonneg hstep hforcing
  apply Metric.tendsto_atTop.mpr
  intro δ hδ
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp hE (δ ^ 2 / 10) (by positivity)
  refine ⟨N, fun n hn => ?_⟩
  have h := hN n hn
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (hnonneg n)] at h
  have hc := hcoerce n
  rw [Real.dist_eq, sub_zero]
  have hsq : (v n) ^ 2 < δ ^ 2 := by nlinarith [sq_nonneg (v (n + 1))]
  exact (sq_lt_sq₀ (abs_nonneg (v n)) hδ.le).mp (by simpa using hsq)

/-- The stable logarithmic recurrence for the tribonacci root.  The cubic
identity is needed upstream to derive the recurrence, not for this stability
lemma. -/
theorem tendsto_zero_of_recurrence {p : ℝ} {v eps : ℕ → ℝ}
    (hp_lower : 1839 / 1000 < p) (hp_upper : p < 46 / 25)
    (hrec : ∀ n, v (n + 2) + (p - 1) * v (n + 1) + (1 / p) * v n = eps n)
    (heps : Tendsto eps atTop (𝓝 0)) : Tendsto v atTop (𝓝 0) := by
  have hp : 0 < p := by linarith
  apply tendsto_zero_of_second_order (a := p - 1) (b := 1 / p)
      (by linarith) (by linarith) _ _ hrec heps
  · apply (le_div_iff₀ hp).mpr
    nlinarith
  · apply (div_le_iff₀ hp).mpr
    nlinarith

end ModAB.RationalOrderStability
