import ModAB.RationalOrderStability
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.FieldSimp

/-!
The asymptotic order of a three-error product recurrence.
This module concerns the recurrence itself.  A separate local-entry argument
must establish that the safeguarded rational solver satisfies this recurrence.
-/
namespace ModAB.RationalOrder
open Filter Topology Set
noncomputable section

def orderPolynomial (x : ℝ) : ℝ := x ^ 3 - x ^ 2 - x - 1

/-- The order polynomial has a zero in a certified rational interval. -/
theorem exists_order_root : ∃ p : ℝ,
    1839 / 1000 < p ∧ p < 46 / 25 ∧ p ^ 3 = p ^ 2 + p + 1 := by
  have hc : Continuous (fun x : ℝ => orderPolynomial x) := by
    unfold orderPolynomial
    fun_prop
  have hlo : orderPolynomial (1839 / 1000) < 0 := by norm_num [orderPolynomial]
  have hhi : 0 < orderPolynomial (46 / 25) := by norm_num [orderPolynomial]
  obtain ⟨p, hp, hp0⟩ := intermediate_value_Icc (by norm_num : (1839 / 1000 : ℝ) ≤ 46 / 25)
    hc.continuousOn ⟨le_of_lt hlo, le_of_lt hhi⟩
  refine ⟨p, ?_, ?_, ?_⟩
  · exact lt_of_le_of_ne hp.1 (by intro he; subst p; linarith)
  · exact lt_of_le_of_ne hp.2 (by intro he; subst p; linarith)
  · dsimp [orderPolynomial] at hp0
    linarith

/-- There is at most one real root of the order polynomial greater than one. -/
theorem order_root_unique {p q : ℝ} (hp : 1 < p) (hq : 1 < q)
    (hep : p ^ 3 = p ^ 2 + p + 1) (heq : q ^ 3 = q ^ 2 + q + 1) : p = q := by
  have hpp : 0 < p * (p - 1) := mul_pos (by linarith) (by linarith)
  have hqq : 0 < q * (q - 1) := mul_pos (by linarith) (by linarith)
  have hprod : 0 < (p - 1) * (q - 1) := mul_pos (by linarith) (by linarith)
  have hpq : 1 < p * q := by nlinarith [hprod]
  have hfactor : 0 < q ^ 2 + p * q + p ^ 2 - q - p - 1 := by nlinarith
  have hz : (q - p) * (q ^ 2 + p * q + p ^ 2 - q - p - 1) = 0 := by
    nlinarith [hep, heq]
  rcases mul_eq_zero.mp hz with h | h
  · linarith
  · linarith

/-- The tribonacci order, chosen from the proved rational isolating interval. -/
def tau : ℝ := Classical.choose exists_order_root

theorem tau_lower : 1839 / 1000 < tau := (Classical.choose_spec exists_order_root).1
theorem tau_upper : tau < 46 / 25 := (Classical.choose_spec exists_order_root).2.1
theorem tau_cubic : tau ^ 3 = tau ^ 2 + tau + 1 :=
  (Classical.choose_spec exists_order_root).2.2
theorem one_lt_tau : 1 < tau := by linarith [tau_lower]
theorem tau_pos : 0 < tau := by linarith [one_lt_tau]

/-- The other factor of the characteristic polynomial is quadratic. -/
theorem characteristic_factorization (x : ℝ) :
    x ^ 3 - x ^ 2 - x - 1 =
      (x - tau) * (x ^ 2 + (tau - 1) * x + 1 / tau) := by
  have ht : tau ≠ 0 := ne_of_gt tau_pos
  field_simp
  linear_combination x * tau_cubic

/-- A nonvanishing three-error recurrence becomes an additive logarithmic recurrence. -/
theorem log_abs_recurrence {e C : ℕ → ℝ}
    (he : ∀ n, e n ≠ 0)
    (hrec : ∀ n, e (n + 3) = C n * e (n + 2) * e (n + 1) * e n) (n : ℕ) :
    Real.log |e (n + 3)| = Real.log |C n| + Real.log |e (n + 2)| +
      Real.log |e (n + 1)| + Real.log |e n| := by
  have hC : C n ≠ 0 := by
    intro hz
    have hz' : e (n + 3) = 0 := by rw [hrec n, hz]; ring
    exact he (n + 3) hz'
  simp only [Real.log_abs]
  rw [hrec n]
  simp [Real.log_mul, hC, he]

/-- Subtracting the growing characteristic mode produces a stable second-order recurrence. -/
theorem transformed_recurrence {u ε : ℕ → ℝ}
    (hu : ∀ n, u (n + 3) = u (n + 2) + u (n + 1) + u n + ε n) (n : ℕ) :
    (u (n + 3) - tau * u (n + 2)) +
      (tau - 1) * (u (n + 2) - tau * u (n + 1)) +
      (1 / tau) * (u (n + 1) - tau * u n) = ε n := by
  have ht : tau ≠ 0 := ne_of_gt tau_pos
  have hc : tau * (tau - 1) - 1 / tau = 1 := by
    field_simp
    nlinarith [tau_cubic]
  have hi : (1 / tau) * tau = 1 := by field_simp
  rw [hu n]
  linear_combination -(u (n + 1)) * hc - u n * hi

/-- The error quotient is the exponential of the first logarithmic difference. -/
theorem ratio_eq_exp_log_difference {e : ℕ → ℝ} (he : ∀ n, e n ≠ 0) (n : ℕ) :
    |e (n + 1)| / |e n| ^ tau =
      Real.exp (Real.log |e (n + 1)| - tau * Real.log |e n|) := by
  rw [Real.exp_sub, Real.exp_log (abs_pos.mpr (he (n + 1))),
    Real.rpow_def_of_pos (abs_pos.mpr (he n))]
  congr 2
  ring

/-- A logarithmic difference limit gives the genuine power-quotient limit. -/
theorem q_ratio_of_log_difference {e : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0)
    (hlog : Tendsto (fun n => Real.log |e (n + 1)| - tau * Real.log |e n|)
      atTop (𝓝 (((tau - 1) / 2) * Real.log |K|))) :
    Tendsto (fun n => |e (n + 1)| / |e n| ^ tau)
      atTop (𝓝 (|K| ^ ((tau - 1) / 2))) := by
  have hh := Real.continuous_exp.continuousAt.tendsto.comp hlog
  simpa only [Function.comp_def, ← ratio_eq_exp_log_difference he,
    Real.rpow_def_of_pos (abs_pos.mpr hK), mul_comm] using hh

/-- Full logarithmic limit, derived from the product recurrence rather than assumed. -/
theorem three_error_recurrence_log_limit {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0)
    (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 3) = C n * e (n + 2) * e (n + 1) * e n) :
    Tendsto (fun n => Real.log |e (n + 1)| - tau * Real.log |e n|)
      atTop (𝓝 (((tau - 1) / 2) * Real.log |K|)) := by
  let u : ℕ → ℝ := fun n => Real.log |e n| + Real.log |K| / 2
  let ε : ℕ → ℝ := fun n => Real.log |C n| - Real.log |K|
  let v : ℕ → ℝ := fun n => u (n + 1) - tau * u n
  have hε : Tendsto ε atTop (𝓝 0) := by
    have hh := (hC.abs.log (abs_ne_zero.mpr hK)).sub_const (Real.log |K|)
    simpa only [sub_self] using hh
  have hu : ∀ n, u (n + 3) = u (n + 2) + u (n + 1) + u n + ε n := by
    intro n
    dsimp [u, ε]
    linear_combination log_abs_recurrence he hrec n
  have hvrec : ∀ n, v (n + 2) + (tau - 1) * v (n + 1) + (1 / tau) * v n = ε n := by
    intro n
    simpa only [v, Nat.add_assoc, Nat.reduceAdd] using transformed_recurrence hu n
  have hv : Tendsto v atTop (𝓝 0) :=
    ModAB.RationalOrderStability.tendsto_zero_of_recurrence tau_lower tau_upper hvrec hε
  have heq : (fun n => Real.log |e (n + 1)| - tau * Real.log |e n|) =
      fun n => v n + ((tau - 1) / 2) * Real.log |K| := by
    funext n
    dsimp [v, u]
    ring
  rw [heq]
  simpa only [zero_add] using hv.add_const (((tau - 1) / 2) * Real.log |K|)

/-- A three-error product recurrence with asymptotically nonzero coefficient
has the asserted power-quotient limit. Convergence of the errors is not needed
for this ratio theorem and is stated separately in the Q-order corollary. -/
theorem three_error_recurrence_q_ratio {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0)
    (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 3) = C n * e (n + 2) * e (n + 1) * e n) :
    Tendsto (fun n => |e (n + 1)| / |e n| ^ tau)
      atTop (𝓝 (|K| ^ ((tau - 1) / 2))) :=
  q_ratio_of_log_difference he hK (three_error_recurrence_log_limit he hK hC hrec)

/-- Q-order includes convergence to zero and a finite strictly positive
asymptotic quotient. -/
def HasQOrder (e : ℕ → ℝ) (p L : ℝ) : Prop :=
  1 < p ∧ 0 < L ∧ Tendsto e atTop (𝓝 0) ∧
    Tendsto (fun n => |e (n + 1)| / |e n| ^ p) atTop (𝓝 L)

/-- The complete genuine Q-order statement for a convergent three-error recurrence. -/
theorem three_error_recurrence_hasQOrder {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0)
    (hlim : Tendsto e atTop (𝓝 0)) (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ n, e (n + 3) = C n * e (n + 2) * e (n + 1) * e n) :
    HasQOrder e tau (|K| ^ ((tau - 1) / 2)) := by
  exact ⟨one_lt_tau, Real.rpow_pos_of_pos (abs_pos.mpr hK) _, hlim,
    three_error_recurrence_q_ratio he hK hC hrec⟩

/-- A finite exceptional prefix does not change the Q-order conclusion. -/
theorem eventual_three_error_recurrence_hasQOrder {e C : ℕ → ℝ} {K : ℝ}
    (he : ∀ n, e n ≠ 0) (hK : K ≠ 0)
    (hlim : Tendsto e atTop (𝓝 0)) (hC : Tendsto C atTop (𝓝 K))
    (hrec : ∀ᶠ n in atTop,
      e (n + 3) = C n * e (n + 2) * e (n + 1) * e n) :
    HasQOrder e tau (|K| ^ ((tau - 1) / 2)) := by
  obtain ⟨N, hN⟩ := eventually_atTop.mp hrec
  have hshift : Tendsto (fun n => C (n + N)) atTop (𝓝 K) :=
    (tendsto_add_atTop_iff_nat N).mpr hC
  have hr : ∀ n, e (n + 3 + N) = C (n + N) *
      e (n + 2 + N) * e (n + 1 + N) * e (n + N) := by
    intro n
    simpa only [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using
      hN (n + N) (Nat.le_add_left N n)
  have hq := three_error_recurrence_q_ratio (e := fun n => e (n + N))
    (C := fun n => C (n + N)) (fun n => he (n + N)) hK hshift hr
  refine ⟨one_lt_tau, Real.rpow_pos_of_pos (abs_pos.mpr hK) _, hlim, ?_⟩
  apply (tendsto_add_atTop_iff_nat N).mp
  simpa only [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hq

end
end ModAB.RationalOrder
