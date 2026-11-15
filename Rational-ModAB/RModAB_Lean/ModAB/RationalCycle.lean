import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Convert

/-!
# The two-step local cycle

This module starts with the two explicit recurrence coefficients supplied by
the local interpolation analysis.  It proves their consequences, rather than
assuming cubic convergence or the two alternating orders.  Establishing entry
into this local cycle from a function and an algorithmic state is a separate
analytic obligation.
-/

namespace ModAB.Rational.Cycle
open Filter Topology
noncomputable section

/-- Coefficient of the first step of a pair, indexed without natural subtraction. -/
def firstCoeff (a b : ℕ → ℝ) (n : ℕ) : ℝ :=
  b (n + 2) / (a (n + 2) * b n * b (n + 1))

/-- Coefficient of the crossing step of the same pair. -/
def secondCoeff (a b : ℕ → ℝ) (n : ℕ) : ℝ :=
  a (n + 3) / (a (n + 2) * b (n + 2) * b n)

/-- The cubic quotient on one of the pair subsequences. -/
def cubicCoeff (b : ℕ → ℝ) (n : ℕ) : ℝ := b (n + 1) / b n ^ 3

/-- Exact elimination of the intermediate point; no limiting argument occurs here. -/
theorem cubicCoeff_identity {a b : ℕ → ℝ}
    (ha : ∀ n, a n ≠ 0) (hb : ∀ n, b n ≠ 0) (n : ℕ) :
    cubicCoeff b (n + 2) =
      firstCoeff a b (n + 1) * secondCoeff a b n / firstCoeff a b n := by
  simp only [cubicCoeff, firstCoeff, secondCoeff, Nat.add_assoc]
  field_simp [ha, hb]

/-- The two signed recurrence coefficients tending to K imply a cubic quotient K. -/
theorem cubicCoeff_tendsto {a b : ℕ → ℝ} {K : ℝ}
    (ha : ∀ n, a n ≠ 0) (hb : ∀ n, b n ≠ 0) (hK : K ≠ 0)
    (hfirst : Tendsto (firstCoeff a b) atTop (𝓝 K))
    (hsecond : Tendsto (secondCoeff a b) atTop (𝓝 K)) :
    Tendsto (cubicCoeff b) atTop (𝓝 K) := by
  have hshift : Tendsto (fun n => firstCoeff a b (n + 1)) atTop (𝓝 K) :=
    (tendsto_add_atTop_iff_nat 1).mpr hfirst
  have hquot := (hshift.mul hsecond).div hfirst hK
  have hlim : Tendsto (fun n => cubicCoeff b (n + 2)) atTop (𝓝 K) := by
    convert hquot using 1
    · ext n
      exact cubicCoeff_identity ha hb n
    · congr 1
      field_simp
  exact (tendsto_add_atTop_iff_nat 2).mp hlim

/-- The absolute-value form of the cubic limit, valid also for negative K. -/
theorem abs_cubicCoeff_tendsto {a b : ℕ → ℝ} {K : ℝ}
    (ha : ∀ n, a n ≠ 0) (hb : ∀ n, b n ≠ 0) (hK : K ≠ 0)
    (hfirst : Tendsto (firstCoeff a b) atTop (𝓝 K))
    (hsecond : Tendsto (secondCoeff a b) atTop (𝓝 K)) :
    Tendsto (fun n => |b (n + 1)| / |b n| ^ 3) atTop (𝓝 |K|) := by
  simpa only [cubicCoeff, abs_div, abs_pow] using
    (cubicCoeff_tendsto ha hb hK hfirst hsecond).abs

/-- A logarithmic normalization identity for positive magnitudes. -/
theorem first_order_identity {A x y z : ℝ}
    (hA : 0 < A) (hx : 0 < x) (hy : 0 < y) (hz : 0 < z) :
    z / A ^ (9 / 5 : ℝ) =
      (z / (A * x * y)) ^ (9 / 5 : ℝ) /
        ((y / x ^ 3) ^ (3 / 5 : ℝ) * (z / y ^ 3) ^ (4 / 5 : ℝ)) := by
  apply Real.log_injOn_pos (Set.mem_Ioi.mpr (by positivity)) (Set.mem_Ioi.mpr (by positivity))
  rw [Real.log_div (by positivity) (by positivity), Real.log_rpow hA]
  rw [Real.log_div (by positivity) (by positivity)]
  rw [Real.log_rpow (by positivity), Real.log_mul (by positivity) (by positivity)]
  rw [Real.log_rpow (by positivity), Real.log_rpow (by positivity)]
  rw [Real.log_div (by positivity) (by positivity),
    Real.log_div (by positivity) (by positivity),
    Real.log_div (by positivity) (by positivity)]
  rw [Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]
  ring

/-- The corresponding identity for the crossing step. -/
theorem second_order_identity {A x y z : ℝ}
    (hA : 0 < A) (hx : 0 < x) (hy : 0 < y) (hz : 0 < z) :
    A / y ^ (5 / 3 : ℝ) =
      (z / y ^ 3) * (y / x ^ 3) ^ (1 / 3 : ℝ) / (z / (A * x * y)) := by
  apply Real.log_injOn_pos (Set.mem_Ioi.mpr (by positivity)) (Set.mem_Ioi.mpr (by positivity))
  rw [Real.log_div (by positivity) (by positivity), Real.log_rpow hy]
  rw [Real.log_div (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_rpow (by positivity)]
  rw [Real.log_div (by positivity) (by positivity),
    Real.log_div (by positivity) (by positivity),
    Real.log_div (by positivity) (by positivity)]
  rw [Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]
  ring

/-- Simplification of the limiting constant for the first step. -/
theorem first_constant {κ : ℝ} (hκ : 0 < κ) :
    κ ^ (9 / 5 : ℝ) / (κ ^ (3 / 5 : ℝ) * κ ^ (4 / 5 : ℝ)) =
      κ ^ (2 / 5 : ℝ) := by
  apply Real.log_injOn_pos (Set.mem_Ioi.mpr (by positivity)) (Set.mem_Ioi.mpr (by positivity))
  rw [Real.log_div (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity)]
  simp only [Real.log_rpow hκ]
  ring

/-- Simplification of the limiting constant for the crossing step. -/
theorem second_constant {κ : ℝ} (hκ : 0 < κ) :
    κ * κ ^ (1 / 3 : ℝ) / κ = κ ^ (1 / 3 : ℝ) := by
  field_simp

/-- Cubic convergence of b and the first-step coefficient determine order 9/5. -/
theorem first_order_tendsto {a b : ℕ → ℝ} {κ : ℝ}
    (ha : ∀ n, 0 < a n) (hb : ∀ n, 0 < b n) (hκ : 0 < κ)
    (hfirst : Tendsto (firstCoeff a b) atTop (𝓝 κ))
    (hcubic : Tendsto (cubicCoeff b) atTop (𝓝 κ)) :
    Tendsto (fun n => b n / (a n) ^ (9 / 5 : ℝ)) atTop
      (𝓝 (κ ^ (2 / 5 : ℝ))) := by
  have hshift : Tendsto (fun n => cubicCoeff b (n + 1)) atTop (𝓝 κ) :=
    (tendsto_add_atTop_iff_nat 1).mpr hcubic
  have hlim := (hfirst.rpow_const (p := (9 / 5 : ℝ)) (Or.inl hκ.ne')).div
    ((hcubic.rpow_const (p := (3 / 5 : ℝ)) (Or.inl hκ.ne')).mul
      (hshift.rpow_const (p := (4 / 5 : ℝ)) (Or.inl hκ.ne')))
    (show κ ^ (3 / 5 : ℝ) * κ ^ (4 / 5 : ℝ) ≠ 0 by positivity)
  have htail : Tendsto (fun n => b (n + 2) / (a (n + 2)) ^ (9 / 5 : ℝ))
      atTop (𝓝 (κ ^ (2 / 5 : ℝ))) := by
    convert hlim using 1
    · ext n
      simpa only [firstCoeff, cubicCoeff, Nat.add_assoc, Nat.reduceAdd, Pi.div_apply] using
        first_order_identity (ha (n + 2)) (hb n) (hb (n + 1)) (hb (n + 2))
    · rw [first_constant hκ]
  exact (tendsto_add_atTop_iff_nat 2).mp htail

/-- The crossing step has order 5/3 and its own asymptotic constant. -/
theorem second_order_tendsto {a b : ℕ → ℝ} {κ : ℝ}
    (ha : ∀ n, 0 < a n) (hb : ∀ n, 0 < b n) (hκ : 0 < κ)
    (hfirst : Tendsto (firstCoeff a b) atTop (𝓝 κ))
    (hcubic : Tendsto (cubicCoeff b) atTop (𝓝 κ)) :
    Tendsto (fun n => a (n + 1) / (b n) ^ (5 / 3 : ℝ)) atTop
      (𝓝 (κ ^ (1 / 3 : ℝ))) := by
  have hc1 : Tendsto (fun n => cubicCoeff b (n + 1)) atTop (𝓝 κ) :=
    (tendsto_add_atTop_iff_nat 1).mpr hcubic
  have hc2 : Tendsto (fun n => cubicCoeff b (n + 2)) atTop (𝓝 κ) :=
    (tendsto_add_atTop_iff_nat 2).mpr hcubic
  have hf1 : Tendsto (fun n => firstCoeff a b (n + 1)) atTop (𝓝 κ) :=
    (tendsto_add_atTop_iff_nat 1).mpr hfirst
  have hlim := (hc2.mul (hc1.rpow_const (p := (1 / 3 : ℝ)) (Or.inl hκ.ne'))).div hf1 hκ.ne'
  have htail : Tendsto (fun n => a (n + 3) / (b (n + 2)) ^ (5 / 3 : ℝ))
      atTop (𝓝 (κ ^ (1 / 3 : ℝ))) := by
    convert hlim using 1
    · ext n
      simpa only [firstCoeff, cubicCoeff, Nat.add_assoc, Nat.reduceAdd, Pi.div_apply] using
        second_order_identity (ha (n + 3)) (hb (n + 1)) (hb (n + 2)) (hb (n + 3))
    · rw [second_constant hκ]
  exact (tendsto_add_atTop_iff_nat 2).mp htail

/-- The two unequal step exponents have product three. -/
theorem alternating_orders_product : (9 / 5 : ℝ) * (5 / 3 : ℝ) = 3 := by norm_num

/-- Their geometric mean is sqrt(3), a cycle-normalized value. -/
theorem alternating_orders_geometric_mean :
    Real.sqrt ((9 / 5 : ℝ) * (5 / 3 : ℝ)) = Real.sqrt 3 := by
  rw [alternating_orders_product]

/-- Passing from signed errors to their magnitudes commutes with the first coefficient. -/
theorem firstCoeff_abs (a b : ℕ → ℝ) (n : ℕ) :
    firstCoeff (fun k => |a k|) (fun k => |b k|) n = |firstCoeff a b n| := by
  simp only [firstCoeff, abs_div, abs_mul]

/-- Passing from signed errors to their magnitudes commutes with the second coefficient. -/
theorem secondCoeff_abs (a b : ℕ → ℝ) (n : ℕ) :
    secondCoeff (fun k => |a k|) (fun k => |b k|) n = |secondCoeff a b n| := by
  simp only [secondCoeff, abs_div, abs_mul]

/-- The full cycle consequence for nonzero signed errors.  In particular, it
applies to K < 0.  The hypotheses are precisely the two recurrence-coefficient
limits; none of the three concluded order limits is assumed. -/
theorem cycle_orders_of_coefficients {a b : ℕ → ℝ} {K : ℝ}
    (ha : ∀ n, a n ≠ 0) (hb : ∀ n, b n ≠ 0) (hK : K ≠ 0)
    (hfirst : Tendsto (firstCoeff a b) atTop (𝓝 K))
    (hsecond : Tendsto (secondCoeff a b) atTop (𝓝 K)) :
    Tendsto (fun n => |b (n + 1)| / |b n| ^ 3) atTop (𝓝 |K|) ∧
    Tendsto (fun n => |b n| / |a n| ^ (9 / 5 : ℝ)) atTop
      (𝓝 (|K| ^ (2 / 5 : ℝ))) ∧
    Tendsto (fun n => |a (n + 1)| / |b n| ^ (5 / 3 : ℝ)) atTop
      (𝓝 (|K| ^ (1 / 3 : ℝ))) := by
  have hκ : 0 < |K| := abs_pos.mpr hK
  have hpa : ∀ n, 0 < |a n| := fun n => abs_pos.mpr (ha n)
  have hpb : ∀ n, 0 < |b n| := fun n => abs_pos.mpr (hb n)
  have hfa : Tendsto (firstCoeff (fun k => |a k|) (fun k => |b k|))
      atTop (𝓝 |K|) := by
    convert hfirst.abs using 1
    ext n
    exact firstCoeff_abs a b n
  have hcb : Tendsto (cubicCoeff (fun k => |b k|)) atTop (𝓝 |K|) :=
    abs_cubicCoeff_tendsto ha hb hK hfirst hsecond
  exact ⟨hcb, first_order_tendsto hpa hpb hκ hfa hcb,
    second_order_tendsto hpa hpb hκ hfa hcb⟩

end
end ModAB.Rational.Cycle
