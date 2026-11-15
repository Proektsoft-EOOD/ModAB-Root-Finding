import ModAB.RationalCycle
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.Linarith

/-!
# Invariant cone for a two-step rational cycle

The interface consists of two explicit step equations with bounded negative
coefficients.  From these equations we prove signs, cone invariance, geometric
contraction, convergence and increasing separation.  Deriving the coefficient
bounds from C3 regularity and the mixed-ordinate interpolation identity remains
a separate analytic task; it is not included as an undocumented conclusion.
-/

namespace ModAB.Rational.CycleEntry
open Filter Topology
noncomputable section

/-- A fresh post-crossing state, expressed in errors relative to the root. -/
structure Cone (θ ε a s d : ℝ) : Prop where
  same : 0 < a * d
  opposite : a * s < 0
  latest_small : |a| ≤ θ * |s|
  opposite_small : |s| ≤ θ * |d|
  history_bound : |d| ≤ ε

theorem Cone.latest_ne_zero {θ ε a s d : ℝ} (h : Cone θ ε a s d) : a ≠ 0 := by
  intro hz
  have hh := h.same
  simp [hz] at hh

theorem Cone.opposite_ne_zero {θ ε a s d : ℝ} (h : Cone θ ε a s d) : s ≠ 0 := by
  intro hz
  have hh := h.opposite
  simp [hz] at hh

theorem Cone.history_ne_zero {θ ε a s d : ℝ} (h : Cone θ ε a s d) : d ≠ 0 := by
  intro hz
  have hh := h.same
  simp [hz] at hh

theorem Cone.ordered_errors {θ ε a s d : ℝ} (h : Cone θ ε a s d)
    (hθ : θ ≤ 1) : |a| ≤ |s| ∧ |s| ≤ |d| := by
  constructor
  · exact h.latest_small.trans (by nlinarith [abs_nonneg s])
  · exact h.opposite_small.trans (by nlinarith [abs_nonneg d])

/-- A uniform bound usable for either step of the pair. -/
theorem product_bound {u a s d C ε : ℝ}
    (hC : 0 ≤ C) (hε : 0 ≤ ε) (hu : |u| ≤ C)
    (ha : |a| ≤ ε) (hd : |d| ≤ ε) :
    |u * a * s * d| ≤ C * ε ^ 2 * |s| := by
  calc
    |u * a * s * d| = |u| * |a| * |d| * |s| := by
      simp only [abs_mul]
      ring
    _ ≤ C * ε * ε * |s| := by gcongr
    _ = C * ε ^ 2 * |s| := by ring

/-- The exact two-step sign change. The first point stays on the latest side;
the second crosses to the opposite side. -/
theorem pair_signs {a s d u v b a' : ℝ}
    (had : 0 < a * d) (has : a * s < 0) (hu : u < 0) (hv : v < 0)
    (hb : b = u * a * s * d) (ha' : a' = v * a * b * d) :
    0 < a * b ∧ 0 < a' * s ∧ a' * b < 0 := by
  have hd : d ≠ 0 := by intro hh; simp [hh] at had
  have hab : 0 < a * b := by
    rw [hb]
    convert mul_pos (mul_pos_of_neg_of_neg hu has) had using 1; ring
  have hbd : 0 < b * d := by
    rw [hb]
    convert mul_pos (mul_pos_of_neg_of_neg hu has) (sq_pos_of_ne_zero hd) using 1; ring
  have hb0 : b ≠ 0 := by intro hh; simp [hh] at hab
  refine ⟨hab, ?_, ?_⟩
  · rw [ha']
    convert mul_pos (mul_pos_of_neg_of_neg hv has) hbd using 1; ring
  · rw [ha']
    convert mul_neg_of_neg_of_pos (mul_neg_of_neg_of_pos hv had)
      (sq_pos_of_ne_zero hb0) using 1; ring

/-- A bounded pair step preserves the cone and reduces its radius by θ. -/
theorem pair_preserves_cone {θ ε C a s d u v b a' : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) (h : Cone θ ε a s d)
    (hu : u < 0) (hv : v < 0) (huC : |u| ≤ C) (hvC : |v| ≤ C)
    (hb : b = u * a * s * d) (ha' : a' = v * a * b * d) :
    Cone θ (θ * ε) a' b s := by
  have hord := h.ordered_errors hθ1.le
  have haε : |a| ≤ ε := hord.1.trans (hord.2.trans h.history_bound)
  have hbb : |b| ≤ θ * |s| := by
    rw [hb]
    exact (product_bound hC hε huC haε h.history_bound).trans
      (mul_le_mul_of_nonneg_right hsmall (abs_nonneg s))
  have haa : |a'| ≤ θ * |b| := by
    rw [ha']
    exact (product_bound hC hε hvC haε h.history_bound).trans
      (mul_le_mul_of_nonneg_right hsmall (abs_nonneg b))
  have hsign := pair_signs h.same h.opposite hu hv hb ha'
  exact ⟨hsign.2.1, hsign.2.2, haa, hbb,
    h.opposite_small.trans (mul_le_mul_of_nonneg_left h.history_bound hθ0.le)⟩

theorem Cone.enlarge_radius {θ ε ε' a s d : ℝ}
    (h : Cone θ ε a s d) (hε : ε ≤ ε') : Cone θ ε' a s d :=
  ⟨h.same, h.opposite, h.latest_small, h.opposite_small, h.history_bound.trans hε⟩

/-- A sequence of pairs with the explicit two-step equations. -/
structure Run (θ ε C : ℝ) where
  a : ℕ → ℝ
  s : ℕ → ℝ
  d : ℕ → ℝ
  u : ℕ → ℝ
  v : ℕ → ℝ
  initial : Cone θ ε (a 0) (s 0) (d 0)
  u_neg : ∀ n, u n < 0
  v_neg : ∀ n, v n < 0
  u_bound : ∀ n, |u n| ≤ C
  v_bound : ∀ n, |v n| ≤ C
  first_step : ∀ n, s (n + 1) = u n * a n * s n * d n
  second_step : ∀ n, a (n + 1) = v n * a n * s (n + 1) * d n
  drop_step : ∀ n, d (n + 1) = s n

variable {θ ε C : ℝ}

/-- Invariance is derived by induction from the step equations and coefficient bounds. -/
theorem Run.all_cones (R : Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) : ∀ n, Cone θ ε (R.a n) (R.s n) (R.d n) := by
  intro n
  induction n with
  | zero => exact R.initial
  | succ n ih =>
    have hn := pair_preserves_cone hθ0 hθ1 hC hε hsmall ih
      (R.u_neg n) (R.v_neg n) (R.u_bound n) (R.v_bound n)
      (R.first_step n) (R.second_step n)
    rw [R.drop_step]
    exact hn.enlarge_radius (by nlinarith)

/-- The selected discarded error decreases geometrically across pairs. -/
theorem Run.history_geometric_bound (R : Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) : ∀ n, |R.d n| ≤ θ ^ n * ε := by
  have hc := R.all_cones hθ0 hθ1 hC hε hsmall
  intro n
  induction n with
  | zero => simpa using R.initial.history_bound
  | succ n ih =>
    rw [R.drop_step]
    calc
      |R.s n| ≤ θ * |R.d n| := (hc n).opposite_small
      _ ≤ θ * (θ ^ n * ε) := mul_le_mul_of_nonneg_left ih hθ0.le
      _ = θ ^ (n + 1) * ε := by rw [pow_succ]; ring

/-- All three error sequences converge to zero; convergence is not a Run field. -/
theorem Run.errors_tendsto (R : Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) :
    Tendsto R.a atTop (𝓝 0) ∧ Tendsto R.s atTop (𝓝 0) ∧ Tendsto R.d atTop (𝓝 0) := by
  have hc := R.all_cones hθ0 hθ1 hC hε hsmall
  have hp : Tendsto (fun n : ℕ => θ ^ n * ε) atTop (𝓝 (0 : ℝ)) := by
    simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one hθ0.le hθ1).mul_const ε
  have hd : Tendsto (fun n => |R.d n|) atTop (𝓝 0) :=
    squeeze_zero (fun n => abs_nonneg _) (R.history_geometric_bound hθ0 hθ1 hC hε hsmall) hp
  have hs : Tendsto (fun n => |R.s n|) atTop (𝓝 0) :=
    squeeze_zero (fun n => abs_nonneg _) (fun n => ((hc n).ordered_errors hθ1.le).2) hd
  have ha : Tendsto (fun n => |R.a n|) atTop (𝓝 0) :=
    squeeze_zero (fun n => abs_nonneg _) (fun n => ((hc n).ordered_errors hθ1.le).1) hs
  exact ⟨(tendsto_zero_iff_abs_tendsto_zero _).mpr ha,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr hs, (tendsto_zero_iff_abs_tendsto_zero _).mpr hd⟩

/-- Separation strengthens along the orbit; a fixed cone is not confused with
the vanishing ratios needed by the asymptotic coefficient calculation. -/
theorem Run.separation_tendsto (R : Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) :
    Tendsto (fun n => R.a n / R.s n) atTop (𝓝 0) ∧
    Tendsto (fun n => R.s n / R.d n) atTop (𝓝 0) := by
  have hc := R.all_cones hθ0 hθ1 hC hε hsmall
  have he := R.errors_tendsto hθ0 hθ1 hC hε hsmall
  have hm : Tendsto (fun n => C * |R.a n| * |R.d n|) atTop (𝓝 (0 : ℝ)) := by
    simpa using ((he.1.abs.const_mul C).mul he.2.2.abs)
  have hfirst : ∀ n, |R.a (n + 1) / R.s (n + 1)| ≤ C * |R.a n| * |R.d n| := by
    intro n
    have hs := (hc (n + 1)).opposite_ne_zero
    have hid : R.a (n + 1) / R.s (n + 1) = R.v n * R.a n * R.d n := by
      rw [R.second_step]
      field_simp
    rw [hid, abs_mul, abs_mul]
    gcongr
    exact R.v_bound n
  have hsecond : ∀ n, |R.s (n + 1) / R.d (n + 1)| ≤ C * |R.a n| * |R.d n| := by
    intro n
    have hs := (hc n).opposite_ne_zero
    have hid : R.s (n + 1) / R.d (n + 1) = R.u n * R.a n * R.d n := by
      rw [R.drop_step, R.first_step]
      field_simp
    rw [hid, abs_mul, abs_mul]
    gcongr
    exact R.u_bound n
  have hf := squeeze_zero (fun n => abs_nonneg _) hfirst hm
  have hs := squeeze_zero (fun n => abs_nonneg _) hsecond hm
  constructor
  · exact (tendsto_add_atTop_iff_nat 1).mp ((tendsto_zero_iff_abs_tendsto_zero _).mpr hf)
  · exact (tendsto_add_atTop_iff_nat 1).mp ((tendsto_zero_iff_abs_tendsto_zero _).mpr hs)

/-- The selected history makes the first coefficient exactly the local step coefficient. -/
theorem Run.firstCoeff_eq (R : Run θ ε C)
    (ha : ∀ n, R.a n ≠ 0) (hs : ∀ n, R.s n ≠ 0) (n : ℕ) :
    Cycle.firstCoeff R.a (fun k => R.s (k + 1)) n = R.u (n + 2) := by
  dsimp only [Cycle.firstCoeff]
  rw [R.first_step]
  have hd : R.d (n + 2) = R.s (n + 1) := R.drop_step (n + 1)
  rw [hd]
  field_simp [ha, hs]

/-- The crossing coefficient is likewise recovered exactly from the state sequence. -/
theorem Run.secondCoeff_eq (R : Run θ ε C)
    (ha : ∀ n, R.a n ≠ 0) (hs : ∀ n, R.s n ≠ 0) (n : ℕ) :
    Cycle.secondCoeff R.a (fun k => R.s (k + 1)) n = R.v (n + 2) := by
  dsimp only [Cycle.secondCoeff]
  rw [show n + 3 = (n + 2) + 1 by omega, R.second_step]
  have hd : R.d (n + 2) = R.s (n + 1) := R.drop_step (n + 1)
  rw [hd]
  field_simp [ha, hs]

/-- The invariant cone and coefficient limits imply the cubic and alternating
order quotients; error convergence was proved separately by `errors_tendsto`. -/
theorem Run.orders (R : Run θ ε C) {K : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) (hK : K ≠ 0)
    (hu : Tendsto R.u atTop (𝓝 K)) (hv : Tendsto R.v atTop (𝓝 K)) :
    Tendsto (fun n => |R.s (n + 2)| / |R.s (n + 1)| ^ 3) atTop (𝓝 |K|) ∧
    Tendsto (fun n => |R.s (n + 1)| / |R.a n| ^ (9 / 5 : ℝ)) atTop
      (𝓝 (|K| ^ (2 / 5 : ℝ))) ∧
    Tendsto (fun n => |R.a (n + 1)| / |R.s (n + 1)| ^ (5 / 3 : ℝ)) atTop
      (𝓝 (|K| ^ (1 / 3 : ℝ))) := by
  have hc := R.all_cones hθ0 hθ1 hC hε hsmall
  have ha : ∀ n, R.a n ≠ 0 := fun n => (hc n).latest_ne_zero
  have hs : ∀ n, R.s n ≠ 0 := fun n => (hc n).opposite_ne_zero
  have hu' : Tendsto (Cycle.firstCoeff R.a (fun k => R.s (k + 1))) atTop (𝓝 K) := by
    convert (tendsto_add_atTop_iff_nat 2).mpr hu using 1
    ext n
    exact R.firstCoeff_eq ha hs n
  have hv' : Tendsto (Cycle.secondCoeff R.a (fun k => R.s (k + 1))) atTop (𝓝 K) := by
    convert (tendsto_add_atTop_iff_nat 2).mpr hv using 1
    ext n
    exact R.secondCoeff_eq ha hs n
  simpa only [Nat.add_assoc, Nat.reduceAdd] using
    Cycle.cycle_orders_of_coefficients ha (fun n => hs (n + 1)) hK hu' hv'

/-- Combining the two unequal step exponents recovers a cubic quotient for a. -/
theorem latest_cubic_identity {A B A' : ℝ}
    (hA : 0 < A) (hB : 0 < B) (hA' : 0 < A') :
    A' / A ^ 3 = (A' / B ^ (5 / 3 : ℝ)) * (B / A ^ (9 / 5 : ℝ)) ^ (5 / 3 : ℝ) := by
  apply Real.log_injOn_pos (Set.mem_Ioi.mpr (by positivity)) (Set.mem_Ioi.mpr (by positivity))
  rw [Real.log_div (by positivity) (by positivity), Real.log_pow,
    Real.log_mul (by positivity) (by positivity),
    Real.log_div (by positivity) (by positivity), Real.log_rpow hB,
    Real.log_rpow (by positivity),
    Real.log_div (by positivity) (by positivity), Real.log_rpow hA]
  ring

theorem latest_cubic_constant {κ : ℝ} (hκ : 0 < κ) :
    κ ^ (1 / 3 : ℝ) * (κ ^ (2 / 5 : ℝ)) ^ (5 / 3 : ℝ) = κ := by
  apply Real.log_injOn_pos (Set.mem_Ioi.mpr (by positivity)) (Set.mem_Ioi.mpr hκ)
  rw [Real.log_mul (by positivity) (by positivity), Real.log_rpow hκ,
    Real.log_rpow (by positivity), Real.log_rpow hκ]
  ring

/-- Both pair subsequences have the same cubic magnitude constant. -/
theorem Run.latest_cubic_order (R : Run θ ε C) {K : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) (hK : K ≠ 0)
    (hu : Tendsto R.u atTop (𝓝 K)) (hv : Tendsto R.v atTop (𝓝 K)) :
    Tendsto (fun n => |R.a (n + 1)| / |R.a n| ^ 3) atTop (𝓝 |K|) := by
  have hc := R.all_cones hθ0 hθ1 hC hε hsmall
  have ho := R.orders hθ0 hθ1 hC hε hsmall hK hu hv
  have hκ : 0 < |K| := abs_pos.mpr hK
  have hlim := ho.2.2.mul (ho.2.1.rpow_const (p := (5 / 3 : ℝ)) (Or.inr (by norm_num)))
  convert hlim using 1
  · ext n
    exact latest_cubic_identity (abs_pos.mpr (hc n).latest_ne_zero)
      (abs_pos.mpr (hc (n + 1)).opposite_ne_zero)
      (abs_pos.mpr (hc (n + 1)).latest_ne_zero)
  · rw [latest_cubic_constant hκ]

end
end ModAB.Rational.CycleEntry
