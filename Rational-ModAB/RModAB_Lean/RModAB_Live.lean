/-
Standalone RModAB paper: exact-real proofs, release 0.8.0 (11 November 2026).
Self-contained source for a Lean 4.34.0 / matching Mathlib environment.
Includes inherited exact-real dependencies and all new rational modules.
It does NOT include the old binary64/FloatLib layer; use the full project
or ModAB_Live_Complete.lean for that layer. This is not an IEEE solver proof.
The file is checked locally with the pinned project; no browser execution is implied.
-/
import Lean.Util.CollectAxioms
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Analytic.IsolatedZeros
import Mathlib.Analysis.Calculus.ContDiff.Deriv
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.DSlope
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.Calculus.LHopital
import Mathlib.Analysis.Calculus.LocalExtr.Rolle
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Basic.Real.Basic
import Mathlib.NumberTheory.Real.GoldenRatio
import Mathlib.Tactic
import Mathlib.Tactic.Convert
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring
import Mathlib.Topology.Instances.Real.Lemmas
import Mathlib.Topology.Order.Basic
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Order.MonotoneConvergence
import Mathlib.Topology.Order.OrderClosed
import Mathlib.Topology.Sequences

/- Original module: ModAB.Switching -/
section RModABBundle0

/-!+# Rounding errors in the ModAB switching criterion

The mathematical statements are over the real numbers. Rounding errors and
range preservation are explicit hypotheses about an evaluation trace. No
claim identifying this trace with Lean Float or a C# execution is assumed
implicitly. All error bounds below are proved from the local hypotheses.
-/

namespace ModAB

noncomputable section

def Within (x y e : ℝ) : Prop := |x - y| ≤ e

lemma within_refl (x : ℝ) : Within x x 0 := by simp [Within]

lemma within_trans {x y z e d : ℝ} (h : Within x y e) (g : Within y z d) :
    Within x z (e + d) := by
  rcases abs_le.mp h with ⟨h₁, h₂⟩
  rcases abs_le.mp g with ⟨g₁, g₂⟩
  apply abs_le.mpr
  constructor <;> linarith

lemma within_add {x x' y y' e d : ℝ} (h : Within x' x e) (g : Within y' y d) :
    Within (x' + y') (x + y) (e + d) := by
  rcases abs_le.mp h with ⟨h₁, h₂⟩
  rcases abs_le.mp g with ⟨g₁, g₂⟩
  apply abs_le.mpr
  constructor <;> linarith

lemma within_sub {x x' y y' e d : ℝ} (h : Within x' x e) (g : Within y' y d) :
    Within (x' - y') (x - y) (e + d) := by
  rcases abs_le.mp h with ⟨h₁, h₂⟩
  rcases abs_le.mp g with ⟨g₁, g₂⟩
  apply abs_le.mpr
  constructor <;> linarith

lemma within_abs {x y e : ℝ} (h : Within x y e) : Within |x| |y| e := by
  exact le_trans (abs_abs_sub_abs_le_abs_sub x y) h

lemma within_half {x y e : ℝ} (h : Within x y e) : Within (x/2) (y/2) (e/2) := by
  rcases abs_le.mp h with ⟨h₁, h₂⟩
  apply abs_le.mpr
  constructor <;> linarith

lemma within_mono {x y e d : ℝ} (h : Within x y e) (g : e ≤ d) : Within x y d :=
  le_trans h g

lemma product_error {x x' y y' A B C D : ℝ}
    (hx : |x'| ≤ A) (hy : Within y' y B)
    (hd : Within x' x C) (hb : |y| ≤ D) :
    Within (x'*y') (x*y) (A*B+C*D) := by
  have hA : 0 ≤ A := le_trans (abs_nonneg _) hx
  have hC : 0 ≤ C := le_trans (abs_nonneg _) hd
  calc
    |x'*y' - x*y| = |x'*(y'-y)+(x'-x)*y| := by congr 1; ring
    _ ≤ |x'*(y'-y)| + |(x'-x)*y| := abs_add_le _ _
    _ = |x'| *|y'-y| + |x'-x| *|y| := by rw [abs_mul, abs_mul]
    _ ≤ A*B+C*D := add_le_add
      (mul_le_mul hx hy (abs_nonneg _) hA)
      (mul_le_mul hd hb (abs_nonneg _) hC)

lemma square_error {x y e : ℝ} (hx : x ∈ Set.Icc (0:ℝ) 1)
    (hy : y ∈ Set.Icc (0:ℝ) 1) (h : Within x y e) :
    Within (x^2) (y^2) (2*e) := by
  have p := product_error (x := y) (x' := x) (y := y) (y' := x)
    (A := 1) (B := e) (C := e) (D := 1)
    (by rw [abs_of_nonneg hx.1]; exact hx.2) h h
    (by rw [abs_of_nonneg hy.1]; exact hy.2)
  simpa [pow_two, mul_two, two_mul] using p

def kappa (r : ℝ) : ℝ := ((1+3*r)/(2*(1+r)))^2

lemma kappa_identity {r : ℝ} (hr : 0 ≤ r) :
    kappa r = (1 - (1-r)/(2*(1+r)))^2 := by
  unfold kappa
  congr 1
  field_simp
  ; ring

lemma kappa_hasDerivAt (r : ℝ) (hr : 0 ≤ r) :
    HasDerivAt kappa ((1+3*r)/(1+r)^3) r := by
  have hd : 2*(1+r) ≠ 0 := by positivity
  have h := (((hasDerivAt_const r (1:ℝ)).add
    ((hasDerivAt_id r).const_mul 3)).div
    (((hasDerivAt_const r (1:ℝ)).add (hasDerivAt_id r)).const_mul 2) hd).pow 2
  convert h using 1 <;> (try funext x) <;> (try dsimp [kappa]) <;> field_simp <;> ring

lemma kappa_derivative_bounds {r : ℝ} (hr : 0 ≤ r) :
    0 < (1+3*r)/(1+r)^3 ∧ (1+3*r)/(1+r)^3 ≤ 1 := by
  have hp : 0 < (1+r)^3 := by positivity
  constructor
  · positivity
  · apply (div_le_iff₀ hp).mpr
    nlinarith [mul_nonneg (sq_nonneg r) (by linarith : 0 ≤ 3+r)]

lemma kappa_lipschitz_ordered {x y : ℝ} (hx : 0 ≤ x) (hxy : x < y) :
    |kappa y-kappa x| ≤ y-x := by
  have hc : ContinuousOn kappa (Set.Icc x y) := fun r hr =>
    (kappa_hasDerivAt r (hx.trans hr.1)).continuousAt.continuousWithinAt
  obtain ⟨c,hc',hslope⟩ := exists_hasDerivAt_eq_slope kappa
    (fun r => (1+3*r)/(1+r)^3) hxy hc
    (fun r hr => kappa_hasDerivAt r (by linarith [hr.1]))
  have hd := kappa_derivative_bounds (show 0≤c by linarith [hc'.1])
  have habs : |(1+3*c)/(1+c)^3| ≤ 1 := by rw [abs_of_pos hd.1]; exact hd.2
  rw [hslope,abs_div,abs_of_pos (sub_pos.mpr hxy)] at habs
  simpa using (div_le_iff₀ (sub_pos.mpr hxy)).mp habs

/-- The manuscript's derivative and mean-value argument, also valid for r≥0. -/
lemma kappa_lipschitz {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) :
    |kappa x-kappa y| ≤ |x-y| := by
  rcases lt_trichotomy x y with h|h|h
  · simpa only [abs_sub_comm (kappa x),abs_sub_comm x,
      abs_of_pos (sub_pos.mpr h)] using kappa_lipschitz_ordered hx h
  · subst y; simp
  · simpa only [abs_of_pos (sub_pos.mpr h)] using kappa_lipschitz_ordered hy h

/-- The separate algebraic proof printed in `thm:formal_switching`. -/
lemma kappa_lipschitz_algebraic {x y : ℝ} (hx : 0≤x) (hy : 0≤y) :
    |kappa x-kappa y|≤|x-y| := by
  let d:=(1+x)^2*(1+y)^2
  let c:=(1+2*x+2*y+3*x*y)/d
  have hd : 0<d:=by dsimp [d]; positivity
  have hc0 : 0≤c:=by dsimp [c]; positivity
  have hc1 : c≤1:=by
    apply (div_le_iff₀ hd).mpr
    dsimp [d]
    have hp : 0≤x*y+x^2+y^2+2*x^2*y+2*x*y^2+x^2*y^2:=by positivity
    nlinarith only [hp]
  have hid : kappa x-kappa y=(x-y)*c:=by
    dsimp [kappa,c,d]
    field_simp
    ring
  rw [hid,abs_mul,abs_of_nonneg hc0]
  nlinarith [abs_nonneg (x-y)]

lemma ratio_range {r : ℝ} (hr : r ∈ Set.Icc (0:ℝ) 1) :
    (1-r)/(2*(1+r)) ∈ Set.Icc (0:ℝ) (1/2) := by
  have hd : 0 < 2*(1+r) := by linarith [hr.1]
  constructor
  · apply div_nonneg <;> linarith [hr.1, hr.2]
  · apply (div_le_iff₀ hd).mpr
    linarith [hr.1]

lemma quotient_error {a A b B e : ℝ} (he : 0 ≤ e)
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : 1 ≤ b) (hB : 1 ≤ B)
    (hae : Within A a e) (hbe : Within B b e) :
    Within ((A/2)/B) (a/(2*b)) e := by
  change |A-a| ≤ e at hae
  have bp : 0 < b := by linarith
  have Bp : 0 < B := by linarith
  have h₁ : Within (A/(2*B)) (a/(2*B)) (e/2) := by
    unfold Within
    rw [← sub_div, abs_div, abs_of_pos (by positivity : 0 < 2*B)]
    apply (div_le_iff₀ (by positivity : 0 < 2*B)).mpr
    have hp := mul_nonneg he (sub_nonneg.mpr hB)
    nlinarith [hae]
  have h₂ : Within (a/(2*B)) (a/(2*b)) (e/2) := by
    have hid : a/(2*B)-a/(2*b) = a*(b-B)/(2*B*b) := by field_simp
    unfold Within
    rw [hid, abs_div, abs_mul, abs_of_nonneg ha.1,
      abs_of_pos (by positivity : 0 < 2*B*b)]
    apply (div_le_iff₀ (by positivity : 0 < 2*B*b)).mpr
    have hab : |b-B| ≤ e := by simpa [Within, abs_sub_comm] using hbe
    have hm : a*|b-B| ≤ e := by
      calc
        a*|b-B| ≤ 1*e := mul_le_mul ha.2 hab (abs_nonneg _) (by norm_num)
        _ = e := one_mul e
    have hprod : 1 ≤ B*b := by nlinarith [mul_nonneg (sub_nonneg.mpr hB) (sub_nonneg.mpr hb)]
    have hep := mul_nonneg he (sub_nonneg.mpr hprod)
    nlinarith
  have h := within_trans h₁ h₂
  convert h using 1 <;> ring

/-- Local rounding contracts for the separately scaled symmetry factor.
The exact halving of `aHat` is represented by its occurrence as `aHat/2`.
The binary64 justification of these local contracts is a separate obligation. -/
structure FactorTrace (e rho : ℝ) where
  (rhoHat aHat bHat vHat rHat kHat : ℝ)
  rho_range : rhoHat ∈ Set.Icc (0:ℝ) 1
  a_range : aHat ∈ Set.Icc (0:ℝ) 1
  b_range : bHat ∈ Set.Icc (1:ℝ) 2
  v_range : vHat ∈ Set.Icc (0:ℝ) (1/2)
  r_range : rHat ∈ Set.Icc (1/2:ℝ) 1
  k_range : kHat ∈ Set.Icc (1/4:ℝ) 1
  rho_error : Within rhoHat rho e
  a_error : Within aHat (1-rhoHat) e
  b_error : Within bHat (1+rhoHat) e
  v_error : Within vHat ((aHat/2)/bHat) e
  r_error : Within rHat (1-vHat) e
  k_error : Within kHat (rHat^2) e

theorem factor_error_from_coefficient {e rho : ℝ} (he : 0≤e)
    (t : FactorTrace e rho) (hp : Within (kappa t.rhoHat) (kappa rho) e) :
    Within t.kHat (kappa rho) (8*e) := by
  let v₀ := (1-t.rhoHat)/(2*(1+t.rhoHat))
  have hv₀ := ratio_range t.rho_range
  have hq := quotient_error he
    (a := 1-t.rhoHat) (A := t.aHat) (b := 1+t.rhoHat) (B := t.bHat)
    (by constructor <;> linarith [t.rho_range.1, t.rho_range.2])
    (by linarith [t.rho_range.1]) t.b_range.1 t.a_error t.b_error
  have hv : Within t.vHat v₀ (2*e) := by
    simpa [v₀, two_mul] using within_trans t.v_error hq
  have hr : Within t.rHat (1-v₀) (3*e) := by
    have hs := within_sub (within_refl (1:ℝ)) hv
    have ht := within_trans t.r_error hs
    convert ht using 1 ; ring
  have hs := square_error
    (x := t.rHat) (y := 1-v₀)
    (by constructor <;> linarith [t.r_range.1, t.r_range.2])
    (by constructor <;> dsimp [v₀] <;> linarith [hv₀.1, hv₀.2]) hr
  have hk : Within t.kHat (kappa t.rhoHat) (7*e) := by
    rw [kappa_identity t.rho_range.1]
    have ht := within_trans t.k_error hs
    convert ht using 1 ; ring
  have ht := within_trans hk hp
  convert ht using 1 ; ring

/-- The first scaled-error proof uses the derivative/mean-value coefficient bound. -/
theorem factor_error {e rho : ℝ} (he : 0≤e)
    (hrho : rho∈Set.Icc (0:ℝ) 1) (t : FactorTrace e rho) :
    Within t.kHat (kappa rho) (8*e) := by
  have hp : Within (kappa t.rhoHat) (kappa rho) e:=
    le_trans (kappa_lipschitz t.rho_range.1 hrho.1) t.rho_error
  exact factor_error_from_coefficient he t hp

/-- The later executable proof uses its own algebraic coefficient argument. -/
theorem factor_error_algebraic {e rho : ℝ} (he : 0≤e)
    (hrho : rho∈Set.Icc (0:ℝ) 1) (t : FactorTrace e rho) :
    Within t.kHat (kappa rho) (8*e) := by
  have hp : Within (kappa t.rhoHat) (kappa rho) e:=
    le_trans (kappa_lipschitz_algebraic t.rho_range.1 hrho.1) t.rho_error
  exact factor_error_from_coefficient he t hp

/-- Local absolute-error contracts for the remaining scaled operations. -/
structure ScaledTrace (e rho p1 p2 p3 : ℝ) where
  factor : FactorTrace e rho
  (p1Hat p2Hat p3Hat sumHat hHat diffHat productH productP rightHat gapHat : ℝ)
  p1_range : |p1Hat| ≤ 1
  p2_range : |p2Hat| ≤ 1
  p3_range : |p3Hat| ≤ 1
  sum_range : |sumHat| ≤ 1
  h_range : |hHat| ≤ 1/2
  p1_error : Within p1Hat p1 e
  p2_error : Within p2Hat p2 e
  p3_error : Within p3Hat p3 e
  sum_error : Within sumHat (p1Hat+p2Hat) e
  h_error : Within hHat (sumHat/2) e
  diff_error : Within diffHat (p3Hat-hHat) e
  productH_error : Within productH (factor.kHat*|hHat|) e
  productP_error : Within productP (factor.kHat*|p3Hat|) e
  right_error : Within rightHat (productH+productP) e
  gap_error : Within gapHat (rightHat-|diffHat|) e

def exactH (p1 p2 : ℝ) : ℝ := (p1+p2)/2
def exactL (p1 p2 p3 : ℝ) : ℝ := |p3-exactH p1 p2|
def exactR (rho p1 p2 p3 : ℝ) : ℝ := kappa rho*(|exactH p1 p2|+|p3|)
def exactG (rho p1 p2 p3 : ℝ) : ℝ := exactR rho p1 p2 p3 - exactL p1 p2 p3

/-- The constants 8, 3, 5, 19 and 25 follow from the local contracts. -/
theorem scaled_errors_from_factor {e rho p1 p2 p3 : ℝ} (he : 0 ≤ e)
    (hh : |exactH p1 p2| ≤ 1/2) (hp3 : |p3| ≤ 1)
    (t : ScaledTrace e rho p1 p2 p3) (hk : Within t.factor.kHat (kappa rho) (8*e)) :
    Within t.factor.kHat (kappa rho) (8*e) ∧
    Within t.hHat (exactH p1 p2) (3*e) ∧
    Within |t.diffHat| (exactL p1 p2 p3) (5*e) ∧
    Within t.rightHat (exactR rho p1 p2 p3) (19*e) ∧
    Within t.gapHat (exactG rho p1 p2 p3) (25*e) := by
  have hsum := within_trans t.sum_error (within_add t.p1_error t.p2_error)
  have hmid := within_trans t.h_error (within_half hsum)
  have hm : Within t.hHat (exactH p1 p2) (3*e) := by
    exact within_mono hmid (by linarith)
  have hl₀ := within_trans t.diff_error (within_sub t.p3_error hm)
  have hl : Within |t.diffHat| (exactL p1 p2 p3) (5*e) := by
    have ha := within_abs hl₀
    convert ha using 1 <;> (try dsimp [exactL]) <;> ring
  have hkb : |t.factor.kHat| ≤ 1 := by
    rw [abs_of_nonneg (by linarith [t.factor.k_range.1] : 0 ≤ t.factor.kHat)]
    exact t.factor.k_range.2
  have pp₀ := product_error hkb (within_abs t.p3_error) hk (by simpa using hp3)
  have pp : Within t.productP (kappa rho*|p3|) (10*e) := by
    have ha := within_trans t.productP_error pp₀
    convert ha using 1 ; ring
  have ph₀ := product_error hkb (within_abs hm) hk (by simpa using hh)
  have ph : Within t.productH (kappa rho*|exactH p1 p2|) (8*e) := by
    have ha := within_trans t.productH_error ph₀
    convert ha using 1 ; ring
  have hright : Within t.rightHat (exactR rho p1 p2 p3) (19*e) := by
    have ha := within_trans t.right_error (within_add ph pp)
    convert ha using 1 <;> (try dsimp [exactR]) <;> ring
  have hgap : Within t.gapHat (exactG rho p1 p2 p3) (25*e) := by
    have ha := within_trans t.gap_error (within_sub hright hl)
    convert ha using 1 <;> (try dsimp [exactG]) <;> ring
  exact ⟨hk, hm, hl, hright, hgap⟩

theorem scaled_errors {e rho p1 p2 p3 : ℝ} (he : 0≤e)
    (hrho : rho∈Set.Icc (0:ℝ) 1) (hh : |exactH p1 p2|≤1/2) (hp3 : |p3|≤1)
    (t : ScaledTrace e rho p1 p2 p3) :
    Within t.factor.kHat (kappa rho) (8*e) ∧
    Within t.hHat (exactH p1 p2) (3*e) ∧
    Within |t.diffHat| (exactL p1 p2 p3) (5*e) ∧
    Within t.rightHat (exactR rho p1 p2 p3) (19*e) ∧
    Within t.gapHat (exactG rho p1 p2 p3) (25*e) :=
  scaled_errors_from_factor he hh hp3 t (factor_error he hrho t.factor)

theorem scaled_errors_algebraic {e rho p1 p2 p3 : ℝ} (he : 0≤e)
    (hrho : rho∈Set.Icc (0:ℝ) 1) (hh : |exactH p1 p2|≤1/2) (hp3 : |p3|≤1)
    (t : ScaledTrace e rho p1 p2 p3) :
    Within t.factor.kHat (kappa rho) (8*e) ∧
    Within t.hHat (exactH p1 p2) (3*e) ∧
    Within |t.diffHat| (exactL p1 p2 p3) (5*e) ∧
    Within t.rightHat (exactR rho p1 p2 p3) (19*e) ∧
    Within t.gapHat (exactG rho p1 p2 p3) (25*e) :=
  scaled_errors_from_factor he hh hp3 t (factor_error_algebraic he hrho t.factor)

theorem guard_positive {e G g : ℝ} (he : 0 < e)
    (herr : Within g G (25*e)) (hguard : 32*e < g) : 0 < G := by
  have h := (abs_le.mp herr).2
  linarith

theorem guard_negative {e G g : ℝ} (he : 0 < e)
    (herr : Within g G (25*e)) (hguard : g < -(32*e)) : G < 0 := by
  have h := (abs_le.mp herr).1
  linarith

theorem complete_positive {e G g : ℝ}
    (herr : Within g G (25*e)) (hmargin : 57*e < G) : 32*e < g := by
  have h := (abs_le.mp herr).1
  linarith

theorem complete_negative {e G g : ℝ}
    (herr : Within g G (25*e)) (hmargin : G < -(57*e)) : g < -(32*e) := by
  have h := (abs_le.mp herr).2
  linarith

theorem certified_switch {e rho p1 p2 p3 : ℝ} (he : 0 < e)
    (hrho : rho ∈ Set.Icc (0:ℝ) 1)
    (hh : |exactH p1 p2| ≤ 1/2) (hp3 : |p3| ≤ 1)
    (t : ScaledTrace e rho p1 p2 p3) (hguard : 32*e < t.gapHat) :
    exactL p1 p2 p3 < exactR rho p1 p2 p3 := by
  have herr := (scaled_errors (le_of_lt he) hrho hh hp3 t).2.2.2.2
  have h := guard_positive he herr hguard
  dsimp [exactG] at h
  linarith

theorem certified_no_switch {e rho p1 p2 p3 : ℝ} (he : 0 < e)
    (hrho : rho ∈ Set.Icc (0:ℝ) 1)
    (hh : |exactH p1 p2| ≤ 1/2) (hp3 : |p3| ≤ 1)
    (t : ScaledTrace e rho p1 p2 p3) (hguard : t.gapHat < -(32*e)) :
    exactR rho p1 p2 p3 < exactL p1 p2 p3 := by
  have herr := (scaled_errors (le_of_lt he) hrho hh hp3 t).2.2.2.2
  have h := guard_negative he herr hguard
  dsimp [exactG] at h
  linarith

def binary64Epsilon : ℝ := 1 / (2:ℝ)^53

theorem binary64Epsilon_pos : 0 < binary64Epsilon := by norm_num [binary64Epsilon]

theorem binary64_guard_value : 32*binary64Epsilon = 1/(2:ℝ)^48 := by
  norm_num [binary64Epsilon]

def originalL (f1 f2 f3 : ℝ) : ℝ := |f3-(f1+f2)/2|
def originalR (k f1 f2 f3 : ℝ) : ℝ := k*(|(f1+f2)/2|+|f3|)

theorem normalization_identity {s : ℝ} (hs : 0 < s) (k f1 f2 f3 : ℝ) :
    k*(|(f1/s+f2/s)/2|+|f3/s|)-|f3/s-(f1/s+f2/s)/2| =
      (originalR k f1 f2 f3-originalL f1 f2 f3)/s := by
  have hm : (f1/s+f2/s)/2 = ((f1+f2)/2)/s := by ring
  have hl : f3/s - ((f1+f2)/2)/s = (f3-(f1+f2)/2)/s := by ring
  rw [hm, hl]
  dsimp [originalR, originalL]
  simp only [abs_div, abs_of_pos hs]
  ring

theorem normalization_sign {s : ℝ} (hs : 0 < s) (k f1 f2 f3 : ℝ) :
    0 < k*(|(f1/s+f2/s)/2|+|f3/s|)-|f3/s-(f1/s+f2/s)/2| ↔
      originalL f1 f2 f3 < originalR k f1 f2 f3 := by
  rw [normalization_identity hs, div_pos_iff_of_pos_right hs, sub_pos]

theorem endpoint_ratio {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    |(-a+b)/2|/(a+b) =
      (1-min a b/max a b)/(2*(1+min a b/max a b)) := by
  rcases le_total a b with h | h
  · rw [min_eq_left h, max_eq_right h, abs_of_nonneg (by linarith : 0 ≤ (-a+b)/2)]
    field_simp
    ; ring
  · rw [min_eq_right h, max_eq_left h, abs_of_nonpos (by linarith : (-a+b)/2 ≤ 0)]
    field_simp
    ; ring

theorem endpoint_rho_range {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    min a b/max a b ∈ Set.Icc (0:ℝ) 1 := by
  have hp : 0 < max a b := lt_of_lt_of_le ha (le_max_left _ _)
  constructor
  · apply div_nonneg
    · exact le_min (le_of_lt ha) (le_of_lt hb)
    · exact le_of_lt hp
  · apply (div_le_iff₀ hp).mpr
    simpa only [one_mul] using le_trans (min_le_left a b) (le_max_left a b)

theorem endpoint_factor {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (1-|(-a+b)/2|/(a+b))^2 = kappa (min a b/max a b) := by
  rw [endpoint_ratio ha hb, kappa_identity (endpoint_rho_range ha hb).1]

/-- Both endpoint orientations are admitted; the midpoint may be zero. -/
theorem normalized_midpoint_bound {p1 p2 : ℝ}
    (h1 : |p1| ≤ 1) (h2 : |p2| ≤ 1)
    (hsign : (p1 ≤ 0 ∧ 0 ≤ p2) ∨ (p2 ≤ 0 ∧ 0 ≤ p1)) :
    |exactH p1 p2| ≤ 1/2 := by
  rcases abs_le.mp h1 with ⟨h1l,h1u⟩
  rcases abs_le.mp h2 with ⟨h2l,h2u⟩
  apply abs_le.mpr
  dsimp [exactH]
  rcases hsign with ⟨ha,hb⟩ | ⟨ha,hb⟩ <;> constructor <;> linarith

theorem scaled_original_switch {e rho f1 f2 f3 s : ℝ} (he : 0 < e) (hs : 0 < s)
    (hrho : rho ∈ Set.Icc (0:ℝ) 1)
    (hh : |exactH (f1/s) (f2/s)| ≤ 1/2) (hp3 : |f3/s| ≤ 1)
    (t : ScaledTrace e rho (f1/s) (f2/s) (f3/s)) (hguard : 32*e < t.gapHat) :
    originalL f1 f2 f3 < originalR (kappa rho) f1 f2 f3 := by
  apply (normalization_sign hs (kappa rho) f1 f2 f3).mp
  have h := certified_switch he hrho hh hp3 t hguard
  dsimp [exactL, exactR, exactH] at h
  linarith

theorem subnormal_exact_scaled_gap :
    exactG (1/2) (-1/2) 1 (1/2) = (13/48:ℝ) := by
  norm_num [exactG, exactR, exactL, exactH, kappa, abs_of_nonneg]

theorem subnormal_margin : (57:ℝ)*binary64Epsilon < 13/48 := by
  norm_num [binary64Epsilon]

theorem symmetric_exact_gap (z : ℝ) : exactG 1 (-1) 1 z = 0 := by
  norm_num [exactG, exactR, exactL, exactH, kappa]

def inputScale (a b z : ℝ) : ℝ := max (max a b) |z|

theorem inputScale_pos {a b z : ℝ} (ha : 0 < a) : 0 < inputScale a b z := by
  exact lt_of_lt_of_le ha (le_trans (le_max_left a b) (le_max_left (max a b) |z|))

theorem normalized_input_bounds {a b z : ℝ} (ha : 0 < a) (hb : 0 < b) :
    |exactH (-a/inputScale a b z) (b/inputScale a b z)| ≤ 1/2 ∧
    |z/inputScale a b z| ≤ 1 := by
  let s := inputScale a b z
  have hs : 0 < s := inputScale_pos ha
  have haS : a ≤ s := le_trans (le_max_left a b) (le_max_left _ _)
  have hbS : b ≤ s := le_trans (le_max_right a b) (le_max_left _ _)
  have hzS : |z| ≤ s := le_max_right _ _
  have h1 : |-a/s| ≤ 1 := by
    rw [abs_div, abs_neg, abs_of_pos ha, abs_of_pos hs]
    exact (div_le_one hs).mpr haS
  have h2 : |b/s| ≤ 1 := by
    rw [abs_div, abs_of_pos hb, abs_of_pos hs]
    exact (div_le_one hs).mpr hbS
  constructor
  · apply normalized_midpoint_bound h1 h2
    left
    constructor
    · exact div_nonpos_of_nonpos_of_nonneg (by linarith) (le_of_lt hs)
    · exact div_nonneg (le_of_lt hb) (le_of_lt hs)
  · rw [abs_div, abs_of_pos hs]
    exact (div_le_one hs).mpr hzS

/-- The guard implies the original condition, with the symmetry factor
defined directly from the two opposite-sign endpoint residuals. -/
theorem certified_from_residuals {e a b z : ℝ} (he : 0 < e)
    (ha : 0 < a) (hb : 0 < b)
    (t : ScaledTrace e (min a b/max a b)
      (-a/inputScale a b z) (b/inputScale a b z) (z/inputScale a b z))
    (hguard : 32*e < t.gapHat) :
    |z-(-a+b)/2| < (1-|(-a+b)/2|/(a+b))^2*(|(-a+b)/2|+|z|) := by
  rw [endpoint_factor ha hb]
  have hbds := normalized_input_bounds (z := z) ha hb
  exact scaled_original_switch he (inputScale_pos ha) (endpoint_rho_range ha hb)
    hbds.1 hbds.2 t hguard

theorem scaled_margin_guarantees {e rho p1 p2 p3 : ℝ} (he : 0 ≤ e)
    (hrho : rho ∈ Set.Icc (0:ℝ) 1)
    (hh : |exactH p1 p2| ≤ 1/2) (hp3 : |p3| ≤ 1)
    (t : ScaledTrace e rho p1 p2 p3) :
    (57*e < exactG rho p1 p2 p3 → 32*e < t.gapHat) ∧
    (exactG rho p1 p2 p3 < -(57*e) → t.gapHat < -(32*e)) := by
  have h := (scaled_errors he hrho hh hp3 t).2.2.2.2
  exact ⟨complete_positive h, complete_negative h⟩

def gamma (j : ℕ) : ℝ := (j:ℝ)*binary64Epsilon/(1-(j:ℝ)*binary64Epsilon)

/-- Exact arithmetic check of the final coefficient in the unscaled proof. -/
theorem unscaled_coefficient :
    (3+binary64Epsilon)/2*gamma 2 + 3*gamma 3 +
      (5*binary64Epsilon+binary64Epsilon^2)/2 < 16*binary64Epsilon := by
  norm_num [gamma, binary64Epsilon]

/-- This result composes already propagated side bounds. Unlike
`scaled_errors`, it does not derive them from individual rounding steps. -/
theorem unscaled_error_budget {L R l r H : ℝ} (hH : 0 ≤ H)
    (hl : Within l L ((2*binary64Epsilon+binary64Epsilon^2/2)*H))
    (hr : Within r R (((3+binary64Epsilon)/2*gamma 2 +
      3*gamma 3 + binary64Epsilon/2)*H)) :
    Within (r-l) (R-L) (16*binary64Epsilon*H) := by
  have h := within_sub hr hl
  apply within_mono h
  have hc := mul_le_mul_of_nonneg_right (le_of_lt unscaled_coefficient) hH
  nlinarith only [hc]

end
end ModAB

end RModABBundle0

/- Original module: ModAB.Phase -/
section RModABBundle1

/-! The width-policy projection of a continuing AB phase. The exception constructor
records the cap; requiring residual progress can only remove permitted executions.
Neither the geometric envelope nor a convergence conclusion is a field of the model. -/
namespace ModAB.Convergence
noncomputable section

def threshold (B S : ℝ) (t : ℕ) : ℝ := 2*B*S/(2:ℝ)^t
def envelopeConstant (B : ℝ) (M : ℕ) : ℝ := (2:ℝ)^(M+1)*B

inductive ActivePhase (B : ℝ) (M : ℕ) (S : ℝ) : ℕ → ℕ → ℝ → Prop
  | start : ActivePhase B M S 0 0 S
  | width {t c : ℕ} {w w' : ℝ} :
      ActivePhase B M S t c w → 0≤w' → w'≤w → w'≤threshold B S (t+1) →
      ActivePhase B M S (t+1) 0 w'
  | exception {t c : ℕ} {w w' : ℝ} :
      ActivePhase B M S t c w → 0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      ActivePhase B M S (t+1) (c+1) w'

theorem active_nonnegative {B S w : ℝ} {M t c : ℕ} (hS : 0≤S)
    (h : ActivePhase B M S t c w) : 0≤w := by
  induction h with
  | start => exact hS
  | width _ hw _ _ _ => exact hw
  | exception _ hw _ _ _ _ => exact hw

theorem active_nonincrease {B S w : ℝ} {M t c : ℕ}
    (h : ActivePhase B M S t c w) : w≤S := by
  induction h with
  | start => exact le_refl _
  | width _ _ hw _ ih => exact hw.trans ih
  | exception _ _ hw _ _ ih => exact hw.trans ih

/-- The manuscript's last successful width test j=t-c, including j=0. -/
theorem last_width_success {B S w : ℝ} {M t c : ℕ} (hB : (1/2:ℝ)≤B) (hS : 0≤S)
    (h : ActivePhase B M S t c w) :
    ∃ j : ℕ, j+c=t ∧ c≤M ∧ w≤threshold B S j := by
  induction h with
  | start =>
      refine ⟨0,rfl,Nat.zero_le _,?_⟩
      simp only [threshold,pow_zero,div_one]
      nlinarith [mul_nonneg hS (show 0≤2*B-1 by linarith)]
  | @width t c w w' old hw hmono htest ih =>
      exact ⟨t+1,by omega,Nat.zero_le _,htest⟩
  | @exception t c w w' old hw hmono htest hcap ih =>
      obtain ⟨j,hj,hc,hbound⟩:=ih
      exact ⟨j,by omega,by omega,hmono.trans hbound⟩

theorem threshold_shift (B S : ℝ) (j c : ℕ) :
    threshold B S j=(2:ℝ)^c*threshold B S (j+c) := by
  unfold threshold
  rw [pow_add]
  field_simp

/-- `lem:adaptive_envelope`, following its last-successful-test argument. -/
theorem active_envelope {B S w : ℝ} {M t c : ℕ} (hB : (1/2:ℝ)≤B) (hS : 0≤S)
    (h : ActivePhase B M S t c w) : w≤envelopeConstant B M*S/(2:ℝ)^t := by
  obtain ⟨j,hj,hc,hw⟩:=last_width_success hB hS h
  have hp : (2:ℝ)^c≤(2:ℝ)^M := pow_le_pow_right₀ (by norm_num) hc
  have ht : 0≤threshold B S t := by unfold threshold; positivity
  calc
    w≤threshold B S j := hw
    _=(2:ℝ)^c*threshold B S t := by rw [threshold_shift B S j c,hj]
    _≤(2:ℝ)^M*threshold B S t := mul_le_mul_of_nonneg_right hp ht
    _=envelopeConstant B M*S/(2:ℝ)^t := by unfold threshold envelopeConstant; rw [pow_succ]; ring

theorem C_at_least_one {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) : 1≤envelopeConstant B M := by
  have hp : (1:ℝ)≤2^M := one_le_pow₀ (by norm_num)
  unfold envelopeConstant
  rw [pow_succ]
  nlinarith

theorem C32 : envelopeConstant 2 3=32 := by norm_num [envelopeConstant]

theorem envelope_to_dyadic {B S w : ℝ} {M t c d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : ActivePhase B M S t c w) (ht : d≤t) : w≤S/(2:ℝ)^(t-d) := by
  have he := active_envelope hB hS h
  have hm := div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_right hC hS)
    (by positivity : 0≤(2:ℝ)^t)
  have ht' : t=d+(t-d) := by omega
  calc
    w≤envelopeConstant B M*S/(2:ℝ)^t := he
    _≤(2:ℝ)^d*S/(2:ℝ)^t := hm
    _=S/(2:ℝ)^(t-d) := by conv_lhs => rw [ht',pow_add]; field_simp

def phaseCredit (d t : ℕ) : ℕ := if t≤d+1 then 1 else t-d

/-- Completed phase: t−1 continuing steps, one last AB update, then bisection. -/
theorem completed_phase {B S u v : ℝ} {M t c d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (ht : 1≤t) (hphase : ActivePhase B M S (t-1) c u) (huv : v≤u) :
    v/2≤S/(2:ℝ)^(phaseCredit d t) ∧ t+1≤(d+2)*phaseCredit d t := by
  by_cases hshort : t≤d+1
  · simp only [phaseCredit,if_pos hshort,pow_one]
    exact ⟨div_le_div_of_nonneg_right (huv.trans (active_nonincrease hphase)) (by norm_num),by omega⟩
  · have hlong : d≤t-1 := by omega
    have hw := envelope_to_dyadic hB hS hC hphase hlong
    have he : t-d=(t-1-d)+1 := by omega
    have hratio : (S/(2:ℝ)^(t-1-d))/2=S/(2:ℝ)^(t-d) := by rw [he,pow_succ]; ring
    simp only [phaseCredit,if_neg hshort]
    constructor
    · exact (div_le_div_of_nonneg_right (huv.trans hw) (by norm_num)).trans_eq hratio
    · have hh : 2≤t-d := by omega
      have hm := Nat.mul_le_mul_left (d+1) (show 1≤t-d by omega)
      have htde : t=d+(t-d) := by omega
      nlinarith

inductive Tail (B : ℝ) (M : ℕ) (S : ℝ) : ℕ → ℝ → Prop
  | active {t c w} : ActivePhase B M S t c w → Tail B M S t w
  | stopped {t c u w} : ActivePhase B M S t c u → 0≤w → w≤u → Tail B M S (t+1) w

theorem tail_nonnegative {B S w : ℝ} {M t : ℕ} (hS : 0≤S)
    (h : Tail B M S t w) : 0≤w := by
  cases h with
  | active h => exact active_nonnegative hS h
  | stopped _ hw _ => exact hw

/-- `lem:tail`: the same four cases (active/stopped, short/long). -/
theorem tail_bound {B S w : ℝ} {M t d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : Tail B M S t w) :
    ∃ k : ℕ, w≤S/(2:ℝ)^k ∧ t≤(d+2)*k+(d+1) := by
  cases h with
  | @active t c w hp =>
      by_cases ht : t≤d
      · exact ⟨0,by simpa using active_nonincrease hp,by omega⟩
      · refine ⟨t-d,envelope_to_dyadic hB hS hC hp (by omega),?_⟩
        have hm := Nat.mul_le_mul_left (d+1) (show 1≤t-d by omega)
        have htde : t=d+(t-d) := by omega
        nlinarith
  | @stopped t c u w hp hw hwu =>
      by_cases ht : t≤d
      · exact ⟨0,by simpa using hwu.trans (active_nonincrease hp),by omega⟩
      · refine ⟨t-d,hwu.trans (envelope_to_dyadic hB hS hC hp (by omega)),?_⟩
        have hm := Nat.zero_le ((d+1)*(t-d))
        have htde : t=d+(t-d) := by omega
        nlinarith

end
end ModAB.Convergence

end RModABBundle1

/- Original module: ModAB.Schedule -/
section RModABBundle2

/-! Exact-arithmetic width histories, their decomposition, and the manuscript's
amortized global estimate. The model permits earlier fallback and an arbitrary
switching test; the residual test only restricts its exception transitions.
-/
namespace ModAB.Convergence
open Filter Topology
noncomputable section

inductive Mode where
  | bisection
  | ab (S : ℝ) (t c : ℕ)

/-- One constructor per permitted width-policy transition. -/
inductive History (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Mode → Prop
  | initial : History B M W0 0 W0 .bisection
  | bisectStay {n w} : History B M W0 n w .bisection →
      History B M W0 (n+1) (w/2) .bisection
  | bisectStart {n w} : History B M W0 n w .bisection →
      History B M W0 (n+1) (w/2) (.ab (w/2) 0 0)
  | abWidth {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → w'≤threshold B S (t+1) →
      History B M W0 (n+1) w' (.ab S (t+1) 0)
  | abException {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      History B M W0 (n+1) w' (.ab S (t+1) (c+1))
  | abStop {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → History B M W0 (n+1) w' .bisection

/-- Completed parts: isolated bisections and AB blocks ending in bisection. -/
inductive CompletedPrefix (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Prop
  | initial : CompletedPrefix B M W0 0 W0
  | bisect {n w} : CompletedPrefix B M W0 n w → CompletedPrefix B M W0 (n+1) (w/2)
  | phase {n w t c u v} : CompletedPrefix B M W0 n w →
      ActivePhase B M w t c u → 0≤v → v≤u →
      CompletedPrefix B M W0 (n+t+2) (v/2)

theorem prefix_nonnegative {B W0 w : ℝ} {M n : ℕ} (hW0 : 0≤W0)
    (h : CompletedPrefix B M W0 n w) : 0≤w := by
  induction h with
  | initial => exact hW0
  | bisect _ ih => positivity
  | phase _ _ hv _ _ => positivity

inductive Decomposed (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Mode → Prop
  | ready {n w} : CompletedPrefix B M W0 n w → Decomposed B M W0 n w .bisection
  | active {n S t c w} : CompletedPrefix B M W0 n S → ActivePhase B M S t c w →
      Decomposed B M W0 (n+t) w (.ab S t c)
  | stopped {n S t c u w} : CompletedPrefix B M W0 n S → ActivePhase B M S t c u →
      0≤w → w≤u → Decomposed B M W0 (n+t+1) w .bisection

/-- The partition in the manuscript is derived from transitions, not assumed. -/
theorem history_decomposes {B W0 w : ℝ} {M n : ℕ} {mode : Mode}
    (h : History B M W0 n w mode) : Decomposed B M W0 n w mode := by
  induction h with
  | initial => exact .ready .initial
  | bisectStay old ih =>
      cases ih with
      | ready hp => exact .ready (.bisect hp)
      | stopped hp ha hw hmono =>
          simpa only [Nat.add_assoc] using Decomposed.ready (CompletedPrefix.phase hp ha hw hmono)
  | bisectStart old ih =>
      cases ih with
      | ready hp => simpa using Decomposed.active (CompletedPrefix.bisect hp) ActivePhase.start
      | stopped hp ha hw hmono =>
          simpa only [Nat.add_zero,Nat.add_assoc] using
            Decomposed.active (CompletedPrefix.phase hp ha hw hmono) ActivePhase.start
  | abWidth old hw hmono ht ih =>
      cases ih with
      | active hp ha =>
          simpa only [Nat.add_assoc] using Decomposed.active hp (ActivePhase.width ha hw hmono ht)
  | abException old hw hmono ht hc ih =>
      cases ih with
      | active hp ha =>
          simpa only [Nat.add_assoc] using Decomposed.active hp (ActivePhase.exception ha hw hmono ht hc)
  | abStop old hw hmono ih =>
      cases ih with
      | active hp ha => exact .stopped hp ha hw hmono

theorem dyadic_compose (W0 : ℝ) (h k : ℕ) :
    (W0/(2:ℝ)^h)/(2:ℝ)^k=W0/(2:ℝ)^(h+k) := by rw [pow_add]; ring

/-- Sum of the credits of completed blocks; same induction as summing their bounds. -/
theorem prefix_credits {B W0 w : ℝ} {M n d : ℕ}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : CompletedPrefix B M W0 n w) :
    ∃ H : ℕ, w≤W0/(2:ℝ)^H ∧ n≤(d+2)*H := by
  induction h with
  | initial => exact ⟨0,by simp,by simp⟩
  | @bisect n w old ih =>
      obtain ⟨H,hw,hn⟩:=ih
      refine ⟨H+1,?_,?_⟩
      · simpa only [pow_one] using
          (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by simpa using dyadic_compose W0 H 1)
      · nlinarith
  | @phase n w t c u v old ha hv hmono ih =>
      obtain ⟨H,hw,hn⟩:=ih
      have hphase : ActivePhase B M w ((t+1)-1) c u := by simpa using ha
      obtain ⟨hb,hcost⟩:=completed_phase hB (prefix_nonnegative hW0 old) hC
        (by omega : 1≤t+1) hphase hmono
      let k:=phaseCredit d (t+1)
      refine ⟨H+k,?_,?_⟩
      · exact hb.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · change n+t+2≤(d+2)*(H+phaseCredit d (t+1))
        nlinarith

theorem history_credits {B W0 w : ℝ} {M n d : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : History B M W0 n w mode) :
    ∃ H : ℕ, w≤W0/(2:ℝ)^H ∧ n≤(d+2)*H+(d+1) := by
  have hd:=history_decomposes h
  cases hd with
  | ready hp =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      exact ⟨H,hw,by omega⟩
  | @active n S t c w hp ha =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      obtain ⟨k,hk,ht⟩:=tail_bound hB (prefix_nonnegative hW0 hp) hC (Tail.active ha)
      refine ⟨H+k,?_,?_⟩
      · exact hk.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · nlinarith
  | @stopped n S t c u w hp ha hw0 hmono =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      obtain ⟨k,hk,ht⟩:=tail_bound hB (prefix_nonnegative hW0 hp) hC (Tail.stopped ha hw0 hmono)
      refine ⟨H+k,?_,?_⟩
      · exact hk.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · nlinarith

/-- Integer step in `thm:quantitative`: n≤qH+q−1 implies floor(n/q)≤H. -/
theorem credit_floor {n q H : ℕ} (hq : 0<q) (hn : n≤q*H+(q-1)) : n/q≤H := by
  apply Nat.le_of_lt_succ
  apply (Nat.div_lt_iff_lt_mul hq).mpr
  have hq' : q-1+1=q := by omega
  nlinarith

theorem dyadic_antitone {W0 : ℝ} (hW0 : 0≤W0) {h k : ℕ} (hh : h≤k) :
    W0/(2:ℝ)^k≤W0/(2:ℝ)^h := by
  exact div_le_div_of_nonneg_left hW0 (by positivity)
    (pow_le_pow_right₀ (by norm_num) hh)

/-- `thm:quantitative`; natural-number division is the displayed floor. -/
theorem global_width_bound {B W0 w : ℝ} {M n d : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : History B M W0 n w mode) : w≤W0/(2:ℝ)^(n/(d+2)) := by
  obtain ⟨H,hw,hn⟩:=history_credits hB hW0 hC h
  have hc : n/(d+2)≤H := credit_floor (by omega) (by simpa using hn)
  exact hw.trans (dyadic_antitone hW0 hc)

def dyadicExponent (C : ℝ) : ℕ := Nat.ceil (Real.logb 2 C)

theorem ceil_log_bound {C : ℝ} (hC : 0<C) : C≤(2:ℝ)^(dyadicExponent C) := by
  have h := (Real.logb_le_iff_le_rpow (by norm_num : (1:ℝ)<2) hC).mp
    (Nat.le_ceil (Real.logb 2 C))
  simpa only [Real.rpow_natCast, dyadicExponent] using h

theorem global_width_bound_log {B W0 w : ℝ} {M n : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (h : History B M W0 n w mode) :
    w≤W0/(2:ℝ)^(n/(dyadicExponent (envelopeConstant B M)+2)) :=
  global_width_bound hB hW0 (ceil_log_bound (lt_of_lt_of_le (by norm_num) (C_at_least_one hB))) h

/-- `cor:C32`, without a numerical approximation to logarithms. -/
theorem C32_global_width {W0 w : ℝ} {n : ℕ} {mode : Mode} (hW0 : 0≤W0)
    (h : History 2 3 W0 n w mode) : w≤W0/(2:ℝ)^(n/7) := by
  exact global_width_bound (d:=5) (by norm_num) hW0 (by norm_num [envelopeConstant]) h

theorem finite_termination {B W0 w τ : ℝ} {M d K : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (hK : W0/(2:ℝ)^K≤τ) (h : History B M W0 ((d+2)*K) w mode) : w≤τ := by
  have hb:=global_width_bound hB hW0 hC h
  have hdiv : ((d+2)*K)/(d+2)=K := Nat.mul_div_right K (by omega)
  rw [hdiv] at hb
  exact hb.trans hK

def toleranceExponent (W0 τ : ℝ) : ℕ := Nat.ceil (Real.logb 2 (W0/τ))

theorem tolerance_exponent_sufficient {W0 τ : ℝ} (hW0 : 0≤W0) (hτ : 0<τ) :
    W0/(2:ℝ)^(toleranceExponent W0 τ)≤τ := by
  by_cases hz : W0=0
  · simp [hz,hτ.le]
  · have hW : 0<W0 := lt_of_le_of_ne hW0 (Ne.symm hz)
    have hp := ceil_log_bound (div_pos hW hτ)
    change W0/τ≤(2:ℝ)^(toleranceExponent W0 τ) at hp
    apply (div_le_iff₀ (by positivity : 0<(2:ℝ)^(toleranceExponent W0 τ))).mpr
    exact (div_le_iff₀ hτ).mp hp |>.trans_eq (mul_comm _ _)

end
end ModAB.Convergence

end RModABBundle2

/- Original module: ModAB.Global -/
section RModABBundle3

/-! Infinite exact-arithmetic runs. The qualitative proof deliberately follows
both bisection cases in `thm:global_simple`, independently of amortized counting. -/
namespace ModAB.Convergence
open Filter Topology
noncomputable section

inductive Step (B : ℝ) (M : ℕ) : ℝ → Mode → ℝ → Mode → Prop
  | bisectStay (w) : Step B M w .bisection (w/2) .bisection
  | bisectStart (w) : Step B M w .bisection (w/2) (.ab (w/2) 0 0)
  | abWidth {w w' S t c} : 0≤w' → w'≤w → w'≤threshold B S (t+1) →
      Step B M w (.ab S t c) w' (.ab S (t+1) 0)
  | abException {w w' S t c} : 0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      Step B M w (.ab S t c) w' (.ab S (t+1) (c+1))
  | abStop {w w' S t c} : 0≤w' → w'≤w → Step B M w (.ab S t c) w' .bisection

structure Run (B : ℝ) (M : ℕ) where
  width : ℕ→ℝ
  mode : ℕ→Mode
  nonneg_initial : 0≤width 0
  initial : mode 0=.bisection
  next : ∀n,Step B M (width n) (mode n) (width (n+1)) (mode (n+1))

theorem step_history {B W0 w w' : ℝ} {M n : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (h : History B M W0 n w m) :
    History B M W0 (n+1) w' m' := by
  cases hs with
  | bisectStay => exact .bisectStay h
  | bisectStart => exact .bisectStart h
  | abWidth h0 hm ht => exact .abWidth h h0 hm ht
  | abException h0 hm ht hc => exact .abException h h0 hm ht hc
  | abStop h0 hm => exact .abStop h h0 hm

theorem step_nonnegative {B w w' : ℝ} {M : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (hw : 0≤w) : 0≤w' := by
  cases hs with
  | bisectStay => positivity
  | bisectStart => positivity
  | abWidth h0 _ _ => exact h0
  | abException h0 _ _ _ => exact h0
  | abStop h0 _ => exact h0

theorem history_of_run {B : ℝ} {M : ℕ} (r : Run B M) (n : ℕ) :
    History B M (r.width 0) n (r.width n) (r.mode n) := by
  induction n with
  | zero => rw [r.initial]; exact .initial
  | succ n ih => exact step_history (r.next n) ih

theorem run_nonnegative {B : ℝ} {M : ℕ} (r : Run B M) : ∀n,0≤r.width n := by
  intro n
  induction n with
  | zero => exact r.nonneg_initial
  | succ n ih => exact step_nonnegative (r.next n) ih

def isBisection : Mode→Bool
  | .bisection=>true
  | .ab _ _ _=>false

def bisectionCount (m : ℕ→Mode) : ℕ→ℕ
  | 0=>0
  | n+1=>bisectionCount m n + if isBisection (m n) then 1 else 0

theorem count_monotone (m : ℕ→Mode) : Monotone (bisectionCount m) := by
  apply monotone_nat_of_le_succ
  intro n
  simp only [bisectionCount]
  omega

theorem step_count_bound {B W0 w w' : ℝ} {M k : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (hw : w≤W0/(2:ℝ)^k) :
    w'≤W0/(2:ℝ)^(k+if isBisection m then 1 else 0) := by
  cases hs with
  | bisectStay =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | bisectStart =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | abWidth _ hm _ => simpa [isBisection] using hm.trans hw
  | abException _ hm _ _ => simpa [isBisection] using hm.trans hw
  | abStop _ hm => simpa [isBisection] using hm.trans hw

/-- Each bisection halves the width; all other steps only decrease it. -/
theorem bisection_count_bound {B : ℝ} {M : ℕ} (r : Run B M) (n : ℕ) :
    r.width n≤r.width 0/(2:ℝ)^(bisectionCount r.mode n) := by
  induction n with
  | zero => simp [bisectionCount]
  | succ n ih => exact step_count_bound (r.next n) ih

theorem count_unbounded {m : ℕ→Mode}
    (h : ∀N,∃n,N≤n ∧ isBisection (m n)=true) : ∀k,∃n,k≤bisectionCount m n := by
  intro k
  induction k with
  | zero => exact ⟨0,by simp [bisectionCount]⟩
  | succ k ih =>
      obtain ⟨n,hn⟩:=ih
      obtain ⟨j,hnj,hj⟩:=h n
      refine ⟨j+1,?_⟩
      have hm:=count_monotone m hnj
      simp only [bisectionCount,hj,ite_true]
      omega

theorem dyadic_tendsto (S : ℝ) : Tendsto (fun n:ℕ=>S/(2:ℝ)^n) atTop (𝓝 0) := by
  have h:= (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num : (0:ℝ)≤1/2)
    (by norm_num : (1/2:ℝ)<1)).const_mul S
  simpa only [mul_zero,div_pow,one_pow,div_eq_mul_inv,one_mul,inv_pow] using h

/-- Case 1 of the manuscript's qualitative proof. -/
theorem infinite_bisections {B : ℝ} {M : ℕ} (r : Run B M)
    (h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  have hc : Tendsto (bisectionCount r.mode) atTop atTop := by
    apply tendsto_atTop_atTop.mpr
    intro k
    obtain ⟨n,hn⟩:=count_unbounded h k
    exact ⟨n,fun j hj=>hn.trans (count_monotone r.mode hj)⟩
  exact squeeze_zero (run_nonnegative r) (bisection_count_bound r) ((dyadic_tendsto _).comp hc)

theorem active_step {B S w w' : ℝ} {M t c : ℕ} {m' : Mode}
    (hs : Step B M w (.ab S t c) w' m') (ha : ActivePhase B M S t c w)
    (hno : isBisection m'=false) :
    ∃c',m'=.ab S (t+1) c' ∧ ActivePhase B M S (t+1) c' w' := by
  cases hs with
  | abWidth h0 hm ht => exact ⟨0,rfl,.width ha h0 hm ht⟩
  | abException h0 hm ht hc => exact ⟨c+1,rfl,.exception ha h0 hm ht hc⟩
  | abStop => simp [isBisection] at hno

theorem active_history {B W0 w : ℝ} {M n : ℕ} {m : Mode}
    (hW0 : 0≤W0) (h : History B M W0 n w m) (hno : isBisection m=false) :
    ∃S t c,m=.ab S t c ∧ 0≤S ∧ ActivePhase B M S t c w := by
  cases history_decomposes h with
  | ready => simp [isBisection] at hno
  | stopped => simp [isBisection] at hno
  | active hp ha => exact ⟨_,_,_,rfl,prefix_nonnegative hW0 hp,ha⟩

/-- In an infinite suffix without bisection, the AB phase cannot end. -/
theorem continuing_suffix {B S : ℝ} {M N t c : ℕ} (r : Run B M)
    (hm : r.mode N=.ab S t c) (ha : ActivePhase B M S t c (r.width N))
    (hno : ∀n,N≤n→isBisection (r.mode n)=false) :
    ∀k,∃c',r.mode (N+k)=.ab S (t+k) c' ∧
      ActivePhase B M S (t+k) c' (r.width (N+k)) := by
  intro k
  induction k with
  | zero => exact ⟨c,by simpa using hm,by simpa using ha⟩
  | succ k ih =>
      obtain ⟨c',hm',ha'⟩:=ih
      have hs:=r.next (N+k)
      rw [hm'] at hs
      obtain ⟨c'',he,hactive⟩:=active_step hs ha' (hno (N+k+1) (by omega))
      exact ⟨c'',by simpa [Nat.add_assoc] using he,by simpa [Nat.add_assoc] using hactive⟩

/-- Case 2: extract the remaining phase from the history and apply its envelope. -/
theorem finite_bisections {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M)
    (h : ¬∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  push_neg at h
  obtain ⟨N,hN⟩:=h
  have hno : ∀n,N≤n→isBisection (r.mode n)=false := by
    intro n hn
    exact Bool.eq_false_iff.mpr (hN n hn)
  obtain ⟨S,t,c,hm,hS,ha⟩:=active_history r.nonneg_initial (history_of_run r N) (hno N (le_refl N))
  have hs:=continuing_suffix r hm ha hno
  have hlim : Tendsto (fun k=>r.width (N+k)) atTop (𝓝 0) := by
    apply squeeze_zero (fun k=>run_nonnegative r (N+k))
      (fun k=>?_) (dyadic_tendsto (envelopeConstant B M*S/(2:ℝ)^t))
    obtain ⟨c',hm',ha'⟩:=hs k
    have hb:=active_envelope hB hS ha'
    simpa only [pow_add,div_mul_eq_div_div] using hb
  exact (tendsto_add_atTop_iff_nat N).mp (by simpa [Nat.add_comm] using hlim)

/-- `thm:global_simple`, width part, with precisely the paper's two cases. -/
theorem width_tendsto_zero {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M) :
    Tendsto r.width atTop (𝓝 0) := by
  by_cases h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true
  · exact infinite_bisections r h
  · exact finite_bisections hB r h

/-- Independent consistency check using the quantitative theorem. -/
theorem width_tendsto_zero_quantitative {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M) :
    Tendsto r.width atTop (𝓝 0) := by
  let d:=dyadicExponent (envelopeConstant B M)
  have he : Tendsto (fun n:ℕ=>n/(d+2)) atTop atTop:=Nat.tendsto_div_const_atTop (by omega)
  exact squeeze_zero (run_nonnegative r)
    (fun n=>global_width_bound_log hB r.nonneg_initial (history_of_run r n))
    ((dyadic_tendsto _).comp he)

end
end ModAB.Convergence

end RModABBundle3

/- Original module: ModAB.Bracketing -/
section RModABBundle4

/-! Exact-arithmetic interval invariants and the final nested-interval argument. -/
namespace ModAB.Convergence
open Set Filter Topology
noncomputable section

/-- The convex-combination identity used in `lem:interiority`. -/
theorem secant_convex {a b u v : ℝ} (hu : u<0) (hv : 0<v) :
    let alpha:=v/(v+|u|)
    0<alpha ∧ alpha<1 ∧ (a*v-b*u)/(v-u)=alpha*a+(1-alpha)*b := by
  dsimp
  rw [abs_of_neg hu]
  have hd : 0<v-u := by linarith
  refine ⟨div_pos hv hd,(div_lt_one hd).mpr (by linarith),?_⟩
  simp only [← sub_eq_add_neg]
  field_simp [ne_of_gt hd]
  <;> ring

theorem secant_interior {a b u v : ℝ} (hab : a<b) (huv : u*v<0) :
    (a*v-b*u)/(v-u) ∈ Ioo a b := by
  rcases mul_neg_iff.mp huv with h|h
  · obtain ⟨halpha0,halpha1,he⟩:=secant_convex (a:=b) (b:=a) h.2 h.1
    have hs : (a*v-b*u)/(v-u)=(b*u-a*v)/(u-v) := by
      have hd : v-u≠0 := by linarith [h.1,h.2]
      have hd' : u-v≠0 := by linarith [h.1,h.2]
      field_simp
      <;> ring
    rw [hs,he]
    constructor <;> nlinarith [mul_pos halpha0 (sub_pos.mpr hab),mul_pos (sub_pos.mpr halpha1) (sub_pos.mpr hab)]
  · obtain ⟨halpha0,halpha1,he⟩:=secant_convex (a:=a) (b:=b) h.1 h.2
    rw [he]
    constructor <;> nlinarith [mul_pos halpha0 (sub_pos.mpr hab),mul_pos (sub_pos.mpr halpha1) (sub_pos.mpr hab)]

/-- A modified ordinate has the sign of its nonzero true residual. -/
def SameSign (u f : ℝ) : Prop := 0<u*f

theorem positive_correction {u f θ : ℝ} (h : SameSign u f) (hθ : 0<θ) :
    SameSign (θ*u) f := by unfold SameSign at *; nlinarith [mul_pos hθ h]

theorem inserted_sign {f : ℝ} (hf : f≠0) : SameSign f f := by
  exact mul_self_pos.mpr hf

theorem working_opposite {u v fa fb : ℝ}
    (hu : SameSign u fa) (hv : SameSign v fb) (hf : fa*fb<0) : u*v<0 := by
  unfold SameSign at hu hv
  have hp : 0<(u*v)*(fa*fb) := by nlinarith [mul_pos hu hv]
  rcases mul_pos_iff.mp hp with h|h
  · linarith [h.2]
  · exact h.1

/-- Nonzero interior evaluation leaves exactly one sign-changing subinterval. -/
theorem sign_update {fa fb fx : ℝ} (hs : fa*fb<0) (hx : fx≠0) :
    (fa*fx<0 ∧ ¬fx*fb<0) ∨ (fx*fb<0 ∧ ¬fa*fx<0) := by
  rcases mul_neg_iff.mp hs with h|h
  · rcases lt_or_gt_of_ne hx with hx|hx
    · exact Or.inl ⟨mul_neg_of_pos_of_neg h.1 hx,not_lt_of_gt (mul_pos_of_neg_of_neg hx h.2)⟩
    · exact Or.inr ⟨mul_neg_of_pos_of_neg hx h.2,not_lt_of_gt (mul_pos h.1 hx)⟩
  · rcases lt_or_gt_of_ne hx with hx|hx
    · exact Or.inr ⟨mul_neg_of_neg_of_pos hx h.2,not_lt_of_gt (mul_pos_of_neg_of_neg h.1 hx)⟩
    · exact Or.inl ⟨mul_neg_of_neg_of_pos h.1 hx,not_lt_of_gt (mul_pos hx h.2)⟩

/-- Both possible updates are nested and do not increase the width. -/
theorem interval_update {a b x : ℝ} (hx : x∈Icc a b) :
    Icc a x ⊆ Icc a b ∧ Icc x b ⊆ Icc a b ∧
    0≤x-a ∧ x-a≤b-a ∧ 0≤b-x ∧ b-x≤b-a := by
  refine ⟨fun y hy=>⟨hy.1,hy.2.trans hx.2⟩,
    fun y hy=>⟨hx.1.trans hy.1,hy.2⟩,?_,?_,?_,?_⟩ <;> linarith [hx.1,hx.2]

theorem bisection_width (a b : ℝ) :
    (a+b)/2-a=(b-a)/2 ∧ b-(a+b)/2=(b-a)/2 := by constructor <;> ring

theorem bracket_contains_zero {f : ℝ→ℝ} {a b : ℝ} (hab : a≤b)
    (hf : ContinuousOn f (Icc a b)) (hs : f a*f b<0) :
    ∃ z∈Icc a b, f z=0 := by
  rcases mul_neg_iff.mp hs with h|h
  · exact intermediate_value_Icc' hab hf ⟨h.2.le,h.1.le⟩
  · exact intermediate_value_Icc hab hf ⟨h.1.le,h.2.le⟩

theorem residual_progress {r old retained η : ℝ} (hold : 0<old)
    (hret : old≤retained) (hη : η<1) (hnew : r<η*old) :
    min retained r=r ∧ r<η*old := by
  have hr : r<retained := by nlinarith
  exact ⟨min_eq_right hr.le,hnew⟩

/-- Repeated substitution; the strict statement requires at least one step. -/
theorem residual_progress_iterated {R : ℕ→ℝ} {η : ℝ} (hη : 0<η)
    (h : ∀ n, R (n+1)<η*R n) {k : ℕ} (hk : 0<k) : R k<η^k*R 0 := by
  cases k with
  | zero => omega
  | succ k =>
      induction k with
      | zero => simpa using h 0
      | succ k ih =>
          have hmul:=mul_lt_mul_of_pos_left (ih (by omega)) hη
          calc R (k+1+1)<η*R (k+1):=h (k+1)
               _<η*(η^(k+1)*R 0):=hmul
               _=η^(k+1+1)*R 0:=by rw [pow_succ]; ring

theorem interval_distance {a b x y : ℝ} (hx : x∈Icc a b) (hy : y∈Icc a b) :
    |x-y|≤b-a := by rw [abs_le]; constructor <;> linarith [hx.1,hx.2,hy.1,hy.2]

/-- Every choice inside shrinking intervals converges to their common point. -/
theorem chosen_points_tendsto {a b : ℕ→ℝ} {x : ℝ} {z : ℕ→ℝ}
    (hx : ∀n,x∈Icc (a n) (b n)) (hz : ∀n,z n∈Icc (a n) (b n))
    (hw : Tendsto (fun n=>b n-a n) atTop (𝓝 0)) : Tendsto z atTop (𝓝 x) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw ε hε
  refine ⟨N,fun n hn=>?_⟩
  rw [Real.dist_eq]
  have hpos : 0≤b n-a n := by linarith [(hx n).1,(hx n).2]
  have hsmall:=hN n hn
  rw [Real.dist_eq,sub_zero,abs_of_nonneg hpos] at hsmall
  exact (interval_distance (hz n) (hx n)).trans_lt hsmall

/-- Last paragraph of `thm:global_simple`, including the chosen-zero sequence. -/
theorem nested_intervals_root {f : ℝ→ℝ} {a b : ℕ→ℝ}
    (hnested : ∀n,Icc (a (n+1)) (b (n+1))⊆Icc (a n) (b n))
    (hzero : ∀n,∃z∈Icc (a n) (b n),f z=0)
    (hf : ContinuousOn f (Icc (a 0) (b 0)))
    (hw : Tendsto (fun n=>b n-a n) atTop (𝓝 0)) :
    ∃x, (⋂n,Icc (a n) (b n))={x} ∧ f x=0 := by
  have hne : ∀n,(Icc (a n) (b n)).Nonempty := fun n=>by
    obtain ⟨z,hz,_⟩:=hzero n; exact ⟨z,hz⟩
  obtain ⟨x,hx⟩:=IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed
    (fun n=>Icc (a n) (b n)) hnested hne isCompact_Icc (fun _=>isClosed_Icc)
  have hxall : ∀n,x∈Icc (a n) (b n) := mem_iInter.mp hx
  have hinter : (⋂n,Icc (a n) (b n))={x} := by
    apply Set.eq_singleton_iff_unique_mem.mpr
    refine ⟨hx,?_⟩
    intro y hy
    have hyall:=mem_iInter.mp hy
    have hlim:=chosen_points_tendsto hxall hyall hw
    exact tendsto_nhds_unique tendsto_const_nhds hlim
  choose z hz hval using hzero
  have hzlim:=chosen_points_tendsto hxall hz hw
  have hmono : Antitone (fun n=>Icc (a n) (b n)):=antitone_nat_of_succ_le hnested
  have hz0 : ∀n,z n∈Icc (a 0) (b 0):=fun n=>hmono (Nat.zero_le n) (hz n)
  have hcont : Tendsto (fun n=>f (z n)) atTop (𝓝 (f x)) :=
    (hf x (hxall 0)).tendsto.comp (tendsto_nhdsWithin_iff.mpr ⟨hzlim,Filter.Eventually.of_forall hz0⟩)
  have heq : f x=0 := by
    have he : (fun n=>f (z n))=(fun _ : ℕ=>0) := funext hval
    rw [he] at hcont
    exact (tendsto_nhds_unique tendsto_const_nhds hcont).symm
  exact ⟨x,hinter,heq⟩

end
end ModAB.Convergence

end RModABBundle4

/- Original module: ModAB.Convergence -/
section RModABBundle5

/-! Assembly of the interval invariants, the width-policy execution, and the
qualitative and quantitative convergence theorems. All arithmetic here is real. -/
namespace ModAB.Convergence
open Set Filter Topology
noncomputable section

/-- The actual adaptive decision after a nonterminating AB update. -/
def adaptiveMode (B S η residual previousMinimum w : ℝ) (M t c : ℕ) : Mode :=
  if w≤threshold B S (t+1) then .ab S (t+1) 0
  else if residual<η*previousMinimum ∧ c<M then .ab S (t+1) (c+1)
  else .bisection

/-- The residual condition restricts the permitted exceptions; no estimate is assumed. -/
theorem adaptive_step {B S η residual previousMinimum w w' : ℝ} {M t c : ℕ}
    (hw' : 0≤w') (hmono : w'≤w) :
    Step B M w (.ab S t c) w' (adaptiveMode B S η residual previousMinimum w' M t c) := by
  unfold adaptiveMode
  split_ifs with ht hr
  · exact .abWidth hw' hmono ht
  · exact .abException hw' hmono (lt_of_not_ge ht) hr.2
  · exact .abStop hw' hmono

/-- Both sign orientations of the secant proposal are admitted. -/
theorem proposal_interior {a b u v fa fb : ℝ} (hab : a<b)
    (hu : SameSign u fa) (hv : SameSign v fb) (hf : fa*fb<0) :
    (a*v-b*u)/(v-u)∈Ioo a b := secant_interior hab (working_opposite hu hv hf)

/-- One nonterminating exact-arithmetic sign update. -/
def Updated (f : ℝ→ℝ) (a b x a' b' : ℝ) : Prop :=
  (a',b') = if f a*f x<0 then (a,x) else (x,b)

theorem updated_invariants {f : ℝ→ℝ} {a b x a' b' : ℝ}
    (hx : x∈Ioo a b) (hs : f a*f b<0) (hz : f x≠0)
    (hu : Updated f a b x a' b') :
    a'<b' ∧ f a'*f b'<0 ∧ Icc a' b'⊆Icc a b := by
  unfold Updated at hu
  have hi:=interval_update ⟨hx.1.le,hx.2.le⟩
  by_cases hax : f a*f x<0
  · rw [if_pos hax,Prod.mk.injEq] at hu
    rcases hu with ⟨rfl,rfl⟩
    exact ⟨hx.1,hax,hi.1⟩
  · rw [if_neg hax,Prod.mk.injEq] at hu
    have hxb : f x*f b<0 := by
      rcases sign_update hs hz with h|h
      · exact False.elim (hax h.1)
      · exact h.1
    rcases hu with ⟨rfl,rfl⟩
    exact ⟨hx.2,hxb,hi.2.1⟩

/-- An infinite nonterminating execution, with the sign update made explicit.
The proposal is any interior point and the modes obey the width policy. This
is the generality used in Remark `rem:proof_independence`; `proposal_interior`
verifies the AB proposal, and a midpoint is interior as well. -/
structure ExactRun (B : ℝ) (M : ℕ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→Mode
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  initial_mode : mode 0=.bisection
  interior : ∀n,point n∈Ioo (a n) (b n)
  nonzero : ∀n,f (point n)≠0
  update : ∀n,Updated f (a n) (b n) (point n) (a (n+1)) (b (n+1))
  width_step : ∀n,Step B M (b n-a n) (mode n) (b (n+1)-a (n+1)) (mode (n+1))

theorem exact_invariants {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) :
    ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 ∧
      Icc (r.a (n+1)) (r.b (n+1))⊆Icc (r.a n) (r.b n) := by
  have hs : ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 := by
    intro n
    induction n with
    | zero => exact ⟨r.initial_order,r.initial_sign⟩
    | succ n ih => exact (updated_invariants (r.interior n) ih.2 (r.nonzero n) (r.update n)).imp_right And.left
  intro n
  exact ⟨(hs n).1,(hs n).2,(updated_invariants (r.interior n) (hs n).2 (r.nonzero n) (r.update n)).2.2⟩

def ExactRun.widthRun {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) : Run B M where
  width n:=r.b n-r.a n
  mode:=r.mode
  nonneg_initial:=by linarith [r.initial_order]
  initial:=r.initial_mode
  next:=r.width_step

theorem exact_zeros {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) :
    ∀n,∃z∈Icc (r.a n) (r.b n),f z=0 := by
  have hi:=exact_invariants r
  have hmono : Antitone (fun n=>Icc (r.a n) (r.b n)):=antitone_nat_of_succ_le (fun n=>(hi n).2.2)
  intro n
  exact bracket_contains_zero (hi n).1.le (r.continuous.mono (hmono (Nat.zero_le n))) (hi n).2.1

/-- Full nonterminating alternative of `thm:global_simple`. Exact-zero termination
is the other alternative and needs no limit argument. -/
theorem global_convergence {B : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (r : ExactRun B M f) :
    Tendsto (fun n=>r.b n-r.a n) atTop (𝓝 0) ∧
    ∃x,(⋂n,Icc (r.a n) (r.b n))={x} ∧ f x=0 := by
  have hw:=width_tendsto_zero hB r.widthRun
  exact ⟨hw,nested_intervals_root (fun n=>(exact_invariants r n).2.2) (exact_zeros r) r.continuous hw⟩

theorem exact_global_width {B : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (r : ExactRun B M f) (n : ℕ) :
    r.b n-r.a n≤(r.b 0-r.a 0)/(2:ℝ)^(n/(dyadicExponent (envelopeConstant B M)+2)) :=
  global_width_bound_log hB r.widthRun.nonneg_initial (history_of_run r.widthRun n)

/-- `cor:finite_termination`: a continued run reaches the width tolerance by qK. -/
theorem exact_finite_tolerance {B τ : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (hτ : 0<τ) (r : ExactRun B M f) :
    let q:=dyadicExponent (envelopeConstant B M)+2
    let K:=toleranceExponent (r.b 0-r.a 0) τ
    r.b (q*K)-r.a (q*K)≤τ := by
  dsimp
  exact finite_termination hB r.widthRun.nonneg_initial
    (ceil_log_bound (lt_of_lt_of_le (by norm_num) (C_at_least_one hB)))
    (tolerance_exponent_sufficient r.widthRun.nonneg_initial hτ)
    (history_of_run r.widthRun _)

/-- `eq:distance_zero_set`, expressed as an explicit nearby zero. -/
theorem distance_to_zero {f : ℝ→ℝ} {a b x τ : ℝ} (hx : x∈Icc a b)
    (hz : ∃z∈Icc a b,f z=0) (hw : b-a≤τ) : ∃z,f z=0 ∧ |x-z|≤τ := by
  obtain ⟨z,hz,hfz⟩:=hz
  exact ⟨z,hfz,(interval_distance hx hz).trans hw⟩

end
end ModAB.Convergence

end RModABBundle5

/- Original module: ModAB.ExamplesBase -/
section RModABBundle6

/-! Exact one-sided AB states and the regularity and limit facts used by both
manuscript examples. Every state transition is an explicit secant and correction.
The sign-update bridge proves when the right endpoint actually remains fixed. -/
namespace ModAB.Examples
open Set Filter Topology
noncomputable section

def gluedQuadratic (A B x : ℝ) : ℝ := if x≤0 then A*x^2 else B*x^2
def gluedSlope (A B x : ℝ) : ℝ := if x≤0 then A*x else B*x

theorem gluedSlope_continuous (A B : ℝ) : Continuous (gluedSlope A B) := by
  unfold gluedSlope
  apply Continuous.if_le (continuous_const.mul continuous_id)
    (continuous_const.mul continuous_id) continuous_id continuous_const
  intro x hx
  change x=0 at hx
  simp [hx]

theorem gluedQuadratic_slope (A B x : ℝ) :
    slope (gluedQuadratic A B) 0 x=gluedSlope A B x := by
  by_cases hx : x=0
  · simp [hx,slope_def_field,gluedQuadratic,gluedSlope]
  · simp only [slope_def_field,gluedQuadratic,gluedSlope,le_refl,if_true,
      zero_pow (by decide : 2≠0),mul_zero,sub_zero]
    split_ifs <;> field_simp <;> ring

theorem gluedQuadratic_hasDerivAt_zero (A B : ℝ) :
    HasDerivAt (gluedQuadratic A B) 0 0 := by
  apply hasDerivAt_iff_tendsto_slope.mpr
  have ht := (gluedSlope_continuous A B).continuousAt.tendsto.mono_left
    (nhdsWithin_le_nhds : 𝓝[≠] (0:ℝ)≤𝓝 0)
  have he : slope (gluedQuadratic A B) 0=gluedSlope A B:=
    funext (gluedQuadratic_slope A B)
  rw [he]
  simpa [gluedSlope] using ht

theorem gluedQuadratic_hasDerivAt (A B x : ℝ) :
    HasDerivAt (gluedQuadratic A B) (2*gluedSlope A B x) x := by
  rcases lt_trichotomy x 0 with hx|hx|hx
  · have h := ((hasDerivAt_id x).pow 2).const_mul A
    have hg : HasDerivAt (fun y : ℝ=>A*y^2) (2*(A*x)) x := by
      convert h using 1 <;> simp <;> ring
    have he := hg.congr_of_eventuallyEq (f₁:=gluedQuadratic A B) (by
      filter_upwards [eventually_lt_nhds hx] with y hy
      simp [gluedQuadratic,hy.le])
    simpa [gluedSlope,hx.le] using he
  · subst x
    simpa [gluedSlope] using gluedQuadratic_hasDerivAt_zero A B
  · have h := ((hasDerivAt_id x).pow 2).const_mul B
    have hg : HasDerivAt (fun y : ℝ=>B*y^2) (2*(B*x)) x := by
      convert h using 1 <;> simp <;> ring
    have he := hg.congr_of_eventuallyEq (f₁:=gluedQuadratic A B) (by
      filter_upwards [eventually_gt_nhds hx] with y hy
      simp [gluedQuadratic,not_le.mpr hy])
    simpa [gluedSlope,not_le.mpr hx] using he

theorem gluedQuadratic_C1 (A B : ℝ) : ContDiff ℝ 1 (gluedQuadratic A B) := by
  apply contDiff_one_iff_deriv.mpr
  constructor
  · exact fun x=>(gluedQuadratic_hasDerivAt A B x).differentiableAt
  · have he : deriv (gluedQuadratic A B)=fun x=>2*gluedSlope A B x :=
      funext (fun x=>(gluedQuadratic_hasDerivAt A B x).deriv)
    rw [he]
    exact continuous_const.mul (gluedSlope_continuous A B)

def oneFunction (x : ℝ) : ℝ := x+gluedQuadratic 0 1 x
def counterFunction (x : ℝ) : ℝ := gluedQuadratic (-1) 7 x

theorem oneFunction_formula (x : ℝ) :
    oneFunction x=if x≤0 then x else x+x^2 := by
  unfold oneFunction gluedQuadratic
  split_ifs <;> ring

theorem counterFunction_formula (x : ℝ) :
    counterFunction x=if x<0 then -x^2 else 7*x^2 := by
  rcases lt_trichotomy x 0 with h|h|h
  · simp [counterFunction,gluedQuadratic,h,h.le]
  · simp [h,counterFunction,gluedQuadratic]
  · simp [counterFunction,gluedQuadratic,not_lt.mpr h.le,not_le.mpr h]

theorem oneFunction_C1 : ContDiff ℝ 1 oneFunction :=
  contDiff_id.add (gluedQuadratic_C1 0 1)

theorem counterFunction_C1 : ContDiff ℝ 1 counterFunction := gluedQuadratic_C1 (-1) 7

theorem oneFunction_left {e : ℝ} (he : 0≤e) : oneFunction (-e)=-e := by
  rw [oneFunction_formula,if_pos (by linarith)]

theorem counterFunction_left {e : ℝ} (he : 0≤e) : counterFunction (-e)=-e^2 := by
  simp [counterFunction,gluedQuadratic,show -e≤0 by linarith]

theorem counterFunction_zero_iff (x : ℝ) : counterFunction x=0 ↔ x=0 := by
  rw [counterFunction_formula]
  split_ifs <;> constructor
  · intro h; nlinarith [sq_nonneg x]
  · intro h; simp [h]
  · intro h; nlinarith [sq_nonneg x]
  · intro h; simp [h]

def leftTrial (p : ℕ) (e v : ℝ) : ℝ := (e*v-e^p)/(v+e^p)

structure LeftState where
  e : ℝ
  v : ℝ
  correct : Bool

def leftStep (p : ℕ) (s : LeftState) : LeftState :=
  let t:=leftTrial p s.e s.v
  ⟨t,if s.correct then s.v*(1-t^p/s.e^p) else s.v,true⟩

def leftOrbit (p : ℕ) (s : LeftState) : ℕ→LeftState
  | 0=>s
  | n+1=>leftStep p (leftOrbit p s n)

theorem left_secant_identity (p : ℕ) (e v : ℝ) :
    -leftTrial p e v=((-e)*v-1*(-e^p))/(v-(-e^p)) := by
  simp only [leftTrial,sub_neg_eq_add,one_mul]
  rw [←neg_div]
  congr 1
  ring

theorem left_updated {f : ℝ→ℝ} {p : ℕ} {e t : ℝ}
    (he : 0<e) (ht : 0<t) (hfa : f (-e)=-e^p) (hfx : f (-t)=-t^p) :
    ModAB.Convergence.Updated f (-e) 1 (-t) (-t) 1 := by
  have hp : 0≤f (-e)*f (-t) := by
    rw [hfa,hfx,neg_mul_neg]
    exact (mul_pos (pow_pos he p) (pow_pos ht p)).le
  simp [ModAB.Convergence.Updated,not_lt.mpr hp]

theorem geometric_bound {e : ℕ→ℝ} {r : ℝ} (hr : 0≤r)
    (hs : ∀n,e (n+1)≤r*e n) : ∀n,e n≤r^n*e 0 := by
  intro n
  induction n with
  | zero=>simp
  | succ n ih=>
    calc
      e (n+1)≤r*e n:=hs n
      _≤r*(r^n*e 0):=mul_le_mul_of_nonneg_left ih hr
      _=r^(n+1)*e 0:=by rw [pow_succ]; ring

theorem geometric_tendsto {e : ℕ→ℝ} {r : ℝ} (he : ∀n,0≤e n)
    (hr : 0≤r) (hr1 : r<1) (hs : ∀n,e (n+1)≤r*e n) :
    Tendsto e atTop (𝓝 0) := by
  have hb:=geometric_bound hr hs
  have ht:= (tendsto_pow_atTop_nhds_zero_of_lt_one hr hr1).mul_const (e 0)
  simp only [zero_mul] at ht
  exact squeeze_zero he hb ht

end
end ModAB.Examples

end RModABBundle6

/- Original module: ModAB.TopEntry -/
section RModABBundle7

/-! Exact-arithmetic bridge for top switching. The switching point is still a
midpoint. Recording its moved side changes the next correction, not the bracket
or the phase width policy. This is not a verification of the C# loop. -/
namespace ModAB.Convergence
noncomputable section

inductive EndpointSide where
  | left | right
  deriving DecidableEq

/-- `GetABFactor`, with its positive fallback. -/
def abFactor (new old : ℝ) : ℝ :=
  if 0 < 1-new/old then 1-new/old else 1/2

theorem abFactor_positive (new old : ℝ) : 0 < abFactor new old := by
  unfold abFactor
  split_ifs with h
  · exact h
  · norm_num

structure OrdinateState where
  leftTrue : ℝ
  rightTrue : ℝ
  leftWork : ℝ
  rightWork : ℝ
  previous : Option EndpointSide

/-- The side is chosen using true function signs outside this ordinate update. -/
def recordUpdate (s : OrdinateState) (moved : EndpointSide) (y : ℝ) : OrdinateState :=
  match moved with
  | .left => ⟨y,s.rightTrue,y,
      (if s.previous = some .left then abFactor y s.leftWork else 1)*s.rightWork,
      some .left⟩
  | .right => ⟨s.leftTrue,y,
      (if s.previous = some .right then abFactor y s.rightWork else 1)*s.leftWork,y,
      some .right⟩

theorem recorded_side (s : OrdinateState) (moved : EndpointSide) (y : ℝ) :
    (recordUpdate s moved y).previous = some moved := by cases moved <;> rfl

/-- Opposite-side first secant: no scaling; the actual new side is recorded. -/
theorem opposite_side_no_correction (a b u v y : ℝ) :
    recordUpdate ⟨a,b,u,v,some .left⟩ .right y = ⟨a,y,u,y,some .right⟩ ∧
    recordUpdate ⟨a,b,u,v,some .right⟩ .left y = ⟨y,b,y,v,some .left⟩ := by
  simp [recordUpdate]

/-- Same-side first secant: exactly the usual positive AB factor is applied. -/
theorem same_side_correction (a b u v y : ℝ) :
    (recordUpdate ⟨a,b,u,v,some .left⟩ .left y).rightWork = abFactor y u*v ∧
    (recordUpdate ⟨a,b,u,v,some .right⟩ .right y).leftWork = abFactor y v*u := by
  simp [recordUpdate]

theorem recordUpdate_working_signs (s : OrdinateState) (moved : EndpointSide) (y : ℝ)
    (hu : SameSign s.leftWork s.leftTrue) (hv : SameSign s.rightWork s.rightTrue)
    (hy : y ≠ 0) :
    SameSign (recordUpdate s moved y).leftWork (recordUpdate s moved y).leftTrue ∧
    SameSign (recordUpdate s moved y).rightWork (recordUpdate s moved y).rightTrue := by
  have hp (old : ℝ) (side : EndpointSide) :
      0 < (if s.previous = some side then abFactor y old else 1) := by
    split_ifs
    · exact abFactor_positive y old
    · norm_num
  cases moved with
  | left => exact ⟨inserted_sign hy,positive_correction hv (hp _ _)⟩
  | right => exact ⟨positive_correction hu (hp _ _),inserted_sign hy⟩

/-- Before the switching bisection, the previous side is cleared and the
working values are reset to the true endpoint residuals. -/
theorem top_entry_ordinates (a b y : ℝ) :
    recordUpdate ⟨a,b,a,b,none⟩ .left y = ⟨y,b,y,b,some .left⟩ ∧
    recordUpdate ⟨a,b,a,b,none⟩ .right y = ⟨a,y,a,y,some .right⟩ := by
  simp [recordUpdate]

/-- At entry the width test succeeds, independently of residuals. Halving
its stored threshold gives exactly the threshold of the first true secant. -/
theorem top_entry_width_policy {B W : ℝ} (hB : (1/2:ℝ) ≤ B) (hW : 0 ≤ W) (M : ℕ) :
    W/2 ≤ B*W ∧ (B*W)/2 = threshold B (W/2) 1 ∧
    Step B M W .bisection (W/2) (.ab (W/2) 0 0) := by
  refine ⟨?_,?_,Step.bisectStart W⟩
  · nlinarith [mul_nonneg (sub_nonneg.mpr hB) hW]
  · simp [threshold]; ring

/-- Concrete common entry bracket for the examples starting from [-3,1]. -/
theorem top_entry_midpoint :
    ((-3:ℝ)+1)/2 = -1 ∧ (1:ℝ)-(-1)=2 := by norm_num

/-- Concrete left-side top initialization used by the example orbit models. -/
theorem top_left_seed (a v : ℝ) :
    recordUpdate ⟨a,v,a,v,none⟩ .left (-1) = ⟨-1,v,-1,v,some .left⟩ :=
  (top_entry_ordinates a v (-1)).1

end
end ModAB.Convergence

namespace ModAB.Examples
noncomputable section
/-- Projection of the recorded top-entry state onto the one-sided orbit model. -/
def topLeftSeed (a v : ℝ) : LeftState :=
  let s := ModAB.Convergence.recordUpdate ⟨a,v,a,v,none⟩ .left (-1)
  ⟨1,s.rightWork,decide (s.previous = some .left)⟩

theorem topLeftSeed_eq (a v : ℝ) : topLeftSeed a v = ⟨1,v,true⟩ := by
  simp [topLeftSeed,ModAB.Convergence.recordUpdate]
end
end ModAB.Examples

end RModABBundle7

/- Original module: ModAB.Rational -/
section RModABBundle8

/-! Independent exact-arithmetic specification of the rational proposal.
The only algorithmic input taken from the new description is its normalized
candidate and secant-relative limit. These results establish algebra and
admissibility; they assert neither an asymptotic order nor an uncapped width
schedule. All bracket decisions must still use true residuals. -/
namespace ModAB.Rational
open Set
noncomputable section

/-- Secant through working ordinates, which preserve the true endpoint signs. -/
def secant (a b ya yb : ℝ) : ℝ := (a*yb-b*ya)/(yb-ya)

/-- Distance of the secant to its nearer endpoint. -/
def radius (a b ya yb : ℝ) : ℝ :=
  min (secant a b ya yb-a) (b-secant a b ya yb)

/-- The description's normalized denominator, q=ya/yb, qp=fp/yb. -/
def normalizedDenominator (w u q qp : ℝ) : ℝ :=
  w*(qp-q)-u*qp*(1-q)

def normalizedNumerator (w u q qp : ℝ) : ℝ := w*u*q*(qp-1)

def normalizedOffset (w u q qp : ℝ) : ℝ :=
  normalizedNumerator w u q qp / normalizedDenominator w u q qp

def denominator (a b ya yb xp fp : ℝ) : ℝ :=
  normalizedDenominator (b-a) (xp-a) (ya/yb) (fp/yb)

def rawPoint (a b ya yb xp fp : ℝ) : ℝ :=
  a+normalizedOffset (b-a) (xp-a) (ya/yb) (fp/yb)

/-- Real-valued projection; machine NaN/infinity handling is a separate obligation. -/
def limited (a b ya yb raw : ℝ) : ℝ :=
  let s := secant a b ya yb
  let d := radius a b ya yb
  max (s-(7/10:ℝ)*d) (min raw (s+(7/10:ℝ)*d))

/-- Guarded exact-real completion, not the literal exceptional path of C#.
Curvature may disable the proposal. A zero rational denominator falls back
to the secant here; C# instead uses IEEE infinity/NaN and its right-endpoint
fallback. Every admitted finite rational value is projected. -/
def candidate (a b ya yb xp fp : ℝ) (allowRational : Bool) : ℝ :=
  if allowRational = true ∧ denominator a b ya yb xp fp ≠ 0
  then limited a b ya yb (rawPoint a b ya yb xp fp)
  else secant a b ya yb

theorem secant_interior {a b ya yb : ℝ} (hab : a<b) (hy : ya*yb<0) :
    secant a b ya yb ∈ Ioo a b :=
  ModAB.Convergence.secant_interior hab hy

theorem working_ratio_negative {ya yb : ℝ} (hy : ya*yb<0) : ya/yb<0 := by
  rcases mul_neg_iff.mp hy with h|h
  · exact div_neg_of_pos_of_neg h.1 h.2
  · exact div_neg_of_neg_of_pos h.1 h.2

theorem working_second_ne_zero {ya yb : ℝ} (hy : ya*yb<0) : yb≠0 := by
  intro h
  simp [h] at hy

theorem secant_normalized {a b ya yb : ℝ} (hy : ya*yb<0) :
    secant a b ya yb = a-(b-a)*(ya/yb)/(1-ya/yb) := by
  have hv := working_second_ne_zero hy
  have hq := working_ratio_negative hy
  have hc : 1-ya/yb≠0 := by linarith
  have hd : yb-ya≠0 := by
    intro h
    have he : yb=ya := sub_eq_zero.mp h
    rw [he] at hy
    nlinarith [sq_nonneg ya]
  unfold secant
  field_simp [hv,hc,hd]
  ring

theorem radius_pos {a b ya yb : ℝ} (hab : a<b) (hy : ya*yb<0) :
    0<radius a b ya yb := by
  have hs := secant_interior hab hy
  unfold radius
  exact lt_min (sub_pos.mpr hs.1) (sub_pos.mpr hs.2)

theorem radius_le_left (a b ya yb : ℝ) :
    radius a b ya yb ≤ secant a b ya yb-a := min_le_left _ _

theorem radius_le_right (a b ya yb : ℝ) :
    radius a b ya yb ≤ b-secant a b ya yb := min_le_right _ _

theorem limited_bounds {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    secant a b ya yb-(7/10:ℝ)*radius a b ya yb ≤ limited a b ya yb raw ∧
    limited a b ya yb raw ≤ secant a b ya yb+(7/10:ℝ)*radius a b ya yb := by
  have hd := radius_pos hab hy
  unfold limited
  dsimp
  exact ⟨le_max_left _ _,max_le (by nlinarith) (min_le_right _ _)⟩

theorem limited_margins {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    a+(3/10:ℝ)*radius a b ya yb ≤ limited a b ya yb raw ∧
    limited a b ya yb raw ≤ b-(3/10:ℝ)*radius a b ya yb := by
  have hb := limited_bounds (raw:=raw) hab hy
  have hl := radius_le_left a b ya yb
  have hr := radius_le_right a b ya yb
  constructor <;> nlinarith [hb.1,hb.2]

theorem limited_interior {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    limited a b ya yb raw ∈ Ioo a b := by
  have hd := radius_pos hab hy
  have hm := limited_margins (raw:=raw) hab hy
  constructor <;> nlinarith [hm.1,hm.2]

theorem limited_eq_self {a b ya yb raw : ℝ}
    (hl : secant a b ya yb-(7/10:ℝ)*radius a b ya yb ≤ raw)
    (hr : raw ≤ secant a b ya yb+(7/10:ℝ)*radius a b ya yb) :
    limited a b ya yb raw=raw := by
  unfold limited
  dsimp
  rw [min_eq_left hr,max_eq_right hl]

theorem candidate_interior {a b ya yb xp fp : ℝ} (allowRational : Bool)
    (hab : a<b) (hy : ya*yb<0) :
    candidate a b ya yb xp fp allowRational ∈ Ioo a b := by
  unfold candidate
  split_ifs
  · exact limited_interior hab hy
  · exact secant_interior hab hy

theorem candidate_fallback_denominator {a b ya yb xp fp : ℝ}
    (allowRational : Bool) (h : denominator a b ya yb xp fp=0) :
    candidate a b ya yb xp fp allowRational=secant a b ya yb := by
  simp [candidate,h]

theorem candidate_fallback_curvature (a b ya yb xp fp : ℝ) :
    candidate a b ya yb xp fp false=secant a b ya yb := by
  simp [candidate]

/-- A conservative exact-real proposal relation for the source branches.
The projected cases cover finite proposals and the two limiting projections of
infinities; the right endpoint admits the source's NaN/rounding fallback.
This is not a theorem about all IEEE operations in the implementation. -/
def SourceProposal (a b ya yb x : ℝ) : Prop :=
  x=secant a b ya yb ∨ (∃raw,x=limited a b ya yb raw) ∨ x=b

theorem sourceProposal_mem_closed {a b ya yb x : ℝ}
    (hab : a<b) (hy : ya*yb<0) (h : SourceProposal a b ya yb x) :
    x ∈ Icc a b := by
  rcases h with rfl|⟨raw,rfl⟩|rfl
  · have hs := secant_interior hab hy
    exact ⟨hs.1.le,hs.2.le⟩
  · have hs := limited_interior (raw:=raw) hab hy
    exact ⟨hs.1.le,hs.2.le⟩
  · exact ⟨hab.le,le_rfl⟩

theorem sourceProposal_interior_of_ne_right {a b ya yb x : ℝ}
    (hab : a<b) (hy : ya*yb<0) (h : SourceProposal a b ya yb x)
    (hne : x≠b) : x ∈ Ioo a b := by
  rcases h with rfl|⟨raw,rfl⟩|rfl
  · exact secant_interior hab hy
  · exact limited_interior hab hy
  · exact False.elim (hne rfl)

theorem candidate_is_sourceProposal (a b ya yb xp fp : ℝ) (allowRational : Bool) :
    SourceProposal a b ya yb (candidate a b ya yb xp fp allowRational) := by
  unfold candidate
  split_ifs
  · exact Or.inr (Or.inl ⟨_,rfl⟩)
  · exact Or.inl rfl

/-- Source-order closing operation, applied AFTER the secant-relative clamp.
The supplied tolerance is computed from the pre-closing point. The source's
second guard prevents a displacement past the opposite endpoint. -/
def closeCandidate (a b tolerance x : ℝ) (mostRecent : ModAB.Convergence.EndpointSide) : ℝ :=
  match mostRecent with
  | .left => if x-a<tolerance/2 ∧ a+tolerance/2<b then a+tolerance/2 else x
  | .right => if b-x<tolerance/2 ∧ a<b-tolerance/2 then b-tolerance/2 else x

theorem closeCandidate_interior {a b tolerance x : ℝ}
    (mostRecent : ModAB.Convergence.EndpointSide)
    (ht : 0<tolerance) (hx : x∈Ioo a b) :
    closeCandidate a b tolerance x mostRecent ∈ Ioo a b := by
  cases mostRecent <;> dsimp only [closeCandidate] <;> split_ifs with h
  · exact ⟨by linarith,h.2⟩
  · exact hx
  · exact ⟨h.2,by linarith⟩
  · exact hx

theorem closeCandidate_mem_closed {a b tolerance x : ℝ}
    (mostRecent : ModAB.Convergence.EndpointSide)
    (ht : 0<tolerance) (hx : x∈Icc a b) :
    closeCandidate a b tolerance x mostRecent ∈ Icc a b := by
  cases mostRecent <;> dsimp only [closeCandidate] <;> split_ifs with h
  · exact ⟨by linarith,h.2.le⟩
  · exact hx
  · exact ⟨h.2.le,by linarith⟩
  · exact hx

theorem closeCandidate_left_separation {a b tolerance x : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    a+tolerance/2 ≤ closeCandidate a b tolerance x .left := by
  have hg : a+tolerance/2<b := by linarith
  simp only [closeCandidate,hg,and_true]
  split_ifs with h <;> linarith

theorem closeCandidate_right_separation {a b tolerance x : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    closeCandidate a b tolerance x .right ≤ b-tolerance/2 := by
  have hg : a<b-tolerance/2 := by linarith
  simp only [closeCandidate,hg,and_true]
  split_ifs with h <;> linarith

theorem closeCandidate_right_fallback_interior {a b tolerance : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    closeCandidate a b tolerance b .right ∈ Ioo a b := by
  have hg : a<b-tolerance/2 := by linarith
  have hz : b-b<tolerance/2 := by linarith
  change (if b-b<tolerance/2 ∧ a<b-tolerance/2 then b-tolerance/2 else b) ∈ Ioo a b
  rw [ite_eq_left ⟨hz,hg⟩]
  exact ⟨hg,by linarith⟩

/-- With an outside history node and opposite endpoint ordinates, the normalized
numerator and denominator cannot vanish together in exact arithmetic. This says
nothing about NaN caused by missing history or floating-point overflow. -/
theorem numerator_ne_zero_of_denominator_zero {w u q qp : ℝ}
    (hw : 0<w) (hu : u<0 ∨ w<u) (hq : q<0)
    (hD : normalizedDenominator w u q qp=0) :
    normalizedNumerator w u q qp≠0 := by
  have hwn : w≠0 := ne_of_gt hw
  have hun : u≠0 := by rcases hu with hu|hu <;> linarith
  have hqn : q≠0 := ne_of_lt hq
  intro hn
  unfold normalizedNumerator at hn
  have hqp : qp-1=0 := (mul_eq_zero.mp hn).resolve_left
    (mul_ne_zero (mul_ne_zero hwn hun) hqn)
  have he : qp=1 := by linarith
  have hwu : w-u≠0 := by rcases hu with hu|hu <;> linarith
  have hc : 1-q≠0 := by linarith
  have hid : normalizedDenominator w u q qp=(w-u)*(1-q) := by
    rw [he]
    unfold normalizedDenominator
    ring
  rw [hid] at hD
  exact (mul_ne_zero hwu hc) hD

/-- Independent coefficient construction: the normalized proposal satisfies
all three cross-multiplied interpolation equations. Genuine function-value
interpolation additionally requires nonzero model denominators at the nodes. -/
theorem normalized_interpolation_equations {w u q qp : ℝ}
    (hw : w≠0) (hq : q≠0) (hD : normalizedDenominator w u q qp≠0) :
    let t := normalizedOffset w u q qp
    let p := 1+t*(1-q)/(q*w)
    let Q := -t/q
    (-t=q*Q ∧ w-t=p*w+Q ∧ u-t=qp*(p*u+Q)) := by
  dsimp [normalizedOffset,normalizedNumerator,normalizedDenominator]
  unfold normalizedDenominator at hD
  constructor
  · field_simp [hw,hq,hD]
  constructor
  · field_simp [hw,hq,hD]
    ring
  · field_simp [hw,hq,hD]
    ring

/-- A history record deliberately keeps discarded TRUE residuals. -/
structure HistoryState where
  left : ℝ
  right : ℝ
  ordinates : ModAB.Convergence.OrdinateState
  dropped : Option (ℝ×ℝ)

/-- Genuine endpoint replacement: this helper records the removed endpoint.
For source executions that reuse the same coordinate, use recordSourceHistory. -/
def recordHistory (s : HistoryState) (moved : ModAB.Convergence.EndpointSide)
    (x y : ℝ) : HistoryState :=
  match moved with
  | .left => ⟨x,s.right,ModAB.Convergence.recordUpdate s.ordinates moved y,
      some (s.left,s.ordinates.leftTrue)⟩
  | .right => ⟨s.left,x,ModAB.Convergence.recordUpdate s.ordinates moved y,
      some (s.right,s.ordinates.rightTrue)⟩

theorem recordHistory_true_values {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    (recordHistory s moved x (f x)).ordinates.leftTrue =
      f (recordHistory s moved x (f x)).left ∧
    (recordHistory s moved x (f x)).ordinates.rightTrue =
      f (recordHistory s moved x (f x)).right := by
  cases moved <;> simp [recordHistory,ModAB.Convergence.recordUpdate,hl,hr]

theorem recordHistory_dropped_true {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    ∃xp,(recordHistory s moved x (f x)).dropped=some (xp,f xp) := by
  cases moved with
  | left => exact ⟨s.left,by simp [recordHistory,hl]⟩
  | right => exact ⟨s.right,by simp [recordHistory,hr]⟩

/-- Truth of an optional history record: an absent history is allowed. -/
def TrueHistory (f : ℝ→ℝ) (history : Option (ℝ×ℝ)) : Prop :=
  ∀xp yp,history=some (xp,yp) → yp=f xp

/-- Source-style history guard. Reusing an endpoint preserves the previous
history, but it does not skip the AB ordinate update. In the source a repeated
same-side reuse can still apply an AB correction before this guard. This local
helper uses the old AB-side update and does not encode the source's independent
lastMoved/sideMoved plateau/reset controller. -/
def recordSourceHistory (s : HistoryState) (moved : ModAB.Convergence.EndpointSide)
    (x y : ℝ) : HistoryState :=
  match moved with
  | .left => ⟨x,s.right,ModAB.Convergence.recordUpdate s.ordinates moved y,
      if x=s.left then s.dropped else some (s.left,s.ordinates.leftTrue)⟩
  | .right => ⟨s.left,x,ModAB.Convergence.recordUpdate s.ordinates moved y,
      if x=s.right then s.dropped else some (s.right,s.ordinates.rightTrue)⟩

theorem recordSourceHistory_true_values {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    (recordSourceHistory s moved x (f x)).ordinates.leftTrue =
      f (recordSourceHistory s moved x (f x)).left ∧
    (recordSourceHistory s moved x (f x)).ordinates.rightTrue =
      f (recordSourceHistory s moved x (f x)).right := by
  cases moved <;> simp [recordSourceHistory,ModAB.Convergence.recordUpdate,hl,hr]

theorem recordSourceHistory_history_true {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right)
    (hd : TrueHistory f s.dropped) :
    TrueHistory f (recordSourceHistory s moved x (f x)).dropped := by
  intro xp yp hh
  cases moved with
  | left =>
      by_cases hx : x=s.left
      · exact hd xp yp (by simpa [recordSourceHistory,hx] using hh)
      · have he : s.left=xp ∧ s.ordinates.leftTrue=yp := by
          simpa [recordSourceHistory,hx] using hh
        calc yp = s.ordinates.leftTrue := he.2.symm
             _ = f s.left := hl
             _ = f xp := congrArg f he.1
  | right =>
      by_cases hx : x=s.right
      · exact hd xp yp (by simpa [recordSourceHistory,hx] using hh)
      · have he : s.right=xp ∧ s.ordinates.rightTrue=yp := by
          simpa [recordSourceHistory,hx] using hh
        calc yp = s.ordinates.rightTrue := he.2.symm
             _ = f s.right := hr
             _ = f xp := congrArg f he.1

theorem recordSourceHistory_reuse_left (s : HistoryState) (y : ℝ) :
    (recordSourceHistory s .left s.left y).dropped=s.dropped := by
  simp [recordSourceHistory]

theorem recordSourceHistory_reuse_right (s : HistoryState) (y : ℝ) :
    (recordSourceHistory s .right s.right y).dropped=s.dropped := by
  simp [recordSourceHistory]

theorem recordSourceHistory_left_eq_recordHistory {s : HistoryState} {x y : ℝ}
    (hx : x≠s.left) :
    recordSourceHistory s .left x y=recordHistory s .left x y := by
  simp [recordSourceHistory,recordHistory,hx]

theorem recordSourceHistory_right_eq_recordHistory {s : HistoryState} {x y : ℝ}
    (hx : x≠s.right) :
    recordSourceHistory s .right x y=recordHistory s .right x y := by
  simp [recordSourceHistory,recordHistory,hx]

/-- The endpoint recorded as most recently evaluated has its true ordinate. -/
def Fresh (s : ModAB.Convergence.OrdinateState) : Prop :=
  match s.previous with
  | some .left => s.leftWork=s.leftTrue
  | some .right => s.rightWork=s.rightTrue
  | none => s.leftWork=s.leftTrue ∧ s.rightWork=s.rightTrue

theorem recordUpdate_fresh (s : ModAB.Convergence.OrdinateState)
    (moved : ModAB.Convergence.EndpointSide) (y : ℝ) :
    Fresh (ModAB.Convergence.recordUpdate s moved y) := by
  cases moved <;> simp [Fresh,ModAB.Convergence.recordUpdate]

theorem recordUpdate_at_most_one_modified (s : ModAB.Convergence.OrdinateState)
    (moved : ModAB.Convergence.EndpointSide) (y : ℝ) :
    (ModAB.Convergence.recordUpdate s moved y).leftWork=
      (ModAB.Convergence.recordUpdate s moved y).leftTrue ∨
    (ModAB.Convergence.recordUpdate s moved y).rightWork=
      (ModAB.Convergence.recordUpdate s moved y).rightTrue := by
  cases moved <;> simp [ModAB.Convergence.recordUpdate]

/-- When a correction is applied on the same side, the old implementation's
working denominator equals the true replaced residual by the freshness invariant. -/
theorem correction_true_denominator (s : ModAB.Convergence.OrdinateState)
    (h : Fresh s) (y : ℝ) :
    (s.previous=some .left → ModAB.Convergence.abFactor y s.leftWork=
      ModAB.Convergence.abFactor y s.leftTrue) ∧
    (s.previous=some .right → ModAB.Convergence.abFactor y s.rightWork=
      ModAB.Convergence.abFactor y s.rightTrue) := by
  constructor
  · intro hp
    have he : s.leftWork=s.leftTrue := by simpa [Fresh,hp] using h
    rw [he]
  · intro hp
    have he : s.rightWork=s.rightTrue := by simpa [Fresh,hp] using h
    rw [he]

/-- Same-sign replacement residuals make the AB multiplier a strict contraction;
together with the old positivity theorem this yields a factor in (0,1). -/
theorem abFactor_lt_one {new old : ℝ} (h : 0<new*old) :
    ModAB.Convergence.abFactor new old<1 := by
  have hr : 0<new/old := by
    rcases mul_pos_iff.mp h with h|h
    · exact div_pos h.1 h.2
    · exact div_pos_of_neg_of_neg h.1 h.2
  unfold ModAB.Convergence.abFactor
  split_ifs <;> linarith

theorem abFactor_mem_Ioo {new old : ℝ} (h : 0<new*old) :
    ModAB.Convergence.abFactor new old ∈ Ioo (0:ℝ) 1 :=
  ⟨ModAB.Convergence.abFactor_positive new old,abFactor_lt_one h⟩

/-- Mixed residuals do not inherit exactness even for f(x)=x: the left ordinate
is multiplied by 1/2, the right ordinate and the dropped value remain true.
The computed point is not the root 0. This is a data-level example, not a claim
that a full run must encounter this state. -/
theorem scaled_linear_candidate_example :
    candidate (-1) 1 (-(1/2:ℝ)) 1 2 2 true=-(1/4:ℝ) := by
  norm_num [candidate,denominator,rawPoint,normalizedOffset,
    normalizedNumerator,normalizedDenominator,limited,radius,secant]

end
end ModAB.Rational

end RModABBundle8

/- Original module: ModAB.RationalBridge -/
section RModABBundle9

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

end RModABBundle9

/- Original module: ModAB.RationalSpacing -/
section RModABBundle10

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

end RModABBundle10

/- Original module: ModAB.RationalResidual -/
section RModABBundle11

/-!
Residual convergence on the original compact interval. These results do not
assume that the bracket width tends to zero, and they do not assume uniqueness
until the final corollary. They therefore apply to the uncapped residual rule
without silently importing the capped width theorem.
-/
namespace ModAB.RationalResidual
open Set Filter Topology
noncomputable section

/-- Every cluster point of a residual-convergent sequence in the interval is a
zero. The membership hypothesis also supplies the domain of continuity. -/
theorem limit_is_zero {f : ℝ → ℝ} {a b z : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hz : z ∈ Icc a b)
    (hx : ∀ n, x n ∈ Icc a b) (hlim : Tendsto x atTop (𝓝 z))
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0)) : f z = 0 := by
  have hfseq : Tendsto (fun n => f (x n)) atTop (𝓝 (f z)) :=
    (hf z hz).tendsto.comp
      (tendsto_nhdsWithin_iff.mpr ⟨hlim, Filter.Eventually.of_forall hx⟩)
  have he : |f z| = 0 := tendsto_nhds_unique hfseq.abs hres
  exact abs_eq_zero.mp he

/-- Residual convergence gives proximity to the zero set, even if the
enclosing intervals have a positive limiting width or the zero is not unique.
The witness is an actual zero in the original interval. -/
theorem eventually_near_zero_set {f : ℝ → ℝ} {a b : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hx : ∀ n, x n ∈ Icc a b)
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0))
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ n in atTop, ∃ z ∈ Icc a b, f z = 0 ∧ |x n - z| < ε := by
  by_contra h
  have hbad : ∃ᶠ n in atTop,
      ¬ (∃ z ∈ Icc a b, f z = 0 ∧ |x n - z| < ε) := not_eventually.mp h
  obtain ⟨ψ, hψ, hψbad⟩ := extraction_of_frequently_atTop hbad
  obtain ⟨z, hz, φ, hφ, ht⟩ :=
    isCompact_Icc.tendsto_subseq (fun n => hx (ψ n))
  have ht' : Tendsto (fun n => x (ψ (φ n))) atTop (𝓝 z) := ht
  have hr' : Tendsto (fun n => |f (x (ψ (φ n)))|) atTop (𝓝 0) :=
    hres.comp (hψ.tendsto_atTop.comp hφ.tendsto_atTop)
  have hz0 : f z = 0 := limit_is_zero hf hz (fun n => hx (ψ (φ n))) ht' hr'
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp ht' ε hε
  have hd : |x (ψ (φ N)) - z| < ε := by
    simpa only [Real.dist_eq] using hN N (le_refl N)
  exact hψbad (φ N) ⟨z, hz, hz0, hd⟩

/-- With a unique zero in the original interval, residual convergence implies
convergence of the evaluated points themselves. No derivative bound is used. -/
theorem tendsto_unique_zero {f : ℝ → ℝ} {a b r : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hx : ∀ n, x n ∈ Icc a b)
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0))
    (hunique : ∀ z ∈ Icc a b, f z = 0 → z = r) :
    Tendsto x atTop (𝓝 r) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  have he := eventually_near_zero_set hf hx hres hε
  obtain ⟨N, hN⟩ := eventually_atTop.mp he
  refine ⟨N, fun n hn => ?_⟩
  obtain ⟨z, hz, hz0, hd⟩ := hN n hn
  have hzr := hunique z hz hz0
  simpa only [hzr, Real.dist_eq] using hd

end
end ModAB.RationalResidual

end RModABBundle11

/- Original module: ModAB.RationalUncapped -/
section RModABBundle12

/-!
Uncapped continuation for a bracket-preserving trial method.

This module does not impose the old cap on residual exceptions.  A continuing
phase advances its clock after every interpolation update and may continue only
by satisfying the geometric width test or the true-residual decrease test.
Switching, curvature checks, and working ordinates are deliberately unrestricted.
The conclusions concern exact arithmetic; no binary64 claim is made here.
-/
namespace ModAB.RationalUncapped
open Set Filter Topology
open ModAB.Convergence
noncomputable section

inductive Mode where
  | bisection
  | interpolation (scale : ℝ) (clock : ℕ)

/-- `m` is the smaller true endpoint residual before the update. -/
inductive Step (B η m : ℝ) : ℝ → ℝ → Mode → ℝ → ℝ → Mode → Prop
  | bisectStay (w R R') : Step B η m w R .bisection (w/2) R' .bisection
  | bisectStart (w R R') :
      Step B η m w R .bisection (w/2) R' (.interpolation (w/2) 0)
  | width {w w' R R' S t} : 0≤w' → w'≤w →
      w'≤threshold B S (t+1) →
      Step B η m w R (.interpolation S t) w' R' (.interpolation S (t+1))
  | residual {w w' R R' S t} : 0≤w' → w'≤w →
      threshold B S (t+1)<w' → R'<η*m →
      Step B η m w R (.interpolation S t) w' R' (.interpolation S (t+1))
  | stop {w w' R R' S t} : 0≤w' → w'≤w →
      Step B η m w R (.interpolation S t) w' R' .bisection

/-- The actual uncapped continuation decision.  Neither successful width tests
nor residual exceptions reset the phase clock. -/
def nextMode (B η S m w' R' : ℝ) (t : ℕ) : Mode :=
  if w'≤threshold B S (t+1) then .interpolation S (t+1)
  else if R'<η*m then .interpolation S (t+1)
  else .bisection

theorem decision_step {B η S m w R w' R' : ℝ} {t : ℕ}
    (h0 : 0≤w') (hmono : w'≤w) :
    Step B η m w R (.interpolation S t) w' R' (nextMode B η S m w' R' t) := by
  unfold nextMode
  split_ifs with hw hr
  · exact .width h0 hmono hw
  · exact .residual h0 hmono (lt_of_not_ge hw) hr
  · exact .stop h0 hmono

theorem forty_two_test (R m : ℝ) : R<(1/42:ℝ)*m ↔ 42*R<m := by
  constructor <;> intro h <;> linarith

/-- Residual `R n` belongs to the point evaluated most recently.  That point is
an endpoint, so the smaller endpoint residual is at most `R n`. -/
structure Run (B η : ℝ) where
  width : ℕ→ℝ
  residual : ℕ→ℝ
  minimum : ℕ→ℝ
  mode : ℕ→Mode
  nonneg_initial : 0≤width 0
  nonneg_residual : ∀n,0≤residual n
  minimum_le : ∀n,minimum n≤residual n
  next : ∀n,Step B η (minimum n) (width n) (residual n) (mode n)
    (width (n+1)) (residual (n+1)) (mode (n+1))

theorem step_nonnegative {B η m w R w' R' : ℝ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : 0≤w) : 0≤w' := by
  cases hs with
  | bisectStay => positivity
  | bisectStart => positivity
  | width h _ _ => exact h
  | residual h _ _ _ => exact h
  | stop h _ => exact h

theorem run_nonnegative {B η : ℝ} (r : Run B η) (n : ℕ) : 0≤r.width n := by
  induction n with
  | zero => exact r.nonneg_initial
  | succ n ih => exact step_nonnegative (r.next n) ih

theorem step_nonincrease {B η m w R w' R' : ℝ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : 0≤w) : w'≤w := by
  cases hs with
  | bisectStay => linarith
  | bisectStart => linarith
  | width _ h _ => exact h
  | residual _ h _ _ => exact h
  | stop _ h => exact h

theorem width_antitone {B η : ℝ} (r : Run B η) : Antitone r.width :=
  antitone_nat_of_succ_le (fun n=>step_nonincrease (r.next n) (run_nonnegative r n))

def isBisection : Mode→Bool
  | .bisection=>true
  | .interpolation _ _=>false

def bisectionCount (m : ℕ→Mode) : ℕ→ℕ
  | 0=>0
  | n+1=>bisectionCount m n + if isBisection (m n) then 1 else 0

theorem count_monotone (m : ℕ→Mode) : Monotone (bisectionCount m) := by
  apply monotone_nat_of_le_succ
  intro n
  simp only [bisectionCount]
  omega

theorem step_count_bound {B η m W0 w R w' R' : ℝ} {k : ℕ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : w≤W0/(2:ℝ)^k) :
    w'≤W0/(2:ℝ)^(k+if isBisection s then 1 else 0) := by
  cases hs with
  | bisectStay =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | bisectStart =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | width _ hm _ => simpa [isBisection] using hm.trans hw
  | residual _ hm _ _ => simpa [isBisection] using hm.trans hw
  | stop _ hm => simpa [isBisection] using hm.trans hw

theorem bisection_count_bound {B η : ℝ} (r : Run B η) (n : ℕ) :
    r.width n≤r.width 0/(2:ℝ)^(bisectionCount r.mode n) := by
  induction n with
  | zero => simp [bisectionCount]
  | succ n ih => exact step_count_bound (r.next n) ih

theorem count_unbounded {m : ℕ→Mode}
    (h : ∀N,∃n,N≤n ∧ isBisection (m n)=true) : ∀k,∃n,k≤bisectionCount m n := by
  intro k
  induction k with
  | zero => exact ⟨0,by simp [bisectionCount]⟩
  | succ k ih =>
      obtain ⟨n,hn⟩:=ih
      obtain ⟨j,hnj,hj⟩:=h n
      refine ⟨j+1,?_⟩
      have hm:=count_monotone m hnj
      simp only [bisectionCount,hj,ite_true]
      omega

theorem infinite_bisections {B η : ℝ} (r : Run B η)
    (h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  have hc : Tendsto (bisectionCount r.mode) atTop atTop := by
    apply tendsto_atTop_atTop.mpr
    intro k
    obtain ⟨n,hn⟩:=count_unbounded h k
    exact ⟨n,fun j hj=>hn.trans (count_monotone r.mode hj)⟩
  exact squeeze_zero (run_nonnegative r) (bisection_count_bound r)
    ((dyadic_tendsto _).comp hc)

/-- There is no reset of a continuing phase's scale or clock. -/
theorem continuing_step {B η m w R w' R' S : ℝ} {t : ℕ} {s' : Mode}
    (hs : Step B η m w R (.interpolation S t) w' R' s')
    (hno : isBisection s'=false) :
    s'=.interpolation S (t+1) ∧
      (w'≤threshold B S (t+1) ∨ R'<η*m) := by
  cases hs with
  | width _ _ ht => exact ⟨rfl,Or.inl ht⟩
  | residual _ _ _ hr => exact ⟨rfl,Or.inr hr⟩
  | stop => simp [isBisection] at hno

theorem continuing_suffix {B η S : ℝ} {N t : ℕ} (r : Run B η)
    (hm : r.mode N=.interpolation S t)
    (hno : ∀n,N≤n→isBisection (r.mode n)=false) :
    ∀k,r.mode (N+k)=.interpolation S (t+k) := by
  intro k
  induction k with
  | zero => simpa using hm
  | succ k ih =>
      have hs:=r.next (N+k)
      rw [ih] at hs
      have he:=(continuing_step hs (hno (N+k+1) (by omega))).1
      simpa [Nat.add_assoc] using he

theorem small_widths_tendsto {B η : ℝ} (r : Run B η)
    (h : ∀ε:ℝ,0<ε→∃n,r.width n<ε) : Tendsto r.width atTop (𝓝 0) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩:=h ε hε
  refine ⟨N,fun n hn=>?_⟩
  rw [Real.dist_eq,sub_zero,abs_of_nonneg (run_nonnegative r n)]
  exact (width_antitone r hn).trans_lt hN

theorem geometric_residual {R : ℕ→ℝ} {η : ℝ} (hη0 : 0≤η) (hη1 : η<1)
    (hR : ∀n,0≤R n) (hstep : ∀n,R (n+1)≤η*R n) :
    Tendsto R atTop (𝓝 0) := by
  have hb : ∀n,R n≤η^n*R 0 := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
      calc R (n+1)≤η*R n:=hstep n
           _≤η*(η^n*R 0):=mul_le_mul_of_nonneg_left ih hη0
           _=η^(n+1)*R 0:=by rw [pow_succ]; ring
  exact squeeze_zero hR hb
    (by simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one hη0 hη1).mul_const (R 0))

/-- The uncapped policy guarantees width convergence or residual convergence.
This is derived from the transition system, not postulated as a run property. -/
theorem width_or_residual {B η : ℝ} (hη0 : 0≤η) (hη1 : η<1) (r : Run B η) :
    Tendsto r.width atTop (𝓝 0) ∨ Tendsto r.residual atTop (𝓝 0) := by
  by_cases hsmall : ∀ε:ℝ,0<ε→∃n,r.width n<ε
  · exact Or.inl (small_widths_tendsto r hsmall)
  push Not at hsmall
  obtain ⟨ε,hε,hfloor⟩:=hsmall
  by_cases hb : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true
  · exact Or.inl (infinite_bisections r hb)
  push Not at hb
  obtain ⟨N,hN⟩:=hb
  have hno : ∀n,N≤n→isBisection (r.mode n)=false := by
    intro n hn
    exact Bool.eq_false_iff.mpr (hN n hn)
  obtain ⟨S,t,hm⟩ : ∃S t,r.mode N=.interpolation S t := by
    cases he : r.mode N with
    | bisection => simpa [he,isBisection] using hno N (le_refl N)
    | interpolation S t => exact ⟨S,t,rfl⟩
  have hmode:=continuing_suffix r hm hno
  have hlim : Tendsto (fun k:ℕ=>threshold B S (t+k+1)) atTop (𝓝 0) := by
    convert dyadic_tendsto (2*B*S/(2:ℝ)^(t+1)) using 1
    funext k
    simp only [threshold,show t+k+1=(t+1)+k by omega,pow_add,div_mul_eq_div_div]
  obtain ⟨K,hK⟩:=Metric.tendsto_atTop.mp hlim ε hε
  have hcontract : ∀k,r.residual (N+K+k+1)≤η*r.residual (N+K+k) := by
    intro k
    have hs:=r.next (N+(K+k))
    rw [hmode (K+k)] at hs
    have hor:=(continuing_step hs (hno (N+(K+k)+1) (by omega))).2
    have ht : threshold B S (t+(K+k)+1)<ε := by
      have hk:=hK (K+k) (by omega)
      rw [Real.dist_eq,sub_zero] at hk
      exact (le_abs_self _).trans_lt hk
    rcases hor with hw|hr
    · have hf:=hfloor (N+(K+k)+1)
      linarith
    · have hh:=hr.le.trans (mul_le_mul_of_nonneg_left (r.minimum_le (N+(K+k))) hη0)
      simpa [Nat.add_assoc] using hh
  have hr : Tendsto (fun k=>r.residual (N+K+k)) atTop (𝓝 0) :=
    geometric_residual hη0 hη1 (fun k=>r.nonneg_residual _) (by simpa [Nat.add_assoc] using hcontract)
  exact Or.inr ((tendsto_add_atTop_iff_nat (N+K)).mp (by simpa [Nat.add_comm] using hr))

/-- Any positive width-or-residual stopping tolerances are attained.  Continuity
is unnecessary for this scheduling statement. -/
theorem finite_combined_tolerance {B η τw τf : ℝ} (hη0 : 0≤η) (hη1 : η<1)
    (hτw : 0<τw) (hτf : 0<τf) (r : Run B η) :
    ∃n,r.width n<τw ∨ r.residual n<τf := by
  rcases width_or_residual hη0 hη1 r with hw|hr
  · obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw τw hτw
    refine ⟨N,Or.inl ?_⟩
    have h:=hN N (le_refl N)
    simpa [Real.dist_eq,abs_of_nonneg (run_nonnegative r N)] using h
  · obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hr τf hτf
    refine ⟨N,Or.inr ?_⟩
    have h:=hN N (le_refl N)
    simpa [Real.dist_eq,abs_of_nonneg (r.nonneg_residual N)] using h

/-- A nonterminating exact run.  `point 0` is an initial endpoint and `point
(n+1)` is the new trial inserted by update `n`. -/
structure ExactRun (B η : ℝ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→Mode
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  initial_point : point 0=a 0 ∨ point 0=b 0
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  candidate_mem : ∀n,point (n+1)∈Icc (a n) (b n)
  nonzero : ∀n,f (point (n+1))≠0
  update : ∀n,Updated f (a n) (b n) (point (n+1)) (a (n+1)) (b (n+1))
  policy : ∀n,Step B η (min |f (a n)| |f (b n)|)
    (b n-a n) |f (point n)| (mode n)
    (b (n+1)-a (n+1)) |f (point (n+1))| (mode (n+1))

theorem exact_invariants {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) :
    ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 ∧
      Icc (r.a (n+1)) (r.b (n+1))⊆Icc (r.a n) (r.b n) := by
  have hs : ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 := by
    intro n
    induction n with
    | zero => exact ⟨r.initial_order,r.initial_sign⟩
    | succ n ih =>
      exact (ModAB.RationalBridge.updated_invariants_closed ih.1 (r.candidate_mem n)
        ih.2 (r.nonzero n) (r.update n)).imp_right And.left
  intro n
  exact ⟨(hs n).1,(hs n).2,
    (ModAB.RationalBridge.updated_invariants_closed (hs n).1 (r.candidate_mem n)
      (hs n).2 (r.nonzero n) (r.update n)).2.2⟩

theorem point_endpoint {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n=r.a n ∨ r.point n=r.b n := by
  cases n with
  | zero => exact r.initial_point
  | succ n =>
    have hu:=r.update n
    unfold Updated at hu
    split_ifs at hu with hs
    · exact Or.inr (congrArg Prod.snd hu).symm
    · exact Or.inl (congrArg Prod.fst hu).symm

theorem point_mem {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n∈Icc (r.a n) (r.b n) := by
  rcases point_endpoint r n with h|h
  · rw [h]; exact ⟨le_refl _,(exact_invariants r n).1.le⟩
  · rw [h]; exact ⟨(exact_invariants r n).1.le,le_refl _⟩

theorem point_mem_initial {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n∈Icc (r.a 0) (r.b 0) := by
  have hm : Antitone (fun n=>Icc (r.a n) (r.b n)) :=
    antitone_nat_of_succ_le (fun n=>(exact_invariants r n).2.2)
  exact hm (Nat.zero_le n) (point_mem r n)

def ExactRun.policyRun {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) : Run B η where
  width n:=r.b n-r.a n
  residual n:=|f (r.point n)|
  minimum n:=min |f (r.a n)| |f (r.b n)|
  mode:=r.mode
  nonneg_initial:=by linarith [r.initial_order]
  nonneg_residual n:=abs_nonneg _
  minimum_le n:=by
    rcases point_endpoint r n with h|h
    · rw [h]; exact min_le_left _ _
    · rw [h]; exact min_le_right _ _
  next:=r.policy

theorem exact_zeros {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) :
    ∀n,∃z∈Icc (r.a n) (r.b n),f z=0 := by
  have hi:=exact_invariants r
  have hm : Antitone (fun n=>Icc (r.a n) (r.b n)):=
    antitone_nat_of_succ_le (fun n=>(hi n).2.2)
  intro n
  exact bracket_contains_zero (hi n).1.le
    (r.continuous.mono (hm (Nat.zero_le n))) (hi n).2.1

theorem residual_of_width {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f)
    (hw : Tendsto (fun n=>r.b n-r.a n) atTop (𝓝 0)) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) := by
  obtain ⟨z,hz,hfz⟩:=nested_intervals_root
    (fun n=>(exact_invariants r n).2.2) (exact_zeros r) r.continuous hw
  have hzall : ∀n,z∈Icc (r.a n) (r.b n) := by
    apply mem_iInter.mp
    rw [hz]
    exact mem_singleton _
  have hp:=chosen_points_tendsto hzall (point_mem r) hw
  have hf : Tendsto (fun n=>f (r.point n)) atTop (𝓝 (f z)) :=
    (r.continuous z (hzall 0)).tendsto.comp
      (tendsto_nhdsWithin_iff.mpr ⟨hp,Filter.Eventually.of_forall (point_mem_initial r)⟩)
  simpa [hfz] using hf.abs

/-- Full residual convergence for continuous sign-bracketed functions and the
actual uncapped width-or-residual continuation schedule. -/
theorem global_residual_convergence {B η : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) := by
  rcases width_or_residual hη0 hη1 r.policyRun with hw|hr
  · exact residual_of_width r hw
  · exact hr

theorem finite_residual_tolerance {B η τ : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (hτ : 0<τ) (r : ExactRun B η f) :
    ∃N,∀n,N≤n→|f (r.point n)|<τ := by
  obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp (global_residual_convergence hη0 hη1 r) τ hτ
  refine ⟨N,fun n hn=>?_⟩
  simpa [Real.dist_eq] using hN n hn

/-- The concrete factor used by the uncapped implementation. -/
theorem eta42_global_residual {B : ℝ} {f : ℝ→ℝ} (r : ExactRun B (1/42) f) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) :=
  global_residual_convergence (by norm_num) (by norm_num) r

theorem eventually_near_zero_set {B η ε : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (hε : 0<ε) (r : ExactRun B η f) :
    ∀ᶠn in atTop,∃z∈Icc (r.a 0) (r.b 0),f z=0 ∧ |r.point n-z|<ε :=
  ModAB.RationalResidual.eventually_near_zero_set r.continuous (point_mem_initial r)
    (global_residual_convergence hη0 hη1 r) hε

theorem unique_root_convergence {B η z : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f)
    (hunique : ∀x∈Icc (r.a 0) (r.b 0),f x=0→x=z) :
    Tendsto r.point atTop (𝓝 z) :=
  ModAB.RationalResidual.tendsto_unique_zero r.continuous (point_mem_initial r)
    (global_residual_convergence hη0 hη1 r) hunique

/-- Conditional width termination from a positive closing distance.  The
spacing premise refers to accepted/evaluated points; it must be checked for
the implementation's closing-step and rounding guards. -/
theorem width_stop_with_closing {B η z δ τ : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f)
    (hunique : ∀x∈Icc (r.a 0) (r.b 0),f x=0→x=z) (hδ : 0<δ)
    (hclose : ∀ᶠn in atTop,τ<r.b n-r.a n→δ≤|r.point (n+1)-r.point n|) :
    ∃n,r.b n-r.a n≤τ :=
  exists_width_stop_of_converges (unique_root_convergence hη0 hη1 r hunique) hδ hclose

end
end ModAB.RationalUncapped

end RModABBundle12

/- Original module: ModAB.RationalAssembly -/
section RModABBundle13

/-!
Operational assembly of the exact-real uncapped method.  Candidate membership
and the scheduling transition are proved from the proposal, closing operation,
true-sign update, and mode decision.  They are not assumed as independent
invariants.  The source proposal relation includes the right-endpoint fallback.
This is an exact-real model of those branches, not a full IEEE refinement.
-/
namespace ModAB.RationalAssembly
open Set Filter Topology
open ModAB.Convergence ModAB.Rational
noncomputable section

def RecentEndpoint (a b latest : ℝ) : EndpointSide→Prop
  | .left => latest=a
  | .right => latest=b

def stepPoint (a b tolerance pre : ℝ) (recent : EndpointSide) :
    ModAB.RationalUncapped.Mode→ℝ
  | .bisection => (a+b)/2
  | .interpolation _ _ => closeCandidate a b tolerance pre recent

/-- The source proposal relation, including its endpoint fallback, supplies
closed membership after the guarded closing operation. -/
theorem source_closed_update {f : ℝ→ℝ} {a b ya yb pre tolerance x a' b' : ℝ}
    (recent : EndpointSide) (hab : a<b) (hs : f a*f b<0)
    (hya : SameSign ya (f a)) (hyb : SameSign yb (f b))
    (hp : SourceProposal a b ya yb pre) (ht : 0<tolerance)
    (hx : x=closeCandidate a b tolerance pre recent)
    (hz : f x≠0) (hu : Updated f a b x a' b') :
    a'<b' ∧ f a'*f b'<0 ∧ Icc a' b'⊆Icc a b := by
  have hpre:=sourceProposal_mem_closed hab (working_opposite hya hyb hs) hp
  have hc : x∈Icc a b := by
    rw [hx]
    exact closeCandidate_mem_closed recent ht hpre
  exact ModAB.RationalBridge.updated_invariants_closed hab hc hs hz hu

/-- Operational data.  The phase clock advances through `nextMode`; neither
bracket membership nor width contraction is postulated.  The working-sign
premises are supplied by positive AB corrections in the ordinate-state layer. -/
structure SourceRun (B η : ℝ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→ModAB.RationalUncapped.Mode
  ya : ℕ→ℝ
  yb : ℕ→ℝ
  preproposal : ℕ→ℝ
  tolerance : ℕ→ℝ
  recent : ℕ→EndpointSide
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  initial_point : point 0=a 0 ∨ point 0=b 0
  initial_mode : mode 0=.bisection
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  working_left : ∀n,SameSign (ya n) (f (a n))
  working_right : ∀n,SameSign (yb n) (f (b n))
  proposal : ∀n,SourceProposal (a n) (b n) (ya n) (yb n) (preproposal n)
  tolerance_pos : ∀n,0<tolerance n
  recent_correct : ∀n S t,mode n=.interpolation S t→
    RecentEndpoint (a n) (b n) (point n) (recent n)
  point_rule : ∀n,point (n+1)=stepPoint (a n) (b n) (tolerance n)
    (preproposal n) (recent n) (mode n)
  nonzero : ∀n,f (point (n+1))≠0
  update : ∀n,Updated f (a n) (b n) (point (n+1)) (a (n+1)) (b (n+1))
  bisect_mode : ∀n,mode n=.bisection→
    mode (n+1)=.bisection ∨ mode (n+1)=.interpolation ((b n-a n)/2) 0
  continue_mode : ∀n S t,mode n=.interpolation S t→
    mode (n+1)=ModAB.RationalUncapped.nextMode B η S
      (min |f (a n)| |f (b n)|) (b (n+1)-a (n+1)) |f (point (n+1))| t

theorem proposal_mem_of_invariants {B η : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (n : ℕ) (hab : s.a n<s.b n)
    (hs : f (s.a n)*f (s.b n)<0) : s.point (n+1)∈Icc (s.a n) (s.b n) := by
  rw [s.point_rule n]
  cases hm : s.mode n with
  | bisection =>
    simp only [stepPoint]
    constructor <;> linarith
  | interpolation S t =>
    simp only [stepPoint]
    exact closeCandidate_mem_closed (s.recent n) (s.tolerance_pos n)
      (sourceProposal_mem_closed hab
        (working_opposite (s.working_left n) (s.working_right n) hs) (s.proposal n))

theorem source_invariants {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) :
    ∀n,s.a n<s.b n ∧ f (s.a n)*f (s.b n)<0 := by
  intro n
  induction n with
  | zero => exact ⟨s.initial_order,s.initial_sign⟩
  | succ n ih =>
    have h:=ModAB.RationalBridge.updated_invariants_closed ih.1
      (proposal_mem_of_invariants s n ih.1 ih.2) ih.2 (s.nonzero n) (s.update n)
    exact ⟨h.1,h.2.1⟩

theorem source_candidate_mem {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    s.point (n+1)∈Icc (s.a n) (s.b n) :=
  proposal_mem_of_invariants s n (source_invariants s n).1 (source_invariants s n).2

theorem source_width_step {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    0≤s.b (n+1)-s.a (n+1) ∧ s.b (n+1)-s.a (n+1)≤s.b n-s.a n := by
  have hx:=source_candidate_mem s n
  have hu:=s.update n
  unfold Updated at hu
  split_ifs at hu with hs
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2]
    constructor <;> linarith [hx.1,hx.2]
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2]
    constructor <;> linarith [hx.1,hx.2]

theorem source_bisection_width {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f)
    (n : ℕ) (hm : s.mode n=.bisection) :
    s.b (n+1)-s.a (n+1)=(s.b n-s.a n)/2 := by
  have hp : s.point (n+1)=(s.a n+s.b n)/2 := by
    simpa [stepPoint,hm] using s.point_rule n
  have hu:=s.update n
  unfold Updated at hu
  split_ifs at hu with hs
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2,hp]
    ring
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2,hp]
    ring

theorem source_policy {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    ModAB.RationalUncapped.Step B η (min |f (s.a n)| |f (s.b n)|)
      (s.b n-s.a n) |f (s.point n)| (s.mode n)
      (s.b (n+1)-s.a (n+1)) |f (s.point (n+1))| (s.mode (n+1)) := by
  cases hm : s.mode n with
  | bisection =>
    have hw:=source_bisection_width s n hm
    rcases s.bisect_mode n hm with hn|hn
    · rw [hw,hn]
      exact .bisectStay _ _ _
    · rw [hw,hn]
      exact .bisectStart _ _ _
  | interpolation S t =>
    rw [s.continue_mode n S t hm]
    exact ModAB.RationalUncapped.decision_step
      (source_width_step s n).1 (source_width_step s n).2

/-- The abstract uncapped run is built from the operational proposal and tests.
Its candidate-membership and policy fields are derived theorems above. -/
def SourceRun.exactRun {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) :
    ModAB.RationalUncapped.ExactRun B η f where
  a:=s.a
  b:=s.b
  point:=s.point
  mode:=s.mode
  initial_order:=s.initial_order
  initial_sign:=s.initial_sign
  initial_point:=s.initial_point
  continuous:=s.continuous
  candidate_mem:=source_candidate_mem s
  nonzero:=s.nonzero
  update:=s.update
  policy:=source_policy s

theorem source_global_residual {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f)
    (hη0 : 0≤η) (hη1 : η<1) :
    Tendsto (fun n=>|f (s.point n)|) atTop (𝓝 0) :=
  ModAB.RationalUncapped.global_residual_convergence hη0 hη1 s.exactRun

theorem source_unique_root_convergence {B η z : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (hη0 : 0≤η) (hη1 : η<1)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z) :
    Tendsto s.point atTop (𝓝 z) :=
  ModAB.RationalUncapped.unique_root_convergence hη0 hη1 s.exactRun hunique

/-- Actual closing-operation equations give the accepted-point spacing. -/
theorem source_closing_separation {B η : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (n : ℕ) {S : ℝ} {t : ℕ}
    (hm : s.mode n=.interpolation S t) (hw : s.tolerance n<s.b n-s.a n) :
    s.tolerance n/2≤|s.point (n+1)-s.point n| := by
  have hp : s.point (n+1)=closeCandidate (s.a n) (s.b n) (s.tolerance n)
      (s.preproposal n) (s.recent n) := by
    simpa [stepPoint,hm] using s.point_rule n
  have hrecent:=s.recent_correct n S t hm
  cases he : s.recent n with
  | left =>
    have hold : s.point n=s.a n := by simpa [RecentEndpoint,he] using hrecent
    have hsep:=closeCandidate_left_separation (x:=s.preproposal n) (s.tolerance_pos n) hw
    rw [hp,he,hold]
    exact (show s.tolerance n/2≤
      closeCandidate (s.a n) (s.b n) (s.tolerance n) (s.preproposal n) .left-s.a n by linarith).trans
      (le_abs_self _)
  | right =>
    have hold : s.point n=s.b n := by simpa [RecentEndpoint,he] using hrecent
    have hsep:=closeCandidate_right_separation (x:=s.preproposal n) (s.tolerance_pos n) hw
    rw [hp,he,hold,abs_sub_comm]
    exact (show s.tolerance n/2≤
      s.b n-closeCandidate (s.a n) (s.b n) (s.tolerance n) (s.preproposal n) .right by linarith).trans
      (le_abs_self _)

/-- Dynamic-tolerance width termination for the actual closing operation in
the exact-real model, under uniqueness and a positive absolute tolerance floor.
There is no assumed spacing inequality in this theorem. -/
theorem source_finite_width_tolerance {B η z τ0 : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (hη0 : 0≤η) (hη1 : η<1)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z)
    (hτ0 : 0<τ0) (hfloor : ∀n,τ0≤s.tolerance n) :
    ∃n,s.b n-s.a n≤s.tolerance n := by
  by_cases hb : ∀N,∃n,N≤n ∧ ModAB.RationalUncapped.isBisection (s.mode n)=true
  · have hw:=ModAB.RationalUncapped.infinite_bisections s.exactRun.policyRun hb
    obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw τ0 hτ0
    have hsmall:=hN N (le_refl N)
    change |s.b N-s.a N-0|<τ0 at hsmall
    have hn : s.b N-s.a N<τ0 := by
      exact (le_abs_self _).trans_lt (by simpa using hsmall)
    exact ⟨N,hn.le.trans (hfloor N)⟩
  push Not at hb
  obtain ⟨N,hN⟩:=hb
  have hp:=source_unique_root_convergence s hη0 hη1 hunique
  by_contra hn
  push Not at hn
  apply ModAB.RationalUncapped.no_convergent_separated_tail hp
    (show 0<τ0/2 by linarith)
  apply eventually_atTop.mpr
  refine ⟨N,fun n hNn=>?_⟩
  cases hm : s.mode n with
  | bisection =>
    have hfalse:=hN n hNn
    simp [hm,ModAB.RationalUncapped.isBisection] at hfalse
  | interpolation S t =>
    exact (show τ0/2≤s.tolerance n/2 by linarith [hfloor n]).trans
      (source_closing_separation s n hm (hn n))

theorem eta42_source_finite_width {B z τ0 : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B (1/42) f)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z)
    (hτ0 : 0<τ0) (hfloor : ∀n,τ0≤s.tolerance n) :
    ∃n,s.b n-s.a n≤s.tolerance n :=
  source_finite_width_tolerance s (by norm_num) (by norm_num) hunique hτ0 hfloor

end
end ModAB.RationalAssembly

end RModABBundle13

/- Original module: ModAB.RationalLocalAlgebra -/
section RModABBundle14

/-!
Exact identities behind the local rational-step analysis.  The nodes and values
are arbitrary real data.  The theorems state every required nonvanishing
condition; they do not assume differentiability or an asymptotic recurrence.
The positivity result at the end is the finite-data crossing mechanism used in
the paper's local K > 0 argument.  Passing from C³ hypotheses to the signs of
divided differences is a separate analytic statement.
-/
namespace ModAB.Rational.LocalAlgebra
noncomputable section

/-- First divided difference, with the function values supplied explicitly. -/
def dd1 (x₁ x₂ v₁ v₂ : ℝ) : ℝ := (v₂-v₁)/(x₂-x₁)

/-- Second divided difference in Newton form. -/
def dd2 (x₁ x₂ x₃ v₁ v₂ v₃ : ℝ) : ℝ :=
  (dd1 x₂ x₃ v₂ v₃-dd1 x₁ x₂ v₁ v₂)/(x₃-x₁)

/-- Barycentric quotient specifying the root of the linear-fractional model.
Its denominator is required nonzero in all theorems that interpret it as a root. -/
def baryRoot (x₁ x₂ x₃ y₁ y₂ y₃ : ℝ) : ℝ :=
  dd2 x₁ x₂ x₃ (x₁/y₁) (x₂/y₂) (x₃/y₃) /
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)

/-- The Newton quadratic through g evaluated at the proposed root r. -/
def poleValue (x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ) : ℝ :=
  g₁-dd1 x₁ x₂ g₁ g₂*(x₁-r)+
    dd2 x₁ x₂ x₃ g₁ g₂ g₃*(x₁-r)*(x₂-r)

theorem dd2_add (x₁ x₂ x₃ a₁ a₂ a₃ b₁ b₂ b₃ : ℝ) :
    dd2 x₁ x₂ x₃ (a₁+b₁) (a₂+b₂) (a₃+b₃) =
      dd2 x₁ x₂ x₃ a₁ a₂ a₃+dd2 x₁ x₂ x₃ b₁ b₂ b₃ := by
  unfold dd2 dd1
  ring

theorem dd2_smul (x₁ x₂ x₃ c a₁ a₂ a₃ : ℝ) :
    dd2 x₁ x₂ x₃ (c*a₁) (c*a₂) (c*a₃) =
      c*dd2 x₁ x₂ x₃ a₁ a₂ a₃ := by
  unfold dd2 dd1
  ring

theorem dd2_barycentric {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ =
      v₁/((x₁-x₂)*(x₁-x₃))+
      v₂/((x₂-x₁)*(x₂-x₃))+v₃/((x₃-x₁)*(x₃-x₂)) := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  have h₁₃ : x₁-x₃≠0 := by intro h; apply h₃₁; linarith
  have h₂₃ : x₂-x₃≠0 := by intro h; apply h₃₂; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂,h₁₃,h₂₃]
  ring

theorem dd2_affine {x₁ x₂ x₃ p q : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) :
    dd2 x₁ x₂ x₃ (p*x₁+q) (p*x₂+q) (p*x₃+q)=0 := by
  have hleft : dd1 x₁ x₂ (p*x₁+q) (p*x₂+q)=p := by
    unfold dd1
    field_simp [h₂₁]
    ring
  have hright : dd1 x₂ x₃ (p*x₂+q) (p*x₃+q)=p := by
    unfold dd1
    field_simp [h₃₂]
    ring
  simp [dd2,hleft,hright]

theorem dd2_swap_last {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₁ x₃ x₂ v₁ v₃ v₂ := by
  have h₂₃ : x₂-x₃≠0 := by intro h; apply h₃₂; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₂₃]
  ring

theorem dd2_swap_first {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₂ x₁ x₃ v₂ v₁ v₃ := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂]
  ring

theorem dd2_cyclic {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₂ x₃ x₁ v₂ v₃ v₁ := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  have h₁₃ : x₁-x₃≠0 := by intro h; apply h₃₁; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂,h₁₃]
  ring

theorem baryRoot_swap_last {x₁ x₂ x₃ y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ y₃ = baryRoot x₁ x₃ x₂ y₁ y₃ y₂ := by
  unfold baryRoot
  rw [dd2_swap_last h₂₁ h₃₂ h₃₁,dd2_swap_last h₂₁ h₃₂ h₃₁]

theorem dd2_change_third {x₁ x₂ x₃ v₁ v₂ v₃ Δ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ (v₃+Δ) =
      dd2 x₁ x₂ x₃ v₁ v₂ v₃+Δ/((x₃-x₁)*(x₃-x₂)) := by
  unfold dd2 dd1
  field_simp [h₃₂,h₃₁]
  ring

/-- This linear identity does not require nonzero ordinate or node denominators. -/
theorem dd2_error_decomposition (x₁ x₂ x₃ y₁ y₂ y₃ r : ℝ) :
    dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) =
      dd2 x₁ x₂ x₃ (x₁/y₁) (x₂/y₂) (x₃/y₃)-
        r*dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃) := by
  unfold dd2 dd1
  ring

theorem baryRoot_error_identity {x₁ x₂ x₃ y₁ y₂ y₃ r : ℝ}
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)≠0) :
    (baryRoot x₁ x₂ x₃ y₁ y₂ y₃-r)*
      dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃) =
        dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) := by
  rw [dd2_error_decomposition]
  unfold baryRoot
  field_simp [hH]

/-- The pole-part identity is exact, not a Taylor approximation. -/
theorem pole_identity {x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0) :
    dd2 x₁ x₂ x₃ (g₁/(x₁-r)) (g₂/(x₂-r)) (g₃/(x₃-r))*
      ((x₁-r)*(x₂-r)*(x₃-r)) = poleValue x₁ x₂ x₃ r g₁ g₂ g₃ := by
  unfold poleValue dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁,h₂,h₃]
  ring

theorem pole_identity_from_values {x₁ x₂ x₃ r y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0) :
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)*((x₁-r)*(x₂-r)*(x₃-r)) =
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) := by
  have h := pole_identity (g₁:=(x₁-r)/y₁) (g₂:=(x₂-r)/y₂)
    (g₃:=(x₃-r)/y₃) h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  have hc (e y : ℝ) (he : e≠0) : (e/y)/e=1/y := by
    calc (e/y)/e = (e*e⁻¹)*y⁻¹ := by ring
         _ = 1/y := by simp [he]
  simpa only [hc (x₁-r) y₁ h₁,hc (x₂-r) y₂ h₂,hc (x₃-r) y₃ h₃] using h

theorem poleValue_margin (x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ) :
    poleValue x₁ x₂ x₃ r g₁ g₂ g₃-
      ((x₁-r)*(x₂-r))*dd2 x₁ x₂ x₃ g₁ g₂ g₃ =
        g₁-dd1 x₁ x₂ g₁ g₂*(x₁-r) := by
  unfold poleValue
  ring

/-- Exact error product for unscaled values.  Analytic convergence of its
coefficient to K is deliberately not an assumption or conclusion here. -/
theorem baryRoot_true_error {x₁ x₂ x₃ r y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ y₃-r =
      dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) /
        poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)*
          ((x₁-r)*(x₂-r)*(x₃-r)) := by
  have hp := pole_identity_from_values (y₁:=y₁) (y₂:=y₂) (y₃:=y₃)
    h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  rw [←hp]
  rw [dd2_error_decomposition]
  unfold baryRoot
  field_simp [hH,h₁,h₂,h₃]

/-- Replacing one ordinate by μ times that ordinate changes its reciprocal
by the exact relative increment δ=1/μ-1. -/
theorem rescale_div {a y μ : ℝ} (hy : y≠0) (hμ : μ≠0) :
    a/(μ*y)=a/y+(1/μ-1)*(a/y) := by
  field_simp [hy,hμ]
  ring

theorem dd2_scaled_third {x₁ x₂ x₃ a₁ a₂ a₃ y₁ y₂ y₃ μ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) (hy₃ : y₃≠0) (hμ : μ≠0) :
    dd2 x₁ x₂ x₃ (a₁/y₁) (a₂/y₂) (a₃/(μ*y₃)) =
      dd2 x₁ x₂ x₃ (a₁/y₁) (a₂/y₂) (a₃/y₃)+
        (1/μ-1)*(a₃/y₃)/((x₃-x₁)*(x₃-x₂)) := by
  rw [rescale_div hy₃ hμ,dd2_change_third h₃₂ h₃₁]

theorem poleValue_scaled_third {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) (hy₃ : y₃≠0) (hμ : μ≠0) :
    poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/(μ*y₃)) =
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)) := by
  unfold poleValue
  rw [dd2_scaled_third h₃₂ h₃₁ hy₃ hμ]
  ring

/-- Exact mixed-ordinate error formula.  This retains the denominator and
therefore remains valid even when two possible leading terms cancel. -/
theorem baryRoot_scaled_error {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃)-r =
      (dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)/((x₃-x₁)*(x₃-x₂))) /
      (poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)))*
      ((x₁-r)*(x₂-r)*(x₃-r)) := by
  rw [baryRoot_true_error h₂₁ h₃₂ h₃₁ h₁ h₂ h₃ hH,
    dd2_scaled_third h₃₂ h₃₁ hy₃ hμ,
    poleValue_scaled_third h₃₂ h₃₁ hy₃ hμ]

/-- A finite-data crossing lemma: the rational error is a strict fraction of
the error at the opposite endpoint, uniformly over nonnegative scaling δ. -/
theorem scaled_crossing_ratio {E G P g δ π : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) :
    E*(G+δ*g/π)/(P+E*(δ*g/π)) ∈ Set.Ioo (0:ℝ) 1 := by
  have ht : 0≤δ*g/π := div_nonneg (mul_nonneg hδ hg.le) hπ.le
  have hP : 0<P := by nlinarith [mul_pos hE hG]
  have hD : 0<P+E*(δ*g/π) := by nlinarith [mul_nonneg hE.le ht]
  constructor
  · exact div_pos (mul_pos hE (by linarith)) hD
  · apply (div_lt_one hD).2
    nlinarith

theorem scaled_pole_pos {E G P g δ π : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) : 0<P+δ*g*E/π := by
  have ht : 0≤δ*g*E/π :=
    div_nonneg (mul_nonneg (mul_nonneg hδ hg.le) hE.le) hπ.le
  nlinarith [mul_pos hE hG]

/-- Under the local positive-coefficient data conditions, the mixed model's
denominator is nonzero; this condition need not be assumed independently. -/
theorem scaled_denominator_ne_zero {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hE : 0<(x₁-r)*(x₂-r))
    (hG : 0<dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃))
    (hg : 0<(x₃-r)/y₃) (hδ : 0≤1/μ-1)
    (hπ : 0<(x₃-x₁)*(x₃-x₂))
    (hmargin : 0<
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)-
        ((x₁-r)*(x₂-r))*
          dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)) :
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0 := by
  have hp := pole_identity_from_values (y₁:=y₁) (y₂:=y₂) (y₃:=μ*y₃)
    h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  rw [poleValue_scaled_third h₃₂ h₃₁ hy₃ hμ] at hp
  have hpos := scaled_pole_pos hE hG hg hδ hπ hmargin
  intro hz
  rw [hz,zero_mul] at hp
  have he : (1/μ-1)*((x₃-r)/y₃)*((x₁-r)*(x₂-r))/((x₃-x₁)*(x₃-x₂)) =
      (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)) := by ring
  rw [he] at hpos
  linarith

theorem scaled_error_positive_side {E G P g δ π es : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) (hes : 0<es) :
    (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) ∈ Set.Ioo (0:ℝ) es := by
  have hr := scaled_crossing_ratio hE hG hg hδ hπ hmargin
  have he : (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) =
      (E*(G+δ*g/π)/(P+E*(δ*g/π)))*es := by ring
  rw [he]
  exact ⟨mul_pos hr.1 hes,by simpa using mul_lt_mul_of_pos_right hr.2 hes⟩

theorem scaled_error_negative_side {E G P g δ π es : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) (hes : es<0) :
    (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) ∈ Set.Ioo es (0:ℝ) := by
  have hr := scaled_crossing_ratio hE hG hg hδ hπ hmargin
  have he : (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) =
      (E*(G+δ*g/π)/(P+E*(δ*g/π)))*es := by ring
  rw [he]
  exact ⟨by simpa using mul_lt_mul_of_neg_right hr.2 hes,mul_neg_of_pos_of_neg hr.1 hes⟩

/-- In a same-side step the scaling cannot grow arbitrarily.  The sign of G
need not be assumed: it follows from the same-side inequality and δ ≥ 0. -/
theorem delta_bound_of_same_side {E G P g δ π : ℝ}
    (hE : 0<E) (hP : 0<P) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π))<0) :
    δ< -π*G/g := by
  have ht : 0≤δ*g/π := div_nonneg (mul_nonneg hδ hg.le) hπ.le
  have hD : 0<P+E*(δ*g/π) := by nlinarith [mul_nonneg hE.le ht]
  have hn : E*(G+δ*g/π)<0 := by
    have h := (div_lt_iff₀ hD).mp hsame
    simpa using h
  have hs : G+δ*g/π<0 := by
    by_contra h
    have hp := mul_nonneg hE.le (le_of_not_gt h)
    linarith
  have hterm : δ*g/π < -G := by linarith
  have hmul := (div_lt_iff₀ hπ).mp hterm
  apply (lt_div_iff₀ hg).2
  nlinarith

theorem negative_coefficient_of_same_side {E G P g δ π : ℝ}
    (hE : 0<E) (hP : 0<P) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π))<0) : G<0 := by
  have h := delta_bound_of_same_side hE hP hg hδ hπ hsame
  have hp : 0< -π*G/g := lt_of_le_of_lt hδ h
  have hn : 0< -π*G := by
    have hh := (lt_div_iff₀ hg).mp hp
    simpa using hh
  nlinarith

/-- An a priori upper bound on δ gives a positive lower bound for the weight. -/
theorem weight_lower_bound {δ M : ℝ} (hδ : 0≤δ) (hM : δ≤M) :
    0<1/(1+M) ∧ 1/(1+M)≤1/(1+δ) := by
  have hδp : 0<1+δ := by linarith
  have hMp : 0<1+M := by linarith
  exact ⟨one_div_pos.mpr hMp,one_div_le_one_div_of_le hδp (by linarith)⟩

theorem weight_from_delta {μ : ℝ} (hμ : μ≠0) : 1/(1+(1/μ-1))=μ := by
  field_simp [hμ]
  ring

/-- The exact mixed-ordinate root crosses toward x₃ and stops strictly before
x₃.  The hypotheses are finite-data signs, rather than an assumed crossing. -/
theorem baryRoot_scaled_between {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0)
    (hE : 0<(x₁-r)*(x₂-r))
    (hG : 0<dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃))
    (hg : 0<(x₃-r)/y₃) (hδ : 0≤1/μ-1)
    (hπ : 0<(x₃-x₁)*(x₃-x₂))
    (hmargin : 0<
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)-
        ((x₁-r)*(x₂-r))*
          dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)) :
    (r<x₃ → baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃) ∈ Set.Ioo r x₃) ∧
    (x₃<r → baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃) ∈ Set.Ioo x₃ r) := by
  let E := (x₁-r)*(x₂-r)
  let G := dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)
  let P := poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)
  let g := (x₃-r)/y₃
  let δ := 1/μ-1
  let π := (x₃-x₁)*(x₃-x₂)
  have he : baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃)-r =
      (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) := by
    rw [baryRoot_scaled_error h₂₁ h₃₂ h₃₁ h₁ h₂ h₃ hy₃ hμ hH]
    dsimp [E,G,P,g,δ,π]
    congr 1
    congr 1
    ring
  constructor
  · intro hs
    have hp := scaled_error_positive_side hE hG hg hδ hπ hmargin (sub_pos.mpr hs)
    change (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) ∈ Set.Ioo 0 (x₃-r) at hp
    rw [←he] at hp
    exact ⟨by linarith [hp.1],by linarith [hp.2]⟩
  · intro hs
    have hp := scaled_error_negative_side hE hG hg hδ hπ hmargin (sub_neg.mpr hs)
    change (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) ∈ Set.Ioo (x₃-r) 0 at hp
    rw [←he] at hp
    exact ⟨by linarith [hp.1],by linarith [hp.2]⟩

/-- The barycentric root agrees with the normalized proposal already used in
the algorithm specification.  Neither quotient is interpreted at a zero
denominator. -/
theorem bary_denominator_eq {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0) :
    dd2 a b p (1/ya) (1/yb) (1/yp) =
      -ModAB.Rational.denominator a b ya yb p yp*yb/
        (ya*yp*(b-a)*(p-b)*(p-a)) := by
  unfold dd2 dd1 ModAB.Rational.denominator ModAB.Rational.normalizedDenominator
  field_simp [hba,hpb,hpa,hya,hyb,hyp]
  ring

theorem bary_denominator_ne_zero {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0) :
    dd2 a b p (1/ya) (1/yb) (1/yp)≠0 := by
  rw [bary_denominator_eq hba hpb hpa hya hyb hyp]
  exact div_ne_zero (mul_ne_zero (neg_ne_zero.mpr hD) hyb)
    (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hya hyp) hba) hpb) hpa)

theorem baryRoot_eq_rawPoint {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0)
    (hH : dd2 a b p (1/ya) (1/yb) (1/yp)≠0) :
    baryRoot a b p ya yb yp = ModAB.Rational.rawPoint a b ya yb p yp := by
  have hD' : ModAB.Rational.normalizedDenominator (b-a) (p-a) (ya/yb) (yp/yb)≠0 := hD
  unfold baryRoot ModAB.Rational.rawPoint ModAB.Rational.normalizedOffset
  field_simp [hH,hD']
  unfold dd2 dd1 ModAB.Rational.normalizedNumerator ModAB.Rational.normalizedDenominator
  field_simp [hba,hpb,hpa,hya,hyb,hyp]
  ring

/-- A source-oriented bridge with just the algorithm's denominator guard. -/
theorem rawPoint_eq_baryRoot {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0) :
    ModAB.Rational.rawPoint a b ya yb p yp = baryRoot a b p ya yb yp := by
  exact (baryRoot_eq_rawPoint hba hpb hpa hya hyb hyp hD
    (bary_denominator_ne_zero hba hpb hpa hya hyb hyp hD)).symm

end
end ModAB.Rational.LocalAlgebra

end RModABBundle14

/- Original module: ModAB.RationalOrderStability -/
section RModABBundle15

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

end RModABBundle15

/- Original module: ModAB.RationalOrder -/
section RModABBundle16

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

end RModABBundle16

/- Original module: ModAB.RationalGoldenOrder -/
section RModABBundle17

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

end RModABBundle17

/- Original module: ModAB.RationalLocalSafeguards -/
section RModABBundle18

/-!
Quantitative local interfaces for the rational interpolation phase.

This file proves the safeguard consequences of explicit local error and
residual bounds. It does not assert that `C³` regularity supplies those bounds,
that an arbitrary global run enters this neighbourhood, or that fixed positive
stopping tolerances are asymptotically inactive. The clamp, zero-tolerance
closing operation and uncapped continuation decision are the existing project
definitions. The recurrence results are an exact-real local model.
-/
namespace ModAB.RationalLocalSafeguards
open Set Filter Topology
open ModAB.Convergence
noncomputable section

/-- Three local errors give a uniform contraction against either endpoint. -/
theorem triple_error_bound {u v d z B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε)
    (hu : |u| ≤ ε) (hv : |v| ≤ ε) (hd : |d| ≤ ε)
    (hz : |z| ≤ B * |u| * |v| * |d|) :
    |z| ≤ (B * ε ^ 2) * min |u| |v| := by
  have hvd : |v| * |d| ≤ ε ^ 2 := by
    calc |v| * |d| ≤ ε * ε := mul_le_mul hv hd (abs_nonneg _) hε
         _ = ε ^ 2 := by ring
  have hud : |u| * |d| ≤ ε ^ 2 := by
    calc |u| * |d| ≤ ε * ε := mul_le_mul hu hd (abs_nonneg _) hε
         _ = ε ^ 2 := by ring
  have hleft : |z| ≤ (B * ε ^ 2) * |u| := by
    calc |z| ≤ (B * |u|) * (|v| * |d|) := by nlinarith [hz]
         _ ≤ (B * |u|) * ε ^ 2 :=
           mul_le_mul_of_nonneg_left hvd (mul_nonneg hB (abs_nonneg _))
         _ = (B * ε ^ 2) * |u| := by ring
  have hright : |z| ≤ (B * ε ^ 2) * |v| := by
    calc |z| ≤ (B * |v|) * (|u| * |d|) := by nlinarith [hz]
         _ ≤ (B * |v|) * ε ^ 2 :=
           mul_le_mul_of_nonneg_left hud (mul_nonneg hB (abs_nonneg _))
         _ = (B * ε ^ 2) * |v| := by ring
  rw [mul_min_of_nonneg |u| |v| (mul_nonneg hB (sq_nonneg ε))]
  exact le_min hleft hright

/-- The usual quadratic secant bound has the matching minimum-error form. -/
theorem pair_error_bound {u v z S ε : ℝ} (hS : 0 ≤ S)
    (hu : |u| ≤ ε) (hv : |v| ≤ ε)
    (hz : |z| ≤ S * |u| * |v|) :
    |z| ≤ (S*ε) * min |u| |v| := by
  have hε : 0 ≤ ε := (abs_nonneg _).trans hu
  have hl := mul_le_mul_of_nonneg_left hv (mul_nonneg hS (abs_nonneg u))
  have hr := mul_le_mul_of_nonneg_left hu (mul_nonneg hS (abs_nonneg v))
  rw [mul_min_of_nonneg |u| |v| (mul_nonneg hS hε)]
  apply le_min <;> nlinarith

/-- The secant error bound implies the lower bound on its actual clamp radius. -/
theorem radius_lower_of_secant_error {a b r ya yb t : ℝ}
    (ha : a ≤ r) (hb : r ≤ b)
    (hs : |ModAB.Rational.secant a b ya yb - r| ≤
      t * min |a-r| |b-r|) :
    (1-t) * min |a-r| |b-r| ≤ ModAB.Rational.radius a b ya yb := by
  have hleft := min_le_left |a-r| |b-r|
  have hright := min_le_right |a-r| |b-r|
  have habsa : |a-r| = r-a := by rw [abs_of_nonpos (sub_nonpos.mpr ha)]; ring
  have habsb : |b-r| = b-r := abs_of_nonneg (sub_nonneg.mpr hb)
  rcases abs_le.mp hs with ⟨hsl, hsr⟩
  unfold ModAB.Rational.radius
  apply le_min
  · nlinarith
  · nlinarith

/-- This is an inactivity theorem for the existing `0.7` clamp, rather than
an independent replacement implementation of a clamp. -/
theorem limited_eq_raw_of_local_estimates {a b r ya yb raw q t : ℝ}
    (ha : a ≤ r) (hb : r ≤ b)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|)
    (hsec : |ModAB.Rational.secant a b ya yb-r| ≤
      t * min |a-r| |b-r|)
    (hsmall : q+t ≤ (7/10:ℝ) * (1-t)) :
    ModAB.Rational.limited a b ya yb raw = raw := by
  have hm : 0 ≤ min |a-r| |b-r| := le_min (abs_nonneg _) (abs_nonneg _)
  have hrad := radius_lower_of_secant_error ha hb hsec
  have hdev : |raw-ModAB.Rational.secant a b ya yb| ≤
      (q+t) * min |a-r| |b-r| := by
    have htri := abs_sub_le raw r (ModAB.Rational.secant a b ya yb)
    rw [abs_sub_comm r (ModAB.Rational.secant a b ya yb)] at htri
    nlinarith
  have hsmallm := mul_le_mul_of_nonneg_right hsmall hm
  have hbound : |raw-ModAB.Rational.secant a b ya yb| ≤
      (7/10:ℝ) * ModAB.Rational.radius a b ya yb := by nlinarith
  rcases abs_le.mp hbound with ⟨hl, hr⟩
  exact ModAB.Rational.limited_eq_self (by linarith) (by linarith)

/-- The local point bound gives the literal factor-42 residual exception. -/
theorem forty_two_residual_of_local_bounds {f : ℝ → ℝ}
    {a b r raw lower upper q : ℝ}
    (hupper : 0 ≤ upper)
    (hm : 0 < min |a-r| |b-r|)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|)
    (hfraw : |f raw| ≤ upper * |raw-r|)
    (hflower : lower * min |a-r| |b-r| ≤ min |f a| |f b|)
    (hsmall : 42 * upper * q < lower) :
    42 * |f raw| < min |f a| |f b| := by
  have hp := mul_le_mul_of_nonneg_left hraw hupper
  have hs := mul_lt_mul_of_pos_right hsmall hm
  nlinarith

/-- The residual exception prevents fallback regardless of the width clock. -/
theorem nextMode_interpolation_of_forty_two {B S m w R : ℝ} (t : ℕ)
    (hR : 42 * R < m) :
    ModAB.RationalUncapped.nextMode B (1/42) S m w R t =
      .interpolation S (t+1) := by
  have hr : R < (1/42:ℝ) * m :=
    (ModAB.RationalUncapped.forty_two_test R m).mpr hR
  unfold ModAB.RationalUncapped.nextMode
  split_ifs <;> simp_all

/-- In the zero-tolerance idealization, closing does not move an interior trial. -/
theorem closeCandidate_zero {a b x : ℝ} (hx : x ∈ Icc a b)
    (side : EndpointSide) : ModAB.Rational.closeCandidate a b 0 x side = x := by
  cases side <;> simp only [ModAB.Rational.closeCandidate] <;> split_ifs with h
  · exfalso; linarith [hx.1, h.1]
  · rfl
  · exfalso; linarith [hx.2, h.1]
  · rfl

/-- A sufficiently small error relative to both endpoints is strictly interior. -/
theorem interior_of_local_error {a b r raw q : ℝ}
    (ha : a < r) (hb : r < b) (hq : q < 1)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|) : raw ∈ Ioo a b := by
  have hm : 0 < min |a-r| |b-r| := by
    apply lt_min <;> apply abs_pos.mpr <;> linarith
  have hsmall : |raw-r| < min |a-r| |b-r| :=
    hraw.trans_lt (by simpa using mul_lt_mul_of_pos_right hq hm)
  rcases abs_lt.mp hsmall with ⟨hl, hr⟩
  have hleft := min_le_left |a-r| |b-r|
  have hright := min_le_right |a-r| |b-r|
  have habsa : |a-r| = r-a := by rw [abs_of_neg (sub_neg.mpr ha)]; ring
  have habsb : |b-r| = b-r := abs_of_pos (sub_pos.mpr hb)
  constructor <;> linarith

/-- Operational bridge with every analytic input displayed: the existing
guarded candidate is the raw rational point, zero-tolerance closing leaves it
unchanged, and the existing factor-42 policy continues interpolation. The
curvature decision and nonzero rational denominator are explicit premises. -/
theorem local_operational_bridge {f : ℝ → ℝ}
    {a b r xp q t lower upper B S w : ℝ} (clock : ℕ) (side : EndpointSide)
    (ha : a < r) (hb : r < b) (hq : q < 1)
    (hden : ModAB.Rational.denominator a b (f a) (f b) xp (f xp) ≠ 0)
    (hraw : |ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)-r| ≤
      q * min |a-r| |b-r|)
    (hsec : |ModAB.Rational.secant a b (f a) (f b)-r| ≤
      t * min |a-r| |b-r|)
    (hclamp : q+t ≤ (7/10:ℝ)*(1-t)) (hupper : 0 ≤ upper)
    (hfraw : |f (ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp))| ≤
      upper * |ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)-r|)
    (hflower : lower * min |a-r| |b-r| ≤ min |f a| |f b|)
    (hres : 42*upper*q < lower) :
    let raw := ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)
    let trial := ModAB.Rational.candidate a b (f a) (f b) xp (f xp) true
    trial = raw ∧ ModAB.Rational.closeCandidate a b 0 trial side = raw ∧
      ModAB.RationalUncapped.nextMode B (1/42) S
        (min |f a| |f b|) w |f trial| clock = .interpolation S (clock+1) := by
  dsimp only
  have hc : ModAB.Rational.candidate a b (f a) (f b) xp (f xp) true =
      ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp) := by
    rw [ModAB.Rational.candidate, ite_eq_left ⟨rfl, hden⟩]
    exact limited_eq_raw_of_local_estimates ha.le hb.le hraw hsec hclamp
  rw [hc]
  have hi := interior_of_local_error ha hb hq hraw
  have hm : 0 < min |a-r| |b-r| := by
    apply lt_min <;> apply abs_pos.mpr <;> linarith
  exact ⟨rfl, closeCandidate_zero ⟨hi.1.le, hi.2.le⟩ side,
    nextMode_interpolation_of_forty_two clock
      (forty_two_residual_of_local_bounds hupper hm hraw hfraw hflower hres)⟩

/-- Positive local coefficient means the new point crosses from the most
recent side to the opposite side, even before any asymptotic limit is used. -/
theorem positive_coefficient_crossing {en es ed ez H : ℝ}
    (hH : 0 < H) (hnd : 0 < en*ed) (hns : en*es < 0)
    (he : ez = H*en*es*ed) : 0 < ez*es ∧ ez*en < 0 := by
  have hs : es ≠ 0 := by intro hs; simp [hs] at hns
  have hh : 0 < H*(en*ed) := mul_pos hH hnd
  rw [he]
  constructor
  · have h := mul_pos hh (sq_pos_of_ne_zero hs)
    convert h using 1
    ring
  · have h := mul_neg_of_pos_of_neg hh hns
    convert h using 1
    ring

/-- Crossing from a fresh left endpoint restores two true ordinates. -/
theorem crossing_left_to_right_restores_true (s : OrdinateState) (y : ℝ)
    (hp : s.previous = some .left) (hf : ModAB.Rational.Fresh s) :
    (recordUpdate s .right y).leftWork = (recordUpdate s .right y).leftTrue ∧
    (recordUpdate s .right y).rightWork = (recordUpdate s .right y).rightTrue := by
  have h : s.leftWork = s.leftTrue := by simpa [ModAB.Rational.Fresh, hp] using hf
  simp [recordUpdate, hp, h]

/-- The symmetric crossing restores two true ordinates as well. -/
theorem crossing_right_to_left_restores_true (s : OrdinateState) (y : ℝ)
    (hp : s.previous = some .right) (hf : ModAB.Rational.Fresh s) :
    (recordUpdate s .left y).leftWork = (recordUpdate s .left y).leftTrue ∧
    (recordUpdate s .left y).rightWork = (recordUpdate s .left y).rightTrue := by
  have h : s.rightWork = s.rightTrue := by simpa [ModAB.Rational.Fresh, hp] using hf
  simp [recordUpdate, hp, h]

/-- A bounded triple-product recurrence preserves a sufficiently small initial
neighbourhood. This is an explicit recurrence interface, not a source-code
refinement or an assumed differentiability theorem. -/
theorem recurrence_local_invariant {e : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 ≤ 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hstep : ∀ n, |e (n+3)| ≤ B * |e (n+2)| * |e (n+1)| * |e n|) :
    ∀ n, |e n| ≤ ε := by
  have hq0 : 0 ≤ B*ε^2 := mul_nonneg hB (sq_nonneg ε)
  have htrip : ∀ n, |e n| ≤ ε ∧ |e (n+1)| ≤ ε ∧ |e (n+2)| ≤ ε := by
    intro n
    induction n with
    | zero => exact ⟨h0, h1, h2⟩
    | succ n ih =>
      have hc := triple_error_bound hB hε ih.2.2 ih.2.1 ih.1 (hstep n)
      have hm : min |e (n+2)| |e (n+1)| ≤ ε := (min_le_left _ _).trans ih.2.2
      have hb := mul_le_mul_of_nonneg_left hm hq0
      have hqe : (B*ε^2)*ε ≤ ε := by
        simpa using mul_le_mul_of_nonneg_right hq hε
      have he : |e (n+3)| ≤ ε := hc.trans (hb.trans hqe)
      exact ⟨ih.2.1, ih.2.2, by simpa [Nat.add_assoc] using he⟩
  exact fun n => (htrip n).1

/-- The same local recurrence entails convergence, not merely boundedness. -/
theorem recurrence_tendsto_zero {e : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hstep : ∀ n, |e (n+3)| ≤ B * |e (n+2)| * |e (n+1)| * |e n|) :
    Tendsto e atTop (𝓝 0) := by
  have hi := recurrence_local_invariant hB hε hq.le h0 h1 h2 hstep
  have hq0 : 0 ≤ B*ε^2 := mul_nonneg hB (sq_nonneg ε)
  have ht : Tendsto (fun n => |e (n+2)|) atTop (𝓝 0) := by
    apply ModAB.RationalUncapped.geometric_residual hq0 hq
    · exact fun _ => abs_nonneg _
    · intro n
      have hc := triple_error_bound hB hε (hi (n+2)) (hi (n+1)) (hi n) (hstep n)
      have hm := mul_le_mul_of_nonneg_left (min_le_left |e (n+2)| |e (n+1)|) hq0
      simpa [Nat.add_assoc] using hc.trans hm
  have habs : Tendsto (fun n => |e n|) atTop (𝓝 0) :=
    (tendsto_add_atTop_iff_nat 2).mp ht
  apply Metric.tendsto_atTop.mpr
  intro δ hδ
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp habs δ hδ
  exact ⟨N, fun n hn => by simpa [Real.dist_eq] using hN n hn⟩

/-- A bounded error coefficient is enough for the local convergence premise
needed by an order theorem; no convergence assumption on the errors is used. -/
theorem bounded_coefficient_recurrence_tendsto {e H : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hH : ∀ n, |H n| ≤ B)
    (he : ∀ n, e (n+3) = H n * e (n+2) * e (n+1) * e n) :
    Tendsto e atTop (𝓝 0) := by
  apply recurrence_tendsto_zero hB hε hq h0 h1 h2
  intro n
  rw [he n, abs_mul, abs_mul, abs_mul]
  gcongr
  exact hH n

end
end ModAB.RationalLocalSafeguards

end RModABBundle18

/- Original module: ModAB.RationalCycle -/
section RModABBundle19

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

end RModABBundle19

/- Original module: ModAB.RationalCycleEntry -/
section RModABBundle20

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

end RModABBundle20

/- Original module: ModAB.RationalDSlope -/
section RModABBundle21

/-! Finite-smoothness division at a removable zero, by explicit derivative
formulas and L'Hopital limits. No analyticity or assumed smooth factor is used. -/
namespace ModAB.Rational.DSlope
open Set Filter Topology
open scoped ContDiff
noncomputable section

theorem power_quotient_limit {N F : ℝ → ℝ} {r : ℝ} {k : ℕ}
    (hk : 0 < k) (hN : ContinuousAt N r) (hzero : N r = 0)
    (hF : ContinuousAt F r)
    (hder : ∀ᶠ x in 𝓝[≠] r, HasDerivAt N ((x-r)^(k-1)*F x) x) :
    Tendsto (fun x => N x/(x-r)^k) (𝓝[≠] r) (𝓝 (F r/k)) := by
  have hkn : (k:ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hk)
  have hx : ∀ᶠ x : ℝ in 𝓝[≠] r, x-r ≠ 0 := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    exact sub_ne_zero.mpr hx
  have hdg : ∀ᶠ x : ℝ in 𝓝[≠] r,
      HasDerivAt (fun y : ℝ => (y-r)^k) ((k:ℝ)*(x-r)^(k-1)) x := by
    exact Filter.Eventually.of_forall (fun x => by
      simpa only [mul_one, id_eq] using! ((hasDerivAt_id x).sub_const r).fun_pow k)
  apply HasDerivAt.lhopital_zero_nhdsNE hder hdg
  · filter_upwards [hx] with x hx
    simpa using mul_ne_zero hkn (pow_ne_zero _ hx)
  · simpa [hzero] using hN.tendsto.mono_left nhdsWithin_le_nhds
  · have hc : ContinuousAt (fun x : ℝ => (x-r)^k) r := by fun_prop
    simpa [Nat.ne_of_gt hk] using hc.tendsto.mono_left nhdsWithin_le_nhds
  · have hlim : Tendsto (fun x => F x/(k:ℝ)) (𝓝[≠] r) (𝓝 (F r/(k:ℝ))) :=
      (hF.tendsto.div_const (k:ℝ)).mono_left nhdsWithin_le_nhds
    apply hlim.congr'
    filter_upwards [hx] with x hx
    field_simp [hx, hkn]

def firstFormula (f : ℝ → ℝ) (r x : ℝ) : ℝ :=
  ((x-r)*deriv f x-(f x-f r))/(x-r)^2

def secondFormula (f : ℝ → ℝ) (r x : ℝ) : ℝ :=
  ((x-r)^2*deriv (deriv f) x-2*(x-r)*deriv f x+2*(f x-f r))/(x-r)^3

def first (f : ℝ → ℝ) (r : ℝ) : ℝ → ℝ :=
  Function.update (firstFormula f r) r (deriv (deriv f) r/2)

def second (f : ℝ → ℝ) (r : ℝ) : ℝ → ℝ :=
  Function.update (secondFormula f r) r (deriv (deriv (deriv f)) r/3)

theorem dslope_hasDerivAt_of_ne {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (hx : x ≠ r) :
    HasDerivAt (dslope f r) (first f r x) x := by
  have hn := (hf.hasDerivAt.sub_const (f r)).div
    ((hasDerivAt_id x).sub_const r) (sub_ne_zero.mpr hx)
  have hs : HasDerivAt (slope f r) (firstFormula f r x) x := by
    have he : slope f r = fun y => (f y-f r)/(y-r) := by
      funext y
      simp [slope_def_field]
    rw [he]
    simpa [firstFormula, mul_comm] using! hn
  simpa [first, Function.update_of_ne hx] using
    hs.congr_of_eventuallyEq (dslope_eventuallyEq_slope_of_ne f hx)

theorem dslope_hasDerivAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 f r) :
    HasDerivAt (dslope f r) (deriv (deriv f) r/2) r := by
  have h1 : ContDiffAt ℝ 1 (deriv f) r := hf.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => f x-f r-deriv f r*(x-r)) r := by
    exact (hf.continuousAt.sub continuousAt_const).sub
      (continuousAt_const.mul (continuousAt_id.sub continuousAt_const))
  have hF : ContinuousAt (dslope (deriv f) r) r :=
    continuousAt_dslope_same.mpr (h1.differentiableAt (by norm_num))
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun x => f x-f r-deriv f r*(x-r))
        ((x-r)^(2-1)*dslope (deriv f) r x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hd := ((hx.differentiableAt (by norm_num)).hasDerivAt.sub_const (f r)).sub
      (((hasDerivAt_id x).sub_const r).const_mul (deriv f r))
    convert hd using 1
    · rfl
    · simpa [smul_eq_mul] using sub_smul_dslope (deriv f) r x
  have hl := power_quotient_limit (k:=2) (by norm_num) hN (by simp) hF hder
  rw [dslope_same] at hl
  apply hasDerivAt_iff_tendsto_slope.mpr
  apply hl.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxr : x-r ≠ 0 := sub_ne_zero.mpr hx
  simp only [slope_def_field, dslope_same, dslope_of_ne f hx]
  field_simp [hxr]

theorem first_numerator_derivative {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (h1 : DifferentiableAt ℝ (deriv f) x) :
    HasDerivAt (fun y => (y-r)*deriv f y-(f y-f r))
      ((x-r)*deriv (deriv f) x) x := by
  convert (((hasDerivAt_id x).sub_const r).mul h1.hasDerivAt).sub
    (hf.hasDerivAt.sub_const (f r)) using 1
  · rfl
  · simp only [id_eq]
    ring

theorem first_continuousAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 f r) : ContinuousAt (first f r) r := by
  have h1 : ContDiffAt ℝ 1 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 0 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => (x-r)*deriv f x-(f x-f r)) r := by
    exact ((continuousAt_id.sub continuousAt_const).mul h1.continuousAt).sub
      (hf.continuousAt.sub continuousAt_const)
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun y => (y-r)*deriv f y-(f y-f r))
        ((x-r)^(2-1)*deriv (deriv f) x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 1 (deriv f) x := hx.derivWithin (by norm_num)
    simpa using first_numerator_derivative
      (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
  exact continuousAt_update_same.mpr
    (power_quotient_limit (k:=2) (by norm_num) hN (by simp) h2.continuousAt hder)

theorem first_hasDerivAt_of_ne {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 2 f x) (hx : x ≠ r) :
    HasDerivAt (first f r) (second f r x) x := by
  have h1 : ContDiffAt ℝ 1 (deriv f) x := hf.derivWithin (by norm_num)
  have hxr := sub_ne_zero.mpr hx
  have hN := first_numerator_derivative
    (r:=r) (hf.differentiableAt (by norm_num)) (h1.differentiableAt (by norm_num))
  have hQ := hN.div (((hasDerivAt_id x).sub_const r).fun_pow 2) (pow_ne_zero 2 hxr)
  have hformula : HasDerivAt (firstFormula f r) (secondFormula f r x) x := by
    convert hQ using 1
    · rfl
    · simp only [secondFormula, id_eq]
      field_simp [hxr]
      ring
  have he : first f r =ᶠ[𝓝 x] firstFormula f r := by
    filter_upwards [eventually_ne_nhds hx] with y hy
    simp [first, hy]
  simpa [second, hx] using hformula.congr_of_eventuallyEq he

theorem first_hasDerivAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) :
    HasDerivAt (first f r) (deriv (deriv (deriv f)) r/3) r := by
  have h1 : ContDiffAt ℝ 2 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have hF : ContinuousAt (dslope (deriv (deriv f)) r) r :=
    continuousAt_dslope_same.mpr (h2.differentiableAt (by norm_num))
  have hN : ContinuousAt (fun x => (x-r)*deriv f x-(f x-f r)-
      (deriv (deriv f) r/2)*(x-r)^2) r := by
    fun_prop
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun x => (x-r)*deriv f x-(f x-f r)-
        (deriv (deriv f) r/2)*(x-r)^2)
      ((x-r)^(3-1)*dslope (deriv (deriv f)) r x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 2 (deriv f) x := hx.derivWithin (by norm_num)
    have hN := first_numerator_derivative
      (r:=r) (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
    have hD := hN.sub ((((hasDerivAt_id x).sub_const r).fun_pow 2).const_mul
      (deriv (deriv f) r/2))
    convert hD using 1
    · rfl
    · have he := sub_smul_dslope (deriv (deriv f)) r x
      simp only [smul_eq_mul] at he
      have he2 := congrArg (fun z : ℝ => (x-r)*z) he
      norm_num [id_eq]
      nlinarith [he2]
  have hl := power_quotient_limit (k:=3) (by norm_num) hN (by simp) hF hder
  rw [dslope_same] at hl
  apply hasDerivAt_iff_tendsto_slope.mpr
  apply hl.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxr : x-r ≠ 0 := sub_ne_zero.mpr hx
  simp only [slope_def_field, first, Function.update_self, Function.update_of_ne hx,
    firstFormula]
  field_simp [hxr]

theorem second_numerator_derivative {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (h1 : DifferentiableAt ℝ (deriv f) x)
    (h2 : DifferentiableAt ℝ (deriv (deriv f)) x) :
    HasDerivAt (fun y => (y-r)^2*deriv (deriv f) y-
      2*(y-r)*deriv f y+2*(f y-f r))
      ((x-r)^2*deriv (deriv (deriv f)) x) x := by
  convert ((((hasDerivAt_id x).sub_const r).fun_pow 2).mul h2.hasDerivAt).sub
    ((((hasDerivAt_id x).sub_const r).const_mul 2).mul h1.hasDerivAt) |>.add
      ((hf.hasDerivAt.sub_const (f r)).const_mul 2) using 1
  · rfl
  · simp only [id_eq]
    ring

theorem second_continuousAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) : ContinuousAt (second f r) r := by
  have h1 : ContDiffAt ℝ 2 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have h3 : ContDiffAt ℝ 0 (deriv (deriv (deriv f))) r := h2.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => (x-r)^2*deriv (deriv f) x-
      2*(x-r)*deriv f x+2*(f x-f r)) r := by
    fun_prop
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun y => (y-r)^2*deriv (deriv f) y-
        2*(y-r)*deriv f y+2*(f y-f r))
        ((x-r)^(3-1)*deriv (deriv (deriv f)) x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 2 (deriv f) x := hx.derivWithin (by norm_num)
    have hddx : ContDiffAt ℝ 1 (deriv (deriv f)) x := hdx.derivWithin (by norm_num)
    simpa using second_numerator_derivative
      (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
      (hddx.differentiableAt (by norm_num))
  exact continuousAt_update_same.mpr
    (power_quotient_limit (k:=3) (by norm_num) hN (by simp) h3.continuousAt hder)

theorem dslope_hasDerivAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 2 f x) : HasDerivAt (dslope f r) (first f r x) x := by
  by_cases hx : x=r
  · subst x
    simpa [first] using dslope_hasDerivAt_same hf
  · exact dslope_hasDerivAt_of_ne (hf.differentiableAt (by norm_num)) hx

theorem first_hasDerivAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 3 f x) : HasDerivAt (first f r) (second f r x) x := by
  by_cases hx : x=r
  · subst x
    simpa [second] using first_hasDerivAt_same hf
  · exact first_hasDerivAt_of_ne (hf.of_le (by norm_num)) hx

theorem second_continuousAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 3 f x) : ContinuousAt (second f r) x := by
  by_cases hx : x=r
  · subst x
    exact second_continuousAt_same hf
  · have h1 : ContDiffAt ℝ 2 (deriv f) x := hf.derivWithin (by norm_num)
    have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) x := h1.derivWithin (by norm_num)
    have hc : ContinuousAt (secondFormula f r) x := by
      unfold secondFormula
      apply ContinuousAt.div
      · fun_prop
      · fun_prop
      · exact pow_ne_zero 3 (sub_ne_zero.mpr hx)
    apply hc.congr_of_eventuallyEq
    filter_upwards [eventually_ne_nhds hx] with y hy
    simp [second, hy]

/-- Division by a simple coordinate factor loses exactly one finite derivative.
This theorem only assumes local `C³` regularity, not analyticity. -/
theorem contDiffAt_dslope {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) : ContDiffAt ℝ 2 (dslope f r) r := by
  have hev := hf.eventually (by norm_num)
  obtain ⟨u, hu, hopen, hr⟩ := mem_nhds_iff.mp hev
  have hgood : ∀ x ∈ u, ContDiffAt ℝ 3 f x := fun x hx => hu hx
  have hder0 : ∀ x ∈ u, HasDerivAt (dslope f r) (first f r x) x :=
    fun x hx => dslope_hasDerivAt ((hgood x hx).of_le (by norm_num))
  have hder1 : ∀ x ∈ u, HasDerivAt (first f r) (second f r x) x :=
    fun x hx => first_hasDerivAt (hgood x hx)
  have hcont2 : ContinuousOn (second f r) u :=
    fun x hx => (second_continuousAt (hgood x hx)).continuousWithinAt
  have hc1 : ContDiffOn ℝ 1 (first f r) u := by
    rw [show (1:ℕ∞ω)=0+1 by rfl, contDiffOn_succ_iff_deriv_of_isOpen hopen]
    refine ⟨fun x hx => (hder1 x hx).differentiableAt.differentiableWithinAt, by norm_num, ?_⟩
    exact contDiffOn_zero.mpr (hcont2.congr (fun x hx => (hder1 x hx).deriv))
  have hc0 : ContDiffOn ℝ 2 (dslope f r) u := by
    rw [show (2:ℕ∞ω)=1+1 by norm_num, contDiffOn_succ_iff_deriv_of_isOpen hopen]
    refine ⟨fun x hx => (hder0 x hx).differentiableAt.differentiableWithinAt, by norm_num, ?_⟩
    exact hc1.congr (fun x hx => (hder0 x hx).deriv)
  exact hc0.contDiffAt (hopen.mem_nhds hr)

/-- First coefficient of the removable divided difference. -/
theorem deriv_dslope_at {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 2 f r) :
    deriv (dslope f r) r = deriv (deriv f) r/2 :=
  (dslope_hasDerivAt_same hf).deriv

/-- Second coefficient of the removable divided difference. -/
theorem second_deriv_dslope_at {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r) :
    deriv (deriv (dslope f r)) r = deriv (deriv (deriv f)) r/3 := by
  have he : deriv (dslope f r) =ᶠ[𝓝 r] first f r := by
    filter_upwards [hf.eventually (by norm_num)] with x hx
    exact (dslope_hasDerivAt (r:=r) (hx.of_le (by norm_num))).deriv
  rw [he.deriv_eq]
  exact (first_hasDerivAt_same hf).deriv

end
end ModAB.Rational.DSlope

end RModABBundle21

/- Original module: ModAB.RationalLocalAnalysis -/
section RModABBundle22

/-!
Analytic bridge for second divided differences.  The mean-value theorem is
proved by applying Rolle's theorem twice to the interpolation remainder.
Derivatives are supplied explicitly, so the statements also apply to local
C² hypotheses after the standard derivative/continuity conversions.
-/
namespace ModAB.Rational.LocalAnalysis
open Set Filter
open scoped Topology
open LocalAlgebra
noncomputable section

/-- The Newton quadratic interpolates its third node. -/
theorem quadratic_at_third {a b c ga gb gc : ℝ}
    (hab : a < b) (hbc : b < c) :
    ga + dd1 a b ga gb * (c-a) +
      dd2 a b c ga gb gc * ((c-a)*(c-b)) = gc := by
  have hab0 : b-a ≠ 0 := ne_of_gt (sub_pos.mpr hab)
  have hbc0 : c-b ≠ 0 := ne_of_gt (sub_pos.mpr hbc)
  have hac0 : c-a ≠ 0 := ne_of_gt (sub_pos.mpr (hab.trans hbc))
  unfold dd2 dd1
  field_simp
  <;> ring

/-- A second divided difference is half a second derivative at an intermediate
point.  Every analytic hypothesis is local to the interval of the three nodes. -/
theorem exists_second_derivative_eq_twice_dd2
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℝ}
    (hab : a < b) (hbc : b < c)
    (hg : ContinuousOn g (Icc a c))
    (hg₁ : ContinuousOn g₁ (Icc a c))
    (hd : ∀ x ∈ Ioo a c, HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ x ∈ Ioo a c, HasDerivAt g₁ (g₂ x) x) :
    ∃ ξ ∈ Ioo a c, g₂ ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  let s := dd1 a b (g a) (g b)
  let D := dd2 a b c (g a) (g b) (g c)
  let q : ℝ → ℝ := fun x => g a+s*(x-a)+D*((x-a)*(x-b))
  let q₁ : ℝ → ℝ := fun x => s+D*(2*x-a-b)
  have hq : ∀ x, HasDerivAt q (q₁ x) x := by
    intro x
    have hxa := (hasDerivAt_id x).sub_const a
    have hxb := (hasDerivAt_id x).sub_const b
    convert ((hasDerivAt_const x (g a)).add (hxa.const_mul s)).add
      ((hxa.mul hxb).const_mul D) using 1 <;> (try funext y) <;> (try dsimp [q, q₁]) <;> ring
  have hq₁ : ∀ x, HasDerivAt q₁ (2*D) x := by
    intro x
    convert (hasDerivAt_const x s).add
      ((((hasDerivAt_id x).const_mul 2).sub_const a).sub_const b |>.const_mul D)
      using 1 <;> (try funext y) <;> (try dsimp [q₁]) <;> ring
  have hqc : Continuous q := continuous_iff_continuousAt.mpr
    (fun x => (hq x).continuousAt)
  have hq₁c : Continuous q₁ := continuous_iff_continuousAt.mpr
    (fun x => (hq₁ x).continuousAt)
  have hqa : q a = g a := by simp [q]
  have hqb : q b = g b := by
    dsimp [q, s]
    unfold dd1
    have h : b-a ≠ 0 := ne_of_gt (sub_pos.mpr hab)
    field_simp
    <;> ring
  have hqcval : q c = g c := quadratic_at_third hab hbc
  let h : ℝ → ℝ := fun x => g x-q x
  let h₁ : ℝ → ℝ := fun x => g₁ x-q₁ x
  have hc : ContinuousOn h (Icc a c) := hg.sub hqc.continuousOn
  have h₁c : ContinuousOn h₁ (Icc a c) := hg₁.sub hq₁c.continuousOn
  have hder : ∀ x ∈ Ioo a c, HasDerivAt h (h₁ x) x :=
    fun x hx => (hd x hx).sub (hq x)
  have hder₁ : ∀ x ∈ Ioo a c, HasDerivAt h₁ (g₂ x-2*D) x :=
    fun x hx => (hd₁ x hx).sub (hq₁ x)
  have ha : h a = 0 := by simp [h, hqa]
  have hb : h b = 0 := by simp [h, hqb]
  have hcv : h c = 0 := by simp [h, hqcval]
  obtain ⟨u, hu, hzu⟩ := exists_hasDerivAt_eq_zero hab
    (hc.mono (Icc_subset_Icc_right hbc.le)) (ha.trans hb.symm)
    (fun x hx => hder x ⟨hx.1, hx.2.trans hbc⟩)
  obtain ⟨v, hv, hzv⟩ := exists_hasDerivAt_eq_zero hbc
    (hc.mono (Icc_subset_Icc_left hab.le)) (hb.trans hcv.symm)
    (fun x hx => hder x ⟨hab.trans hx.1, hx.2⟩)
  have huv : u < v := hu.2.trans hv.1
  have huva : Icc u v ⊆ Icc a c :=
    Icc_subset_Icc hu.1.le hv.2.le
  obtain ⟨ξ, hξ, hz⟩ := exists_hasDerivAt_eq_zero huv
    (h₁c.mono huva) (hzu.trans hzv.symm)
    (fun x hx => hder₁ x ⟨hu.1.trans hx.1, hx.2.trans hv.2⟩)
  exact ⟨ξ, ⟨hu.1.trans hξ.1, hξ.2.trans hv.2⟩, sub_eq_zero.mp hz⟩


/-- First divided differences converge when both distinct ordered nodes coalesce. -/
theorem tendsto_dd1_of_derivatives
    {g g₁ : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hab : ∀ n, a n < b n)
    (hg : ∀ n, ContinuousOn g (Icc (a n) (b n)))
    (hd : ∀ n x, x ∈ Ioo (a n) (b n) → HasDerivAt g (g₁ x) x)
    (hg₁ : ContinuousAt g₁ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (g₁ r)) := by
  choose ξ hξ heq using fun n =>
    exists_hasDerivAt_eq_slope g g₁ (hab n) (hg n) (hd n)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hb
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq, dd1] using hg₁.tendsto.comp hlim

/-- Second divided differences converge without a node-separation assumption. -/
theorem tendsto_dd2_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (hg : ∀ n, ContinuousOn g (Icc (a n) (c n)))
    (hg₁ : ∀ n, ContinuousOn g₁ (Icc (a n) (c n)))
    (hd : ∀ n x, x ∈ Ioo (a n) (c n) → HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ n x, x ∈ Ioo (a n) (c n) → HasDerivAt g₁ (g₂ x) x)
    (hg₂ : ContinuousAt g₂ r)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g₂ r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_eq_twice_dd2
    (hab n) (hbc n) (hg n) (hg₁ n) (hd n) (hd₁ n)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hc
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hv := (hg₂.tendsto.comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- C² specialization of the divided-difference limit. -/
theorem tendsto_dd2_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  exact tendsto_dd2_of_derivatives hab hbc
    (fun _ => hg.continuous.continuousOn)
    (fun _ => hg'.continuous.continuousOn)
    (fun _ x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (fun _ x _ => (hg'.differentiable (by norm_num) x).hasDerivAt)
    (hg'.continuous_deriv (by norm_num)).continuousAt ha hc

/-- The interpolation polynomial evaluated at the limit root tends to g(r). -/
theorem tendsto_poleValue_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) (c n) r
      (g (a n)) (g (b n)) (g (c n))) atTop (𝓝 (g r)) := by
  have hb : Tendsto b atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hc
      (fun n => (hab n).le) (fun n => (hbc n).le)
  have h₁ := tendsto_dd1_of_derivatives hab
    (fun _ => hg.continuous.continuousOn)
    (fun _ x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (hg.continuous_deriv (by norm_num)).continuousAt ha hb
  have h₂ := tendsto_dd2_of_contDiff_two hg hab hbc ha hc
  have hv := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero,
    poleValue] using hv

/-- The exact true-data error coefficient has its analytic limit. -/
theorem tendsto_error_coefficient_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g) (hgr : g r ≠ 0)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) :=
  (tendsto_dd2_of_contDiff_two hg hab hbc ha hc).div
    (tendsto_poleValue_of_contDiff_two hg hab hbc ha hc) hgr

/-- The second-derivative mean-value statement is independent of node order. -/
theorem exists_second_derivative_dd2_arbitrary_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℝ}
    (hg : ContinuousOn g (Icc (min a (min b c)) (max a (max b c))))
    (hg₁ : ContinuousOn g₁ (Icc (min a (min b c)) (max a (max b c))))
    (hd : ∀ x ∈ Ioo (min a (min b c)) (max a (max b c)), HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ x ∈ Ioo (min a (min b c)) (max a (max b c)), HasDerivAt g₁ (g₂ x) x)
    (hab : a ≠ b) (hbc : b ≠ c) (hac : a ≠ c) :
    ∃ ξ ∈ Ioo (min a (min b c)) (max a (max b c)),
      g₂ ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  let L := min a (min b c)
  let U := max a (max b c)
  have hs : ∀ u v w : ℝ, L ≤ u → w ≤ U → u < v → v < w →
      ∃ ξ ∈ Ioo u w, g₂ ξ = 2 * dd2 u v w (g u) (g v) (g w) := by
    intro u v w hLu hwU huv hvw
    exact exists_second_derivative_eq_twice_dd2 huv hvw
      (hg.mono (Icc_subset_Icc hLu hwU))
      (hg₁.mono (Icc_subset_Icc hLu hwU))
      (fun x hx => hd x ⟨hLu.trans_lt hx.1, hx.2.trans_le hwU⟩)
      (fun x hx => hd₁ x ⟨hLu.trans_lt hx.1, hx.2.trans_le hwU⟩)
  have hLa : L ≤ a := min_le_left _ _
  have hLb : L ≤ b := (min_le_right _ _).trans (min_le_left _ _)
  have hLc : L ≤ c := (min_le_right _ _).trans (min_le_right _ _)
  have haU : a ≤ U := le_max_left _ _
  have hbU : b ≤ U := (le_max_left _ _).trans (le_max_right _ _)
  have hcU : c ≤ U := (le_max_right _ _).trans (le_max_right _ _)
  have hba0 : b-a ≠ 0 := sub_ne_zero.mpr hab.symm
  have hcb0 : c-b ≠ 0 := sub_ne_zero.mpr hbc.symm
  have hca0 : c-a ≠ 0 := sub_ne_zero.mpr hac.symm
  have hab0 : a-b ≠ 0 := sub_ne_zero.mpr hab
  have hbc0 : b-c ≠ 0 := sub_ne_zero.mpr hbc
  have hac0 : a-c ≠ 0 := sub_ne_zero.mpr hac
  have hswap := dd2_swap_last (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hfirst := dd2_swap_first (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hcyc := dd2_cyclic (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hcyc₂ := dd2_cyclic (v₁ := g b) (v₂ := g c) (v₃ := g a) hcb0 hac0 hab0
  have hrev := dd2_swap_last (v₁ := g c) (v₂ := g a) (v₃ := g b) hac0 hba0 hbc0
  have lift {u w D : ℝ} (hLu : L ≤ u) (hwU : w ≤ U)
      (ht : ∃ ξ ∈ Ioo u w, g₂ ξ = 2*D)
      (hD : dd2 a b c (g a) (g b) (g c) = D) :
      ∃ ξ ∈ Ioo L U, g₂ ξ = 2*dd2 a b c (g a) (g b) (g c) := by
    obtain ⟨ξ,hξ,heq⟩ := ht
    exact ⟨ξ, ⟨hLu.trans_lt hξ.1, hξ.2.trans_le hwU⟩, by rw [hD]; exact heq⟩
  rcases lt_or_gt_of_ne hab with hab' | hba'
  · rcases lt_or_gt_of_ne hbc with hbc' | hcb'
    · exact lift hLa hcU (hs a b c hLa hcU hab' hbc') rfl
    · rcases lt_or_gt_of_ne hac with hac' | hca'
      · exact lift hLa hbU (hs a c b hLa hbU hac' hcb') hswap
      · exact lift hLc hbU (hs c a b hLc hbU hca' hab') (hcyc.trans hcyc₂)
  · rcases lt_or_gt_of_ne hac with hac' | hca'
    · exact lift hLb hcU (hs b a c hLb hcU hba' hac') hfirst
    · rcases lt_or_gt_of_ne hbc with hbc' | hcb'
      · exact lift hLb haU (hs b c a hLb haU hbc' hca') hcyc
      · exact lift hLc haU (hs c b a hLc haU hcb' hba') ((hcyc.trans hcyc₂).trans hrev)

/-- C² specialization of the unordered mean-value theorem. -/
theorem exists_second_derivative_dd2_arbitrary
    {g : ℝ → ℝ} {a b c : ℝ} (hg : ContDiff ℝ 2 g)
    (hab : a ≠ b) (hbc : b ≠ c) (hac : a ≠ c) :
    ∃ ξ ∈ Ioo (min a (min b c)) (max a (max b c)),
      deriv (deriv g) ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  exact exists_second_derivative_dd2_arbitrary_of_derivatives
    hg.continuous.continuousOn hg'.continuous.continuousOn
    (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (fun x _ => (hg'.differentiable (by norm_num) x).hasDerivAt) hab hbc hac

/-- C² gives the second divided-difference limit for every distinct-node ordering. -/
theorem tendsto_dd2_arbitrary_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_dd2_arbitrary
    hg (hab n) (hbc n) (hac n)
  have hl : Tendsto (fun n => min (a n) (min (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.min (hb.min hc)
  have hu : Tendsto (fun n => max (a n) (max (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.max (hb.max hc)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le hl hu
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  have hv := (((hg'.continuous_deriv (by norm_num)).continuousAt.tendsto).comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- Symmetry of the first divided difference. -/
theorem dd1_comm (a b ga gb : ℝ) : dd1 a b ga gb = dd1 b a gb ga := by
  unfold dd1
  have hn : ga-gb = -(gb-ga) := by ring
  have hd : a-b = -(b-a) := by ring
  rw [hn, hd, neg_div_neg_eq]

/-- First divided differences converge for either ordering of the nodes. -/
theorem tendsto_dd1_arbitrary_of_contDiff_one
    {g : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n ≠ b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (deriv g r)) := by
  have hmvt : ∀ n, ∃ ξ ∈ Ioo (min (a n) (b n)) (max (a n) (b n)),
      deriv g ξ = dd1 (a n) (b n) (g (a n)) (g (b n)) := by
    intro n
    rcases lt_or_gt_of_ne (hab n) with h | h
    · simpa only [min_eq_left h.le, max_eq_right h.le, dd1] using
        exists_hasDerivAt_eq_slope g (deriv g) h hg.continuous.continuousOn
          (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    · rw [min_eq_right h.le, max_eq_left h.le, dd1_comm (a n) (b n)]
      exact exists_hasDerivAt_eq_slope g (deriv g) h hg.continuous.continuousOn
        (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
  choose ξ hξ heq using hmvt
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le (by simpa using ha.min hb)
      (by simpa using ha.max hb) (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq] using
    (hg.continuous_deriv_one.continuousAt.tendsto).comp hlim

/-- The exact coefficient limit is independent of node ordering and relative spacing. -/
theorem tendsto_error_coefficient_arbitrary_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g) (hgr : g r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) := by
  have h₁ := tendsto_dd1_arbitrary_of_contDiff_one (hg.of_le (by norm_num)) hab ha hb
  have h₂ := tendsto_dd2_arbitrary_of_contDiff_two hg hab hbc hac ha hb hc
  have hv := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  have hp : Tendsto (fun n => poleValue (a n) (b n) (c n) r
      (g (a n)) (g (b n)) (g (c n))) atTop (𝓝 (g r)) := by
    simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero, poleValue] using hv
  exact h₂.div hp hgr

/-- Unordered coalescing-node limit from explicit local derivatives. -/
theorem tendsto_dd2_arbitrary_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ∀ n, ContinuousOn g (Icc (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n)))))
    (hg₁ : ∀ n, ContinuousOn g₁ (Icc (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n)))))
    (hd : ∀ n x, x ∈ Ioo (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n))) → HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ n x, x ∈ Ioo (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n))) → HasDerivAt g₁ (g₂ x) x)
    (hg₂ : ContinuousAt g₂ r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g₂ r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_dd2_arbitrary_of_derivatives
    (hg n) (hg₁ n) (hd n) (hd₁ n) (hab n) (hbc n) (hac n)
  have hl : Tendsto (fun n => min (a n) (min (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.min (hb.min hc)
  have hu : Tendsto (fun n => max (a n) (max (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.max (hb.max hc)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le hl hu
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hv := (hg₂.tendsto.comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- Local derivative form of the first divided-difference limit. -/
theorem tendsto_dd1_arbitrary_of_derivatives
    {g g₁ : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hg : ∀ n, ContinuousOn g (Icc (min (a n) (b n)) (max (a n) (b n))))
    (hd : ∀ n x, x ∈ Ioo (min (a n) (b n)) (max (a n) (b n)) → HasDerivAt g (g₁ x) x)
    (hg₁ : ContinuousAt g₁ r) (hab : ∀ n, a n ≠ b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (g₁ r)) := by
  have hmvt : ∀ n, ∃ ξ ∈ Ioo (min (a n) (b n)) (max (a n) (b n)),
      g₁ ξ = dd1 (a n) (b n) (g (a n)) (g (b n)) := by
    intro n
    have hgn := hg n
    have hdn := hd n
    rcases lt_or_gt_of_ne (hab n) with h | h
    · simp only [min_eq_left h.le, max_eq_right h.le] at hgn hdn ⊢
      exact exists_hasDerivAt_eq_slope g g₁ h hgn hdn
    · simp only [min_eq_right h.le, max_eq_left h.le] at hgn hdn ⊢
      rw [dd1_comm (a n) (b n)]
      exact exists_hasDerivAt_eq_slope g g₁ h hgn hdn
  choose ξ hξ heq using hmvt
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le (by simpa using ha.min hb)
      (by simpa using ha.max hb) (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq] using hg₁.tendsto.comp hlim

/-- Only C² regularity at the coalescence point is required; arbitrary finite
prefixes and arbitrary node orderings are allowed. -/
theorem tendsto_divided_differences_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n))) atTop (𝓝 (deriv g r)) ∧
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp (hg.eventually (by simp))
  let J := Ioo (r-δ) (r+δ)
  have hrJ : J ∈ 𝓝 r := Ioo_mem_nhds (by linarith) (by linarith)
  have hJ : ∀ x ∈ J, ContDiffAt ℝ 2 g x := by
    intro x hx
    apply hball
    rw [Metric.mem_ball, Real.dist_eq]
    apply abs_lt.mpr
    dsimp [J] at hx
    constructor <;> linarith [hx.1,hx.2]
  obtain ⟨N,hN⟩ := eventually_atTop.mp ((ha.eventually hrJ).and ((hb.eventually hrJ).and (hc.eventually hrJ)))
  have hm (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hreg₃ (n : ℕ) (x : ℝ)
      (hx : x ∈ Icc (min (a (n+N)) (min (b (n+N)) (c (n+N))))
        (max (a (n+N)) (max (b (n+N)) (c (n+N))))) : ContDiffAt ℝ 2 g x := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (lt_min (hm n).2.1.1 (hm n).2.2.1)).trans_le hx.1,
      hx.2.trans_lt (max_lt (hm n).1.2 (max_lt (hm n).2.1.2 (hm n).2.2.2))⟩
  have hreg₂ (n : ℕ) (x : ℝ)
      (hx : x ∈ Icc (min (a (n+N)) (b (n+N))) (max (a (n+N)) (b (n+N)))) :
      ContDiffAt ℝ 2 g x := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (hm n).2.1.1).trans_le hx.1,
      hx.2.trans_lt (max_lt (hm n).1.2 (hm n).2.1.2)⟩
  have hg' : ContDiffAt ℝ 1 (deriv g) r := hg.derivWithin (by norm_num)
  have hg'' : ContDiffAt ℝ 0 (deriv (deriv g)) r := hg'.derivWithin (by norm_num)
  constructor
  · apply (tendsto_add_atTop_iff_nat N).mp
    exact tendsto_dd1_arbitrary_of_derivatives
      (fun n x hx => (hreg₂ n x hx).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₂ n x ⟨hx.1.le,hx.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
      hg'.continuousAt (fun n => hab (n+N))
      (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
  · apply (tendsto_add_atTop_iff_nat N).mp
    exact tendsto_dd2_arbitrary_of_derivatives
      (fun n x hx => (hreg₃ n x hx).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₃ n x hx).derivWithin (m:=1) (by norm_num)).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₃ n x ⟨hx.1.le,hx.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
      (fun n x hx => (((hreg₃ n x ⟨hx.1.le,hx.2.le⟩).derivWithin (m:=1) (by norm_num)).differentiableAt (by norm_num)).hasDerivAt)
      hg''.continuousAt (fun n => hab (n+N)) (fun n => hbc (n+N)) (fun n => hac (n+N))
      (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
      (hc.comp (tendsto_add_atTop_nat N))

/-- The Newton pole value limit under local C² regularity. -/
theorem tendsto_poleValue_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g r)) := by
  obtain ⟨h₁,h₂⟩ := tendsto_divided_differences_of_contDiffAt hg hab hbc hac ha hb hc
  have hv := ((hg.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero, poleValue] using hv

/-- The exact coefficient limit under local C² regularity. -/
theorem tendsto_error_coefficient_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) :=
  (tendsto_divided_differences_of_contDiffAt hg hab hbc hac ha hb hc).2.div
    (tendsto_poleValue_of_contDiffAt hg hab hbc hac ha hb hc) hgr

end
end ModAB.Rational.LocalAnalysis

end RModABBundle22

/- Original module: ModAB.RationalRegularity -/
section RModABBundle23

/-!
Regularized reciprocal at a simple root.

The finite `C³` to `C²` division theorem is proved in `RationalDSlope`; the
wrappers below use it to establish the paper's removable-quotient regularity
and derivative coefficient under its actual finite-smoothness assumptions.
An independent analytic route and explicit `C²` factor interfaces are retained.
All derivatives and regularity statements are Mathlib's actual analytic
notions, not uninterpreted coefficient fields.
-/
namespace ModAB.Rational.Regularity
open Set Filter Topology
noncomputable section

/-- `(x-r)/f(x)`, with the reciprocal derivative at the removable singularity. -/
def regularizedReciprocal (f : ℝ → ℝ) (r : ℝ) (x : ℝ) : ℝ :=
  (dslope f r x)⁻¹

theorem regularizedReciprocal_at (f : ℝ → ℝ) (r : ℝ) :
    regularizedReciprocal f r r = 1 / deriv f r := by
  simp [regularizedReciprocal, dslope_same]

theorem regularizedReciprocal_off_root {f : ℝ → ℝ} {r x : ℝ}
    (hf : f r = 0) (hx : x ≠ r) :
    regularizedReciprocal f r x = (x-r)/f x := by
  simp [regularizedReciprocal, dslope_of_ne f hx, slope_def_field, hf, div_eq_mul_inv]

/-- Continuity at the filled-in singularity needs only differentiability. -/
theorem regularizedReciprocal_continuousAt {f : ℝ → ℝ} {r : ℝ}
    (hf : DifferentiableAt ℝ f r) (hd : deriv f r ≠ 0) :
    ContinuousAt (regularizedReciprocal f r) r := by
  exact (continuousAt_dslope_same.mpr hf).inv₀ (by simpa using hd)

/-- An explicit factor-regularity interface for finite smoothness. -/
theorem regularizedReciprocal_contDiffAt_of_dslope {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 (dslope f r) r) (hd : deriv f r ≠ 0) :
    ContDiffAt ℝ 2 (regularizedReciprocal f r) r := by
  exact hf.inv (by simpa using hd)

/-- The power series of the divided difference is obtained by deleting its
constant coefficient. No smooth-factor assumption is needed in this theorem. -/
theorem analyticAt_dslope {f : ℝ → ℝ} {r : ℝ} (hf : AnalyticAt ℝ f r) :
    AnalyticAt ℝ (dslope f r) r := by
  obtain ⟨p, hp⟩ := hf
  exact ⟨p.fslope, hp.has_fpower_series_dslope_fslope⟩

theorem regularizedReciprocal_analyticAt {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hd : deriv f r ≠ 0) :
    AnalyticAt ℝ (regularizedReciprocal f r) r := by
  exact (analyticAt_dslope hf).inv (by simpa using hd)

/-- In particular the coalescing-divided-difference theorem can use the actual
`C²` regularity of this removable quotient for every analytic simple root. -/
theorem regularizedReciprocal_contDiffAt {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hd : deriv f r ≠ 0) :
    ContDiffAt ℝ 2 (regularizedReciprocal f r) r :=
  (regularizedReciprocal_analyticAt hf hd).contDiffAt

theorem reciprocal_second_derivative {φ : ℝ → ℝ} {r : ℝ}
    (hφ : ContDiffAt ℝ 2 φ r) (hne : φ r ≠ 0) :
    iteratedDeriv 2 (fun x => (φ x)⁻¹) r =
      (2 * (deriv φ r)^2 - φ r * iteratedDeriv 2 φ r) / (φ r)^3 := by
  let g : ℝ → ℝ := fun x => (φ x)⁻¹
  have hg : ContDiffAt ℝ 2 g r := hφ.inv hne
  have hfirst : deriv g r = -deriv φ r / (φ r)^2 :=
    deriv_fun_inv'' (hφ.differentiableAt (by norm_num)) hne
  have hprod : (fun x => φ x * g x) =ᶠ[𝓝 r] (fun _ => (1:ℝ)) := by
    filter_upwards [hφ.continuousAt.eventually_ne hne] with x hx
    exact mul_inv_cancel₀ hx
  have hp := hprod.iteratedDeriv_eq 2
  have hm := iteratedDeriv_fun_mul hφ hg
  norm_num [Finset.sum_range_succ, iteratedDeriv_zero, iteratedDeriv_one,
    iteratedDeriv_const, g, hfirst] at hp hm
  rw [hm] at hp
  apply (eq_div_iff (pow_ne_zero 3 hne)).mpr
  field_simp [hne] at hp
  simp only [one_div] at hp
  nlinarith

theorem iteratedDeriv_shift_at (r : ℝ) (n : ℕ) :
    iteratedDeriv n (fun x : ℝ => x-r) r = if n=1 then 1 else 0 := by
  have hd : deriv (fun x : ℝ => x-r) = fun _ => 1 := by
    funext x
    exact ((hasDerivAt_id x).sub_const r).deriv
  cases n with
  | zero => simp
  | succ n =>
    rw [iteratedDeriv_succ', hd]
    cases n <;> simp [iteratedDeriv_const]

theorem iteratedDeriv_shift_at_point (r x : ℝ) (n : ℕ) :
    iteratedDeriv n (fun y : ℝ => y-r) x =
      if n=0 then x-r else if n=1 then 1 else 0 := by
  have hd : deriv (fun y : ℝ => y-r) = fun _ => 1 := by
    funext y
    exact ((hasDerivAt_id y).sub_const r).deriv
  cases n with
  | zero => simp
  | succ n =>
    rw [iteratedDeriv_succ', hd]
    cases n <;> simp [iteratedDeriv_const]

/-- Multiplication by the vanishing linear factor only requires continuity of
the other factor at the point, rather than its differentiability. -/
theorem hasDerivAt_vanishing_factor {ψ : ℝ → ℝ} {r : ℝ}
    (hψ : ContinuousAt ψ r) : HasDerivAt (fun x => (x-r)*ψ x) (ψ r) r := by
  apply hasDerivAt_iff_tendsto_slope.mpr
  have heq : slope (fun x => (x-r)*ψ x) r =ᶠ[𝓝[≠] r] ψ := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    simpa only [smul_eq_mul] using slope_sub_smul ψ (show r ≠ x from Ne.symm hx)
  exact (hψ.tendsto.mono_left nhdsWithin_le_nhds).congr' heq.symm

theorem factor_second_derivative {φ : ℝ → ℝ} {r x : ℝ}
    (hφ : ContDiffAt ℝ 2 φ x) :
    iteratedDeriv 2 (fun y => (y-r)*φ y) x =
      2 * deriv φ x + (x-r)*iteratedDeriv 2 φ x := by
  have hlin : ContDiffAt ℝ 2 (fun y : ℝ => y-r) x :=
    contDiffAt_id.sub contDiffAt_const
  rw [iteratedDeriv_fun_mul (n:=2) hlin hφ]
  norm_num [Finset.sum_range_succ, iteratedDeriv_shift_at_point,
    iteratedDeriv_zero, iteratedDeriv_one]
  ring

/-- At the zero of the linear factor, two continuous derivatives of `φ`
already suffice for the third derivative of `(x-r)φ(x)`. -/
theorem factor_third_derivative_at {φ : ℝ → ℝ} {r : ℝ}
    (hφ : ContDiffAt ℝ 2 φ r) :
    iteratedDeriv 3 (fun x => (x-r)*φ x) r = 3 * iteratedDeriv 2 φ r := by
  have hφ1 : ContDiffAt ℝ 1 (deriv φ) r := hφ.derivWithin (by norm_num)
  have hφ2 : ContDiffAt ℝ 0 (deriv (deriv φ)) r := hφ1.derivWithin (by norm_num)
  have hcont : ContinuousAt (iteratedDeriv 2 φ) r := by
    simpa only [show (2:ℕ)=1+1 from rfl, iteratedDeriv_succ, iteratedDeriv_zero]
      using hφ2.continuousAt
  have heq : iteratedDeriv 2 (fun x => (x-r)*φ x) =ᶠ[𝓝 r]
      (fun x => 2 * deriv φ x + (x-r)*iteratedDeriv 2 φ x) := by
    filter_upwards [hφ.eventually (by norm_num)] with x hx
    exact factor_second_derivative hx
  have hfirst := (hφ1.differentiableAt (by norm_num)).hasDerivAt.const_mul 2
  have hsum := hfirst.add (hasDerivAt_vanishing_factor hcont)
  have hsum' : HasDerivAt (fun x => 2 * deriv φ x + (x-r)*iteratedDeriv 2 φ x)
      (3 * iteratedDeriv 2 φ r) r := by
    convert hsum using 1
    simp only [show (2:ℕ)=1+1 from rfl, iteratedDeriv_succ, iteratedDeriv_zero]
    ring
  have hd := (hsum'.congr_of_eventuallyEq heq).deriv
  simpa only [show (3:ℕ)=2+1 from rfl, iteratedDeriv_succ] using hd

/-- Finite-smoothness factor interface: it does not assume analyticity. -/
theorem dslope_derivative_coefficients_of_contDiff {f : ℝ → ℝ} {r : ℝ}
    (hφ : ContDiffAt ℝ 2 (dslope f r) r) (hzero : f r = 0) :
    iteratedDeriv 2 f r = 2 * deriv (dslope f r) r ∧
    iteratedDeriv 3 f r = 3 * iteratedDeriv 2 (dslope f r) r := by
  let φ := dslope f r
  have heq : f = fun x => (x-r)*φ x := by
    funext x
    simpa only [smul_eq_mul, φ] using (sub_smul_dslope_of_zero hzero x).symm
  constructor
  · have h := factor_second_derivative (r:=r) hφ
    change iteratedDeriv 2 (fun x => (x-r)*φ x) r = _ at h
    rw [←heq] at h
    simpa using h
  · have h := factor_third_derivative_at hφ
    change iteratedDeriv 3 (fun x => (x-r)*φ x) r = _ at h
    rw [←heq] at h
    exact h

/-- The analytic factor is an actual removable divided difference; its first
two derivatives determine the second and third derivatives of `f`. -/
theorem dslope_derivative_coefficients {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) :
    iteratedDeriv 2 f r = 2 * deriv (dslope f r) r ∧
    iteratedDeriv 3 f r = 3 * iteratedDeriv 2 (dslope f r) r := by
  let φ := dslope f r
  have hφ : AnalyticAt ℝ φ r := analyticAt_dslope hf
  have hlin2 : ContDiffAt ℝ 2 (fun x : ℝ => x-r) r :=
    contDiffAt_id.sub contDiffAt_const
  have hlin3 : ContDiffAt ℝ 3 (fun x : ℝ => x-r) r :=
    contDiffAt_id.sub contDiffAt_const
  have heq : f = fun x => (x-r)*φ x := by
    funext x
    simpa only [smul_eq_mul, φ] using (sub_smul_dslope_of_zero hzero x).symm
  constructor
  · change iteratedDeriv 2 f r = 2 * deriv φ r
    calc iteratedDeriv 2 f r = iteratedDeriv 2 (fun x => (x-r)*φ x) r := by rw [heq]
         _ = 2 * deriv φ r := by
           rw [iteratedDeriv_fun_mul (n:=2) hlin2 hφ.contDiffAt]
           norm_num [Finset.sum_range_succ, iteratedDeriv_shift_at,
             iteratedDeriv_zero, iteratedDeriv_one]
  · change iteratedDeriv 3 f r = 3 * iteratedDeriv 2 φ r
    calc iteratedDeriv 3 f r = iteratedDeriv 3 (fun x => (x-r)*φ x) r := by rw [heq]
         _ = 3 * iteratedDeriv 2 φ r := by
           rw [iteratedDeriv_fun_mul (n:=3) hlin3 hφ.contDiffAt]
           norm_num [Finset.sum_range_succ, iteratedDeriv_shift_at,
             iteratedDeriv_zero, iteratedDeriv_one]

/-- Actual derivative identity from the explicitly supplied `C²` divided
factor. No assumption on formal Taylor-coefficient symbols is used. -/
theorem regularizedReciprocal_second_coefficient_of_dslope {f : ℝ → ℝ} {r : ℝ}
    (hφ : ContDiffAt ℝ 2 (dslope f r) r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    iteratedDeriv 2 (regularizedReciprocal f r) r / 2 =
      ((iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r)) / deriv f r := by
  have hn : dslope f r r ≠ 0 := by simpa using hd
  obtain ⟨h2, h3⟩ := dslope_derivative_coefficients_of_contDiff hφ hzero
  change iteratedDeriv 2 (fun x => (dslope f r x)⁻¹) r / 2 = _
  rw [reciprocal_second_derivative hφ hn, h2, h3, dslope_same]
  field_simp [hd]
  ring

/-- The preceding factor interface is fully discharged for analytic simple
roots by the proved power-series division theorem. -/
theorem regularizedReciprocal_second_coefficient {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    iteratedDeriv 2 (regularizedReciprocal f r) r / 2 =
      ((iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r)) / deriv f r :=
  regularizedReciprocal_second_coefficient_of_dslope
    (analyticAt_dslope hf).contDiffAt hzero hd

theorem regularizedReciprocal_error_constant_of_dslope {f : ℝ → ℝ} {r : ℝ}
    (hφ : ContDiffAt ℝ 2 (dslope f r) r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    (iteratedDeriv 2 (regularizedReciprocal f r) r / 2) /
        regularizedReciprocal f r r =
      (iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r) := by
  rw [regularizedReciprocal_second_coefficient_of_dslope hφ hzero hd,
    regularizedReciprocal_at]
  field_simp [hd]

/-- This is exactly the limiting coefficient `g''(r)/(2g(r))` used in the
three-point error formula. -/
theorem regularizedReciprocal_error_constant {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    (iteratedDeriv 2 (regularizedReciprocal f r) r / 2) /
        regularizedReciprocal f r r =
      (iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r) := by
  rw [regularizedReciprocal_second_coefficient hf hzero hd, regularizedReciprocal_at]
  field_simp [hd]

/-- A genuine closed neighbourhood with `C²` regularity is available for the
divided-difference mean-value argument. -/
theorem regularizedReciprocal_exists_contDiffOn_interval {f : ℝ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ContDiffOn ℝ 2 (regularizedReciprocal f r) (Icc (r-ε) (r+ε)) := by
  obtain ⟨u, hu, hcont⟩ :=
    (regularizedReciprocal_contDiffAt hf hd).contDiffOn le_rfl (by simp)
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp hu
  refine ⟨δ/2, by positivity, hcont.mono ?_⟩
  intro x hx
  apply hball
  rw [Metric.mem_ball, Real.dist_eq]
  apply abs_lt.mpr
  constructor <;> linarith [hx.1, hx.2]

/-- The paper's finite `C³` hypothesis supplies the removable quotient's `C²`
regularity. Division by the coordinate factor has been proved, not assumed. -/
theorem regularizedReciprocal_contDiffAt_of_contDiff {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hd : deriv f r ≠ 0) :
    ContDiffAt ℝ 2 (regularizedReciprocal f r) r :=
  regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd

theorem regularizedReciprocal_second_coefficient_of_contDiff {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    iteratedDeriv 2 (regularizedReciprocal f r) r / 2 =
      ((iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r)) / deriv f r :=
  regularizedReciprocal_second_coefficient_of_dslope
    (DSlope.contDiffAt_dslope hf) hzero hd

/-- The fully derived limiting error coefficient `K=c₂²-c₃`, under local
`C³` regularity and a simple root. -/
theorem regularizedReciprocal_error_constant_of_contDiff {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    (iteratedDeriv 2 (regularizedReciprocal f r) r / 2) /
        regularizedReciprocal f r r =
      (iteratedDeriv 2 f r / (2 * deriv f r))^2 -
        iteratedDeriv 3 f r / (6 * deriv f r) :=
  regularizedReciprocal_error_constant_of_dslope
    (DSlope.contDiffAt_dslope hf) hzero hd

theorem regularizedReciprocal_exists_contDiffOn_interval_of_contDiff
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ContDiffOn ℝ 2 (regularizedReciprocal f r) (Icc (r-ε) (r+ε)) := by
  obtain ⟨u, hu, hcont⟩ :=
    (regularizedReciprocal_contDiffAt_of_contDiff hf hd).contDiffOn le_rfl (by simp)
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp hu
  refine ⟨δ/2, by positivity, hcont.mono ?_⟩
  intro x hx
  apply hball
  rw [Metric.mem_ball, Real.dist_eq]
  apply abs_lt.mpr
  constructor <;> linarith [hx.1, hx.2]

end
end ModAB.Rational.Regularity

end RModABBundle23

/- Original module: ModAB.RationalAnalyticBridge -/
section RModABBundle24

/-!
From an analytic simple root to the actual three-point rational error constant.
The argument uses local regularity of the removable quotient, the proved
divided-difference mean-value theorem, and the exact rational-step identity.
No error recurrence or divided-difference limit is assumed.
-/
namespace ModAB.Rational.AnalyticBridge
open Set Filter
open scoped Topology
open LocalAlgebra LocalAnalysis Regularity
noncomputable section

/-- The derivative coefficient appearing in the local rational-step error. -/
def errorConstant (f : ℝ → ℝ) (r : ℝ) : ℝ :=
  (iteratedDeriv 2 f r / (2 * deriv f r))^2 -
    iteratedDeriv 3 f r / (6 * deriv f r)

/-- The coefficient in the exact divided-difference formula tends to K. -/
theorem tendsto_analytic_error_coefficient
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (regularizedReciprocal f r (a n))
        (regularizedReciprocal f r (b n)) (regularizedReciprocal f r (c n)) /
      poleValue (a n) (b n) (c n) r (regularizedReciprocal f r (a n))
        (regularizedReciprocal f r (b n)) (regularizedReciprocal f r (c n)))
      atTop (𝓝 (errorConstant f r)) := by
  have hg := regularizedReciprocal_contDiffAt hf hd
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hlim := tendsto_error_coefficient_of_contDiffAt hg hgr hab hbc hac ha hb hc
  have hk := regularizedReciprocal_error_constant hf hzero hd
  have hsecond : iteratedDeriv 2 (regularizedReciprocal f r) r =
      deriv (deriv (regularizedReciprocal f r)) r := by
    rw [show 2 = 1+1 from rfl, iteratedDeriv_succ', iteratedDeriv_one]
  rw [hsecond] at hk
  simpa only [hk, errorConstant] using hlim

/-- At analytic simple roots the interpolation denominator is eventually
nonzero, rather than being added as an independent local hypothesis. -/
theorem eventually_analytic_denominator_ne_zero
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    ∀ᶠ n in atTop, dd2 (a n) (b n) (c n) (1/f (a n)) (1/f (b n)) (1/f (c n)) ≠ 0 := by
  have hg := regularizedReciprocal_contDiffAt hf hd
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hp := tendsto_poleValue_of_contDiffAt hg hab hbc hac ha hb hc
  filter_upwards [hp.eventually_ne hgr] with n hn
  have h₁ := sub_ne_zero.mpr (har n)
  have h₂ := sub_ne_zero.mpr (hbr n)
  have h₃ := sub_ne_zero.mpr (hcr n)
  have hi := pole_identity_from_values (y₁ := f (a n)) (y₂ := f (b n)) (y₃ := f (c n))
    (sub_ne_zero.mpr (hab n).symm) (sub_ne_zero.mpr (hbc n).symm)
    (sub_ne_zero.mpr (hac n).symm) h₁ h₂ h₃
  rw [← regularizedReciprocal_off_root hzero (har n),
    ← regularizedReciprocal_off_root hzero (hbr n),
    ← regularizedReciprocal_off_root hzero (hcr n)] at hi
  intro h
  rw [h,zero_mul] at hi
  exact hn hi.symm

/-- Full analytic local-error theorem for the actual rational interpolant:
the error divided by the product of the three node errors tends to K.
The nodes may have any ordering and may coalesce at arbitrarily different rates. -/
theorem tendsto_normalized_rational_error
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      (baryRoot (a n) (b n) (c n) (f (a n)) (f (b n)) (f (c n))-r) /
        ((a n-r)*(b n-r)*(c n-r))) atTop (𝓝 (errorConstant f r)) := by
  have hlim := tendsto_analytic_error_coefficient hf hzero hd hab hbc hac ha hb hc
  apply hlim.congr'
  filter_upwards [eventually_analytic_denominator_ne_zero hf hzero hd hab hbc hac
    har hbr hcr ha hb hc] with n hn
  have h₁ := sub_ne_zero.mpr (har n)
  have h₂ := sub_ne_zero.mpr (hbr n)
  have h₃ := sub_ne_zero.mpr (hcr n)
  have he := baryRoot_true_error (y₁ := f (a n)) (y₂ := f (b n)) (y₃ := f (c n))
    (sub_ne_zero.mpr (hab n).symm) (sub_ne_zero.mpr (hbc n).symm)
    (sub_ne_zero.mpr (hac n).symm) h₁ h₂ h₃ hn
  rw [he, regularizedReciprocal_off_root hzero (har n),
    regularizedReciprocal_off_root hzero (hbr n),
    regularizedReciprocal_off_root hzero (hcr n)]
  field_simp [h₁,h₂,h₃]


/-- The actual implementation data and normalized denominator are eventually
nonzero near an analytic simple root. -/
theorem eventually_raw_denominator_ne_zero
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    ∀ᶠ n in atTop, f (a n) ≠ 0 ∧ f (b n) ≠ 0 ∧ f (c n) ≠ 0 ∧
      ModAB.Rational.denominator (a n) (b n) (f (a n)) (f (b n)) (c n) (f (c n)) ≠ 0 := by
  have hg := (regularizedReciprocal_contDiffAt hf hd).continuousAt
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have get_ne {u : ℕ → ℝ} (hu : Tendsto u atTop (𝓝 r)) (hur : ∀ n, u n ≠ r) :
      ∀ᶠ n in atTop, f (u n) ≠ 0 := by
    filter_upwards [(hg.tendsto.comp hu).eventually_ne hgr] with n hn
    intro hz
    apply hn
    simp only [Function.comp_def]
    rw [regularizedReciprocal_off_root hzero (hur n),hz,div_zero]
  filter_upwards [get_ne ha har, get_ne hb hbr, get_ne hc hcr,
    eventually_analytic_denominator_ne_zero hf hzero hd hab hbc hac har hbr hcr ha hb hc]
    with n hfa hfb hfc hH
  refine ⟨hfa,hfb,hfc,?_⟩
  intro hz
  apply hH
  rw [bary_denominator_eq (sub_ne_zero.mpr (hab n).symm)
    (sub_ne_zero.mpr (hbc n).symm) (sub_ne_zero.mpr (hac n).symm) hfa hfb hfc,
    hz,neg_zero,zero_mul,zero_div]

/-- The implemented normalized rational proposal has the same local error
constant K. This theorem concerns raw rational steps before clamping/closing. -/
theorem tendsto_normalized_rawPoint_error
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : AnalyticAt ℝ f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      (ModAB.Rational.rawPoint (a n) (b n) (f (a n)) (f (b n)) (c n) (f (c n))-r) /
        ((a n-r)*(b n-r)*(c n-r))) atTop (𝓝 (errorConstant f r)) := by
  have hlim := tendsto_normalized_rational_error hf hzero hd hab hbc hac har hbr hcr ha hb hc
  apply hlim.congr'
  filter_upwards [eventually_raw_denominator_ne_zero hf hzero hd hab hbc hac
    har hbr hcr ha hb hc] with n hn
  rw [rawPoint_eq_baryRoot (sub_ne_zero.mpr (hab n).symm)
    (sub_ne_zero.mpr (hbc n).symm) (sub_ne_zero.mpr (hac n).symm)
    hn.1 hn.2.1 hn.2.2.1 hn.2.2.2]

end
end ModAB.Rational.AnalyticBridge

end RModABBundle24

/- Original module: ModAB.RationalSmoothBridge -/
section RModABBundle25

/-!
The complete finite-smoothness local error theorem for the raw rational proposal.
Only local C³ regularity, a simple root, and three distinct off-root node
sequences converging to it are assumed. The removable quotient's C² regularity,
its derivative coefficient, and eventual denominator guards are all proved.
-/
namespace ModAB.Rational.SmoothBridge
open Set Filter
open scoped Topology
open LocalAlgebra LocalAnalysis Regularity AnalyticBridge
noncomputable section

/-- The coefficient in the exact divided-difference formula tends to K. -/
theorem tendsto_smooth_error_coefficient
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (regularizedReciprocal f r (a n))
        (regularizedReciprocal f r (b n)) (regularizedReciprocal f r (c n)) /
      poleValue (a n) (b n) (c n) r (regularizedReciprocal f r (a n))
        (regularizedReciprocal f r (b n)) (regularizedReciprocal f r (c n)))
      atTop (𝓝 (errorConstant f r)) := by
  have hg := regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hlim := tendsto_error_coefficient_of_contDiffAt hg hgr hab hbc hac ha hb hc
  have hk := regularizedReciprocal_error_constant_of_dslope (DSlope.contDiffAt_dslope hf) hzero hd
  have hsecond : iteratedDeriv 2 (regularizedReciprocal f r) r =
      deriv (deriv (regularizedReciprocal f r)) r := by
    rw [show 2 = 1+1 from rfl, iteratedDeriv_succ', iteratedDeriv_one]
  rw [hsecond] at hk
  simpa only [hk, errorConstant] using hlim

/-- At C³ simple roots the interpolation denominator is eventually
nonzero, rather than being added as an independent local hypothesis. -/
theorem eventually_smooth_denominator_ne_zero
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    ∀ᶠ n in atTop, dd2 (a n) (b n) (c n) (1/f (a n)) (1/f (b n)) (1/f (c n)) ≠ 0 := by
  have hg := regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hp := tendsto_poleValue_of_contDiffAt hg hab hbc hac ha hb hc
  filter_upwards [hp.eventually_ne hgr] with n hn
  have h₁ := sub_ne_zero.mpr (har n)
  have h₂ := sub_ne_zero.mpr (hbr n)
  have h₃ := sub_ne_zero.mpr (hcr n)
  have hi := pole_identity_from_values (y₁ := f (a n)) (y₂ := f (b n)) (y₃ := f (c n))
    (sub_ne_zero.mpr (hab n).symm) (sub_ne_zero.mpr (hbc n).symm)
    (sub_ne_zero.mpr (hac n).symm) h₁ h₂ h₃
  rw [← regularizedReciprocal_off_root hzero (har n),
    ← regularizedReciprocal_off_root hzero (hbr n),
    ← regularizedReciprocal_off_root hzero (hcr n)] at hi
  intro h
  rw [h,zero_mul] at hi
  exact hn hi.symm

/-- Full finite-smoothness local-error theorem for the actual rational interpolant:
the error divided by the product of the three node errors tends to K.
The nodes may have any ordering and may coalesce at arbitrarily different rates. -/
theorem tendsto_normalized_rational_error
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      (baryRoot (a n) (b n) (c n) (f (a n)) (f (b n)) (f (c n))-r) /
        ((a n-r)*(b n-r)*(c n-r))) atTop (𝓝 (errorConstant f r)) := by
  have hlim := tendsto_smooth_error_coefficient hf hzero hd hab hbc hac ha hb hc
  apply hlim.congr'
  filter_upwards [eventually_smooth_denominator_ne_zero hf hzero hd hab hbc hac
    har hbr hcr ha hb hc] with n hn
  have h₁ := sub_ne_zero.mpr (har n)
  have h₂ := sub_ne_zero.mpr (hbr n)
  have h₃ := sub_ne_zero.mpr (hcr n)
  have he := baryRoot_true_error (y₁ := f (a n)) (y₂ := f (b n)) (y₃ := f (c n))
    (sub_ne_zero.mpr (hab n).symm) (sub_ne_zero.mpr (hbc n).symm)
    (sub_ne_zero.mpr (hac n).symm) h₁ h₂ h₃ hn
  rw [he, regularizedReciprocal_off_root hzero (har n),
    regularizedReciprocal_off_root hzero (hbr n),
    regularizedReciprocal_off_root hzero (hcr n)]
  field_simp [h₁,h₂,h₃]


/-- The actual implementation data and normalized denominator are eventually
nonzero near an C³ simple root. -/
theorem eventually_raw_denominator_ne_zero
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    ∀ᶠ n in atTop, f (a n) ≠ 0 ∧ f (b n) ≠ 0 ∧ f (c n) ≠ 0 ∧
      ModAB.Rational.denominator (a n) (b n) (f (a n)) (f (b n)) (c n) (f (c n)) ≠ 0 := by
  have hg := (regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd).continuousAt
  have hgr : regularizedReciprocal f r r ≠ 0 := by
    rw [regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have get_ne {u : ℕ → ℝ} (hu : Tendsto u atTop (𝓝 r)) (hur : ∀ n, u n ≠ r) :
      ∀ᶠ n in atTop, f (u n) ≠ 0 := by
    filter_upwards [(hg.tendsto.comp hu).eventually_ne hgr] with n hn
    intro hz
    apply hn
    simp only [Function.comp_def]
    rw [regularizedReciprocal_off_root hzero (hur n),hz,div_zero]
  filter_upwards [get_ne ha har, get_ne hb hbr, get_ne hc hcr,
    eventually_smooth_denominator_ne_zero hf hzero hd hab hbc hac har hbr hcr ha hb hc]
    with n hfa hfb hfc hH
  refine ⟨hfa,hfb,hfc,?_⟩
  intro hz
  apply hH
  rw [bary_denominator_eq (sub_ne_zero.mpr (hab n).symm)
    (sub_ne_zero.mpr (hbc n).symm) (sub_ne_zero.mpr (hac n).symm) hfa hfb hfc,
    hz,neg_zero,zero_mul,zero_div]

/-- The implemented normalized rational proposal has the same local error
constant K. This theorem concerns raw rational steps before clamping/closing. -/
theorem tendsto_normalized_rawPoint_error
    {f : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (har : ∀ n, a n ≠ r) (hbr : ∀ n, b n ≠ r) (hcr : ∀ n, c n ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      (ModAB.Rational.rawPoint (a n) (b n) (f (a n)) (f (b n)) (c n) (f (c n))-r) /
        ((a n-r)*(b n-r)*(c n-r))) atTop (𝓝 (errorConstant f r)) := by
  have hlim := tendsto_normalized_rational_error hf hzero hd hab hbc hac har hbr hcr ha hb hc
  apply hlim.congr'
  filter_upwards [eventually_raw_denominator_ne_zero hf hzero hd hab hbc hac
    har hbr hcr ha hb hc] with n hn
  rw [rawPoint_eq_baryRoot (sub_ne_zero.mpr (hab n).symm)
    (sub_ne_zero.mpr (hbc n).symm) (sub_ne_zero.mpr (hac n).symm)
    hn.1 hn.2.1 hn.2.2.1 hn.2.2.2]

end
end ModAB.Rational.SmoothBridge

end RModABBundle25

/- Original module: ModAB.RationalOneSided -/
section RModABBundle26

/-!
The one-sided local regime: nonvanishing of the accumulated weight is obtained
from a finite bound on δ, not from an assumed superlinear rate.  The coefficient
limit is proved for actual divided differences with a fixed opposite endpoint.
The entry into this regime and exclusion of the exceptional zero coefficient
are separate hypotheses in the paper's local-order result.
-/
namespace ModAB.Rational.OneSided
open Filter Set
open scoped Topology
open ModAB.Rational.LocalAlgebra ModAB.Rational.LocalAnalysis
noncomputable section

/-- The accumulated weight is decreasing when each AB multiplier is in (0,1]. -/
theorem weight_antitone_of_factors {μ α : ℕ→ℝ}
    (hμ : ∀ n, 0≤μ n) (hα : ∀ n, α n≤1)
    (hstep : ∀ n, μ (n+1)=μ n*α n) : Antitone μ := by
  apply antitone_nat_of_succ_le
  intro n
  rw [hstep]
  nlinarith [hμ n,hα n]

theorem delta_monotone_of_weight_antitone {μ : ℕ→ℝ}
    (hμ : ∀ n, 0<μ n) (hanti : Antitone μ) :
    Monotone (fun n => 1/μ n-1) := by
  intro i j hij
  have h := one_div_le_one_div_of_le (hμ j) (hanti hij)
  linarith

/-- A bounded monotone δ sequence has a finite nonnegative limit, giving a
strictly positive limiting weight without any hypothesis about convergence order. -/
theorem exists_delta_weight_limit {δ : ℕ→ℝ} {M : ℝ}
    (hδ : ∀ n, 0≤δ n) (hmono : Monotone δ) (hM : ∀ n, δ n≤M) :
    ∃ L : ℝ, 0≤L ∧ L≤M ∧ Tendsto δ atTop (𝓝 L) ∧
      Tendsto (fun n => 1/(1+δ n)) atTop (𝓝 (1/(1+L))) ∧
      0<1/(1+L) ∧ 1/(1+M)≤1/(1+L) := by
  have hbdd : BddAbove (range δ) := by
    refine ⟨M,?_⟩
    rintro _ ⟨n,rfl⟩
    exact hM n
  let L : ℝ := ⨆ n, δ n
  have hlim : Tendsto δ atTop (𝓝 L) := tendsto_atTop_ciSup hmono hbdd
  have hL0 : 0≤L := ge_of_tendsto hlim (Eventually.of_forall hδ)
  have hLM : L≤M := le_of_tendsto hlim (Eventually.of_forall hM)
  have hL : 0<1+L := by linarith
  have hw : Tendsto (fun n => 1/(1+δ n)) atTop (𝓝 (1/(1+L))) :=
    tendsto_const_nhds.div (tendsto_const_nhds.add hlim) (ne_of_gt hL)
  exact ⟨L,hL0,hLM,hlim,hw,one_div_pos.mpr hL,(weight_lower_bound hL0 hLM).2⟩

theorem exists_positive_weight_limit {μ : ℕ→ℝ} {M : ℝ}
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) :
    ∃ L : ℝ, 0<L ∧ 1/(1+M)≤L ∧ Tendsto μ atTop (𝓝 L) := by
  have hδ : ∀ n, 0≤1/μ n-1 := by
    intro n
    have h : 1≤1/μ n := by
      simpa only [one_div_one] using one_div_le_one_div_of_le (hμ n) (hμone n)
    linarith
  obtain ⟨L,hL0,hLM,hlim,hw,hpos,hlower⟩ :=
    exists_delta_weight_limit hδ (delta_monotone_of_weight_antitone hμ hanti) hbound
  refine ⟨1/(1+L),hpos,hlower,?_⟩
  simpa only [weight_from_delta (ne_of_gt (hμ _))] using hw

/-- The no-crossing sign test supplies the bound needed for a positive weight
limit.  No convergence rate or summability of residual ratios is assumed. -/
theorem exists_positive_weight_limit_of_same_side
    {μ E G P g π : ℕ→ℝ} {M : ℝ}
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hE : ∀ n, 0<E n) (hP : ∀ n, 0<P n)
    (hg : ∀ n, 0<g n) (hπ : ∀ n, 0<π n)
    (hsame : ∀ n, E n*(G n+(1/μ n-1)*g n/π n)/
      (P n+E n*((1/μ n-1)*g n/π n))<0)
    (hcap : ∀ n, -π n*G n/g n≤M) :
    ∃ L : ℝ, 0<L ∧ 1/(1+M)≤L ∧ Tendsto μ atTop (𝓝 L) := by
  apply exists_positive_weight_limit hμ hμone hanti
  intro n
  have hδ : 0≤1/μ n-1 := by
    have h : 1≤1/μ n := by
      simpa only [one_div_one] using one_div_le_one_div_of_le (hμ n) (hμone n)
    linarith
  exact (delta_bound_of_same_side (hE n) (hP n) (hg n) hδ (hπ n) (hsame n)).le.trans (hcap n)

/-- The second divided difference with two coincident limiting nodes and one
fixed node.  It is defined by the derivative and a secant, not by division by zero. -/
def confluentDD2 (g : ℝ→ℝ) (r s : ℝ) : ℝ :=
  ((g s-g r)/(s-r)-deriv g r)/(s-r)

theorem tendsto_dd2_fixed_endpoint {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (confluentDD2 g r s)) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hga := hg.continuous.continuousAt.tendsto.comp ha
  have hgb := hg.continuous.continuousAt.tendsto.comp hb
  have hfirst := tendsto_dd1_arbitrary_of_contDiff_one hg hab ha hb
  have hsec : Tendsto (fun n => (g s-g (b n))/(s-b n)) atTop
      (𝓝 ((g s-g r)/(s-r))) :=
    (tendsto_const_nhds.sub hgb).div (tendsto_const_nhds.sub hb) hsr
  convert (hsec.sub hfirst).div (tendsto_const_nhds.sub ha) hsr using 1
  · funext n
    rfl
  · rfl

theorem tendsto_poleValue_fixed_endpoint {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (g r)) := by
  have hfirst := tendsto_dd1_arbitrary_of_contDiff_one hg hab ha hb
  have hsecond := tendsto_dd2_fixed_endpoint hg hab hs ha hb
  have h := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (hfirst.mul (ha.sub_const r))).add
      ((hsecond.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [poleValue,Function.comp_def,sub_self,mul_zero,sub_zero,add_zero] using h

/-- Coefficient multiplying the two fresh errors in the exact mixed-ordinate
formula when the opposite endpoint s is fixed. -/
def mixedCoefficient (g : ℝ→ℝ) (r s δ a b : ℝ) : ℝ :=
  (dd2 a b s (g a) (g b) (g s)+δ*g s/((s-a)*(s-b))) /
    (poleValue a b s r (g a) (g b) (g s)+
      δ*g s*(a-r)*(b-r)/((s-a)*(s-b)))*(s-r)

/-- The coefficient whose limit is proved below is exactly the coefficient of
the new rational error, provided g(x)=(x-r)/f(x) at the three supplied nodes. -/
theorem baryRoot_error_eq_mixedCoefficient {g : ℝ→ℝ} {a b s r ya yb ys μ : ℝ}
    (hba : b-a≠0) (hsb : s-b≠0) (hsa : s-a≠0)
    (har : a-r≠0) (hbr : b-r≠0) (hsr : s-r≠0)
    (hys : ys≠0) (hμ : μ≠0)
    (hH : dd2 a b s (1/ya) (1/yb) (1/(μ*ys))≠0)
    (hga : (a-r)/ya=g a) (hgb : (b-r)/yb=g b) (hgs : (s-r)/ys=g s) :
    baryRoot a b s ya yb (μ*ys)-r =
      mixedCoefficient g r s (1/μ-1) a b*((a-r)*(b-r)) := by
  rw [baryRoot_scaled_error hba hsb hsa har hbr hsr hys hμ hH,hga,hgb,hgs]
  unfold mixedCoefficient
  ring

theorem tendsto_mixedCoefficient {g : ℝ→ℝ} {a b δ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiff ℝ 1 g) (hgr : g r≠0) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hδ : Tendsto δ atTop (𝓝 L)) :
    Tendsto (fun n => mixedCoefficient g r s (δ n) (a n) (b n)) atTop
      (𝓝 ((confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r))) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hπ : Tendsto (fun n => (s-a n)*(s-b n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)
  have hπne : (s-r)*(s-r)≠0 := mul_ne_zero hsr hsr
  have hsecond := tendsto_dd2_fixed_endpoint hg hab hs ha hb
  have hpole := tendsto_poleValue_fixed_endpoint hg hab hs ha hb
  have hcorr := (hδ.mul_const (g s)).div hπ hπne
  have hsmall := (((hδ.mul_const (g s)).mul (ha.sub_const r)).mul
    (hb.sub_const r)).div hπ hπne
  have hden : Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s)+
      δ n*g s*(a n-r)*(b n-r)/((s-a n)*(s-b n))) atTop (𝓝 (g r)) := by
    simpa only [sub_self,mul_zero,zero_div,add_zero,Pi.div_apply] using hpole.add hsmall
  simpa only [mixedCoefficient,Pi.div_apply] using ((hsecond.add hcorr).div hden hgr).mul_const (s-r)

/-- Only C¹ regularity near the limiting root is needed for the two moving
nodes; the function may have poles elsewhere. -/
theorem tendsto_dd1_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n))) atTop (𝓝 (deriv g r)) := by
  obtain ⟨δ,hδ,hball⟩ := Metric.mem_nhds_iff.mp (hg.eventually (by simp))
  let J := Ioo (r-δ) (r+δ)
  have hrJ : J ∈ 𝓝 r := Ioo_mem_nhds (by linarith) (by linarith)
  have hJ : ∀ z ∈ J, ContDiffAt ℝ 1 g z := by
    intro z hz
    apply hball
    rw [Metric.mem_ball,Real.dist_eq]
    apply abs_lt.mpr
    dsimp [J] at hz
    constructor <;> linarith [hz.1,hz.2]
  obtain ⟨N,hN⟩ := eventually_atTop.mp ((ha.eventually hrJ).and (hb.eventually hrJ))
  have hm (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hreg (n : ℕ) (z : ℝ)
      (hz : z ∈ Icc (min (a (n+N)) (b (n+N))) (max (a (n+N)) (b (n+N)))) :
      ContDiffAt ℝ 1 g z := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (hm n).2.1).trans_le hz.1,
      hz.2.trans_lt (max_lt (hm n).1.2 (hm n).2.2)⟩
  have hg' : ContDiffAt ℝ 0 (deriv g) r := hg.derivWithin (by norm_num)
  apply (tendsto_add_atTop_iff_nat N).mp
  exact tendsto_dd1_arbitrary_of_derivatives
    (fun n z hz => (hreg n z hz).continuousAt.continuousWithinAt)
    (fun n z hz => ((hreg n z ⟨hz.1.le,hz.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
    hg'.continuousAt (fun n => hab (n+N))
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))


/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_dd2_fixed_endpoint_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (confluentDD2 g r s)) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hga := hg.continuousAt.tendsto.comp ha
  have hgb := hg.continuousAt.tendsto.comp hb
  have hfirst := tendsto_dd1_of_contDiffAt hg hab ha hb
  have hsec : Tendsto (fun n => (g s-g (b n))/(s-b n)) atTop
      (𝓝 ((g s-g r)/(s-r))) :=
    (tendsto_const_nhds.sub hgb).div (tendsto_const_nhds.sub hb) hsr
  convert (hsec.sub hfirst).div (tendsto_const_nhds.sub ha) hsr using 1
  · funext n
    rfl
  · rfl

/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_poleValue_fixed_endpoint_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (g r)) := by
  have hfirst := tendsto_dd1_of_contDiffAt hg hab ha hb
  have hsecond := tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have h := ((hg.continuousAt.tendsto.comp ha).sub
    (hfirst.mul (ha.sub_const r))).add
      ((hsecond.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [poleValue,Function.comp_def,sub_self,mul_zero,sub_zero,add_zero] using h

/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_mixedCoefficient_of_contDiffAt {g : ℝ→ℝ} {a b δ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hδ : Tendsto δ atTop (𝓝 L)) :
    Tendsto (fun n => mixedCoefficient g r s (δ n) (a n) (b n)) atTop
      (𝓝 ((confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r))) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hπ : Tendsto (fun n => (s-a n)*(s-b n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)
  have hπne : (s-r)*(s-r)≠0 := mul_ne_zero hsr hsr
  have hsecond := tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have hpole := tendsto_poleValue_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have hcorr := (hδ.mul_const (g s)).div hπ hπne
  have hsmall := (((hδ.mul_const (g s)).mul (ha.sub_const r)).mul
    (hb.sub_const r)).div hπ hπne
  have hden : Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s)+
      δ n*g s*(a n-r)*(b n-r)/((s-a n)*(s-b n))) atTop (𝓝 (g r)) := by
    simpa only [sub_self,mul_zero,zero_div,add_zero,Pi.div_apply] using hpole.add hsmall
  simpa only [mixedCoefficient,Pi.div_apply] using ((hsecond.add hcorr).div hden hgr).mul_const (s-r)

end
end ModAB.Rational.OneSided

end RModABBundle26

/- Original module: ModAB.RationalOneSidedAssembly -/
section RModABBundle27

/-!
Connected one-sided local-order theorems.  The hypotheses describe a convergent
nonvanishing rational tail with one fixed opposite endpoint.  Its recurrence
coefficient and golden-ratio Q-order are conclusions, not assumed rates.
Entry of the safeguarded program into such a tail is not claimed here.
-/
namespace ModAB.Rational.OneSidedAssembly
open Filter Set
open scoped Topology
open ModAB.Rational.LocalAlgebra ModAB.Rational.OneSided
noncomputable section

def limitingCoefficient (g : ℝ→ℝ) (r s L : ℝ) : ℝ :=
  (confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r)

/-- Actual barycentric rational updates with a fixed opposite endpoint have
golden-ratio Q-order when their derived limiting coefficient is nonzero. -/
theorem barycentric_one_sided_hasQOrder
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, μ n≠0) (hfs : f s≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hδ : Tendsto (fun n => 1/μ n-1) atTop (𝓝 L))
    (hH : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0)
    (hstep : ∀ᶠ n in atTop,
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s))
    (hC : limitingCoefficient g r s L≠0) :
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s L|^(Real.goldenRatio-1)) := by
  let C : ℕ→ℝ := fun n => mixedCoefficient g r s (1/μ n-1) (x (n+1)) (x n)
  have hxnext : Tendsto (fun n => x (n+1)) atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 1).mpr hx
  have hClim : Tendsto C atTop (𝓝 (limitingCoefficient g r s L)) :=
    tendsto_mixedCoefficient_of_contDiffAt hg hgr hadj hs hxnext hx hδ
  have hrec : ∀ᶠ n in atTop,
      x (n+2)-r=C n*(x (n+1)-r)*(x n-r) := by
    filter_upwards [hH,hstep] with n hHn hstepn
    rw [hstepn]
    have h := baryRoot_error_eq_mixedCoefficient
      (sub_ne_zero.mpr (hadj n).symm) (sub_ne_zero.mpr (hxs n).symm)
      (sub_ne_zero.mpr (hxs (n+1)).symm)
      (sub_ne_zero.mpr (hxr (n+1))) (sub_ne_zero.mpr (hxr n)) (sub_ne_zero.mpr hs)
      hfs (hμ n) hHn (hquot (n+1)) (hquot n) hquots
    simpa only [C,mul_assoc] using h
  exact ModAB.RationalGoldenOrder.eventual_two_error_recurrence_hasQOrder
    (fun n => sub_ne_zero.mpr (hxr n)) hC (by simpa using hx.sub_const r) hClim hrec

/-- Bounded decreasing positive weights supply their own nonzero limiting
weight and δ-limit.  The limiting coefficient is evaluated at sup δ. -/
theorem barycentric_one_sided_hasQOrder_of_bounded_scaling
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s M : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) (hfs : f s≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hH : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0)
    (hstep : ∀ᶠ n in atTop,
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1)≠0) :
    (∃ w : ℝ, 0<w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  have hbdd : BddAbove (range (fun n => 1/μ n-1)) := by
    refine ⟨M,?_⟩
    rintro _ ⟨n,rfl⟩
    exact hbound n
  have hδ := tendsto_atTop_ciSup (delta_monotone_of_weight_antitone hμ hanti) hbdd
  obtain ⟨w,hw,_,hwl⟩ := exists_positive_weight_limit hμ hμone hanti hbound
  exact ⟨⟨w,hw,hwl⟩,barycentric_one_sided_hasQOrder hg hgr hs hx hxr hadj hxs
    (fun n => ne_of_gt (hμ n)) hfs hquot hquots hδ hH hstep hC⟩

/-- Reorder the source's endpoint--opposite--history data into the two-fresh-
nodes--scaled-opposite order used by the local identities. -/
theorem rawPoint_barycentric_bridge {a b s ya yb ys μ : ℝ}
    (hba : b-a≠0) (hsb : s-b≠0) (hsa : s-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hys : ys≠0) (hμ : μ≠0)
    (hD : ModAB.Rational.denominator a s ya (μ*ys) b yb≠0) :
    dd2 a b s (1/ya) (1/yb) (1/(μ*ys))≠0 ∧
    ModAB.Rational.rawPoint a s ya (μ*ys) b yb=baryRoot a b s ya yb (μ*ys) := by
  have hbs : b-s≠0 := by intro h; apply hsb; linarith
  have hh := bary_denominator_ne_zero hsa hbs hba hya (mul_ne_zero hμ hys) hyb hD
  have hswap := dd2_swap_last (v₁:=1/ya) (v₂:=1/yb) (v₃:=1/(μ*ys)) hba hsb hsa
  constructor
  · rw [hswap]
    exact hh
  · rw [rawPoint_eq_baryRoot hsa hbs hba hya (mul_ne_zero hμ hys) hyb hD]
    exact (baryRoot_swap_last hba hsb hsa).symm

/-- Connected source-proposal result: the normalized rawPoint iteration,
bounded accumulated scaling, and the analytic hypotheses imply golden Q-order.
The clamp/closing/fallback must already be inactive on the stated tail. -/
theorem rawPoint_one_sided_hasQOrder_of_bounded_scaling
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s M : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) (hfs : f s≠0) (hfx : ∀ n, f (x n)≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n))≠0)
    (hstep : ∀ᶠ n in atTop, x (n+2)=
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1)≠0) :
    (∃ w : ℝ, 0<w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  have hbridge : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0 ∧
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s) := by
    filter_upwards [hD,hstep] with n hDn hstepn
    have h := rawPoint_barycentric_bridge (sub_ne_zero.mpr (hadj n).symm)
      (sub_ne_zero.mpr (hxs n).symm) (sub_ne_zero.mpr (hxs (n+1)).symm)
      (hfx (n+1)) (hfx n) hfs (ne_of_gt (hμ n)) hDn
    exact ⟨h.1,hstepn.trans h.2⟩
  exact barycentric_one_sided_hasQOrder_of_bounded_scaling hg hgr hs hx hxr hadj hxs
    hμ hμone hanti hbound hfs hquot hquots (hbridge.mono fun _ h => h.1)
    (hbridge.mono fun _ h => h.2) hC

end
end ModAB.Rational.OneSidedAssembly

end RModABBundle27

/- Original module: ModAB.RationalLocalAssembly -/
section RModABBundle28

/-!
Local convergence and the tribonacci Q-order in one theorem. The errors are
not assumed to converge: bounded coefficients and three sufficiently small
seeds imply convergence, and the coefficient limit determines the Q-order.
The displayed recurrence remains an explicit interface to interpolation.
-/
namespace ModAB.RationalLocalAssembly
open Filter Topology
noncomputable section

/-- Three small seeds, a uniformly bounded coefficient and its nonzero limit
suffice for genuine tribonacci Q-order. No error-limit premise is present. -/
theorem bounded_coefficient_recurrence_hasQOrder {e H : ℕ → ℝ} {B ε K : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B * ε ^ 2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hH : ∀ n, |H n| ≤ B) (hK : K ≠ 0)
    (hHlim : Tendsto H atTop (𝓝 K)) (hnz : ∀ n, e n ≠ 0)
    (hrec : ∀ n, e (n + 3) = H n * e (n + 2) * e (n + 1) * e n) :
    ModAB.RationalOrder.HasQOrder e ModAB.RationalOrder.tau
      (|K| ^ ((ModAB.RationalOrder.tau - 1) / 2)) := by
  have hlim := ModAB.RationalLocalSafeguards.bounded_coefficient_recurrence_tendsto
    hB hε hq h0 h1 h2 hH hrec
  exact ModAB.RationalOrder.three_error_recurrence_hasQOrder hnz hK hlim hHlim hrec

/-- The invariant two-step cone and the two coefficient limits imply actual
Q-order three for both subsequences, as well as the unequal consecutive-step
orders. Neither error convergence nor a quotient limit is assumed. -/
theorem cycle_cone_hasQOrders {θ ε C K : ℝ}
    (R : ModAB.Rational.CycleEntry.Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) (hK : K ≠ 0)
    (hu : Tendsto R.u atTop (𝓝 K)) (hv : Tendsto R.v atTop (𝓝 K)) :
    ModAB.RationalOrder.HasQOrder R.a 3 |K| ∧
    ModAB.RationalOrder.HasQOrder (fun n => R.s (n + 1)) 3 |K| ∧
    Tendsto (fun n => |R.s (n + 1)| / |R.a n| ^ (9 / 5 : ℝ)) atTop
      (𝓝 (|K| ^ (2 / 5 : ℝ))) ∧
    Tendsto (fun n => |R.a (n + 1)| / |R.s (n + 1)| ^ (5 / 3 : ℝ)) atTop
      (𝓝 (|K| ^ (1 / 3 : ℝ))) := by
  have he := R.errors_tendsto hθ0 hθ1 hC hε hsmall
  have ho := R.orders hθ0 hθ1 hC hε hsmall hK hu hv
  have hκ : 0 < |K| := abs_pos.mpr hK
  have haq := R.latest_cubic_order hθ0 hθ1 hC hε hsmall hK hu hv
  refine ⟨⟨by norm_num, hκ, he.1, ?_⟩,
    ⟨by norm_num, hκ, (tendsto_add_atTop_iff_nat 1).mpr he.2.1, ?_⟩, ho.2⟩
  · simpa only [Real.rpow_ofNat] using haq
  · simpa only [Nat.add_assoc, Nat.reduceAdd, Real.rpow_ofNat] using ho.1

end
end ModAB.RationalLocalAssembly

end RModABBundle28

/- Original module: ModAB.RationalLocalNeighborhood -/
section RModABBundle29

/-!
Uniform local bounds are consequences of C³ regularity and the simple-root
condition.  The sequential criterion below converts the already proved,
spacing-independent error limit to a single neighbourhood valid for every
admissible triple.  No orbit or coefficient bounds are hypotheses here.
-/
namespace ModAB.Rational.LocalNeighborhood
open Set Filter
open scoped Topology
noncomputable section

def tripleRadius (r a b c : ℝ) : ℝ := max |a-r| (max |b-r| |c-r|)

def Admissible (r a b c : ℝ) : Prop :=
  a ≠ b ∧ b ≠ c ∧ a ≠ c ∧ a ≠ r ∧ b ≠ r ∧ c ≠ r

def rawCoefficient (f : ℝ → ℝ) (r a b c : ℝ) : ℝ :=
  (rawPoint a b (f a) (f b) c (f c)-r) / ((a-r)*(b-r)*(c-r))

theorem nodes_tendsto_of_small_radius {a : ℕ → ℝ} {r : ℝ}
    (ha : ∀ n, |a n-r| ≤ 1 / ((n:ℝ)+1)) : Tendsto a atTop (𝓝 r) := by
  have hlim : Tendsto (fun n : ℕ => (1:ℝ)/((n:ℝ)+1)) atTop (𝓝 0) :=
    tendsto_one_div_add_atTop_nhds_zero_nat
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩ := Metric.tendsto_atTop.mp hlim ε hε
  refine ⟨N,fun n hn => ?_⟩
  rw [Real.dist_eq]
  have hp : 0 ≤ (1:ℝ)/((n:ℝ)+1) := by positivity
  exact (ha n).trans_lt (by simpa only [Real.dist_eq,sub_zero,abs_of_nonneg hp] using hN n hn)

/-- A sequential uniformity criterion on any class of three-node data. -/
theorem uniform_three_of_sequential {P Q : ℝ → ℝ → ℝ → Prop} {r : ℝ}
    (hseq : ∀ a b c : ℕ → ℝ, (∀ n, P (a n) (b n) (c n)) →
      Tendsto a atTop (𝓝 r) → Tendsto b atTop (𝓝 r) → Tendsto c atTop (𝓝 r) →
      ∀ᶠ n in atTop, Q (a n) (b n) (c n)) :
    ∃ ε > 0, ∀ a b c, P a b c → tripleRadius r a b c ≤ ε → Q a b c := by
  by_contra h
  push Not at h
  have hbad : ∀ n : ℕ, ∃ a b c, P a b c ∧
      tripleRadius r a b c ≤ 1 / ((n:ℝ)+1) ∧ ¬ Q a b c := by
    intro n
    exact h (1 / ((n:ℝ)+1)) (by positivity)
  choose a b c hP hsmall hnot using hbad
  have ha := nodes_tendsto_of_small_radius (r:=r) (a:=a) (fun n =>
    (le_max_left _ _).trans (hsmall n))
  have hb := nodes_tendsto_of_small_radius (r:=r) (a:=b) (fun n =>
    ((le_max_left _ _).trans (le_max_right _ _)).trans (hsmall n))
  have hc := nodes_tendsto_of_small_radius (r:=r) (a:=c) (fun n =>
    ((le_max_right _ _).trans (le_max_right _ _)).trans (hsmall n))
  obtain ⟨N,hN⟩ := eventually_atTop.mp (hseq a b c hP ha hb hc)
  exact hnot N (hN N le_rfl)

/-- One common radius gives the raw implementation denominator guard and
an arbitrary prescribed accuracy for its normalized error coefficient. -/
theorem exists_uniform_raw_coefficient {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    {η : ℝ} (hη : 0 < η) :
    ∃ ε > 0, ∀ a b c, Admissible r a b c → tripleRadius r a b c ≤ ε →
      f a ≠ 0 ∧ f b ≠ 0 ∧ f c ≠ 0 ∧
      denominator a b (f a) (f b) c (f c) ≠ 0 ∧
      |rawCoefficient f r a b c - AnalyticBridge.errorConstant f r| < η := by
  apply uniform_three_of_sequential
  intro a b c hP ha hb hc
  have hab := fun n => (hP n).1
  have hbc := fun n => (hP n).2.1
  have hac := fun n => (hP n).2.2.1
  have har := fun n => (hP n).2.2.2.1
  have hbr := fun n => (hP n).2.2.2.2.1
  have hcr := fun n => (hP n).2.2.2.2.2
  have hden := SmoothBridge.eventually_raw_denominator_ne_zero hf hzero hd
    hab hbc hac har hbr hcr ha hb hc
  have hlim := SmoothBridge.tendsto_normalized_rawPoint_error hf hzero hd
    hab hbc hac har hbr hcr ha hb hc
  have hclose := hlim.eventually (Metric.ball_mem_nhds _ hη)
  filter_upwards [hden,hclose] with n hn hc
  exact ⟨hn.1,hn.2.1,hn.2.2.1,hn.2.2.2,by simpa [rawCoefficient,Metric.mem_ball,Real.dist_eq] using hc⟩

/-- Positive K supplies a uniform strictly positive, bounded coefficient,
including all actual denominator guards. -/
theorem exists_positive_raw_bounds {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hK : 0 < AnalyticBridge.errorConstant f r) :
    ∃ ε > 0, ∃ B > 0, ∀ a b c, Admissible r a b c →
      tripleRadius r a b c ≤ ε →
      f a ≠ 0 ∧ f b ≠ 0 ∧ f c ≠ 0 ∧
      denominator a b (f a) (f b) c (f c) ≠ 0 ∧
      0 < rawCoefficient f r a b c ∧ |rawCoefficient f r a b c| ≤ B := by
  obtain ⟨ε,hε,hεbounds⟩ := exists_uniform_raw_coefficient hf hzero hd (η:=AnalyticBridge.errorConstant f r / 2) (by positivity)
  refine ⟨ε,hε,2*AnalyticBridge.errorConstant f r,by positivity,?_⟩
  intro a b c hp hs
  obtain ⟨ha,hb,hc,hden,he⟩ := hεbounds a b c hp hs
  rcases abs_lt.mp he with ⟨hl,hu⟩
  have hpos : 0 < rawCoefficient f r a b c := by linarith
  exact ⟨ha,hb,hc,hden,hpos,by rw [abs_of_pos hpos]; linarith⟩

/-- Exact recurrence identity for the normalized raw coefficient. -/
theorem raw_error_eq_coefficient {f : ℝ → ℝ} {r a b c : ℝ}
    (ha : a ≠ r) (hb : b ≠ r) (hc : c ≠ r) :
    rawPoint a b (f a) (f b) c (f c)-r =
      rawCoefficient f r a b c * (a-r)*(b-r)*(c-r) := by
  unfold rawCoefficient
  field_simp [sub_ne_zero.mpr ha,sub_ne_zero.mpr hb,sub_ne_zero.mpr hc]

end
end ModAB.Rational.LocalNeighborhood

end RModABBundle29

/- Original module: ModAB.RationalLocalState -/
section RModABBundle30

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

end RModABBundle30

/- Original module: ModAB.RationalAlternating -/
section RModABBundle31

/-! Inductive local trajectories of the actual guarded interpolation update.
No convergence or error recurrence is assumed. The two discarded-side records
are refreshed by the true-sign transition and the coefficient limit follows
from the finite-C³ raw-step theorem. -/
namespace ModAB.Rational.Alternating
open Set Filter
open scoped Topology
open LocalState LocalNeighborhood
noncomputable section

def orbit (f : ℝ → ℝ) (allow : State → Bool) (initial : State) : ℕ → State
  | 0 => initial
  | n+1 => step f allow (orbit f allow initial n)

@[simp] theorem orbit_zero (f : ℝ → ℝ) (allow : State → Bool) (s : State) :
    orbit f allow s 0 = s := rfl

@[simp] theorem orbit_succ (f : ℝ → ℝ) (allow : State → Bool) (s : State) (n : ℕ) :
    orbit f allow s (n+1) = step f allow (orbit f allow s n) := rfl

theorem valid_orbit {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    ∀ n, Valid f r ε (orbit f allow s n) := by
  intro n
  induction n with
  | zero => exact h
  | succ n ih => exact (step_preserves N ih 2).1

theorem continue_orbit {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (B : ℝ) (n : ℕ) :
    nextMode B f allow (orbit f allow s n) =
      .interpolation (orbit f allow s n).scale ((orbit f allow s n).clock+1) :=
  (step_preserves N (valid_orbit N h n) B).2.2.2.1

theorem opposite_shift {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    opposite (orbit f allow s (n+1)) = latest (orbit f allow s n) :=
  (step_preserves N (valid_orbit N h n) 2).2.2.2.2.1

theorem history_two_shifts {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    selectedHistory (orbit f allow s (n+2)) = latest (orbit f allow s n) := by
  have hh := (step_preserves N (valid_orbit N h (n+1)) 2).2.2.2.2.2.1
  simpa only [← orbit_succ, Nat.add_assoc, Nat.reduceAdd, opposite_shift N h] using hh

theorem latest_ne_root {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : latest s ≠ r := by
  cases hs : s.recent <;> simp only [latest,hs]
  · exact ne_of_lt h.left
  · exact ne_of_gt h.right

/-- The source-state iteration is geometrically bounded and converges. -/
theorem latest_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => latest (orbit f allow s n)) atTop (𝓝 r) := by
  have habs : Tendsto (fun n => |latest (orbit f allow s n)-r|) atTop (𝓝 0) :=
    ModAB.RationalUncapped.geometric_residual
      (mul_nonneg N.bound_nonneg (sq_nonneg ε)) N.contraction
      (fun _ => abs_nonneg _) (fun n => (step_preserves N (valid_orbit N h n) 2).2.2.2.2.2.2)
  apply Metric.tendsto_atTop.mpr
  intro δ hδ
  obtain ⟨n,hn⟩ := Metric.tendsto_atTop.mp habs δ hδ
  exact ⟨n,fun k hk => by simpa [Real.dist_eq] using hn k hk⟩

theorem opposite_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => opposite (orbit f allow s n)) atTop (𝓝 r) := by
  apply (tendsto_add_atTop_iff_nat 1).mp
  simpa only [opposite_shift N h] using latest_tendsto N h

theorem history_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => selectedHistory (orbit f allow s n)) atTop (𝓝 r) := by
  apply (tendsto_add_atTop_iff_nat 2).mp
  simpa only [history_two_shifts N h] using latest_tendsto N h

theorem endpoints_min_max {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) :
    min (latest s) (opposite s) = s.a ∧ max (latest s) (opposite s) = s.b := by
  have hab : s.a ≤ s.b := le_of_lt (h.left.trans h.right)
  cases hs : s.recent <;> simp [latest,opposite,hs,min_eq_left hab,
    min_eq_right hab,max_eq_left hab,max_eq_right hab]

theorem endpoints_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => (orbit f allow s n).a) atTop (𝓝 r) ∧
    Tendsto (fun n => (orbit f allow s n).b) atTop (𝓝 r) := by
  have hl := latest_tendsto N h
  have ho := opposite_tendsto N h
  constructor
  · have heq : (fun n => (orbit f allow s n).a) =
        fun n => min (latest (orbit f allow s n)) (opposite (orbit f allow s n)) :=
      funext (fun n => (endpoints_min_max (valid_orbit N h n)).1.symm)
    rw [heq]
    simpa only [min_self] using hl.min ho
  · have heq : (fun n => (orbit f allow s n).b) =
        fun n => max (latest (orbit f allow s n)) (opposite (orbit f allow s n)) :=
      funext (fun n => (endpoints_min_max (valid_orbit N h n)).2.symm)
    rw [heq]
    simpa only [max_self] using hl.max ho

theorem endpoint_product (s : State) (r : ℝ) :
    (s.a-r)*(s.b-r) = (latest s-r)*(opposite s-r) := by
  cases hs : s.recent <;> simp [latest,opposite,hs,mul_comm]

/-- The three-iterate recurrence is a consequence of history updates, not an
orbit hypothesis. It starts once both discarded-side records have refreshed. -/
theorem raw_recurrence {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    latest (orbit f allow s (n+3))-r =
      rawCoefficient f r (orbit f allow s (n+2)).a (orbit f allow s (n+2)).b
        (selectedHistory (orbit f allow s (n+2))) *
      (latest (orbit f allow s (n+2))-r) *
      (latest (orbit f allow s (n+1))-r) *
      (latest (orbit f allow s n)-r) := by
  have hv := valid_orbit N h (n+2)
  have hp := step_preserves N hv 2
  have hx : latest (orbit f allow s (n+3)) =
      rawPoint (orbit f allow s (n+2)).a (orbit f allow s (n+2)).b
        (f (orbit f allow s (n+2)).a) (f (orbit f allow s (n+2)).b)
        (selectedHistory (orbit f allow s (n+2)))
        (f (selectedHistory (orbit f allow s (n+2)))) := by
    rw [show n+3 = (n+2)+1 by omega, orbit_succ,step,latest_update]
    exact hp.2.2.1.trans hp.2.1
  rw [hx,raw_error_eq_coefficient hv.admissible.2.2.2.1 hv.admissible.2.2.2.2.1
    hv.admissible.2.2.2.2.2]
  have he := endpoint_product (orbit f allow s (n+2)) r
  have ho := opposite_shift N h (n+1)
  simp only [Nat.add_assoc,Nat.reduceAdd] at ho
  rw [ho] at he
  rw [history_two_shifts N h]
  linear_combination rawCoefficient f r (orbit f allow s (n+2)).a
    (orbit f allow s (n+2)).b (latest (orbit f allow s n)) *
    (latest (orbit f allow s n)-r) * he

/-- Complete trajectory-level Q-order theorem. Local estimates ensure all
operational tests, state invariance and convergence; finite C³ regularity
supplies the coefficient limit for the actual rational proposal. -/
theorem hasQOrder_of_neighborhood {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s)
    (hf : ContDiffAt ℝ 3 f r) (hz : f r=0) (hd : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r ≠ 0) :
    ModAB.RationalOrder.HasQOrder (fun n => latest (orbit f allow s n)-r)
      ModAB.RationalOrder.tau
      (|AnalyticBridge.errorConstant f r| ^ ((ModAB.RationalOrder.tau-1)/2)) := by
  let aa := fun n => (orbit f allow s (n+2)).a
  let bb := fun n => (orbit f allow s (n+2)).b
  let dd := fun n => selectedHistory (orbit f allow s (n+2))
  have had : ∀ n, Admissible r (aa n) (bb n) (dd n) :=
    fun n => (valid_orbit N h (n+2)).admissible
  have ha : Tendsto aa atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (endpoints_tendsto N h).1
  have hb : Tendsto bb atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (endpoints_tendsto N h).2
  have hh : Tendsto dd atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (history_tendsto N h)
  have hcoeff := SmoothBridge.tendsto_normalized_rawPoint_error hf hz hd
    (fun n => (had n).1) (fun n => (had n).2.1) (fun n => (had n).2.2.1)
    (fun n => (had n).2.2.2.1) (fun n => (had n).2.2.2.2.1)
    (fun n => (had n).2.2.2.2.2) ha hb hh
  apply ModAB.RationalOrder.three_error_recurrence_hasQOrder
    (fun n => sub_ne_zero.mpr (latest_ne_root (valid_orbit N h n))) hK
    (by simpa using (latest_tendsto N h).sub_const r) hcoeff
  exact raw_recurrence N h

end
end ModAB.Rational.Alternating

end RModABBundle31

/- Original module: ModAB.RationalCurvature -/
section RModABBundle32

/-! The source curvature selector uses true residuals, including both discarded
histories. These theorems derive its value from the ordinary second derivative.
No curvature decision is assumed as a separate algorithmic premise. -/
namespace ModAB.Rational.Curvature
open Set Filter Topology
open LocalAlgebra LocalAnalysis
noncomputable section

def curvatureNumerator (f : ℝ → ℝ) (a b c : ℝ) : ℝ :=
  (f c - f b) * (b-a) - (f b-f a) * (c-b)

def curvatureChangesSign (f : ℝ → ℝ) (l a b h : ℝ) : Bool :=
  if (curvatureNumerator f l a b < 0 ∧ 0 < curvatureNumerator f a b h) ∨
      (curvatureNumerator f a b h < 0 ∧ 0 < curvatureNumerator f l a b)
  then true else false

theorem curvatureNumerator_eq (f : ℝ → ℝ) {a b c : ℝ}
    (hab : a < b) (hbc : b < c) :
    curvatureNumerator f a b c =
      dd2 a b c (f a) (f b) (f c) * ((b-a)*(c-b)*(c-a)) := by
  have h₁ : b-a ≠ 0 := ne_of_gt (sub_pos.mpr hab)
  have h₂ : c-b ≠ 0 := ne_of_gt (sub_pos.mpr hbc)
  have h₃ : c-a ≠ 0 := ne_of_gt (sub_pos.mpr (hab.trans hbc))
  unfold curvatureNumerator dd2 dd1
  field_simp

theorem curvatureNumerator_sign {f : ℝ → ℝ} {a b c k : ℝ}
    (hab : a < b) (hbc : b < c)
    (hf : ∀ x ∈ Icc a c, ContDiffAt ℝ 2 f x)
    (hk : ∀ x ∈ Ioo a c, 0 < deriv (deriv f) x * k) :
    0 < curvatureNumerator f a b c * k := by
  obtain ⟨ξ,hξ,heq⟩ := exists_second_derivative_eq_twice_dd2 hab hbc
    (fun x hx => (hf x hx).continuousAt.continuousWithinAt)
    (fun x hx => ((hf x hx).derivWithin (m:=1) (by norm_num)).continuousAt.continuousWithinAt)
    (fun x hx => ((hf x ⟨hx.1.le,hx.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
    (fun x hx => (((hf x ⟨hx.1.le,hx.2.le⟩).derivWithin (m:=1)
      (by norm_num)).differentiableAt (by norm_num)).hasDerivAt)
  have hD : 0 < dd2 a b c (f a) (f b) (f c) * k := by
    have hh := hk ξ hξ
    rw [heq] at hh
    nlinarith
  rw [curvatureNumerator_eq f hab hbc]
  have hp : 0 < (b-a)*(c-b)*(c-a) :=
    mul_pos (mul_pos (sub_pos.mpr hab) (sub_pos.mpr hbc))
      (sub_pos.mpr (hab.trans hbc))
  nlinarith [mul_pos hD hp]

theorem curvatureChangesSign_eq_false {f : ℝ → ℝ} {l a b h k : ℝ}
    (hla : l < a) (hab : a < b) (hbh : b < h)
    (hf : ∀ x ∈ Icc l h, ContDiffAt ℝ 2 f x)
    (hk : ∀ x ∈ Ioo l h, 0 < deriv (deriv f) x * k) :
    curvatureChangesSign f l a b h = false := by
  have hl := curvatureNumerator_sign hla hab
    (fun x hx => hf x ⟨hx.1,hx.2.trans hbh.le⟩)
    (fun x hx => hk x ⟨hx.1,hx.2.trans hbh⟩)
  have hr := curvatureNumerator_sign hab hbh
    (fun x hx => hf x ⟨hla.le.trans hx.1,hx.2⟩)
    (fun x hx => hk x ⟨hla.trans hx.1,hx.2⟩)
  unfold curvatureChangesSign
  have hn : ¬ ((curvatureNumerator f l a b < 0 ∧ 0 < curvatureNumerator f a b h) ∨
      (curvatureNumerator f a b h < 0 ∧ 0 < curvatureNumerator f l a b)) := by
    intro hc
    rcases hc with hc|hc
    · have hkneg : k < 0 := by
        by_contra hh
        have hkn : 0 ≤ k := le_of_not_gt hh
        nlinarith [mul_nonpos_of_nonpos_of_nonneg hc.1.le hkn]
      nlinarith [mul_nonpos_of_nonneg_of_nonpos hc.2.le hkneg.le]
    · have hkneg : k < 0 := by
        by_contra hh
        have hkn : 0 ≤ k := le_of_not_gt hh
        nlinarith [mul_nonpos_of_nonpos_of_nonneg hc.1.le hkn]
      nlinarith [mul_nonpos_of_nonneg_of_nonpos hc.2.le hkneg.le]
  simp [hn]

theorem exists_radius_no_curvature_change {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 f r) (h₂ : deriv (deriv f) r ≠ 0) :
    ∃ ε > 0, ∀ l a b h : ℝ,
      r-ε < l → l < a → a < b → b < h → h < r+ε →
      curvatureChangesSign f l a b h = false := by
  have hf' : ContDiffAt ℝ 1 (deriv f) r := hf.derivWithin (by norm_num)
  have hf'' : ContDiffAt ℝ 0 (deriv (deriv f)) r := hf'.derivWithin (by norm_num)
  have hc : ContinuousAt (fun x => deriv (deriv f) x * deriv (deriv f) r) r :=
    hf''.continuousAt.mul_const _
  have hs : ∀ᶠ x in 𝓝 r, 0 < deriv (deriv f) x * deriv (deriv f) r :=
    hc.eventually (Ioi_mem_nhds (by simpa only [sq] using sq_pos_of_ne_zero h₂))
  obtain ⟨ε,hε,hball⟩ := Metric.mem_nhds_iff.mp ((hf.eventually (by simp)).and hs)
  refine ⟨ε,hε,?_⟩
  intro l a b h hl hla hab hbh hh
  have hlocal (x : ℝ) (hx : x ∈ Icc l h) :
      ContDiffAt ℝ 2 f x ∧ 0 < deriv (deriv f) x * deriv (deriv f) r := by
    apply hball
    rw [Metric.mem_ball, Real.dist_eq]
    apply abs_lt.mpr
    constructor <;> linarith [hx.1,hx.2]
  exact curvatureChangesSign_eq_false hla hab hbh
    (fun x hx => (hlocal x hx).1)
    (fun x hx => (hlocal x ⟨hx.1.le,hx.2.le⟩).2)

end
end ModAB.Rational.Curvature

end RModABBundle32

/- Original module: ModAB.RationalCycleUniform -/
section RModABBundle33

/-!
Uniform local coefficients for the negative-coefficient two-step cycle.
The coefficients below are the exact divided-difference expressions. Their
limits are proved from local C² regularity of the removable reciprocal, and
uniform small-cone bounds follow by a sequential contradiction. Thus the
sign and magnitude bounds used by `CycleEntry` are derived, not assumed.
-/
namespace ModAB.Rational.CycleUniform
open Filter Set
open scoped Topology
open LocalAlgebra LocalAnalysis
noncomputable section

def localK (g : ℝ → ℝ) (r : ℝ) : ℝ := (deriv (deriv g) r / 2) / g r

def firstPole (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  poleValue (r+a) (r+s) (r+d) r (g (r+a)) (g (r+s)) (g (r+d))

def firstCoefficient (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  dd2 (r+a) (r+s) (r+d) (g (r+a)) (g (r+s)) (g (r+d)) /
    firstPole g r a s d

def firstError (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  firstCoefficient g r a s d * a * s * d

def ordinateRatio (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  firstCoefficient g r a s d * s * d * (g (r+a) / g (r+firstError g r a s d))

def normalizedGap (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  (1-firstError g r a s d / s) * (1-a/s)

def secondPole (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  let b := firstError g r a s d
  let t := ordinateRatio g r a s d
  poleValue (r+b) (r+a) (r+s) r (g (r+b)) (g (r+a)) (g (r+s)) +
    (t/(1-t)) * g (r+s) * (a/s) * (b/s) / normalizedGap g r a s d

def secondCoefficient (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  let b := firstError g r a s d
  let u := firstCoefficient g r a s d
  let t := ordinateRatio g r a s d
  ((s/d) * dd2 (r+b) (r+a) (r+s) (g (r+b)) (g (r+a)) (g (r+s)) +
    u * g (r+s) * (g (r+a)/g (r+b)) /
      (normalizedGap g r a s d * (1-t))) / secondPole g r a s d

theorem tendsto_firstCoefficient
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstCoefficient g r (a n) (s n) (d n)) atTop (𝓝 (localK g r)) := by
  have h := tendsto_error_coefficient_of_contDiffAt (a:=fun n => r+a n)
    (b:=fun n => r+s n) (c:=fun n => r+d n) hg hgr
    (fun n => by simpa using has n) (fun n => by simpa using hsd n)
    (fun n => by simpa using had n)
    (by simpa using ha.const_add r) (by simpa using hs.const_add r)
    (by simpa using hd.const_add r)
  exact h

theorem tendsto_firstPole
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) := by
  exact tendsto_poleValue_of_contDiffAt hg
    (fun n => by simpa using has n) (fun n => by simpa using hsd n)
    (fun n => by simpa using had n)
    (by simpa using ha.const_add r) (by simpa using hs.const_add r)
    (by simpa using hd.const_add r)

/-- Arbitrarily many coincidences in a finite initial prefix do not affect
the local divided-difference limits. -/
theorem tendsto_divided_differences_eventually
    {g : ℝ → ℝ} {r : ℝ} {a b c : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hdistinct : ∀ᶠ n in atTop, a n ≠ b n ∧ b n ≠ c n ∧ a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) ∧
    Tendsto (fun n => poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g r)) := by
  obtain ⟨N,hN⟩ := eventually_atTop.mp hdistinct
  have hn (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hdd := (tendsto_divided_differences_of_contDiffAt hg
    (fun n => (hn n).1) (fun n => (hn n).2.1) (fun n => (hn n).2.2)
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
    (hc.comp (tendsto_add_atTop_nat N))).2
  have hp := tendsto_poleValue_of_contDiffAt hg
    (fun n => (hn n).1) (fun n => (hn n).2.1) (fun n => (hn n).2.2)
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
    (hc.comp (tendsto_add_atTop_nat N))
  exact ⟨(tendsto_add_atTop_iff_nat N).mp hdd, (tendsto_add_atTop_iff_nat N).mp hp⟩

theorem eventually_ne_of_ratio_zero {a b : ℕ → ℝ}
    (hb : ∀ n, b n ≠ 0) (hlim : Tendsto (fun n => a n / b n) atTop (𝓝 0)) :
    ∀ᶠ n in atTop, a n ≠ b n := by
  filter_upwards [hlim.eventually_ne (by norm_num : (0:ℝ) ≠ 1)] with n hn
  intro h
  apply hn
  simp [h,hb n]

/-- The true first step creates a smaller fresh node and nonzero ordinates
without adding these properties as hypotheses. -/
theorem first_step_limits
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha0 : ∀ n, a n ≠ 0) (hs0 : ∀ n, s n ≠ 0) (hd0 : ∀ n, d n ≠ 0)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstError g r (a n) (s n) (d n)) atTop (𝓝 0) ∧
    Tendsto (fun n => firstError g r (a n) (s n) (d n) / a n) atTop (𝓝 0) ∧
    Tendsto (fun n => firstError g r (a n) (s n) (d n) / s n) atTop (𝓝 0) ∧
    (∀ᶠ n in atTop, firstError g r (a n) (s n) (d n) ≠ 0 ∧
      firstError g r (a n) (s n) (d n) ≠ a n ∧
      firstError g r (a n) (s n) (d n) ≠ s n) := by
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  have hb : Tendsto (fun n => firstError g r (a n) (s n) (d n)) atTop (𝓝 0) := by
    simpa only [firstError,mul_zero] using ((hu.mul ha).mul hs).mul hd
  have hba : Tendsto (fun n => firstError g r (a n) (s n) (d n)/a n) atTop (𝓝 0) := by
    have he (n : ℕ) : firstError g r (a n) (s n) (d n)/a n =
        firstCoefficient g r (a n) (s n) (d n)*s n*d n := by
      unfold firstError
      field_simp [ha0 n]
    simp_rw [he]
    simpa only [mul_zero] using (hu.mul hs).mul hd
  have hbs : Tendsto (fun n => firstError g r (a n) (s n) (d n)/s n) atTop (𝓝 0) := by
    have he (n : ℕ) : firstError g r (a n) (s n) (d n)/s n =
        firstCoefficient g r (a n) (s n) (d n)*a n*d n := by
      unfold firstError
      field_simp [hs0 n]
    simp_rw [he]
    simpa only [mul_zero] using (hu.mul ha).mul hd
  refine ⟨hb,hba,hbs,?_⟩
  filter_upwards [hu.eventually_ne hK, eventually_ne_of_ratio_zero ha0 hba,
    eventually_ne_of_ratio_zero hs0 hbs] with n hn hna hns
  exact ⟨mul_ne_zero (mul_ne_zero (mul_ne_zero hn (ha0 n)) (hs0 n)) (hd0 n),hna,hns⟩

/-- The scaled second coefficient has limit K when the two separation ratios
and the radius tend to zero. There is no assumed second error recurrence. -/
theorem second_step_limits
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha0 : ∀ n, a n ≠ 0) (hs0 : ∀ n, s n ≠ 0) (hd0 : ∀ n, d n ≠ 0)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0))
    (haslim : Tendsto (fun n => a n/s n) atTop (𝓝 0))
    (hsdlim : Tendsto (fun n => s n/d n) atTop (𝓝 0)) :
    Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 0) ∧
    Tendsto (fun n => normalizedGap g r (a n) (s n) (d n)) atTop (𝓝 1) ∧
    Tendsto (fun n => secondPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) ∧
    Tendsto (fun n => secondCoefficient g r (a n) (s n) (d n)) atTop (𝓝 (localK g r)) := by
  let b := fun n => firstError g r (a n) (s n) (d n)
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  obtain ⟨hb,hba,hbs,hbn⟩ := first_step_limits hg hgr hK has hsd had ha0 hs0 hd0 ha hs hd
  have ha' : Tendsto (fun n => r+a n) atTop (𝓝 r) := by simpa using ha.const_add r
  have hs' : Tendsto (fun n => r+s n) atTop (𝓝 r) := by simpa using hs.const_add r
  have hb' : Tendsto (fun n => r+b n) atTop (𝓝 r) := by simpa using hb.const_add r
  have hga := hg.continuousAt.tendsto.comp ha'
  have hgs := hg.continuousAt.tendsto.comp hs'
  have hgb := hg.continuousAt.tendsto.comp hb'
  have hJ : Tendsto (fun n => g (r+a n)/g (r+b n)) atTop (𝓝 1) := by
    have hh := hga.div hgb hgr
    rw [div_self hgr] at hh
    convert hh using 1
    ext n
    rfl
  have ht : Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 0) := by
    simpa only [ordinateRatio,Function.comp_def,b,mul_zero,zero_mul] using ((hu.mul hs).mul hd).mul hJ
  have hq : Tendsto (fun n => normalizedGap g r (a n) (s n) (d n)) atTop (𝓝 1) := by
    simpa only [normalizedGap,sub_zero,mul_one] using
      ((tendsto_const_nhds (x:=(1:ℝ))).sub hbs).mul ((tendsto_const_nhds (x:=(1:ℝ))).sub haslim)
  have hdist : ∀ᶠ n in atTop, r+b n ≠ r+a n ∧ r+a n ≠ r+s n ∧ r+b n ≠ r+s n := by
    filter_upwards [hbn] with n hn
    exact ⟨by simpa using hn.2.1,by simpa using has n,by simpa using hn.2.2⟩
  obtain ⟨hG,hP⟩ := tendsto_divided_differences_eventually hg hdist hb' ha' hs'
  have hμ : Tendsto (fun n => 1-ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 1) := by
    simpa using tendsto_const_nhds.sub ht
  have hδ : Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n) /
      (1-ordinateRatio g r (a n) (s n) (d n))) atTop (𝓝 0) := by
    have hh := ht.div hμ (by norm_num : (1:ℝ) ≠ 0)
    rw [zero_div] at hh
    convert hh using 1
  have hcorr := (((hδ.mul hgs).mul haslim).mul hbs).div hq (by norm_num : (1:ℝ) ≠ 0)
  have hP2 : Tendsto (fun n => secondPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) := by
    simpa only [secondPole,b,Function.comp_def,Pi.div_apply,mul_zero,zero_mul,zero_div,add_zero] using hP.add hcorr
  have hdom := ((hu.mul hgs).mul hJ).div (hq.mul hμ) (by norm_num : (1:ℝ)*1 ≠ 0)
  have hV := ((hsdlim.mul hG).add hdom).div hP2 hgr
  refine ⟨ht,hq,hP2,?_⟩
  simp only [mul_one,zero_mul,zero_add,div_one,mul_div_cancel_right₀ _ hgr] at hV
  convert hV using 1
  ext n
  rfl


/-- A concrete positive null sequence, used only to extract uniform cones. -/
def smallRadius (n : ℕ) : ℝ := 1 / ((n : ℝ)+2)

theorem smallRadius_pos (n : ℕ) : 0 < smallRadius n := by
  unfold smallRadius
  positivity

theorem smallRadius_lt_one (n : ℕ) : smallRadius n < 1 := by
  unfold smallRadius
  apply (div_lt_one (by positivity : (0:ℝ) < (n:ℝ)+2)).mpr
  nlinarith [Nat.cast_nonneg (α:=ℝ) n]

theorem smallRadius_tendsto : Tendsto smallRadius atTop (𝓝 0) := by
  have h : Tendsto (fun n : ℕ => (n:ℝ)+2) atTop atTop :=
    tendsto_natCast_atTop_atTop.atTop_add tendsto_const_nhds
  change Tendsto (fun n : ℕ => 1/((n:ℝ)+2)) atTop (𝓝 0)
  simpa only [one_div,Function.comp_def] using tendsto_inv_atTop_zero.comp h

/-- Sequential extraction of a single uniform neighbourhood; the witness
triples are selected from the failure of uniformity, not supplied externally. -/
theorem exists_uniform_cone_of_sequential {P : ℝ → ℝ → ℝ → Prop}
    (hseq : ∀ a s d : ℕ → ℝ,
      (∀ n, CycleEntry.Cone (smallRadius n) (smallRadius n) (a n) (s n) (d n)) →
      ∀ᶠ n in atTop, P (a n) (s n) (d n)) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d → P a s d := by
  by_contra h
  push Not at h
  choose a s d hc hp using fun n => h (smallRadius n) (smallRadius n)
    (smallRadius_pos n) (smallRadius_lt_one n) (smallRadius_pos n)
  obtain ⟨n,hn⟩ := (hseq a s d hc).exists
  exact hp n hn

theorem Cone.distinct {θ ε a s d : ℝ} (h : CycleEntry.Cone θ ε a s d)
    (hθ : θ < 1) : a ≠ s ∧ s ≠ d ∧ a ≠ d := by
  have has : a ≠ s := by
    intro he
    have hh := h.opposite
    rw [he] at hh
    nlinarith [sq_nonneg s]
  have hsd : s ≠ d := by
    intro he
    have hh := h.opposite
    rw [he] at hh
    linarith [h.same]
  have had : a ≠ d := by
    intro he
    have horder := h.ordered_errors hθ.le
    have hsmall := h.opposite_small
    rw [he] at horder
    have hp := abs_pos.mpr h.history_ne_zero
    nlinarith
  exact ⟨has,hsd,had⟩

/-- A shrinking cone forces both independent separation ratios to zero. -/
theorem shrinking_cone_limits {a s d : ℕ → ℝ}
    (hc : ∀ n, CycleEntry.Cone (smallRadius n) (smallRadius n) (a n) (s n) (d n)) :
    Tendsto a atTop (𝓝 0) ∧ Tendsto s atTop (𝓝 0) ∧ Tendsto d atTop (𝓝 0) ∧
    Tendsto (fun n => a n/s n) atTop (𝓝 0) ∧
    Tendsto (fun n => s n/d n) atTop (𝓝 0) := by
  have hord (n : ℕ) := (hc n).ordered_errors (smallRadius_lt_one n).le
  have hdabs := squeeze_zero (fun n => abs_nonneg (d n))
    (fun n => (hc n).history_bound) smallRadius_tendsto
  have hsabs := squeeze_zero (fun n => abs_nonneg (s n)) (fun n => (hord n).2) hdabs
  have haabs := squeeze_zero (fun n => abs_nonneg (a n)) (fun n => (hord n).1) hsabs
  have hratio₁ : ∀ n, |a n/s n| ≤ smallRadius n := by
    intro n
    rw [abs_div]
    exact (div_le_iff₀ (abs_pos.mpr (hc n).opposite_ne_zero)).mpr (hc n).latest_small
  have hratio₂ : ∀ n, |s n/d n| ≤ smallRadius n := by
    intro n
    rw [abs_div]
    exact (div_le_iff₀ (abs_pos.mpr (hc n).history_ne_zero)).mpr (hc n).opposite_small
  exact ⟨(tendsto_zero_iff_abs_tendsto_zero _).mpr haabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr hsabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr hdabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr
      (squeeze_zero (fun n => abs_nonneg _) hratio₁ smallRadius_tendsto),
    (tendsto_zero_iff_abs_tendsto_zero _).mpr
      (squeeze_zero (fun n => abs_nonneg _) hratio₂ smallRadius_tendsto)⟩

/-- The uniform signed coefficient bounds needed for the invariant-cycle
proof follow from local regularity. The same neighbourhood excludes every
singular denominator in the normalized pair formulas. -/
theorem exists_uniform_negative_coefficients
    {g : ℝ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r < 0) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        firstCoefficient g r a s d ∈ Ioo (2*localK g r) (localK g r/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*localK g r) (localK g r/2) ∧
        firstPole g r a s d ≠ 0 ∧ secondPole g r a s d ≠ 0 ∧
        g (r+a) ≠ 0 ∧ g (r+s) ≠ 0 ∧ g (r+d) ≠ 0 ∧
        g (r+firstError g r a s d) ≠ 0 ∧
        firstError g r a s d ≠ 0 ∧ firstError g r a s d ≠ a ∧
        firstError g r a s d ≠ s ∧
        0 < normalizedGap g r a s d ∧ ordinateRatio g r a s d < 1 ∧
        0 < g (r+a) / g (r+firstError g r a s d) := by
  apply exists_uniform_cone_of_sequential
  intro a s d hc
  obtain ⟨ha,hs,hd,haslim,hsdlim⟩ := shrinking_cone_limits hc
  have hdist (n : ℕ) := Cone.distinct (hc n) (smallRadius_lt_one n)
  have has := fun n => (hdist n).1
  have hsd := fun n => (hdist n).2.1
  have had := fun n => (hdist n).2.2
  have ha0 := fun n => (hc n).latest_ne_zero
  have hs0 := fun n => (hc n).opposite_ne_zero
  have hd0 := fun n => (hc n).history_ne_zero
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  have hp := tendsto_firstPole hg has hsd had ha hs hd
  obtain ⟨hb,hba,hbs,hbn⟩ := first_step_limits hg hgr hK.ne has hsd had ha0 hs0 hd0 ha hs hd
  obtain ⟨ht,hq,hp₂,hv⟩ := second_step_limits hg hgr hK.ne has hsd had
    ha0 hs0 hd0 ha hs hd haslim hsdlim
  have hga : Tendsto (fun n => g (r+a n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using ha.const_add r)
  have hgs : Tendsto (fun n => g (r+s n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hs.const_add r)
  have hgd : Tendsto (fun n => g (r+d n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hd.const_add r)
  have hgb : Tendsto (fun n => g (r+firstError g r (a n) (s n) (d n))) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hb.const_add r)
  have hJ : Tendsto (fun n => g (r+a n)/g (r+firstError g r (a n) (s n) (d n))) atTop (𝓝 1) := by
    have hh := hga.div hgb hgr
    rw [div_self hgr] at hh
    convert hh using 1
  have hI : Ioo (2*localK g r) (localK g r/2) ∈ 𝓝 (localK g r) :=
    Ioo_mem_nhds (by linarith) (by linarith)
  filter_upwards [hu.eventually hI,hv.eventually hI,hp.eventually_ne hgr,
    hp₂.eventually_ne hgr,hga.eventually_ne hgr,hgs.eventually_ne hgr,
    hgd.eventually_ne hgr,hgb.eventually_ne hgr,hbn,
    (tendsto_order.1 hq).1 0 (by norm_num),
    (tendsto_order.1 ht).2 1 (by norm_num),
    (tendsto_order.1 hJ).1 0 (by norm_num)] with n hu hv hp hp₂ hga hgs hgd hgb hb hq ht hJ
  exact ⟨hu,hv,hp,hp₂,hga,hgs,hgd,hgb,hb.1,hb.2.1,hb.2.2,hq,ht,hJ⟩


/-- The coefficient in the uniform cone theorem is exactly the paper's
C³ coefficient, not a separately chosen constant. -/
theorem localK_regularizedReciprocal
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    localK (Regularity.regularizedReciprocal f r) r = AnalyticBridge.errorConstant f r := by
  have h := Regularity.regularizedReciprocal_error_constant_of_contDiff hf hzero hd
  have hder : iteratedDeriv 2 (Regularity.regularizedReciprocal f r) r =
      deriv (deriv (Regularity.regularizedReciprocal f r)) r := by
    rw [show 2 = 1+1 from rfl, iteratedDeriv_succ', iteratedDeriv_one]
  rw [hder] at h
  simpa only [localK,AnalyticBridge.errorConstant] using h

/-- Full finite-C³ uniform-cone theorem for the actual removable reciprocal. -/
theorem smooth_exists_uniform_negative_coefficients
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r < 0) :
    let g := Regularity.regularizedReciprocal f r
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        firstCoefficient g r a s d ∈ Ioo (2*AnalyticBridge.errorConstant f r) (AnalyticBridge.errorConstant f r/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*AnalyticBridge.errorConstant f r) (AnalyticBridge.errorConstant f r/2) ∧
        firstPole g r a s d ≠ 0 ∧ secondPole g r a s d ≠ 0 ∧
        g (r+a) ≠ 0 ∧ g (r+s) ≠ 0 ∧ g (r+d) ≠ 0 ∧
        g (r+firstError g r a s d) ≠ 0 ∧
        firstError g r a s d ≠ 0 ∧ firstError g r a s d ≠ a ∧
        firstError g r a s d ≠ s ∧
        0 < normalizedGap g r a s d ∧ ordinateRatio g r a s d < 1 ∧
        0 < g (r+a) / g (r+firstError g r a s d) := by
  dsimp only
  have hg := Regularity.regularizedReciprocal_contDiffAt_of_contDiff hf hd
  have hgr : Regularity.regularizedReciprocal f r r ≠ 0 := by
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hk := localK_regularizedReciprocal hf hzero hd
  have result := exists_uniform_negative_coefficients hg hgr (by rw [hk]; exact hK)
  simpa only [hk] using result

end
end ModAB.Rational.CycleUniform

end RModABBundle33

/- Original module: ModAB.RationalCycleState -/
section RModABBundle34

/-!
The negative-coefficient local pair orbit is constructed from the interpolation
coefficients themselves. Cone invariance, convergence and coefficient limits
are derived, rather than supplied as fields of a hypothetical orbit.
-/
namespace ModAB.Rational.CycleState
open Filter Set
open scoped Topology
open CycleUniform
noncomputable section

structure PairState where
  a : ℝ
  s : ℝ
  d : ℝ

def nextPair (g : ℝ → ℝ) (r : ℝ) (p : PairState) : PairState :=
  let b := firstError g r p.a p.s p.d
  ⟨secondCoefficient g r p.a p.s p.d * p.a * b * p.d, b, p.s⟩

def orbit (g : ℝ → ℝ) (r : ℝ) (p : PairState) : ℕ → PairState
  | 0 => p
  | n+1 => nextPair g r (orbit g r p n)

structure Bounds (g : ℝ → ℝ) (r θ ε C : ℝ) : Prop where
  theta_pos : 0 < θ
  theta_lt : θ < 1
  radius_pos : 0 < ε
  bound_pos : 0 < C
  small : C * ε^2 ≤ θ
  coefficients : ∀ a s d, CycleEntry.Cone θ ε a s d →
    firstCoefficient g r a s d < 0 ∧ secondCoefficient g r a s d < 0 ∧
    |firstCoefficient g r a s d| ≤ C ∧ |secondCoefficient g r a s d| ≤ C

theorem cone_enlarge_theta {θ θ' ε a s d : ℝ}
    (h : CycleEntry.Cone θ ε a s d) (hθ : θ ≤ θ') :
    CycleEntry.Cone θ' ε a s d := by
  exact ⟨h.same,h.opposite,
    h.latest_small.trans (mul_le_mul_of_nonneg_right hθ (abs_nonneg s)),
    h.opposite_small.trans (mul_le_mul_of_nonneg_right hθ (abs_nonneg d)),h.history_bound⟩

/-- Analytic guard neighbourhoods can be intersected with the coefficient
cone without losing the derived coefficient bounds. -/
theorem Bounds.restrict {g : ℝ → ℝ} {r θ ε C θ' ε' : ℝ}
    (N : Bounds g r θ ε C) (hθ0 : 0 < θ') (hθ : θ' ≤ θ)
    (hε0 : 0 < ε') (hε : ε' ≤ ε) (hsmall : C*ε'^2 ≤ θ') :
    Bounds g r θ' ε' C := by
  refine ⟨hθ0,hθ.trans_lt N.theta_lt,hε0,N.bound_pos,hsmall,?_⟩
  intro a s d hc
  exact N.coefficients _ _ _ ((cone_enlarge_theta hc hθ).enlarge_radius hε)

theorem exists_bounds {g : ℝ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r < 0) :
    ∃ θ ε C, Bounds g r θ ε C := by
  obtain ⟨θ,δ,hθ0,hθ1,hδ,hall⟩ := exists_uniform_negative_coefficients hg hgr hK
  let C := -2*localK g r
  have hC : 0 < C := by dsimp [C]; linarith
  let ε := min δ (min 1 (θ/C))
  have hε : 0 < ε := lt_min hδ (lt_min (by norm_num) (div_pos hθ0 hC))
  have hεδ : ε ≤ δ := min_le_left _ _
  have hε1 : ε ≤ 1 := (min_le_right _ _).trans (min_le_left _ _)
  have hεC : ε ≤ θ/C := (min_le_right _ _).trans (min_le_right _ _)
  have hCε : C*ε ≤ θ := by
    have hh := (le_div_iff₀ hC).mp hεC
    nlinarith
  refine ⟨θ,ε,C,⟨hθ0,hθ1,hε,hC,?_,?_⟩⟩
  · have hsq : ε^2 ≤ ε := by nlinarith
    exact (mul_le_mul_of_nonneg_left hsq hC.le).trans hCε
  · intro a s d hc
    have hh := hall a s d (hc.enlarge_radius hεδ)
    have hu : firstCoefficient g r a s d < 0 := by linarith [hh.1.2]
    have hv : secondCoefficient g r a s d < 0 := by linarith [hh.2.1.2]
    refine ⟨hu,hv,?_,?_⟩
    · rw [abs_of_neg hu]
      dsimp [C]
      linarith [hh.1.1]
    · rw [abs_of_neg hv]
      dsimp [C]
      linarith [hh.2.1.1]

theorem orbit_cone {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    ∀ n, CycleEntry.Cone θ ε (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d := by
  intro n
  induction n with
  | zero => exact hp
  | succ n ih =>
    obtain ⟨hu,hv,huC,hvC⟩ := N.coefficients _ _ _ ih
    have hc := CycleEntry.pair_preserves_cone N.theta_pos N.theta_lt
      N.bound_pos.le N.radius_pos.le N.small ih hu hv huC hvC
      (show firstError g r (orbit g r p n).a (orbit g r p n).s
        (orbit g r p n).d = _ by rfl) (show _ = _ by rfl)
    have he : θ*ε ≤ ε := by nlinarith [N.theta_lt,N.radius_pos]
    exact hc.enlarge_radius he

/-- The established abstract recurrence interface is now populated by a
concrete recursively defined interpolation orbit. -/
def toRun {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    CycleEntry.Run θ ε C where
  a n := (orbit g r p n).a
  s n := (orbit g r p n).s
  d n := (orbit g r p n).d
  u n := firstCoefficient g r (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d
  v n := secondCoefficient g r (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d
  initial := hp
  u_neg n := (N.coefficients _ _ _ (orbit_cone N hp n)).1
  v_neg n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.1
  u_bound n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.2.1
  v_bound n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.2.2
  first_step _n := rfl
  second_step _n := rfl
  drop_step _n := rfl

/-- Local regularity and the derived cone imply both coefficient limits;
neither convergence nor asymptotic separation is assumed. -/
theorem coefficient_limits {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    Tendsto (toRun N hp).u atTop (𝓝 (localK g r)) ∧
    Tendsto (toRun N hp).v atTop (𝓝 (localK g r)) := by
  let R := toRun N hp
  have hc := R.all_cones N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have he := R.errors_tendsto N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have hs := R.separation_tendsto N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have hd (n : ℕ) := CycleUniform.Cone.distinct (hc n) N.theta_lt
  have hu := tendsto_firstCoefficient hg hgr (fun n => (hd n).1)
    (fun n => (hd n).2.1) (fun n => (hd n).2.2) he.1 he.2.1 he.2.2
  have hv := second_step_limits hg hgr hK (fun n => (hd n).1)
    (fun n => (hd n).2.1) (fun n => (hd n).2.2)
    (fun n => (hc n).latest_ne_zero) (fun n => (hc n).opposite_ne_zero)
    (fun n => (hc n).history_ne_zero) he.1 he.2.1 he.2.2 hs.1 hs.2
  exact ⟨hu,hv.2.2.2⟩

theorem orbit_hasQOrders {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    ModAB.RationalOrder.HasQOrder (fun n => (orbit g r p n).a) 3 |localK g r| ∧
    ModAB.RationalOrder.HasQOrder (fun n => (orbit g r p (n+1)).s) 3 |localK g r| ∧
    Tendsto (fun n => |(orbit g r p (n+1)).s| /
      |(orbit g r p n).a| ^ (9/5:ℝ)) atTop (𝓝 (|localK g r| ^ (2/5:ℝ))) ∧
    Tendsto (fun n => |(orbit g r p (n+1)).a| /
      |(orbit g r p (n+1)).s| ^ (5/3:ℝ)) atTop (𝓝 (|localK g r| ^ (1/3:ℝ))) := by
  have hh := coefficient_limits hg hgr hK N hp
  exact ModAB.RationalLocalAssembly.cycle_cone_hasQOrders (toRun N hp)
    N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small hK hh.1 hh.2

end
end ModAB.Rational.CycleState

end RModABBundle34

/- Original module: ModAB.RationalCycleUniformAlgebra -/
section RModABBundle35

/-!
Exact identities for the AB-corrected second step of a local rational pair.
The normalization removes all divisions by a vanishing node from the limit
calculation except the two separation ratios. No asymptotic bound is assumed.
-/
namespace ModAB.Rational.CycleUniformAlgebra
open LocalAlgebra
noncomputable section

/-- The denominator in the AB multiplier is retained exactly. -/
theorem delta_one_sub {t : ℝ} (ht : 1-t ≠ 0) :
    1/(1-t)-1=t/(1-t) := by
  field_simp [ht]
  ring

/-- Residual ratio recovered from the actual first-step product identity. -/
theorem residual_ratio_from_step {a b s d u ya yb : ℝ}
    (_ha : a ≠ 0) (hb : b ≠ 0) (hya : ya ≠ 0) (hyb : yb ≠ 0)
    (hstep : b=u*a*s*d) :
    yb/ya=u*s*d*((a/ya)/(b/yb)) := by
  field_simp [hb,hya,hyb]
  rw [hstep]
  ring

/-- The relative node denominator has a finite limit in a separated cone. -/
def relativeDenominator (a b s : ℝ) : ℝ := (1-b/s)*(1-a/s)

theorem relativeDenominator_eq {a b s : ℝ} (hs : s ≠ 0) :
    relativeDenominator a b s = ((s-b)*(s-a))/s^2 := by
  unfold relativeDenominator
  field_simp [hs]

theorem relativeDenominator_ne_zero {a b s : ℝ}
    (hs : s ≠ 0) (hsb : s-b ≠ 0) (hsa : s-a ≠ 0) :
    relativeDenominator a b s ≠ 0 := by
  rw [relativeDenominator_eq hs]
  exact div_ne_zero (mul_ne_zero hsb hsa) (pow_ne_zero 2 hs)

/-- Continuous expression for the second coefficient, once G and P are
extended by their local divided-difference limits. J is g(a)/g(b). -/
def normalizedSecond (G P gs u J a b s d : ℝ) : ℝ :=
  ((s/d)*G+u*gs*J/(relativeDenominator a b s*(1-u*s*d*J))) /
    (P+((u*s*d*J)/(1-u*s*d*J))*gs*(a/s)*(b/s)/relativeDenominator a b s)

/-- The actual mixed pole denominator equals its nonsingular relative form. -/
theorem mixed_denominator_identity {P gs a b s t : ℝ}
    (hs : s ≠ 0) (hsb : s-b ≠ 0) (hsa : s-a ≠ 0) (ht : 1-t ≠ 0) :
    P+(1/(1-t)-1)*gs*b*a/((s-b)*(s-a)) =
      P+(t/(1-t))*gs*(a/s)*(b/s)/relativeDenominator a b s := by
  rw [delta_one_sub ht,relativeDenominator_eq hs]
  field_simp [hs,hsb,hsa]

/-- The normalized pole form controls the barycentric denominator exactly. -/
theorem bary_mixed_pole_identity {r a b s ya yb ys t : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hys : ys ≠ 0) (ht : 1-t ≠ 0) :
    dd2 (r+b) (r+a) (r+s) (1/yb) (1/ya) (1/((1-t)*ys))*(b*a*s) =
      poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys) +
        (t/(1-t))*(s/ys)*(a/s)*(b/s)/relativeDenominator a b s := by
  have hba : (r+a)-(r+b) ≠ 0 := by simpa using hab
  have hpa : (r+s)-(r+b) ≠ 0 := by simpa using hsb
  have hpb : (r+s)-(r+a) ≠ 0 := by simpa using hsa
  have h := pole_identity_from_values (r:=r) (x₁:=r+b) (x₂:=r+a) (x₃:=r+s)
    (y₁:=yb) (y₂:=ya) (y₃:=(1-t)*ys)
    hba hpb hpa (by simpa using hb) (by simpa using ha) (by simpa using hs)
  rw [poleValue_scaled_third hpb hpa hys ht] at h
  simp only [add_sub_cancel_left,add_sub_add_left_eq_sub] at h
  rw [h]
  exact mixed_denominator_identity hs hsb hsa ht

theorem bary_mixed_denominator_ne_zero {r a b s ya yb ys t : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hys : ys ≠ 0) (ht : 1-t ≠ 0)
    (hP : poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys) +
      (t/(1-t))*(s/ys)*(a/s)*(b/s)/relativeDenominator a b s ≠ 0) :
    dd2 (r+b) (r+a) (r+s) (1/yb) (1/ya) (1/((1-t)*ys)) ≠ 0 := by
  intro hz
  apply hP
  rw [←bary_mixed_pole_identity ha hb hs hab hsa hsb hys ht,hz,zero_mul]

/-- Nonvanishing of the nonsingular mixed pole implies the implementation's
denominator guard in actual bracket order, so no independent guard is needed. -/
theorem source_denominator_ne_zero_of_mixed_pole {r a b s ya yb ys t : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hys : ys ≠ 0) (ht : 1-t ≠ 0)
    (hP : poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys) +
      (t/(1-t))*(s/ys)*(a/s)*(b/s)/relativeDenominator a b s ≠ 0) :
    ModAB.Rational.denominator (r+b) (r+s) yb ((1-t)*ys) (r+a) ya ≠ 0 := by
  have hba : (r+a)-(r+b) ≠ 0 := by simpa using hab
  have hpa : (r+s)-(r+b) ≠ 0 := by simpa using hsb
  have hpb : (r+s)-(r+a) ≠ 0 := by simpa using hsa
  have hap : (r+a)-(r+s) ≠ 0 := by intro hz; apply hpb; linarith
  have hH := bary_mixed_denominator_ne_zero ha hb hs hab hsa hsb hys ht hP
  rw [dd2_swap_last hba hpb hpa] at hH
  intro hz
  apply hH
  rw [bary_denominator_eq hpa hap hba hyb (mul_ne_zero ht hys) hya,hz]
  simp

/-- Exact finite-data normalized coefficient, including the case of a zero
mixed pole value (where the totalized quotients on both sides vanish). -/
theorem normalized_second_identity {G P gs u J a b s d t : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hsb : s-b ≠ 0) (hsa : s-a ≠ 0) (ht : 1-t ≠ 0)
    (hres : t=u*s*d*J) :
    ((G+(1/(1-t)-1)*gs/((s-b)*(s-a))) /
      (P+(1/(1-t)-1)*gs*b*a/((s-b)*(s-a)))*(b*a*s))/(a*b*d) =
      normalizedSecond G P gs u J a b s d := by
  rw [delta_one_sub ht]
  have hq := relativeDenominator_ne_zero hs hsb hsa
  have hden : P+t/(1-t)*gs*b*a/((s-b)*(s-a)) =
      P+(t/(1-t))*gs*(a/s)*(b/s)/relativeDenominator a b s := by
    rw [relativeDenominator_eq hs]
    field_simp [hs,hsb,hsa]
  rw [hden]
  unfold normalizedSecond
  rw [←hres]
  by_cases hD : P+t/(1-t)*gs*(a/s)*(b/s)/relativeDenominator a b s=0
  · simp only [hD,div_zero,zero_mul,zero_div]
  · field_simp [ha,hb,hs,hd,hsb,hsa,ht,hq,hD]
    rw [hres] at *
    rw [relativeDenominator_eq hs]
    field_simp [hs]

/-- The same identity before normalization of a*b*d. -/
theorem second_error_eq {G P gs u J a b s d t : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hsb : s-b ≠ 0) (hsa : s-a ≠ 0) (ht : 1-t ≠ 0)
    (hres : t=u*s*d*J) :
    (G+(1/(1-t)-1)*gs/((s-b)*(s-a))) /
      (P+(1/(1-t)-1)*gs*b*a/((s-b)*(s-a)))*(b*a*s) =
      normalizedSecond G P gs u J a b s d*(a*b*d) := by
  have h := normalized_second_identity (G:=G) (P:=P) (gs:=gs) ha hb hs hd hsb hsa ht hres
  exact (div_eq_iff (mul_ne_zero (mul_ne_zero ha hb) hd)).mp h

/-- Actual mixed barycentric proposal, with arbitrary root and supplied data.
The hypotheses are precisely nonvanishing guards and the first step equation. -/
theorem baryRoot_second_error {r a b s d u ya yb ys : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hys : ys ≠ 0)
    (ht : 1-yb/ya ≠ 0) (hstep : b=u*a*s*d)
    (hH : dd2 (r+b) (r+a) (r+s) (1/yb) (1/ya)
      (1/((1-yb/ya)*ys)) ≠ 0) :
    baryRoot (r+b) (r+a) (r+s) yb ya ((1-yb/ya)*ys)-r =
      normalizedSecond
        (dd2 (r+b) (r+a) (r+s) (b/yb) (a/ya) (s/ys))
        (poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys))
        (s/ys) u ((a/ya)/(b/yb)) a b s d*(a*b*d) := by
  have h := baryRoot_scaled_error (r:=r) (x₁:=r+b) (x₂:=r+a) (x₃:=r+s)
    (y₁:=yb) (y₂:=ya) (y₃:=ys) (μ:=1-yb/ya)
    (by simpa using hab) (by simpa using hsa) (by simpa using hsb)
    (by simpa using hb) (by simpa using ha) (by simpa using hs) hys ht hH
  simp only [add_sub_cancel_left,add_sub_add_left_eq_sub] at h
  rw [h]
  exact second_error_eq ha hb hs hd hsb hsa ht
    (residual_ratio_from_step ha hb hya hyb hstep)

/-- Source normalized rational proposal has the same second-step identity. -/
theorem rawPoint_second_error {r a b s d u ya yb ys : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hys : ys ≠ 0)
    (ht : 1-yb/ya ≠ 0) (hstep : b=u*a*s*d)
    (hD : ModAB.Rational.denominator (r+b) (r+a) yb ya
      (r+s) ((1-yb/ya)*ys) ≠ 0) :
    ModAB.Rational.rawPoint (r+b) (r+a) yb ya (r+s) ((1-yb/ya)*ys)-r =
      normalizedSecond
        (dd2 (r+b) (r+a) (r+s) (b/yb) (a/ya) (s/ys))
        (poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys))
        (s/ys) u ((a/ya)/(b/yb)) a b s d*(a*b*d) := by
  have hba : (r+a)-(r+b) ≠ 0 := by simpa using hab
  have hpa : (r+s)-(r+b) ≠ 0 := by simpa using hsb
  have hpb : (r+s)-(r+a) ≠ 0 := by simpa using hsa
  have hyp := mul_ne_zero ht hys
  rw [rawPoint_eq_baryRoot hba hpb hpa hyb hya hyp hD]
  exact baryRoot_second_error ha hb hs hd hab hsa hsb hya hyb hys ht hstep
    (bary_denominator_ne_zero hba hpb hpa hyb hya hyp hD)

/-- The same identity in implementation order: b and s are the current bracket
endpoints; the dropped point a is supplied with its true value. -/
theorem source_rawPoint_second_error {r a b s d u ya yb ys : ℝ}
    (ha : a ≠ 0) (hb : b ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hab : a-b ≠ 0) (hsa : s-a ≠ 0) (hsb : s-b ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hys : ys ≠ 0)
    (ht : 1-yb/ya ≠ 0) (hstep : b=u*a*s*d)
    (hD : ModAB.Rational.denominator (r+b) (r+s) yb ((1-yb/ya)*ys)
      (r+a) ya ≠ 0) :
    ModAB.Rational.rawPoint (r+b) (r+s) yb ((1-yb/ya)*ys) (r+a) ya-r =
      normalizedSecond
        (dd2 (r+b) (r+a) (r+s) (b/yb) (a/ya) (s/ys))
        (poleValue (r+b) (r+a) (r+s) r (b/yb) (a/ya) (s/ys))
        (s/ys) u ((a/ya)/(b/yb)) a b s d*(a*b*d) := by
  have hba : (r+a)-(r+b) ≠ 0 := by simpa using hab
  have hpa : (r+s)-(r+b) ≠ 0 := by simpa using hsb
  have hpb : (r+s)-(r+a) ≠ 0 := by simpa using hsa
  have hap : (r+a)-(r+s) ≠ 0 := by intro hz; apply hpb; linarith
  have hyp := mul_ne_zero ht hys
  rw [rawPoint_eq_baryRoot hpa hap hba hyb hyp hya hD]
  rw [←baryRoot_swap_last hba hpb hpa]
  apply baryRoot_second_error ha hb hs hd hab hsa hsb hya hyb hys ht hstep
  rw [dd2_swap_last hba hpb hpa]
  exact bary_denominator_ne_zero hpa hap hba hyb hyp hya hD

/-- Exchanging the two current endpoints preserves the raw denominator guard. -/
theorem denominator_swap_ne_zero {a b p ya yb yp : ℝ}
    (hba : b-a ≠ 0) (hpb : p-b ≠ 0) (hpa : p-a ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hyp : yp ≠ 0)
    (hD : ModAB.Rational.denominator a b ya yb p yp ≠ 0) :
    ModAB.Rational.denominator b a yb ya p yp ≠ 0 := by
  have hab : a-b ≠ 0 := by intro hz; apply hba; linarith
  have hH := bary_denominator_ne_zero hba hpb hpa hya hyb hyp hD
  rw [dd2_swap_first hba hpb hpa] at hH
  intro hz
  apply hH
  rw [bary_denominator_eq hab hpa hpb hyb hya hyp,hz]
  simp

/-- Endpoint ordering changes the normalization, but not the interpolating zero. -/
theorem rawPoint_swap {a b p ya yb yp : ℝ}
    (hba : b-a ≠ 0) (hpb : p-b ≠ 0) (hpa : p-a ≠ 0)
    (hya : ya ≠ 0) (hyb : yb ≠ 0) (hyp : yp ≠ 0)
    (hD : ModAB.Rational.denominator a b ya yb p yp ≠ 0) :
    ModAB.Rational.rawPoint a b ya yb p yp =
      ModAB.Rational.rawPoint b a yb ya p yp := by
  have hab : a-b ≠ 0 := by intro hz; apply hba; linarith
  have hD' := denominator_swap_ne_zero hba hpb hpa hya hyb hyp hD
  rw [rawPoint_eq_baryRoot hba hpb hpa hya hyb hyp hD,
    rawPoint_eq_baryRoot hab hpa hpb hyb hya hyp hD']
  unfold baryRoot
  rw [dd2_swap_first hba hpb hpa,dd2_swap_first hba hpb hpa]

end
end ModAB.Rational.CycleUniformAlgebra

end RModABBundle35

/- Original module: ModAB.RationalCycleUniformBridge -/
section RModABBundle36

/-!
Connect the uniformly controlled coefficients to the actual function values
and the source-normalized raw proposal. Nonvanishing guards follow from the
normalized pole values and regularized reciprocal values.
-/
namespace ModAB.Rational.CycleUniformBridge
open Set
open LocalAlgebra Regularity CycleUniform CycleUniformAlgebra
noncomputable section

theorem reciprocal_value_at_error {f : ℝ → ℝ} {r e : ℝ}
    (hzero : f r = 0) (he : e ≠ 0) :
    e / f (r+e) = regularizedReciprocal f r (r+e) := by
  symm
  simpa only [add_sub_cancel_left] using
    regularizedReciprocal_off_root hzero (show r+e ≠ r by simpa using he)

theorem ordinate_ne_zero_of_reciprocal {f : ℝ → ℝ} {r e : ℝ}
    (hzero : f r = 0) (he : e ≠ 0)
    (hg : regularizedReciprocal f r (r+e) ≠ 0) : f (r+e) ≠ 0 := by
  intro hz
  apply hg
  rw [←reciprocal_value_at_error hzero he,hz,div_zero]

/-- Uniform control of the first pole and reciprocal values supplies the
actual raw-step denominator guard. -/
theorem first_source_denominator_ne_zero {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (has : a ≠ s) (hsd : s ≠ d) (had : a ≠ d)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgd : regularizedReciprocal f r (r+d) ≠ 0)
    (hP : firstPole (regularizedReciprocal f r) r a s d ≠ 0) :
    ModAB.Rational.denominator (r+a) (r+s) (f (r+a)) (f (r+s))
      (r+d) (f (r+d)) ≠ 0 := by
  have hsa : (r+s)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr has.symm
  have hds : (r+d)-(r+s) ≠ 0 := by simpa using sub_ne_zero.mpr hsd.symm
  have hda : (r+d)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr had.symm
  have hp := pole_identity_from_values (r:=r) (x₁:=r+a) (x₂:=r+s) (x₃:=r+d)
    (y₁:=f (r+a)) (y₂:=f (r+s)) (y₃:=f (r+d)) hsa hds hda
    (by simpa using ha) (by simpa using hs) (by simpa using hd)
  simp only [add_sub_cancel_left] at hp
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hs,
    reciprocal_value_at_error hzero hd] at hp
  have hH : dd2 (r+a) (r+s) (r+d) (1/f (r+a)) (1/f (r+s)) (1/f (r+d)) ≠ 0 := by
    intro hz
    apply hP
    change poleValue (r+a) (r+s) (r+d) r
      (regularizedReciprocal f r (r+a)) (regularizedReciprocal f r (r+s))
      (regularizedReciprocal f r (r+d)) = 0
    rw [←hp,hz,zero_mul]
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfd := ordinate_ne_zero_of_reciprocal hzero hd hgd
  intro hz
  apply hH
  rw [bary_denominator_eq hsa hds hda hfa hfs hfd,hz]
  simp

/-- The first coefficient is the coefficient of the implemented raw proposal,
not merely a separate divided-difference model. -/
theorem first_rawPoint_error {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (has : a ≠ s) (hsd : s ≠ d) (had : a ≠ d)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgd : regularizedReciprocal f r (r+d) ≠ 0)
    (hP : firstPole (regularizedReciprocal f r) r a s d ≠ 0) :
    ModAB.Rational.rawPoint (r+a) (r+s) (f (r+a)) (f (r+s))
      (r+d) (f (r+d))-r = firstError (regularizedReciprocal f r) r a s d := by
  have hD := first_source_denominator_ne_zero hzero ha hs hd has hsd had hga hgs hgd hP
  have hsa : (r+s)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr has.symm
  have hds : (r+d)-(r+s) ≠ 0 := by simpa using sub_ne_zero.mpr hsd.symm
  have hda : (r+d)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr had.symm
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfd := ordinate_ne_zero_of_reciprocal hzero hd hgd
  rw [rawPoint_eq_baryRoot hsa hds hda hfa hfs hfd hD]
  rw [baryRoot_true_error hsa hds hda (by simpa using ha) (by simpa using hs)
    (by simpa using hd) (bary_denominator_ne_zero hsa hds hda hfa hfs hfd hD)]
  simp only [add_sub_cancel_left]
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hs,
    reciprocal_value_at_error hzero hd]
  unfold firstError firstCoefficient firstPole
  ring

/-- The exact residual ratio used by AB equals the normalized local expression. -/
theorem residual_ratio_eq_ordinateRatio {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0) :
    f (r+firstError (regularizedReciprocal f r) r a s d) / f (r+a) =
      ordinateRatio (regularizedReciprocal f r) r a s d := by
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  have h := residual_ratio_from_step ha hb hfa hfb
    (show firstError (regularizedReciprocal f r) r a s d =
      firstCoefficient (regularizedReciprocal f r) r a s d*a*s*d from rfl)
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb] at h
  exact h

/-- The second normalized pole is exactly the pole in the implementation's
AB-corrected data. -/
theorem second_source_denominator_ne_zero {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hba : firstError (regularizedReciprocal f r) r a s d ≠ a)
    (hbs : firstError (regularizedReciprocal f r) r a s d ≠ s) (has : a ≠ s)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0)
    (ht : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0)
    (hP : secondPole (regularizedReciprocal f r) r a s d ≠ 0) :
    let b := firstError (regularizedReciprocal f r) r a s d
    ModAB.Rational.denominator (r+b) (r+s) (f (r+b))
      ((1-f (r+b)/f (r+a))*f (r+s)) (r+a) (f (r+a)) ≠ 0 := by
  dsimp only
  let g := regularizedReciprocal f r
  let b := firstError g r a s d
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  apply source_denominator_ne_zero_of_mixed_pole ha hb hs
    (sub_ne_zero.mpr hba.symm) (sub_ne_zero.mpr has.symm) (sub_ne_zero.mpr hbs.symm)
    hfa hfb hfs (by rw [hratio]; exact ht)
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb,
    reciprocal_value_at_error hzero hs,hratio]
  exact hP

/-- Actual second raw proposal equals the uniformly controlled coefficient
times a*b*d. The AB factor is formed from the true residuals. -/
theorem second_rawPoint_error {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hba : firstError (regularizedReciprocal f r) r a s d ≠ a)
    (hbs : firstError (regularizedReciprocal f r) r a s d ≠ s) (has : a ≠ s)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0)
    (ht : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0)
    (hP : secondPole (regularizedReciprocal f r) r a s d ≠ 0) :
    let g := regularizedReciprocal f r
    let b := firstError g r a s d
    ModAB.Rational.rawPoint (r+b) (r+s) (f (r+b))
      ((1-f (r+b)/f (r+a))*f (r+s)) (r+a) (f (r+a))-r =
      secondCoefficient g r a s d*(a*b*d) := by
  dsimp only
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  have hD := second_source_denominator_ne_zero hzero ha hs hb hba hbs has hga hgs hgb ht hP
  have h := source_rawPoint_second_error ha hb hs hd (sub_ne_zero.mpr hba.symm)
    (sub_ne_zero.mpr has.symm) (sub_ne_zero.mpr hbs.symm) hfa hfb hfs
    (by rw [hratio]; exact ht)
    (show firstError (regularizedReciprocal f r) r a s d =
      firstCoefficient (regularizedReciprocal f r) r a s d*a*s*d from rfl) hD
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb,
    reciprocal_value_at_error hzero hs] at h
  exact h

/-- Finite-C³ hypotheses alone give a cone on which both actual rational
proposals have bounded negative coefficients, the AB multiplier lies strictly
between zero and one, and both implementation denominator guards hold.
No divided-difference estimate or asymptotic recurrence is a premise. -/
theorem smooth_exists_uniform_raw_pair
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hder : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r < 0) :
    let g := regularizedReciprocal f r
    let K := AnalyticBridge.errorConstant f r
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        let b := firstError g r a s d
        let μ := 1-f (r+b)/f (r+a)
        firstCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
        ModAB.Rational.denominator (r+a) (r+s) (f (r+a)) (f (r+s))
          (r+d) (f (r+d)) ≠ 0 ∧
        ModAB.Rational.rawPoint (r+a) (r+s) (f (r+a)) (f (r+s))
          (r+d) (f (r+d))-r = b ∧
        μ ∈ Ioo (0:ℝ) 1 ∧
        ModAB.Rational.denominator (r+b) (r+s) (f (r+b))
          (μ*f (r+s)) (r+a) (f (r+a)) ≠ 0 ∧
        ModAB.Rational.rawPoint (r+b) (r+s) (f (r+b))
          (μ*f (r+s)) (r+a) (f (r+a))-r =
          secondCoefficient g r a s d*(a*b*d) ∧
        b ≠ 0 ∧ b ≠ a ∧ b ≠ s ∧ 0 < a*b := by
  dsimp only
  obtain ⟨θ,ε,hθ0,hθ1,hε,huniform⟩ :=
    smooth_exists_uniform_negative_coefficients hf hzero hder hK
  refine ⟨θ,ε,hθ0,hθ1,hε,?_⟩
  intro a s d hc
  obtain ⟨hU,hV,hP₁,hP₂,hga,hgs,hgd,hgb,hb,hba,hbs,hq,ht,hJ⟩ :=
    huniform a s d hc
  have ha := hc.latest_ne_zero
  have hs := hc.opposite_ne_zero
  have hd := hc.history_ne_zero
  obtain ⟨has,hsd,had⟩ := CycleUniform.Cone.distinct hc hθ1
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hmu : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0 := by
    linarith
  have hu : firstCoefficient (regularizedReciprocal f r) r a s d < 0 := by
    linarith [hU.2]
  have hsdneg : s*d < 0 := by
    have hprod : a^2*(s*d) < 0 := by
      convert mul_neg_of_neg_of_pos hc.opposite hc.same using 1; ring
    by_contra hh
    have hn := mul_nonneg (sq_nonneg a) (le_of_not_gt hh)
    linarith
  have hupos : 0 < firstCoefficient (regularizedReciprocal f r) r a s d*s*d := by
    convert mul_pos_of_neg_of_neg hu hsdneg using 1; ring
  have htpos : 0 < ordinateRatio (regularizedReciprocal f r) r a s d :=
    mul_pos hupos hJ
  have hmubounds : 1-f (r+firstError (regularizedReciprocal f r) r a s d)/f (r+a)
      ∈ Ioo (0:ℝ) 1 := by
    rw [hratio]
    exact ⟨by linarith,by linarith⟩
  have habpos : 0 < a*firstError (regularizedReciprocal f r) r a s d := by
    unfold firstError
    convert mul_pos (sq_pos_of_ne_zero ha) hupos using 1; ring
  exact ⟨hU,hV,
    first_source_denominator_ne_zero hzero ha hs hd has hsd had hga hgs hgd hP₁,
    first_rawPoint_error hzero ha hs hd has hsd had hga hgs hgd hP₁,
    hmubounds,
    second_source_denominator_ne_zero hzero ha hs hb hba hbs has hga hgs hgb hmu hP₂,
    second_rawPoint_error hzero ha hs hd hb hba hbs has hga hgs hgb hmu hP₂,
    hb,hba,hbs,habpos⟩

end
end ModAB.Rational.CycleUniformBridge

end RModABBundle36

/- Original module: ModAB.RationalSecantLocal -/
section RModABBundle37

/-! Uniform simple-root estimates and the secant part of the actual clamp.
The estimates follow from differentiability, not from a separately assumed
secant-error or residual-comparability hypothesis. -/
namespace ModAB.Rational.SecantLocal
open Set Filter Topology
noncomputable section

def normalizedSlope (f : ℝ → ℝ) (r x : ℝ) : ℝ :=
  dslope f r x / deriv f r

theorem normalizedSlope_at {f : ℝ → ℝ} {r : ℝ} (hd : deriv f r ≠ 0) :
    normalizedSlope f r r = 1 := by simp [normalizedSlope, dslope_same, hd]

theorem normalizedSlope_continuousAt {f : ℝ → ℝ} {r : ℝ}
    (hf : DifferentiableAt ℝ f r) : ContinuousAt (normalizedSlope f r) r :=
  (continuousAt_dslope_same.mpr hf).div_const _

theorem function_eq_normalizedSlope {f : ℝ → ℝ} {r x : ℝ}
    (hf : f r = 0) (hd : deriv f r ≠ 0) :
    f x = deriv f r * (x-r) * normalizedSlope f r x := by
  by_cases hx : x=r
  · simp [hx,hf]
  · simp only [normalizedSlope, dslope_of_ne f hx, slope_def_field, hf, sub_zero]
    field_simp

theorem exists_radius_normalizedSlope {f : ℝ → ℝ} {r δ : ℝ}
    (hf : DifferentiableAt ℝ f r) (hd : deriv f r ≠ 0) (hδ : 0 < δ) :
    ∃ ε > 0, ∀ x : ℝ, |x-r| < ε → |normalizedSlope f r x - 1| < δ := by
  obtain ⟨ε,hε,hball⟩ := Metric.continuousAt_iff.mp
    (normalizedSlope_continuousAt hf) δ hδ
  refine ⟨ε,hε,fun x hx => ?_⟩
  simpa only [Real.dist_eq, normalizedSlope_at hd] using hball (by simpa [Real.dist_eq] using hx)

theorem residual_bounds_of_normalizedSlope {f : ℝ → ℝ} {r x : ℝ}
    (hf : f r = 0) (hd : deriv f r ≠ 0)
    (hq : |normalizedSlope f r x - 1| ≤ (1/2:ℝ)) :
    |deriv f r|/2 * |x-r| ≤ |f x| ∧
      |f x| ≤ (2*|deriv f r|) * |x-r| := by
  rcases abs_le.mp hq with ⟨hl,hu⟩
  have hqp : 0 ≤ normalizedSlope f r x := by linarith
  rw [function_eq_normalizedSlope hf hd, abs_mul, abs_mul, abs_of_nonneg hqp]
  have hp : 0 ≤ |deriv f r| * |x-r| := mul_nonneg (abs_nonneg _) (abs_nonneg _)
  constructor
  · nlinarith [mul_nonneg hp (show 0 ≤ normalizedSlope f r x - 1/2 by linarith)]
  · nlinarith [mul_nonneg hp (show 0 ≤ 2 - normalizedSlope f r x by linarith)]

theorem sign_product_of_normalizedSlope {f : ℝ → ℝ} {r x y : ℝ}
    (hf : f r = 0) (hd : deriv f r ≠ 0)
    (hx : 0 < normalizedSlope f r x) (hy : 0 < normalizedSlope f r y) :
    (0 < f x * f y ↔ 0 < (x-r)*(y-r)) ∧
      (f x * f y < 0 ↔ (x-r)*(y-r) < 0) := by
  have hL : 0 < (deriv f r)^2 * (normalizedSlope f r x * normalizedSlope f r y) :=
    mul_pos (sq_pos_of_ne_zero hd) (mul_pos hx hy)
  have he : f x * f y =
      ((deriv f r)^2 * (normalizedSlope f r x * normalizedSlope f r y)) * ((x-r)*(y-r)) := by
    rw [function_eq_normalizedSlope hf hd (x:=x), function_eq_normalizedSlope hf hd (x:=y)]
    ring
  rw [he]
  constructor
  · exact mul_pos_iff_of_pos_left hL
  · constructor
    · intro hh
      by_contra he
      have hp := mul_nonneg hL.le (le_of_not_gt he)
      linarith
    · intro hh
      exact mul_neg_of_pos_of_neg hL hh

theorem exists_radius_simple_root_bounds {f : ℝ → ℝ} {r : ℝ}
    (hf : DifferentiableAt ℝ f r) (hr : f r = 0) (hd : deriv f r ≠ 0) :
    ∃ ε > 0, ∀ x : ℝ, |x-r| < ε →
      0 < normalizedSlope f r x ∧
      |deriv f r|/2 * |x-r| ≤ |f x| ∧
      |f x| ≤ (2*|deriv f r|)*|x-r| := by
  obtain ⟨ε,hε,h⟩ := exists_radius_normalizedSlope hf hd (by norm_num : (0:ℝ)<1/2)
  refine ⟨ε,hε,fun x hx => ?_⟩
  have hq := h x hx
  refine ⟨by linarith [(abs_lt.mp hq).1], residual_bounds_of_normalizedSlope hr hd hq.le⟩

/-- Algebraic secant estimate valid for any pair of positive normalized slopes.
It also applies to a mildly rescaled retained ordinate. -/
theorem secant_error_of_normalized_values {a b r L qa qb δ : ℝ}
    (ha : a < r) (hb : r < b) (hL : L ≠ 0)
    (hδ : 0 ≤ δ) (hδh : δ ≤ (1/2:ℝ))
    (hqa : |qa-1| ≤ δ) (hqb : |qb-1| ≤ δ) :
    |secant a b (L*(a-r)*qa) (L*(b-r)*qb)-r| ≤
      (4*δ) * min |a-r| |b-r| := by
  let u := r-a
  let v := b-r
  have hu : 0 < u := sub_pos.mpr ha
  have hv : 0 < v := sub_pos.mpr hb
  have hqa' := abs_le.mp hqa
  have hqb' := abs_le.mp hqb
  have hqap : (1/2:ℝ) ≤ qa := by linarith [hqa'.1]
  have hqbp : (1/2:ℝ) ≤ qb := by linarith [hqb'.1]
  have hden : 0 < v*qb+u*qa := by positivity
  have hdenE : L*(b-r)*qb-L*(a-r)*qa = L*(v*qb+u*qa) := by dsimp [u,v]; ring
  have hdenL : L*(b-r)*qb-L*(a-r)*qa ≠ 0 := by
    rw [hdenE]
    exact mul_ne_zero hL hden.ne'
  have he : secant a b (L*(a-r)*qa) (L*(b-r)*qb)-r =
      u*v*(qa-qb)/(v*qb+u*qa) := by
    unfold secant
    rw [hdenE]
    field_simp [hL, hden.ne']
    dsimp [u,v]
    ring
  have habsa : |a-r| = u := by rw [abs_of_neg (sub_neg.mpr ha)]; dsimp [u]; ring
  have habsb : |b-r| = v := abs_of_pos (sub_pos.mpr hb)
  rw [he, habsa, habsb,abs_div,abs_mul,abs_mul,abs_of_pos hu,abs_of_pos hv,abs_of_pos hden]
  have hdiff : |qa-qb| ≤ 2*δ := by
    apply abs_le.mpr
    constructor <;> linarith [hqa'.1,hqa'.2,hqb'.1,hqb'.2]
  have hdlo : (u+v)/2 ≤ v*qb+u*qa := by
    nlinarith [mul_le_mul_of_nonneg_left hqap hu.le,
      mul_le_mul_of_nonneg_left hqbp hv.le]
  apply (div_le_iff₀ hden).mpr
  have hdifmul := mul_le_mul_of_nonneg_left hdiff (mul_nonneg hu.le hv.le)
  rcases le_total u v with huv|hvu
  · rw [min_eq_left huv]
    have hscale := mul_le_mul_of_nonneg_left hdlo (show 0 ≤ 4*δ*u by positivity)
    nlinarith [mul_nonneg hδ (sq_nonneg u)]
  · rw [min_eq_right hvu]
    have hscale := mul_le_mul_of_nonneg_left hdlo (show 0 ≤ 4*δ*v by positivity)
    nlinarith [mul_nonneg hδ (sq_nonneg v)]

theorem exists_radius_secant_error {f : ℝ → ℝ} {r t : ℝ}
    (hf : DifferentiableAt ℝ f r) (hr : f r = 0) (hd : deriv f r ≠ 0)
    (ht : 0 < t) :
    ∃ ε > 0, ∀ a b : ℝ, a < r → r < b → |a-r| < ε → |b-r| < ε →
      |secant a b (f a) (f b)-r| ≤ t * min |a-r| |b-r| := by
  let δ := min (t/4) (1/2:ℝ)
  have hδ : 0 < δ := lt_min (by positivity) (by norm_num)
  obtain ⟨ε,hε,hq⟩ := exists_radius_normalizedSlope hf hd hδ
  refine ⟨ε,hε,?_⟩
  intro a b ha hb hεa hεb
  rw [function_eq_normalizedSlope hr hd (x:=a), function_eq_normalizedSlope hr hd (x:=b)]
  have hh := secant_error_of_normalized_values ha hb hd hδ.le
    (min_le_right (t/4) (1/2:ℝ)) (hq a hεa).le (hq b hεb).le
  refine hh.trans ?_
  apply mul_le_mul_of_nonneg_right _ (le_min (abs_nonneg _) (abs_nonneg _))
  have hδt : δ ≤ t/4 := min_le_left _ _
  linarith

/-- A factor near one preserves the normalized slope estimate with an explicit
constant. This is used for the AB-modified second step of the pair. -/
theorem normalized_product_bound {μ q δ : ℝ}
    (_hδ : 0 ≤ δ) (hδ1 : δ ≤ 1)
    (hμ : |μ-1| ≤ δ) (hq : |q-1| ≤ δ) :
    |μ*q-1| ≤ 3*δ := by
  have hμ2 : |μ| ≤ 2 := by
    apply abs_le.mpr
    rcases abs_le.mp hμ with ⟨hl,hu⟩
    constructor <;> linarith
  have he : μ*q-1 = μ*(q-1)+(μ-1) := by ring
  rw [he]
  calc
    |μ*(q-1)+(μ-1)| ≤ |μ*(q-1)|+|μ-1| := abs_add_le _ _
    _ = |μ| * |q-1|+|μ-1| := by rw [abs_mul]
    _ ≤ 2*δ+δ := add_le_add (mul_le_mul hμ2 hq (abs_nonneg _) (by norm_num)) hμ
    _ = 3*δ := by ring

/-- Uniform secant estimate allowing either endpoint ordinate to be scaled by
a factor sufficiently close to one. The permitted factor size is returned
explicitly for the pair's AB-ratio bound. -/
theorem exists_radius_scaled_secant_error {f : ℝ → ℝ} {r t : ℝ}
    (hf : DifferentiableAt ℝ f r) (hr : f r = 0) (hd : deriv f r ≠ 0)
    (ht : 0 < t) :
    ∃ δ > 0, ∃ ε > 0, ∀ a b μa μb : ℝ,
      a < r → r < b → |a-r| < ε → |b-r| < ε →
      |μa-1| ≤ δ → |μb-1| ≤ δ →
      |secant a b (μa*f a) (μb*f b)-r| ≤ t * min |a-r| |b-r| := by
  let δ := min (t/12) (1/6:ℝ)
  have hδ : 0 < δ := lt_min (by positivity) (by norm_num)
  have hδ6 : δ ≤ (1/6:ℝ) := min_le_right _ _
  obtain ⟨ε,hε,hq⟩ := exists_radius_normalizedSlope hf hd hδ
  refine ⟨δ,hδ,ε,hε,?_⟩
  intro a b μa μb ha hb hεa hεb hμa hμb
  have hqa := normalized_product_bound hδ.le (by linarith) hμa (hq a hεa).le
  have hqb := normalized_product_bound hδ.le (by linarith) hμb (hq b hεb).le
  have hfa : μa*f a = deriv f r*(a-r)*(μa*normalizedSlope f r a) := by
    rw [function_eq_normalizedSlope hr hd (x:=a)]
    ring
  have hfb : μb*f b = deriv f r*(b-r)*(μb*normalizedSlope f r b) := by
    rw [function_eq_normalizedSlope hr hd (x:=b)]
    ring
  rw [hfa,hfb]
  have hh := secant_error_of_normalized_values ha hb hd (by positivity : 0 ≤ 3*δ)
    (by linarith) hqa hqb
  refine hh.trans ?_
  apply mul_le_mul_of_nonneg_right _ (le_min (abs_nonneg _) (abs_nonneg _))
  have hδt : δ ≤ t/12 := min_le_left _ _
  linarith

end
end ModAB.Rational.SecantLocal

end RModABBundle37

/- Original module: ModAB.RationalCycleSafeguards -/
section RModABBundle38

/-! Uniform safeguards for both steps of the negative-coefficient local pair.
The hypotheses on raw error and denominator are supplied by the interpolation
analysis. Secant, curvature, closing and residual decisions are derived here.
The same neighbourhood works at every smaller radius. -/
namespace ModAB.Rational.CycleSafeguards
open Set Filter Topology
open SecantLocal Curvature ModAB.RationalLocalSafeguards ModAB.Convergence
noncomputable section

theorem residual_ratio_bound {f : ℝ → ℝ} {r x y D q : ℝ}
    (hD : 0 < D) (hx : x ≠ r) (hq : 0 ≤ q)
    (hflower : D/2*|x-r| ≤ |f x|)
    (hfyupper : |f y| ≤ (2*D)*|y-r|)
    (hstep : |y-r| ≤ q*|x-r|) :
    |f y / f x| ≤ 4*q := by
  have hxp : 0 < |x-r| := abs_pos.mpr (sub_ne_zero.mpr hx)
  have hfx : 0 < |f x| := (mul_pos (by linarith) hxp).trans_le hflower
  rw [abs_div]
  apply (div_le_iff₀ hfx).mpr
  have hy := mul_le_mul_of_nonneg_left hstep (show 0 ≤ 2*D by positivity)
  have hl := mul_le_mul_of_nonneg_left hflower (show 0 ≤ 4*q by positivity)
  nlinarith

/-- The literal AB factor selects its ordinary branch and stays near one.
This does not assume either conclusion as a scaling hypothesis. -/
theorem abFactor_near_one {f : ℝ → ℝ} {r x y D q : ℝ}
    (hD : 0 < D) (hx : x ≠ r) (hq : 0 ≤ q) (hsmall : 4*q < 1)
    (hflower : D/2*|x-r| ≤ |f x|)
    (hfyupper : |f y| ≤ (2*D)*|y-r|)
    (hstep : |y-r| ≤ q*|x-r|) :
    abFactor (f y) (f x) = 1-f y/f x ∧
      |abFactor (f y) (f x)-1| ≤ 4*q := by
  have hbound := residual_ratio_bound hD hx hq hflower hfyupper hstep
  have hh : 0 < 1-f y/f x := by linarith [le_abs_self (f y/f x)]
  have he : abFactor (f y) (f x) = 1-f y/f x := by simp [abFactor,hh]
  refine ⟨he,?_⟩
  rw [he]
  simpa only [sub_sub_cancel_left,abs_neg] using hbound

structure ScaledNeighborhood (f : ℝ → ℝ) (r C ε δ : ℝ) : Prop where
  radius_pos : 0 < ε
  factor_pos : 0 < δ
  coefficient_nonneg : 0 ≤ C
  deriv_nonzero : deriv f r ≠ 0
  small : C*ε^2 < (1/1000:ℝ)
  factor_small : 4*C*ε^2 ≤ δ
  secant : ∀ a b μa μb : ℝ, a < r → r < b → |a-r| ≤ ε → |b-r| ≤ ε →
    |μa-1| ≤ δ → |μb-1| ≤ δ →
    |ModAB.Rational.secant a b (μa*f a) (μb*f b)-r| ≤
      (1/10:ℝ)*min |a-r| |b-r|
  compare : ∀ x, |x-r| ≤ ε →
    |deriv f r|/2*|x-r| ≤ |f x| ∧ |f x| ≤ (2*|deriv f r|)*|x-r|
  curvature : ∀ l a b h : ℝ, |l-r| ≤ ε → |h-r| ≤ ε →
    l < a → a < b → b < h → curvatureChangesSign f l a b h = false

theorem exists_scaled_neighborhood {f : ℝ → ℝ} {r C : ℝ}
    (hf : ContDiffAt ℝ 2 f r) (hr : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hC : 0 ≤ C) :
    ∃ δ > 0, ∃ ε > 0, ScaledNeighborhood f r C ε δ := by
  obtain ⟨δ,hδ,es,hes,hs⟩ := exists_radius_scaled_secant_error
    (hf.differentiableAt (by norm_num)) hr hd (by norm_num : (0:ℝ)<1/10)
  obtain ⟨ef,hef,hfbound⟩ := exists_radius_simple_root_bounds
    (hf.differentiableAt (by norm_num)) hr hd
  obtain ⟨ec,hec,hcurv⟩ := exists_radius_no_curvature_change hf h₂
  let Q := min (1/1000:ℝ) (δ/4)
  have hQ : 0 < Q := lt_min (by norm_num) (by positivity)
  have hc : ContinuousAt (fun x : ℝ => C*x^2) 0 := by fun_prop
  obtain ⟨eq,heq,hq⟩ := Metric.continuousAt_iff.mp hc Q hQ
  let ε := min (min es ef) (min ec eq) / 2
  have hε : 0 < ε := by dsimp [ε]; positivity
  have hεs : ε < es := by
    have h := min_le_left (min es ef) (min ec eq)
    have h' := min_le_left es ef
    dsimp [ε] at *
    linarith
  have hεf : ε < ef := by
    have h := min_le_left (min es ef) (min ec eq)
    have h' := min_le_right es ef
    dsimp [ε] at *
    linarith
  have hεc : ε < ec := by
    have h := min_le_right (min es ef) (min ec eq)
    have h' := min_le_left ec eq
    dsimp [ε] at *
    linarith
  have hεq : ε < eq := by
    have h := min_le_right (min es ef) (min ec eq)
    have h' := min_le_right ec eq
    dsimp [ε] at *
    linarith
  have hsmall : C*ε^2 < Q := by
    have h := hq (x:=ε) (by simpa [Real.dist_eq,abs_of_pos hε] using hεq)
    simpa [Real.dist_eq,abs_of_nonneg hC,abs_of_nonneg (mul_nonneg hC (sq_nonneg ε))] using h
  refine ⟨δ,hδ,ε,hε,⟨hε,hδ,hC,hd,hsmall.trans_le (min_le_left _ _),?_,?_,?_,?_⟩⟩
  · have hδq := hsmall.trans_le (min_le_right (1/1000:ℝ) (δ/4))
    linarith
  · intro a b μa μb ha hb hna hnb hμa hμb
    exact hs a b μa μb ha hb (hna.trans_lt hεs) (hnb.trans_lt hεs) hμa hμb
  · intro x hx
    exact (hfbound x (hx.trans_lt hεf)).2
  · intro l a b h hnl hnh hla hab hbh
    exact hcurv l a b h (by linarith [(abs_le.mp hnl).1]) hla hab hbh
      (by linarith [(abs_le.mp hnh).2])

theorem ScaledNeighborhood.shrink {f : ℝ → ℝ} {r C ε ε' δ : ℝ}
    (N : ScaledNeighborhood f r C ε δ) (hε' : 0 < ε') (hε : ε' ≤ ε) :
    ScaledNeighborhood f r C ε' δ := by
  have hq : C*ε'^2 ≤ C*ε^2 :=
    mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hε'.le hε 2) N.coefficient_nonneg
  refine ⟨hε',N.factor_pos,N.coefficient_nonneg,N.deriv_nonzero,
    hq.trans_lt N.small,by nlinarith [N.factor_small],?_,?_,?_⟩
  · intro a b μa μb ha hb hna hnb hμa hμb
    exact N.secant a b μa μb ha hb (hna.trans hε) (hnb.trans hε) hμa hμb
  · exact fun x hx => N.compare x (hx.trans hε)
  · exact fun l a b h hl hh hla hab hbh =>
      N.curvature l a b h (hl.trans hε) (hh.trans hε) hla hab hbh

/-- Uniform actual candidate, zero closing and continuation consequences.
The ordinates may both be scaled, and the selector uses only true values. -/
theorem ScaledNeighborhood.guards {f : ℝ → ℝ} {r C ε δ : ℝ}
    (N : ScaledNeighborhood f r C ε δ)
    {l a b h xp fp μa μb : ℝ} (side : EndpointSide)
    (ha : a < r) (hb : r < b) (hla : l < a) (hbh : b < h)
    (hnl : |l-r| ≤ ε) (hnh : |h-r| ≤ ε)
    (hna : |a-r| ≤ ε) (hnb : |b-r| ≤ ε)
    (hμa : |μa-1| ≤ δ) (hμb : |μb-1| ≤ δ)
    (hden : denominator a b (μa*f a) (μb*f b) xp fp ≠ 0)
    (hraw : |rawPoint a b (μa*f a) (μb*f b) xp fp-r| ≤
      (C*ε^2)*min |a-r| |b-r|) :
    let raw := rawPoint a b (μa*f a) (μb*f b) xp fp
    let candidate' := candidate a b (μa*f a) (μb*f b) xp fp
      (!curvatureChangesSign f l a b h)
    candidate' = raw ∧ closeCandidate a b 0 candidate' side = raw ∧
      42*|f raw| < min |f a| |f b| := by
  dsimp only
  have hq : C*ε^2 < 1 := N.small.trans (by norm_num)
  have hc := N.curvature l a b h hnl hnh hla (ha.trans hb) hbh
  have hs := N.secant a b μa μb ha hb hna hnb hμa hμb
  have hclamp : C*ε^2+(1/10:ℝ) ≤ (7/10:ℝ)*(1-(1/10:ℝ)) := by linarith [N.small]
  have he : candidate a b (μa*f a) (μb*f b) xp fp
      (!curvatureChangesSign f l a b h) = rawPoint a b (μa*f a) (μb*f b) xp fp := by
    unfold candidate
    rw [ite_eq_left ⟨by simp [hc],hden⟩]
    exact limited_eq_raw_of_local_estimates ha.le hb.le hraw hs hclamp
  rw [he]
  have hi := interior_of_local_error ha hb hq hraw
  have hnraw : |rawPoint a b (μa*f a) (μb*f b) xp fp-r| ≤ ε := by
    rcases abs_le.mp hna with ⟨hal,har⟩
    rcases abs_le.mp hnb with ⟨hbl,hbr⟩
    exact abs_le.mpr ⟨by linarith [hi.1],by linarith [hi.2]⟩
  have hm : 0 < min |a-r| |b-r| := by
    apply lt_min <;> apply abs_pos.mpr <;> linarith
  have hD : 0 < |deriv f r| := abs_pos.mpr N.deriv_nonzero
  have hl : |deriv f r|/2*min |a-r| |b-r| ≤ min |f a| |f b| := by
    rw [mul_min_of_nonneg _ _ (by positivity)]
    exact min_le_min (N.compare a hna).1 (N.compare b hnb).1
  have hres : 42*(2*|deriv f r|)*(C*ε^2) < |deriv f r|/2 := by
    nlinarith [mul_lt_mul_of_pos_left N.small hD]
  exact ⟨rfl,closeCandidate_zero ⟨hi.1.le,hi.2.le⟩ side,
    forty_two_residual_of_local_bounds (by positivity) hm hraw
      (N.compare _ hnraw).2 hl hres⟩

theorem ScaledNeighborhood.ab_factor {f : ℝ → ℝ} {r C ε δ x y : ℝ}
    (N : ScaledNeighborhood f r C ε δ)
    (hx : x ≠ r) (hnx : |x-r| ≤ ε) (hny : |y-r| ≤ ε)
    (hstep : |y-r| ≤ (C*ε^2)*|x-r|) :
    abFactor (f y) (f x) = 1-f y/f x ∧
      |abFactor (f y) (f x)-1| ≤ δ := by
  have hsmall : 4*(C*ε^2) < 1 := by linarith [N.small]
  have hq := mul_nonneg N.coefficient_nonneg (sq_nonneg ε)
  have hh := abFactor_near_one (abs_pos.mpr N.deriv_nonzero) hx hq hsmall
    (N.compare x hnx).1 (N.compare y hny).2 hstep
  exact ⟨hh.1,hh.2.trans (by nlinarith [N.factor_small])⟩

end
end ModAB.Rational.CycleSafeguards

end RModABBundle38

/- Original module: ModAB.RationalCycleNeighborhood -/
section RModABBundle39

/-! One common neighbourhood for the exact negative-K pair and every source
safeguard. All fields are derived from C³, a simple root and nonzero curvature.
In particular neither future coefficient bounds nor guard decisions are input.
-/
namespace ModAB.Rational.CycleNeighborhood
open Set Filter
open scoped Topology
open Regularity CycleUniform CycleState CycleSafeguards SecantLocal
noncomputable section

structure PairNeighborhood (f : ℝ → ℝ) (r θ ε C δ : ℝ) : Prop where
  bounds : Bounds (regularizedReciprocal f r) r θ ε C
  guards : ScaledNeighborhood f r C ε δ
  same_sign : ∀ x y, |x-r| ≤ ε → |y-r| ≤ ε →
    (0 < f x*f y ↔ 0 < (x-r)*(y-r))
  raw_pair : ∀ a s d, CycleEntry.Cone θ ε a s d →
    let g := regularizedReciprocal f r
    let K := AnalyticBridge.errorConstant f r
    let b := firstError g r a s d
    let μ := 1-f (r+b)/f (r+a)
    firstCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
    secondCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
    denominator (r+a) (r+s) (f (r+a)) (f (r+s)) (r+d) (f (r+d)) ≠ 0 ∧
    rawPoint (r+a) (r+s) (f (r+a)) (f (r+s)) (r+d) (f (r+d))-r = b ∧
    μ ∈ Ioo (0:ℝ) 1 ∧
    denominator (r+b) (r+s) (f (r+b)) (μ*f (r+s)) (r+a) (f (r+a)) ≠ 0 ∧
    rawPoint (r+b) (r+s) (f (r+b)) (μ*f (r+s)) (r+a) (f (r+a))-r =
      secondCoefficient g r a s d*(a*b*d) ∧
    b ≠ 0 ∧ b ≠ a ∧ b ≠ s ∧ 0 < a*b

/-- Simultaneous constants for the interpolation coefficients, local cone,
actual AB update, scaled secant, curvature and true-residual tests. -/
theorem exists_pair_neighborhood {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : AnalyticBridge.errorConstant f r < 0) :
    ∃ θ ε C δ, PairNeighborhood f r θ ε C δ := by
  obtain ⟨θ,ε₀,hθ0,hθ1,hε₀,hraw⟩ :=
    CycleUniformBridge.smooth_exists_uniform_raw_pair hf hz hd hK
  let C := -2*AnalyticBridge.errorConstant f r
  have hC : 0 < C := by dsimp [C]; linarith
  obtain ⟨δ,hδ,ε₁,hε₁,hguards⟩ := exists_scaled_neighborhood
    (hf.of_le (by norm_num)) hz hd h₂ hC.le
  obtain ⟨ε₂,hε₂,hcompare⟩ := exists_radius_simple_root_bounds
    (hf.differentiableAt (by norm_num)) hz hd
  let ε := min ε₀ (min ε₁ (min (ε₂/2) (min 1 (θ/C))))
  have hε : 0 < ε := by dsimp [ε]; positivity
  have e₀ : ε ≤ ε₀ := by simp [ε,min_le_iff]
  have e₁ : ε ≤ ε₁ := by simp [ε,min_le_iff]
  have e₂ : ε < ε₂ := by
    have hh : ε ≤ ε₂/2 := by simp [ε,min_le_iff]
    linarith
  have eone : ε ≤ 1 := by simp [ε,min_le_iff]
  have eC : ε ≤ θ/C := by simp [ε,min_le_iff]
  have hsmall : C*ε^2 ≤ θ := by
    have hCe : ε*C ≤ θ := (le_div_iff₀ hC).mp eC
    have hs := mul_nonneg (mul_nonneg hC.le hε.le) (sub_nonneg.mpr eone)
    nlinarith
  have hrawε := fun a s d (hc : CycleEntry.Cone θ ε a s d) =>
    hraw a s d (hc.enlarge_radius e₀)
  refine ⟨θ,ε,C,δ,{
    bounds := ⟨hθ0,hθ1,hε,hC,hsmall,?_⟩
    guards := hguards.shrink hε e₁
    same_sign := ?_
    raw_pair := hrawε
  }⟩
  · intro a s d hc
    have hu := (hrawε a s d hc).1
    have hv := (hrawε a s d hc).2.1
    have hU : firstCoefficient (regularizedReciprocal f r) r a s d < 0 := by linarith [hu.2]
    have hV : secondCoefficient (regularizedReciprocal f r) r a s d < 0 := by linarith [hv.2]
    refine ⟨hU,hV,?_,?_⟩
    · rw [abs_of_neg hU]
      dsimp [C]
      linarith [hu.1]
    · rw [abs_of_neg hV]
      dsimp [C]
      linarith [hv.1]
  · intro x y hx hy
    exact (sign_product_of_normalizedSlope hz hd
      (hcompare x (hx.trans_lt e₂)).1 (hcompare y (hy.trans_lt e₂)).1).1

end
end ModAB.Rational.CycleNeighborhood

end RModABBundle39

/- Original module: ModAB.RationalCycleIntermediate -/
section RModABBundle40

/-!
The intermediate same-side update in a rational pair. All formulas follow
from the literal source ordinate/history update; the AB multiplier and the
selected discarded point are not assumed as a recurrence interface.
-/
namespace ModAB.Rational.CycleIntermediate
open Set Filter
open scoped Topology
open ModAB.Convergence LocalState
noncomputable section

def latestWork (s : State) : ℝ :=
  match s.recent with | .left => s.ordinates.leftWork | .right => s.ordinates.rightWork

def oppositeWork (s : State) : ℝ :=
  match s.recent with | .left => s.ordinates.rightWork | .right => s.ordinates.leftWork

theorem same_side_update_left {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .left) (hx : x ≠ s.a)
    (hbranch : 0 < f s.a * f x) :
    update f s x =
      ⟨x,s.b,s.a,s.rightHistory,
        ⟨f x,f s.b,f x,abFactor (f x) (f s.a)*f s.b,some .left⟩,
        .left,s.scale,s.clock+1⟩ := by
  simp [update,hbranch,hx,h.fresh,hs,freshOrdinates,recordUpdate]

theorem same_side_update_right {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .right) (hx : x ≠ s.b)
    (hbranch : ¬ 0 < f s.a * f x) :
    update f s x =
      ⟨s.a,x,s.leftHistory,s.b,
        ⟨f s.a,f x,abFactor (f x) (f s.b)*f s.a,f x,some .right⟩,
        .right,s.scale,s.clock+1⟩ := by
  simp [update,hbranch,hx,h.fresh,hs,freshOrdinates,recordUpdate]

/-- Orientation-independent source update. In particular, the model's next
third node is exactly the previous latest endpoint with its true ordinate. -/
theorem same_side_update {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    (update f s x).recent = s.recent ∧
    latest (update f s x) = x ∧
    opposite (update f s x) = opposite s ∧
    selectedHistory (update f s x) = latest s ∧
    latestWork (update f s x) = f x ∧
    oppositeWork (update f s x) = abFactor (f x) (f (latest s))*f (opposite s) ∧
    (update f s x).ordinates.previous = some s.recent ∧
    (update f s x).scale = s.scale ∧
    (update f s x).clock = s.clock+1 := by
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame
    have hbranch := (hsgn _ _ h.near_left hxn).mpr hsame
    rw [same_side_update_left h hs (ne_of_gt hxi.1) hbranch]
    simp [latest,opposite,selectedHistory,latestWork,oppositeWork,hs]
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame
    have hrx : r < x := by nlinarith [h.right]
    have hbranch : ¬ 0 < f s.a*f x := by
      rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg
        (by linarith [h.left]) (by linarith))
    rw [same_side_update_right h hs (ne_of_lt hxi.2) hbranch]
    simp [latest,opposite,selectedHistory,latestWork,oppositeWork,hs]

theorem abFactor_eq_one_sub {new old : ℝ} (h : new/old < 1) :
    abFactor new old = 1-new/old := by
  simp [abFactor,show 0 < 1-new/old by linarith]

/-- For the local pair the executed AB branch is the untruncated factor
1-f(new)/f(old), rather than its fallback value 1/2. -/
theorem same_side_working_ordinates {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r)))
    (hratio : f x/f (latest s) < 1) :
    latestWork (update f s x) = f x ∧
    oppositeWork (update f s x) = (1-f x/f (latest s))*f (opposite s) ∧
    selectedHistory (update f s x) = latest s := by
  have hu := same_side_update h hxi hsame hsgn
  refine ⟨hu.2.2.2.2.1,?_,hu.2.2.2.1⟩
  rw [hu.2.2.2.2.2.1,abFactor_eq_one_sub hratio]


/-- The same-side intermediate state keeps the strict bracket and both
outside histories in the local neighbourhood. Freshness is deliberately
not included: its retained working ordinate has just been rescaled. -/
theorem same_side_geometry {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x : ℝ}
    (h : Valid f r ε s) (hxi : x ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    let t := update f s x
    t.a < r ∧ r < t.b ∧ t.leftHistory < t.a ∧ t.b < t.rightHistory ∧
    |t.a-r| ≤ ε ∧ |t.b-r| ≤ ε ∧
    |t.leftHistory-r| ≤ ε ∧ |t.rightHistory-r| ≤ ε := by
  dsimp only
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame
    have hxr : x < r := by nlinarith [h.left]
    have hbranch := (hsgn _ _ h.near_left hxn).mpr hsame
    rw [same_side_update_left h hs (ne_of_gt hxi.1) hbranch]
    exact ⟨hxr,h.right,hxi.1,h.right_history,hxn,h.near_right,
      h.near_left,h.near_right_history⟩
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame
    have hrx : r < x := by nlinarith [h.right]
    have hbranch : ¬ 0 < f s.a*f x := by
      rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg
        (by linarith [h.left]) (by linarith))
    rw [same_side_update_right h hs (ne_of_lt hxi.2) hbranch]
    exact ⟨h.left,hrx,h.left_history,hxi.2,h.near_left,hxn,
      h.near_left_history,h.near_right⟩

end
end ModAB.Rational.CycleIntermediate

end RModABBundle40

/- Original module: ModAB.RationalLocalComplete -/
section RModABBundle41

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

end RModABBundle41

/- Original module: ModAB.RationalCycleStateGuards -/
section RModABBundle42

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

end RModABBundle42

/- Original module: ModAB.RationalCycleStateUpdate -/
section RModABBundle43

/-!
Two successive genuine source updates in the negative-coefficient cycle.
The first update applies the actual Anderson--Björck factor to the retained
ordinate. The crossing update restores true ordinates and the correct two
discarded-side histories; no state-invariance premise is used.
-/
namespace ModAB.Rational.CycleStateUpdate
open Set Filter
open scoped Topology
open ModAB.Convergence LocalState
noncomputable section

theorem pair_update_left {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .left)
    (hax : s.a < x) (hxr : x < r) (hry : r < y) (hyb : y < s.b)
    (hfirst : 0 < f s.a * f x) (hsecond : ¬ 0 < f x * f y) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxe : x ≠ s.a := ne_of_gt hax
  have hye : y ≠ s.b := ne_of_lt hyb
  have hxnear := interior_near h.left h.right h.near_left h.near_right
    (show x ∈ Ioo s.a s.b from ⟨hax,by linarith [h.right]⟩)
  have hynear := interior_near h.left h.right h.near_left h.near_right
    (show y ∈ Ioo s.a s.b from ⟨by linarith [h.left],hyb⟩)
  have hu : update f (update f s x) y =
      ⟨x,y,s.a,s.b,freshOrdinates f x y .right,.right,s.scale,s.clock+2⟩ := by
    simp [update,hfirst,hsecond,hxe,hye,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨hxr,hry,hax,hyb,hxnear,hynear,h.near_left,h.near_right,rfl⟩,
    rfl,rfl,?_⟩
  simp [opposite,hs,selectedHistory]

theorem pair_update_right {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s) (hs : s.recent = .right)
    (hrx : r < x) (hxb : x < s.b) (hay : s.a < y) (hyr : y < r)
    (hfirst : ¬ 0 < f s.a * f x) (hsecond : 0 < f s.a * f y) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxe : x ≠ s.b := ne_of_lt hxb
  have hye : y ≠ s.a := ne_of_gt hay
  have hxnear := interior_near h.left h.right h.near_left h.near_right
    (show x ∈ Ioo s.a s.b from ⟨by linarith [h.left],hxb⟩)
  have hynear := interior_near h.left h.right h.near_left h.near_right
    (show y ∈ Ioo s.a s.b from ⟨hay,by linarith [h.right]⟩)
  have hu : update f (update f s x) y =
      ⟨y,x,s.a,s.b,freshOrdinates f y x .left,.left,s.scale,s.clock+2⟩ := by
    simp [update,hfirst,hsecond,hxe,hye,h.fresh,hs,freshOrdinates,recordUpdate]
  rw [hu]
  refine ⟨⟨hyr,hrx,hay,hxb,hynear,hxnear,h.near_left,h.near_right,rfl⟩,
    rfl,rfl,?_⟩
  simp [opposite,hs,selectedHistory]

/-- The pair update needs only strict interiority and the two proved step
signs. True signs select the source branches. Both history records are
updated by the actual guarded assignments. -/
theorem pair_update_preserves {f : ℝ → ℝ} {r ε : ℝ} {s : State} {x y : ℝ}
    (h : Valid f r ε s)
    (hxi : x ∈ Ioo s.a s.b) (hyi : y ∈ Ioo s.a s.b)
    (hsame : 0 < (latest s-r)*(x-r))
    (hcross : (latest s-r)*(y-r) < 0)
    (hsgn : ∀ u v, |u-r| ≤ ε → |v-r| ≤ ε →
      (0 < f u*f v ↔ 0 < (u-r)*(v-r))) :
    Valid f r ε (update f (update f s x) y) ∧
    latest (update f (update f s x) y) = y ∧
    opposite (update f (update f s x) y) = x ∧
    selectedHistory (update f (update f s x) y) = opposite s := by
  have hxn := interior_near h.left h.right h.near_left h.near_right hxi
  have hyn := interior_near h.left h.right h.near_left h.near_right hyi
  cases hs : s.recent with
  | left =>
    have hl : latest s = s.a := by simp [latest,hs]
    rw [hl] at hsame hcross
    have hxr : x < r := by nlinarith [h.left]
    have hry : r < y := by nlinarith [h.left]
    apply pair_update_left h hs hxi.1 hxr hry hyi.2
    · exact (hsgn _ _ h.near_left hxn).mpr hsame
    · rw [hsgn _ _ hxn hyn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg (by linarith) (by linarith))
  | right =>
    have hl : latest s = s.b := by simp [latest,hs]
    rw [hl] at hsame hcross
    have hrx : r < x := by nlinarith [h.right]
    have hyr : y < r := by nlinarith [h.right]
    apply pair_update_right h hs hrx hxi.2 hyi.1 hyr
    · rw [hsgn _ _ h.near_left hxn]
      exact not_lt.mpr (mul_nonpos_of_nonpos_of_nonneg (by linarith [h.left]) (by linarith))
    · exact (hsgn _ _ h.near_left hyn).mpr
        (mul_pos_of_neg_of_neg (sub_neg.mpr h.left) (sub_neg.mpr hyr))

end
end ModAB.Rational.CycleStateUpdate

end RModABBundle43

/- Original module: ModAB.RationalCycleTransport -/
section RModABBundle44

/-!
Transport from a pointwise two-step source law to the entire safeguarded
trajectory. The law is about one local input state; invariance, iteration
identities, convergence and Q-order of the source trajectory are derived.
-/
namespace ModAB.Rational.CycleTransport
open Set Filter
open scoped Topology
open LocalState CycleState CycleUniform
noncomputable section

def view (r : ℝ) (s : State) : PairState :=
  ⟨latest s-r, opposite s-r, selectedHistory s-r⟩

def pairStep (f : ℝ → ℝ) (allow : State → Bool) (s : State) : State :=
  step f allow (step f allow s)

structure PairLaw (f g : ℝ → ℝ) (r θ ε : ℝ) (allow : State → Bool) : Prop where
  advance : ∀ s, Valid f r ε s →
    CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d →
    Valid f r ε (pairStep f allow s) ∧
    view r (pairStep f allow s) = nextPair g r (view r s) ∧
    latest (step f allow s)-r = (nextPair g r (view r s)).s ∧
    ∀ B, nextMode B f allow s = .interpolation s.scale (s.clock+1) ∧
      nextMode B f allow (step f allow s) =
        .interpolation (step f allow s).scale ((step f allow s).clock+1)

theorem orbit_even_succ (f : ℝ → ℝ) (allow : State → Bool) (s : State) (n : ℕ) :
    Alternating.orbit f allow s (2*(n+1)) =
      pairStep f allow (Alternating.orbit f allow s (2*n)) := by
  rw [show 2*(n+1) = (2*n+1)+1 by omega]
  rfl

theorem even_orbit {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d) :
    ∀ n, Valid f r ε (Alternating.orbit f allow s (2*n)) ∧
      view r (Alternating.orbit f allow s (2*n)) = CycleState.orbit g r (view r s) n := by
  intro n
  induction n with
  | zero => exact ⟨h,rfl⟩
  | succ n ih =>
    have hc : CycleEntry.Cone θ ε
        (view r (Alternating.orbit f allow s (2*n))).a
        (view r (Alternating.orbit f allow s (2*n))).s
        (view r (Alternating.orbit f allow s (2*n))).d := by
      rw [ih.2]
      exact orbit_cone N hp n
    have hh := L.advance _ ih.1 hc
    rw [orbit_even_succ]
    refine ⟨hh.1,?_⟩
    rw [hh.2.1,ih.2]
    rfl

theorem odd_latest {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d)
    (n : ℕ) :
    latest (Alternating.orbit f allow s (2*n+1))-r =
      (CycleState.orbit g r (view r s) (n+1)).s := by
  have he := even_orbit N L h hp n
  have hc : CycleEntry.Cone θ ε
      (view r (Alternating.orbit f allow s (2*n))).a
      (view r (Alternating.orbit f allow s (2*n))).s
      (view r (Alternating.orbit f allow s (2*n))).d := by
    rw [he.2]
    exact orbit_cone N hp n
  have hh := (L.advance _ he.1 hc).2.2.1
  rw [Alternating.orbit_succ,hh,he.2]
  rfl

/-- Both actual continuation decisions pass at every pair, even though the
odd state has a modified retained ordinate and need not be fresh. -/
theorem continue_pair {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d)
    (B : ℝ) (n : ℕ) :
    nextMode B f allow (Alternating.orbit f allow s (2*n)) =
      .interpolation (Alternating.orbit f allow s (2*n)).scale
        ((Alternating.orbit f allow s (2*n)).clock+1) ∧
    nextMode B f allow (Alternating.orbit f allow s (2*n+1)) =
      .interpolation (Alternating.orbit f allow s (2*n+1)).scale
        ((Alternating.orbit f allow s (2*n+1)).clock+1) := by
  have he := even_orbit N L h hp n
  have hc : CycleEntry.Cone θ ε
      (view r (Alternating.orbit f allow s (2*n))).a
      (view r (Alternating.orbit f allow s (2*n))).s
      (view r (Alternating.orbit f allow s (2*n))).d := by
    rw [he.2]
    exact orbit_cone N hp n
  exact (L.advance _ he.1 hc).2.2.2 B

/-- The actual evaluated subsequences have the proved orders. There are no
error-convergence or coefficient-limit assumptions in this transport. -/
theorem source_hasQOrders {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d) :
    ModAB.RationalOrder.HasQOrder
      (fun n => latest (Alternating.orbit f allow s (2*n))-r) 3 |localK g r| ∧
    ModAB.RationalOrder.HasQOrder
      (fun n => latest (Alternating.orbit f allow s (2*n+1))-r) 3 |localK g r| ∧
    Tendsto (fun n => |latest (Alternating.orbit f allow s (2*n+1))-r| /
      |latest (Alternating.orbit f allow s (2*n))-r| ^ (9/5:ℝ)) atTop
      (𝓝 (|localK g r| ^ (2/5:ℝ))) ∧
    Tendsto (fun n => |latest (Alternating.orbit f allow s (2*(n+1)))-r| /
      |latest (Alternating.orbit f allow s (2*n+1))-r| ^ (5/3:ℝ)) atTop
      (𝓝 (|localK g r| ^ (1/3:ℝ))) := by
  have hev (n : ℕ) : latest (Alternating.orbit f allow s (2*n))-r =
      (CycleState.orbit g r (view r s) n).a := by
    exact congrArg PairState.a (even_orbit N L h hp n).2
  have hod := odd_latest N L h hp
  simp_rw [hev,hod]
  exact orbit_hasQOrders hg hgr hK N hp

end
end ModAB.Rational.CycleTransport

end RModABBundle44

/- Original module: ModAB.RationalCycleComplete -/
section RModABBundle45

/-! Complete negative-coefficient local theorem for the guarded source-state
iteration. The raw coefficient, AB factor, both accepted trials and both
continuation decisions are derived from local C³ data. -/
namespace ModAB.Rational.CycleComplete
open Set Filter
open scoped Topology
open LocalState LocalComplete CycleUniform CycleNeighborhood CycleStateGuards
open CycleIntermediate CycleTransport CycleStateUpdate
noncomputable section

theorem endpoint_min (s : State) (r : ℝ) :
    min |latest s-r| |opposite s-r| = min |s.a-r| |s.b-r| := by
  cases hs : s.recent <;> simp [latest,opposite,hs,min_comm]

theorem valid_latest_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |latest s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [latest,hs]
  · exact h.near_left
  · exact h.near_right

theorem valid_opposite_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |opposite s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [opposite,hs]
  · exact h.near_right
  · exact h.near_left

theorem valid_history_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |selectedHistory s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [selectedHistory,hs]
  · exact h.near_left_history
  · exact h.near_right_history

theorem valid_working {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) :
    latestWork s = f (latest s) ∧ oppositeWork s = f (opposite s) := by
  cases hs : s.recent <;> simp [latestWork,oppositeWork,latest,opposite,hs,h.fresh,freshOrdinates]

/-- One fresh local cone state executes a genuine same-side update followed
by a genuine crossing, including the two actual safeguarded proposals. -/
theorem pair_law_of_neighborhood {f : ℝ → ℝ} {r θ ε C δ : ℝ}
    (N : PairNeighborhood f r θ ε C δ) :
    PairLaw f (Regularity.regularizedReciprocal f r) r θ ε (allowRational f) := by
  refine ⟨?_⟩
  intro st hs hc
  let g := Regularity.regularizedReciprocal f r
  let a := latest st-r
  let s := opposite st-r
  let d := selectedHistory st-r
  let b := firstError g r a s d
  let v := secondCoefficient g r a s d
  let c := v*a*b*d
  let x := r+b
  let y := r+c
  let μ := 1-f x/f (latest st)
  have hcone : CycleEntry.Cone θ ε a s d := hc
  have hh := N.raw_pair a s d hcone
  dsimp only at hh
  have hea : r+a = latest st := by dsimp [a]; ring
  have hes : r+s = opposite st := by dsimp [s]; ring
  have hed : r+d = selectedHistory st := by dsimp [d]; ring
  obtain ⟨hU,hV,hD₁,hR₁,hμ,hD₂,hR₂,hb0,hba,hbs,hab⟩ := hh
  have hD₁' : denominator (latest st) (opposite st) (f (latest st))
      (f (opposite st)) (selectedHistory st) (f (selectedHistory st)) ≠ 0 := by
    simpa only [hea,hes,hed] using hD₁
  have hR₁' : rawPoint (latest st) (opposite st) (f (latest st))
      (f (opposite st)) (selectedHistory st) (f (selectedHistory st))-r = b := by
    simpa only [hea,hes,hed] using hR₁
  have hμ' : μ ∈ Ioo (0:ℝ) 1 := by
    simpa only [hea] using hμ
  have hD₂' : denominator x (opposite st) (f x) (μ*f (opposite st))
      (latest st) (f (latest st)) ≠ 0 := by
    simpa only [hea,hes] using hD₂
  have hR₂' : rawPoint x (opposite st) (f x) (μ*f (opposite st))
      (latest st) (f (latest st))-r = c := by
    simpa only [x,μ,c,v,b,g,hea,hes,mul_assoc] using hR₂
  obtain ⟨hu,hv,huC,hvC⟩ := N.bounds.coefficients a s d hcone
  have hq0 : 0 ≤ C*ε^2 := mul_nonneg N.bounds.bound_pos.le (sq_nonneg ε)
  have hq1 : C*ε^2 < 1 := N.guards.small.trans (by norm_num)
  have haε : |a| ≤ ε := valid_latest_near hs
  have hsε : |s| ≤ ε := valid_opposite_near hs
  have hdε : |d| ≤ ε := valid_history_near hs
  have hbprod : |b| ≤ C*|a| *|s| *|d| := by
    dsimp [b,firstError]
    simp only [abs_mul]
    gcongr
  have hbmin : |b| ≤ (C*ε^2)*min |a| |s| :=
    ModAB.RationalLocalSafeguards.triple_error_bound N.bounds.bound_pos.le
      N.bounds.radius_pos.le haε hsε hdε hbprod
  have hraw₁ : |rawPoint (latest st) (opposite st) (f (latest st))
      ((1:ℝ)*f (opposite st)) (selectedHistory st) (f (selectedHistory st))-r| ≤
      (C*ε^2)*min |latest st-r| |opposite st-r| := by
    simpa only [one_mul,hR₁'] using hbmin
  have hw := valid_working hs
  have guard₁ := guards_latest N.guards (Geometry.of_valid hs) hw.1
    (by simpa using hw.2) (by norm_num : (1:ℝ) ≠ 0)
    (by simpa using N.guards.factor_pos.le)
    (by simpa using hD₁') hraw₁
  have htrial₁ : trial f (allowRational f) st = x := by
    have hh := guard₁.1
    rw [one_mul] at hh
    dsimp only [x]
    linarith only [hR₁',hh]
  have hbminOrd : |x-r| ≤ (C*ε^2)*min |st.a-r| |st.b-r| := by
    simpa only [x,add_sub_cancel_left,←endpoint_min st r] using hbmin
  have hxi := ModAB.RationalLocalSafeguards.interior_of_local_error hs.left hs.right hq1 hbminOrd
  have hxn := interior_near hs.left hs.right hs.near_left hs.near_right hxi
  have hsame : 0 < (latest st-r)*(x-r) := by
    simpa only [x,add_sub_cancel_left] using hab
  have hi := same_side_update hs hxi hsame N.same_sign
  let mid := update f st x
  have hstep₁ : step f (allowRational f) st = mid := by
    simp only [step,htrial₁]; rfl
  have hgeom := Geometry.of_same_side hs hxi hsame N.same_sign
  have hmidlatest : latest mid = x := hi.2.1
  have hmidopp : opposite mid = opposite st := hi.2.2.1
  have hmidhist : selectedHistory mid = latest st := hi.2.2.2.1
  have hbLatest : |x-r| ≤ (C*ε^2)*|latest st-r| := by
    exact (show |x-r| ≤ (C*ε^2)*min |latest st-r| |opposite st-r| by
      simpa only [x,add_sub_cancel_left] using hbmin).trans
        (mul_le_mul_of_nonneg_left (min_le_left _ _) hq0)
  have habfactor := N.guards.ab_factor (Alternating.latest_ne_root hs)
    (valid_latest_near hs) hxn hbLatest
  have hmidwork : latestWork mid = f (latest mid) := by rw [hmidlatest]; exact hi.2.2.2.2.1
  have hmidoppwork : oppositeWork mid = μ*f (opposite mid) := by
    rw [hmidopp,hi.2.2.2.2.2.1,habfactor.1]
  have hμnear : |μ-1| ≤ δ := by simpa only [habfactor.1] using habfactor.2
  have hb_le_min : |b| ≤ min |a| |s| := hbmin.trans
    (mul_le_of_le_one_left (le_min (abs_nonneg a) (abs_nonneg s)) hq1.le)
  have hb_le_s : |b| ≤ |s| := hb_le_min.trans (min_le_right _ _)
  have hcBound : |c| ≤ (C*ε^2)*|b| := by
    exact CycleEntry.product_bound N.bounds.bound_pos.le N.bounds.radius_pos.le hvC haε hdε
  have hraw₂ : |rawPoint (latest mid) (opposite mid) (f (latest mid))
      (μ*f (opposite mid)) (selectedHistory mid) (f (selectedHistory mid))-r| ≤
      (C*ε^2)*min |latest mid-r| |opposite mid-r| := by
    rw [hmidlatest,hmidopp,hmidhist,hR₂']
    have hh : |b| ≤ |opposite st-r| := hb_le_s
    simpa only [x,add_sub_cancel_left,min_eq_left hh] using hcBound
  have hden₂ : denominator (latest mid) (opposite mid) (f (latest mid))
      (μ*f (opposite mid)) (selectedHistory mid) (f (selectedHistory mid)) ≠ 0 := by
    simpa only [hmidlatest,hmidopp,hmidhist] using hD₂'
  have guard₂ := guards_latest N.guards hgeom hmidwork hmidoppwork
    hμ'.1.ne' hμnear hden₂ hraw₂
  have htrial₂ : trial f (allowRational f) mid = y := by
    have hh := guard₂.1
    rw [hmidlatest,hmidopp,hmidhist] at hh
    dsimp only [y]
    linarith only [hR₂',hh]
  have hpair : pairStep f (allowRational f) st = update f mid y := by
    unfold pairStep
    rw [hstep₁]
    change update f mid (trial f (allowRational f) mid) = _
    rw [htrial₂]
  have hcyOrd : |y-r| ≤ (C*ε^2)*min |st.a-r| |st.b-r| := by
    have hcmin := hcBound.trans (mul_le_mul_of_nonneg_left hb_le_min hq0)
    simpa only [y,add_sub_cancel_left,←endpoint_min st r] using hcmin
  have hyi := ModAB.RationalLocalSafeguards.interior_of_local_error hs.left hs.right hq1 hcyOrd
  have hcross : (latest st-r)*(y-r) < 0 := by
    have hneg := mul_neg_of_neg_of_pos (mul_neg_of_neg_of_pos hv hcone.same) hab
    change a*((r+c)-r) < 0
    dsimp only [c]
    convert hneg using 1; ring
  have hout := pair_update_preserves hs hxi hyi hsame hcross N.same_sign
  rw [hpair]
  refine ⟨hout.1,?_,?_,?_⟩
  · dsimp only [view,CycleState.nextPair]
    rw [hout.2.1,hout.2.2.1,hout.2.2.2]
    simp only [y,x,add_sub_cancel_left]
    rfl
  · rw [hstep₁,hmidlatest]
    dsimp [x,CycleState.nextPair,view,a,s,d,b,g]
    ring
  · intro B
    refine ⟨guard₁.2 B,?_⟩
    rw [hstep₁]
    exact guard₂.2 B


/-- Propositions 10--11 in operational form. Finite C³ smoothness and the
initial geometric cone choose one common neighbourhood. Every later pair
executes both actual guarded source steps, stays in the cone and has the
proved two distinct consecutive-step orders. No future-state, guard,
recurrence, coefficient-limit or convergence premise is present. -/
theorem local_pair_complete {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : AnalyticBridge.errorConstant f r < 0) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ st : State, Valid f r ε st →
        CycleEntry.Cone θ ε (view r st).a (view r st).s (view r st).d →
        (∀ n, Valid f r ε (Alternating.orbit f (allowRational f) st (2*n)) ∧
          CycleEntry.Cone θ ε
            (view r (Alternating.orbit f (allowRational f) st (2*n))).a
            (view r (Alternating.orbit f (allowRational f) st (2*n))).s
            (view r (Alternating.orbit f (allowRational f) st (2*n))).d) ∧
        (∀ B n,
          nextMode B f (allowRational f) (Alternating.orbit f (allowRational f) st (2*n)) =
            .interpolation (Alternating.orbit f (allowRational f) st (2*n)).scale
              ((Alternating.orbit f (allowRational f) st (2*n)).clock+1) ∧
          nextMode B f (allowRational f) (Alternating.orbit f (allowRational f) st (2*n+1)) =
            .interpolation (Alternating.orbit f (allowRational f) st (2*n+1)).scale
              ((Alternating.orbit f (allowRational f) st (2*n+1)).clock+1)) ∧
        ModAB.RationalOrder.HasQOrder
          (fun n => latest (Alternating.orbit f (allowRational f) st (2*n))-r)
          3 |AnalyticBridge.errorConstant f r| ∧
        ModAB.RationalOrder.HasQOrder
          (fun n => latest (Alternating.orbit f (allowRational f) st (2*n+1))-r)
          3 |AnalyticBridge.errorConstant f r| ∧
        Tendsto (fun n => |latest (Alternating.orbit f (allowRational f) st (2*n+1))-r| /
          |latest (Alternating.orbit f (allowRational f) st (2*n))-r| ^ (9/5:ℝ)) atTop
          (𝓝 (|AnalyticBridge.errorConstant f r| ^ (2/5:ℝ))) ∧
        Tendsto (fun n => |latest (Alternating.orbit f (allowRational f) st (2*(n+1)))-r| /
          |latest (Alternating.orbit f (allowRational f) st (2*n+1))-r| ^ (5/3:ℝ)) atTop
          (𝓝 (|AnalyticBridge.errorConstant f r| ^ (1/3:ℝ))) := by
  obtain ⟨θ,ε,C,δ,N⟩ := exists_pair_neighborhood hf hz hd h₂ hK
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 2 g r :=
    Regularity.regularizedReciprocal_contDiffAt_of_contDiff hf hd
  have hgr : g r ≠ 0 := by
    dsimp only [g]
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hk : localK g r = AnalyticBridge.errorConstant f r :=
    localK_regularizedReciprocal hf hz hd
  have L := pair_law_of_neighborhood N
  refine ⟨θ,ε,N.bounds.theta_pos,N.bounds.theta_lt,N.bounds.radius_pos,?_⟩
  intro st hs hc
  refine ⟨?_,?_,?_⟩
  · intro n
    have he := even_orbit N.bounds L hs hc n
    refine ⟨he.1,?_⟩
    rw [he.2]
    exact CycleState.orbit_cone N.bounds hc n
  · exact fun B n => continue_pair N.bounds L hs hc B n
  · have hh := source_hasQOrders hg hgr (by rw [hk]; exact hK.ne)
      N.bounds L hs hc
    simpa only [hk] using hh

end
end ModAB.Rational.CycleComplete

end RModABBundle45

/- Original module: ModAB.RationalOneSidedComplete -/
section RModABBundle46

/-!
The bounded-scaling premise of the earlier connected one-sided theorem is
derived here from the actual raw rational update and its no-crossing sign.
The bounding cap has a finite limit by local C¹ regularity. No convergence
rate, summability of residual ratios, or lower bound on the accumulated weight
is assumed. The hypothesis of an unmodified, convergent, nondegenerate
one-sided tail is retained, exactly as in Proposition 12.
-/
namespace ModAB.Rational.OneSidedComplete
open Filter Set
open scoped Topology
open LocalAlgebra OneSided OneSidedAssembly
noncomputable section

/-- An increasing quantity bounded on a tail is bounded everywhere by the
same constant. This removes any finite-prefix issue in the weight bound. -/
theorem monotone_bound_of_eventually {d : ℕ → ℝ} {M : ℝ}
    (hd : Monotone d) (hM : ∀ᶠ n in atTop, d n ≤ M) : ∀ n, d n ≤ M := by
  obtain ⟨N, hN⟩ := eventually_atTop.mp hM
  intro n
  exact (hd (le_max_left n N)).trans (hN (max n N) (le_max_right n N))

/-- The upper bound extracted from the no-crossing sign has a finite limit
for two coalescing fresh nodes and a fixed opposite endpoint. -/
theorem tendsto_no_cross_cap
    {g : ℝ → ℝ} {a b : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n ≠ b n)
    (hs : s ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => -((s-a n)*(s-b n)) *
      dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s) / g s)
      atTop (𝓝 (-((s-r)*(s-r))*confluentDD2 g r s/g s)) := by
  exact (((tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)).neg.mul
    (tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb)).div_const (g s)

/-- Sign-invariant version of the scaling bound: the pole value and fixed
reciprocal value need only have the same nonzero sign. -/
theorem delta_bound_of_same_side_same_sign {E G P g δ π : ℝ}
    (hE : 0 < E) (hPg : 0 < P*g) (hδ : 0 ≤ δ) (hπ : 0 < π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π)) < 0) :
    δ < -π*G/g := by
  rcases mul_pos_iff.mp hPg with hpos | hneg
  · exact delta_bound_of_same_side hE hpos.1 hpos.2 hδ hπ hsame
  · have hsame' : E*(-G+δ*(-g)/π)/(-P+E*(δ*(-g)/π)) < 0 := by
      have hnum : -G+δ*(-g)/π = -(G+δ*g/π) := by ring
      have hden : -P+E*(δ*(-g)/π) = -(P+E*(δ*g/π)) := by ring
      rw [hnum,hden,mul_neg,neg_div_neg_eq]
      exact hsame
    have h := delta_bound_of_same_side hE (neg_pos.mpr hneg.1)
      (neg_pos.mpr hneg.2) hδ hπ hsame'
    convert h using 1
    ring

/-- Exact no-crossing data give the pointwise scaling bound. The sign of the
divided difference is a consequence, not an assumption. -/
theorem delta_bound_of_barycentric_no_cross
    {f g : ℝ → ℝ} {a b s r μ y : ℝ}
    (hba : b-a ≠ 0) (hsb : s-b ≠ 0) (hsa : s-a ≠ 0)
    (har : a-r ≠ 0) (hbr : b-r ≠ 0) (hsr : s-r ≠ 0)
    (hfs : f s ≠ 0) (hμ : 0 < μ) (hμone : μ ≤ 1)
    (hquotA : (a-r)/f a = g a) (hquotB : (b-r)/f b = g b)
    (hquotS : (s-r)/f s = g s)
    (hH : dd2 a b s (1/f a) (1/f b) (1/(μ*f s)) ≠ 0)
    (hstep : y = baryRoot a b s (f a) (f b) (μ*f s))
    (hE : 0 < (a-r)*(b-r))
    (hP : 0 < poleValue a b s r (g a) (g b) (g s) * g s)
    (hπ : 0 < (s-a)*(s-b))
    (hsame : (y-r)*(s-r) < 0) :
    1/μ-1 < -((s-a)*(s-b))*dd2 a b s (g a) (g b) (g s)/g s := by
  have hδ : 0 ≤ 1/μ-1 := by
    have h : 1 ≤ 1/μ := by
      simpa only [one_div_one] using one_div_le_one_div_of_le hμ hμone
    linarith
  let E := (a-r)*(b-r)
  let G := dd2 a b s (g a) (g b) (g s)
  let P := poleValue a b s r (g a) (g b) (g s)
  let π := (s-a)*(s-b)
  let δ := 1/μ-1
  let z := E*(G+δ*g s/π)/(P+E*(δ*g s/π))
  have hid : y-r = z*(s-r) := by
    rw [hstep, baryRoot_error_eq_mixedCoefficient hba hsb hsa har hbr hsr
      hfs (ne_of_gt hμ) hH hquotA hquotB hquotS]
    dsimp [mixedCoefficient,z,E,G,P,π,δ]
    have heq : (1/μ-1)*g s*(a-r)*(b-r)/((s-a)*(s-b)) =
        ((a-r)*(b-r))*((1/μ-1)*g s/((s-a)*(s-b))) := by ring
    rw [heq]
    ring
  have hz : z < 0 := by
    rw [hid] at hsame
    by_contra hn
    have hn0 : 0 ≤ z := le_of_not_gt hn
    have hnonneg := mul_nonneg hn0 (mul_self_nonneg (s-r))
    nlinarith
  exact delta_bound_of_same_side_same_sign hE hP hδ hπ hz

/-- The actual source-oriented rawPoint tail supplies a global finite bound
on accumulated scaling. No bound on scaling is assumed. -/
theorem exists_scaling_bound_of_rawPoint_one_sided
    {f g : ℝ → ℝ} {x μ : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgrgs : 0 < g r * g s) (hs : s ≠ r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n ≠ r)
    (hadj : ∀ n, x (n+1) ≠ x n) (hxs : ∀ n, x n ≠ s)
    (hμ : ∀ n, 0 < μ n) (hμone : ∀ n, μ n ≤ 1) (hanti : Antitone μ)
    (hfs : f s ≠ 0) (hfx : ∀ n, f (x n) ≠ 0)
    (hquot : ∀ n, (x n-r)/f (x n) = g (x n))
    (hquots : (s-r)/f s = g s)
    (hside : ∀ n, (x n-r)*(s-r) < 0)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0)
    (hstep : ∀ᶠ n in atTop, x (n+2) =
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n))) :
    ∃ M : ℝ, ∀ n, 1/μ n-1 ≤ M := by
  have hxnext : Tendsto (fun n => x (n+1)) atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 1).mpr hx
  have hpole := tendsto_poleValue_fixed_endpoint_of_contDiffAt hg hadj hs hxnext hx
  have hP : ∀ᶠ n in atTop,
      0 < poleValue (x (n+1)) (x n) s r (g (x (n+1))) (g (x n)) (g s) * g s :=
    (hpole.mul_const (g s)).eventually (Ioi_mem_nhds hgrgs)
  have hπlim : Tendsto (fun n => (s-x (n+1))*(s-x n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub hxnext).mul (tendsto_const_nhds.sub hx)
  have hsr : s-r ≠ 0 := sub_ne_zero.mpr hs
  have hπ : ∀ᶠ n in atTop, 0 < (s-x (n+1))*(s-x n) :=
    hπlim.eventually (Ioi_mem_nhds (mul_self_pos.mpr hsr))
  let M := -((s-r)*(s-r))*confluentDD2 g r s/g s + 1
  have hcap := tendsto_no_cross_cap hg hadj hs hxnext hx
  have hcapM : ∀ᶠ n in atTop, -((s-x (n+1))*(s-x n))*
      dd2 (x (n+1)) (x n) s (g (x (n+1))) (g (x n)) (g s) / g s < M :=
    hcap.eventually (Iio_mem_nhds (show
      -((s-r)*(s-r))*confluentDD2 g r s/g s < M by dsimp [M]; linarith))
  refine ⟨M, monotone_bound_of_eventually
    (delta_monotone_of_weight_antitone hμ hanti) ?_⟩
  filter_upwards [hP,hπ,hcapM,hD,hstep] with n hPn hπn hcapn hDn hstepn
  have hbridge := rawPoint_barycentric_bridge (sub_ne_zero.mpr (hadj n).symm)
    (sub_ne_zero.mpr (hxs n).symm) (sub_ne_zero.mpr (hxs (n+1)).symm)
    (hfx (n+1)) (hfx n) hfs (ne_of_gt (hμ n)) hDn
  have hE : 0 < (x (n+1)-r)*(x n-r) := by
    have hm := mul_pos_of_neg_of_neg (hside (n+1)) (hside n)
    by_contra hn
    have hn0 : (x (n+1)-r)*(x n-r) ≤ 0 := le_of_not_gt hn
    have hprod := mul_nonpos_of_nonpos_of_nonneg hn0 (mul_self_nonneg (s-r))
    nlinarith
  have hb := delta_bound_of_barycentric_no_cross
    (sub_ne_zero.mpr (hadj n).symm) (sub_ne_zero.mpr (hxs n).symm)
    (sub_ne_zero.mpr (hxs (n+1)).symm)
    (sub_ne_zero.mpr (hxr (n+1))) (sub_ne_zero.mpr (hxr n)) hsr
    hfs (hμ n) (hμone n) (hquot (n+1)) (hquot n) hquots hbridge.1
    (hstepn.trans hbridge.2) hE hPn hπn (hside (n+2))
  exact hb.le.trans hcapn.le

/-- Complete nondegenerate one-sided order theorem for actual raw proposals.
The positive limiting accumulated weight and golden Q-order are derived from
the tail's noncrossing signs and local regularity, without bounded-scaling
premises. The unmodified-tail and nonzero-leading-coefficient hypotheses are
the intended hypotheses of Proposition 12, not a claim of universal entry. -/
theorem rawPoint_one_sided_hasQOrder
    {f g : ℝ → ℝ} {x μ : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgrgs : 0 < g r * g s) (hs : s ≠ r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n ≠ r)
    (hadj : ∀ n, x (n+1) ≠ x n) (hxs : ∀ n, x n ≠ s)
    (hμ : ∀ n, 0 < μ n) (hμone : ∀ n, μ n ≤ 1) (hanti : Antitone μ)
    (hfs : f s ≠ 0) (hfx : ∀ n, f (x n) ≠ 0)
    (hquot : ∀ n, (x n-r)/f (x n) = g (x n))
    (hquots : (s-r)/f s = g s)
    (hside : ∀ n, (x n-r)*(s-r) < 0)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0)
    (hstep : ∀ᶠ n in atTop, x (n+2) =
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1) ≠ 0) :
    (∃ w : ℝ, 0 < w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  obtain ⟨M,hM⟩ := exists_scaling_bound_of_rawPoint_one_sided hg hgrgs hs hx hxr
    hadj hxs hμ hμone hanti hfs hfx hquot hquots hside hD hstep
  have hgr : g r ≠ 0 := by intro h; rw [h,zero_mul] at hgrgs; exact lt_irrefl _ hgrgs
  exact rawPoint_one_sided_hasQOrder_of_bounded_scaling hg hgr hs hx hxr
    hadj hxs hμ hμone hanti hM hfs hfx hquot hquots hD hstep hC

/-- The actual AB correction is at most one whenever the two fresh residuals
have the same sign. Its positive fallback is included. -/
theorem abFactor_le_one {new old : ℝ} (h : 0 < new/old) :
    ModAB.Convergence.abFactor new old ≤ 1 := by
  unfold ModAB.Convergence.abFactor
  split_ifs <;> linarith

/-- The weight hypotheses used above follow from the source's actual AB
factors, including its positive fallback. -/
theorem weights_of_AB_updates {μ new old : ℕ → ℝ}
    (hμ0 : 0 < μ 0) (hμ01 : μ 0 ≤ 1)
    (hratio : ∀ n, 0 < new n / old n)
    (hstep : ∀ n, μ (n+1) = μ n * ModAB.Convergence.abFactor (new n) (old n)) :
    (∀ n, 0 < μ n) ∧ (∀ n, μ n ≤ 1) ∧ Antitone μ := by
  have hμ : ∀ n, 0 < μ n := by
    intro n
    induction n with
    | zero => exact hμ0
    | succ n ih => rw [hstep]; exact mul_pos ih (ModAB.Convergence.abFactor_positive _ _)
  have hanti : Antitone μ := weight_antitone_of_factors (fun n => (hμ n).le)
    (fun n => abFactor_le_one (hratio n)) hstep
  exact ⟨hμ,fun n => (hanti (Nat.zero_le n)).trans hμ01,hanti⟩

/-- At every C³ simple root, the regularized reciprocal has the same nonzero
sign throughout some neighbourhood. No normalization of the function is
needed: negative derivatives are covered as well. -/
theorem exists_regularized_sign_neighbourhood {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ∀ z : ℝ, |z-r| < ε →
      0 < Regularity.regularizedReciprocal f r r *
        Regularity.regularizedReciprocal f r z := by
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 2 g r :=
    Regularity.regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd
  have hgr : g r ≠ 0 := by
    dsimp [g]
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hcont : ContinuousAt (fun z => g r * g z) r := continuousAt_const.mul hg.continuousAt
  have hevent : ∀ᶠ z in 𝓝 r, 0 < g r * g z :=
    hcont.tendsto.eventually (Ioi_mem_nhds (mul_self_pos.mpr hgr))
  obtain ⟨ε,hε,hball⟩ := Metric.mem_nhds_iff.mp hevent
  refine ⟨ε,hε,?_⟩
  intro z hz
  exact hball (by simpa only [Metric.mem_ball,Real.dist_eq] using hz)

/-- Algebraic conversion from a common nonzero sign relative to a reference
to a common sign between the two factors themselves. -/
theorem common_sign_product {q a b : ℝ}
    (ha : 0 < q*a) (hb : 0 < q*b) : 0 < a*b := by
  rcases mul_pos_iff.mp ha with ha | ha <;>
    rcases mul_pos_iff.mp hb with hb | hb
  · exact mul_pos ha.2 hb.2
  · linarith [ha.1,hb.1]
  · linarith [ha.1,hb.1]
  · exact mul_pos_of_neg_of_neg ha.2 hb.2

/-- Removing the reciprocal factors translates the geometric opposite-side
condition into the corresponding true-residual sign condition. -/
theorem opposite_true_residuals {a s r fa fs ga gs : ℝ}
    (hfa : fa ≠ 0) (hfs : fs ≠ 0)
    (hqa : (a-r)/fa = ga) (hqs : (s-r)/fs = gs)
    (hg : 0 < ga*gs) (he : (a-r)*(s-r) < 0) : fa*fs < 0 := by
  have hea : a-r = ga*fa := (div_eq_iff hfa).mp hqa
  have hes : s-r = gs*fs := (div_eq_iff hfs).mp hqs
  rw [hea,hes] at he
  by_contra hn
  have hn0 : 0 ≤ fa*fs := le_of_not_gt hn
  have hnonneg := mul_nonneg hg.le hn0
  nlinarith

/-- C³ version of Proposition 12. The neighbourhood, signs of the removable
reciprocal, positivity/monotonicity of the actual AB products, and their
uniform scaling bound are all derived. Only the explicitly stated unmodified
convergent tail and nonzero limiting coefficient remain assumptions. -/
theorem smooth_rawPoint_one_sided_hasQOrder
    {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ∀ (s : ℝ) (x μ : ℕ → ℝ),
      s ≠ r → |s-r| < ε →
      (∀ n, |x n-r| < ε) → Tendsto x atTop (𝓝 r) →
      (∀ n, x n ≠ r) → (∀ n, x (n+1) ≠ x n) → (∀ n, x n ≠ s) →
      (∀ n, (x n-r)*(s-r) < 0) →
      0 < μ 0 → μ 0 ≤ 1 →
      (∀ n, μ (n+1) = μ n *
        ModAB.Convergence.abFactor (f (x (n+2))) (f (x (n+1)))) →
      (∀ᶠ n in atTop, ModAB.Rational.denominator (x (n+1)) s
        (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0) →
      (∀ᶠ n in atTop, x (n+2) = ModAB.Rational.rawPoint (x (n+1)) s
        (f (x (n+1))) (μ n*f s) (x n) (f (x n))) →
      limitingCoefficient (Regularity.regularizedReciprocal f r) r s
        (⨆ n, 1/μ n-1) ≠ 0 →
      (∃ w : ℝ, 0 < w ∧ Tendsto μ atTop (𝓝 w)) ∧
      ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
        (|limitingCoefficient (Regularity.regularizedReciprocal f r) r s
          (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  obtain ⟨ε,hε,hsign⟩ := exists_regularized_sign_neighbourhood hf hd
  refine ⟨ε,hε,?_⟩
  intro s x μ hs hsmallS hsmallX hx hxr hadj hxs hside hμ0 hμ01 hμstep hD hstep hC
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 1 g r :=
    (Regularity.regularizedReciprocal_contDiffAt_of_dslope
      (DSlope.contDiffAt_dslope hf) hd).of_le (by norm_num)
  have hgs : 0 < g r * g s := hsign s hsmallS
  have hgx : ∀ n, 0 < g r * g (x n) := fun n => hsign (x n) (hsmallX n)
  have hquot : ∀ n, (x n-r)/f (x n) = g (x n) :=
    fun n => (Regularity.regularizedReciprocal_off_root hzero (hxr n)).symm
  have hquots : (s-r)/f s = g s :=
    (Regularity.regularizedReciprocal_off_root hzero hs).symm
  have hfs : f s ≠ 0 := by
    intro h
    have hg0 : g s = 0 := by rw [← hquots,h,div_zero]
    rw [hg0,mul_zero] at hgs
    exact lt_irrefl _ hgs
  have hfx : ∀ n, f (x n) ≠ 0 := by
    intro n h
    have hg0 : g (x n) = 0 := by rw [← hquot n,h,div_zero]
    have hgn := hgx n
    rw [hg0,mul_zero] at hgn
    exact lt_irrefl _ hgn
  have hres : ∀ n, f (x n)*f s < 0 := by
    intro n
    exact opposite_true_residuals (hfx n) hfs (hquot n) hquots
      (common_sign_product (hgx n) hgs) (hside n)
  have hratio : ∀ n, 0 < f (x (n+2))/f (x (n+1)) := by
    intro n
    have hm : 0 < f (x (n+2))*f (x (n+1)) :=
      common_sign_product (q := -f s) (by nlinarith [hres (n+2)])
        (by nlinarith [hres (n+1)])
    exact div_pos_iff.mpr (mul_pos_iff.mp hm)
  obtain ⟨hμ,hμone,hanti⟩ := weights_of_AB_updates hμ0 hμ01 hratio hμstep
  exact rawPoint_one_sided_hasQOrder hg hgs hs hx hxr hadj hxs hμ hμone hanti
    hfs hfx hquot hquots hside hD hstep hC

end
end ModAB.Rational.OneSidedComplete

end RModABBundle46

open Lean Elab Command

run_cmd do
  let declarations : Array Name := #[
    ``ModAB.Rational.secant_interior,
    ``ModAB.Rational.working_ratio_negative,
    ``ModAB.Rational.working_second_ne_zero,
    ``ModAB.Rational.secant_normalized,
    ``ModAB.Rational.radius_pos,
    ``ModAB.Rational.radius_le_left,
    ``ModAB.Rational.radius_le_right,
    ``ModAB.Rational.limited_bounds,
    ``ModAB.Rational.limited_margins,
    ``ModAB.Rational.limited_interior,
    ``ModAB.Rational.limited_eq_self,
    ``ModAB.Rational.candidate_interior,
    ``ModAB.Rational.candidate_fallback_denominator,
    ``ModAB.Rational.candidate_fallback_curvature,
    ``ModAB.Rational.sourceProposal_mem_closed,
    ``ModAB.Rational.sourceProposal_interior_of_ne_right,
    ``ModAB.Rational.candidate_is_sourceProposal,
    ``ModAB.Rational.closeCandidate_interior,
    ``ModAB.Rational.closeCandidate_mem_closed,
    ``ModAB.Rational.closeCandidate_left_separation,
    ``ModAB.Rational.closeCandidate_right_separation,
    ``ModAB.Rational.closeCandidate_right_fallback_interior,
    ``ModAB.Rational.numerator_ne_zero_of_denominator_zero,
    ``ModAB.Rational.normalized_interpolation_equations,
    ``ModAB.Rational.recordHistory_true_values,
    ``ModAB.Rational.recordHistory_dropped_true,
    ``ModAB.Rational.recordSourceHistory_true_values,
    ``ModAB.Rational.recordSourceHistory_history_true,
    ``ModAB.Rational.recordSourceHistory_reuse_left,
    ``ModAB.Rational.recordSourceHistory_reuse_right,
    ``ModAB.Rational.recordSourceHistory_left_eq_recordHistory,
    ``ModAB.Rational.recordSourceHistory_right_eq_recordHistory,
    ``ModAB.Rational.recordUpdate_fresh,
    ``ModAB.Rational.recordUpdate_at_most_one_modified,
    ``ModAB.Rational.correction_true_denominator,
    ``ModAB.Rational.abFactor_lt_one,
    ``ModAB.Rational.abFactor_mem_Ioo,
    ``ModAB.Rational.scaled_linear_candidate_example,
    ``ModAB.RationalBridge.updated_invariants_closed,
    ``ModAB.RationalBridge.closed_invariants,
    ``ModAB.RationalBridge.closed_zeros,
    ``ModAB.RationalBridge.capped_global_convergence,
    ``ModAB.RationalBridge.capped_global_width,
    ``ModAB.RationalBridge.capped_width_seven,
    ``ModAB.RationalBridge.capped_finite_tolerance,
    ``ModAB.RationalUncapped.eventually_short_steps,
    ``ModAB.RationalUncapped.no_convergent_separated_tail,
    ``ModAB.RationalUncapped.exists_width_stop_of_converges,
    ``ModAB.RationalResidual.limit_is_zero,
    ``ModAB.RationalResidual.eventually_near_zero_set,
    ``ModAB.RationalResidual.tendsto_unique_zero,
    ``ModAB.RationalUncapped.decision_step,
    ``ModAB.RationalUncapped.forty_two_test,
    ``ModAB.RationalUncapped.step_nonnegative,
    ``ModAB.RationalUncapped.run_nonnegative,
    ``ModAB.RationalUncapped.step_nonincrease,
    ``ModAB.RationalUncapped.width_antitone,
    ``ModAB.RationalUncapped.count_monotone,
    ``ModAB.RationalUncapped.step_count_bound,
    ``ModAB.RationalUncapped.bisection_count_bound,
    ``ModAB.RationalUncapped.count_unbounded,
    ``ModAB.RationalUncapped.infinite_bisections,
    ``ModAB.RationalUncapped.continuing_step,
    ``ModAB.RationalUncapped.continuing_suffix,
    ``ModAB.RationalUncapped.small_widths_tendsto,
    ``ModAB.RationalUncapped.geometric_residual,
    ``ModAB.RationalUncapped.width_or_residual,
    ``ModAB.RationalUncapped.finite_combined_tolerance,
    ``ModAB.RationalUncapped.exact_invariants,
    ``ModAB.RationalUncapped.point_endpoint,
    ``ModAB.RationalUncapped.point_mem,
    ``ModAB.RationalUncapped.point_mem_initial,
    ``ModAB.RationalUncapped.exact_zeros,
    ``ModAB.RationalUncapped.residual_of_width,
    ``ModAB.RationalUncapped.global_residual_convergence,
    ``ModAB.RationalUncapped.finite_residual_tolerance,
    ``ModAB.RationalUncapped.eta42_global_residual,
    ``ModAB.RationalUncapped.eventually_near_zero_set,
    ``ModAB.RationalUncapped.unique_root_convergence,
    ``ModAB.RationalUncapped.width_stop_with_closing,
    ``ModAB.RationalAssembly.source_closed_update,
    ``ModAB.RationalAssembly.proposal_mem_of_invariants,
    ``ModAB.RationalAssembly.source_invariants,
    ``ModAB.RationalAssembly.source_candidate_mem,
    ``ModAB.RationalAssembly.source_width_step,
    ``ModAB.RationalAssembly.source_bisection_width,
    ``ModAB.RationalAssembly.source_policy,
    ``ModAB.RationalAssembly.source_global_residual,
    ``ModAB.RationalAssembly.source_unique_root_convergence,
    ``ModAB.RationalAssembly.source_closing_separation,
    ``ModAB.RationalAssembly.source_finite_width_tolerance,
    ``ModAB.RationalAssembly.eta42_source_finite_width,
    ``ModAB.Rational.LocalAlgebra.dd2_add,
    ``ModAB.Rational.LocalAlgebra.dd2_smul,
    ``ModAB.Rational.LocalAlgebra.dd2_barycentric,
    ``ModAB.Rational.LocalAlgebra.dd2_affine,
    ``ModAB.Rational.LocalAlgebra.dd2_swap_last,
    ``ModAB.Rational.LocalAlgebra.dd2_swap_first,
    ``ModAB.Rational.LocalAlgebra.dd2_cyclic,
    ``ModAB.Rational.LocalAlgebra.baryRoot_swap_last,
    ``ModAB.Rational.LocalAlgebra.dd2_change_third,
    ``ModAB.Rational.LocalAlgebra.dd2_error_decomposition,
    ``ModAB.Rational.LocalAlgebra.baryRoot_error_identity,
    ``ModAB.Rational.LocalAlgebra.pole_identity,
    ``ModAB.Rational.LocalAlgebra.pole_identity_from_values,
    ``ModAB.Rational.LocalAlgebra.poleValue_margin,
    ``ModAB.Rational.LocalAlgebra.baryRoot_true_error,
    ``ModAB.Rational.LocalAlgebra.rescale_div,
    ``ModAB.Rational.LocalAlgebra.dd2_scaled_third,
    ``ModAB.Rational.LocalAlgebra.poleValue_scaled_third,
    ``ModAB.Rational.LocalAlgebra.baryRoot_scaled_error,
    ``ModAB.Rational.LocalAlgebra.scaled_crossing_ratio,
    ``ModAB.Rational.LocalAlgebra.scaled_pole_pos,
    ``ModAB.Rational.LocalAlgebra.scaled_denominator_ne_zero,
    ``ModAB.Rational.LocalAlgebra.scaled_error_positive_side,
    ``ModAB.Rational.LocalAlgebra.scaled_error_negative_side,
    ``ModAB.Rational.LocalAlgebra.delta_bound_of_same_side,
    ``ModAB.Rational.LocalAlgebra.negative_coefficient_of_same_side,
    ``ModAB.Rational.LocalAlgebra.weight_lower_bound,
    ``ModAB.Rational.LocalAlgebra.weight_from_delta,
    ``ModAB.Rational.LocalAlgebra.baryRoot_scaled_between,
    ``ModAB.Rational.LocalAlgebra.bary_denominator_eq,
    ``ModAB.Rational.LocalAlgebra.bary_denominator_ne_zero,
    ``ModAB.Rational.LocalAlgebra.baryRoot_eq_rawPoint,
    ``ModAB.Rational.LocalAlgebra.rawPoint_eq_baryRoot,
    ``ModAB.RationalOrderStability.tendsto_zero_of_contraction,
    ``ModAB.RationalOrderStability.energy_coercive,
    ``ModAB.RationalOrderStability.forcing_sq_le_energy,
    ``ModAB.RationalOrderStability.energy_update,
    ``ModAB.RationalOrderStability.energy_contracts,
    ``ModAB.RationalOrderStability.tendsto_zero_of_second_order,
    ``ModAB.RationalOrderStability.tendsto_zero_of_recurrence,
    ``ModAB.RationalOrder.exists_order_root,
    ``ModAB.RationalOrder.order_root_unique,
    ``ModAB.RationalOrder.tau_lower,
    ``ModAB.RationalOrder.tau_upper,
    ``ModAB.RationalOrder.tau_cubic,
    ``ModAB.RationalOrder.one_lt_tau,
    ``ModAB.RationalOrder.tau_pos,
    ``ModAB.RationalOrder.characteristic_factorization,
    ``ModAB.RationalOrder.log_abs_recurrence,
    ``ModAB.RationalOrder.transformed_recurrence,
    ``ModAB.RationalOrder.ratio_eq_exp_log_difference,
    ``ModAB.RationalOrder.q_ratio_of_log_difference,
    ``ModAB.RationalOrder.three_error_recurrence_log_limit,
    ``ModAB.RationalOrder.three_error_recurrence_q_ratio,
    ``ModAB.RationalOrder.three_error_recurrence_hasQOrder,
    ``ModAB.RationalOrder.eventual_three_error_recurrence_hasQOrder,
    ``ModAB.RationalGoldenOrder.tendsto_zero_of_first_order,
    ``ModAB.RationalGoldenOrder.log_abs_recurrence,
    ``ModAB.RationalGoldenOrder.two_error_recurrence_log_limit,
    ``ModAB.RationalGoldenOrder.two_error_recurrence_q_ratio,
    ``ModAB.RationalGoldenOrder.two_error_recurrence_hasQOrder,
    ``ModAB.RationalGoldenOrder.eventual_two_error_recurrence_hasQOrder,
    ``ModAB.RationalLocalSafeguards.triple_error_bound,
    ``ModAB.RationalLocalSafeguards.pair_error_bound,
    ``ModAB.RationalLocalSafeguards.radius_lower_of_secant_error,
    ``ModAB.RationalLocalSafeguards.limited_eq_raw_of_local_estimates,
    ``ModAB.RationalLocalSafeguards.forty_two_residual_of_local_bounds,
    ``ModAB.RationalLocalSafeguards.nextMode_interpolation_of_forty_two,
    ``ModAB.RationalLocalSafeguards.closeCandidate_zero,
    ``ModAB.RationalLocalSafeguards.interior_of_local_error,
    ``ModAB.RationalLocalSafeguards.local_operational_bridge,
    ``ModAB.RationalLocalSafeguards.positive_coefficient_crossing,
    ``ModAB.RationalLocalSafeguards.crossing_left_to_right_restores_true,
    ``ModAB.RationalLocalSafeguards.crossing_right_to_left_restores_true,
    ``ModAB.RationalLocalSafeguards.recurrence_local_invariant,
    ``ModAB.RationalLocalSafeguards.recurrence_tendsto_zero,
    ``ModAB.RationalLocalSafeguards.bounded_coefficient_recurrence_tendsto,
    ``ModAB.Rational.Cycle.cubicCoeff_identity,
    ``ModAB.Rational.Cycle.cubicCoeff_tendsto,
    ``ModAB.Rational.Cycle.abs_cubicCoeff_tendsto,
    ``ModAB.Rational.Cycle.first_order_identity,
    ``ModAB.Rational.Cycle.second_order_identity,
    ``ModAB.Rational.Cycle.first_constant,
    ``ModAB.Rational.Cycle.second_constant,
    ``ModAB.Rational.Cycle.first_order_tendsto,
    ``ModAB.Rational.Cycle.second_order_tendsto,
    ``ModAB.Rational.Cycle.alternating_orders_product,
    ``ModAB.Rational.Cycle.alternating_orders_geometric_mean,
    ``ModAB.Rational.Cycle.firstCoeff_abs,
    ``ModAB.Rational.Cycle.secondCoeff_abs,
    ``ModAB.Rational.Cycle.cycle_orders_of_coefficients,
    ``ModAB.Rational.CycleEntry.Cone.latest_ne_zero,
    ``ModAB.Rational.CycleEntry.Cone.opposite_ne_zero,
    ``ModAB.Rational.CycleEntry.Cone.history_ne_zero,
    ``ModAB.Rational.CycleEntry.Cone.ordered_errors,
    ``ModAB.Rational.CycleEntry.product_bound,
    ``ModAB.Rational.CycleEntry.pair_signs,
    ``ModAB.Rational.CycleEntry.pair_preserves_cone,
    ``ModAB.Rational.CycleEntry.Cone.enlarge_radius,
    ``ModAB.Rational.CycleEntry.Run.all_cones,
    ``ModAB.Rational.CycleEntry.Run.history_geometric_bound,
    ``ModAB.Rational.CycleEntry.Run.errors_tendsto,
    ``ModAB.Rational.CycleEntry.Run.separation_tendsto,
    ``ModAB.Rational.CycleEntry.Run.firstCoeff_eq,
    ``ModAB.Rational.CycleEntry.Run.secondCoeff_eq,
    ``ModAB.Rational.CycleEntry.Run.orders,
    ``ModAB.Rational.CycleEntry.latest_cubic_identity,
    ``ModAB.Rational.CycleEntry.latest_cubic_constant,
    ``ModAB.Rational.CycleEntry.Run.latest_cubic_order,
    ``ModAB.Rational.DSlope.power_quotient_limit,
    ``ModAB.Rational.DSlope.dslope_hasDerivAt_of_ne,
    ``ModAB.Rational.DSlope.dslope_hasDerivAt_same,
    ``ModAB.Rational.DSlope.first_numerator_derivative,
    ``ModAB.Rational.DSlope.first_continuousAt_same,
    ``ModAB.Rational.DSlope.first_hasDerivAt_of_ne,
    ``ModAB.Rational.DSlope.first_hasDerivAt_same,
    ``ModAB.Rational.DSlope.second_numerator_derivative,
    ``ModAB.Rational.DSlope.second_continuousAt_same,
    ``ModAB.Rational.DSlope.dslope_hasDerivAt,
    ``ModAB.Rational.DSlope.first_hasDerivAt,
    ``ModAB.Rational.DSlope.second_continuousAt,
    ``ModAB.Rational.DSlope.contDiffAt_dslope,
    ``ModAB.Rational.DSlope.deriv_dslope_at,
    ``ModAB.Rational.DSlope.second_deriv_dslope_at,
    ``ModAB.Rational.LocalAnalysis.quadratic_at_third,
    ``ModAB.Rational.LocalAnalysis.exists_second_derivative_eq_twice_dd2,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd1_of_derivatives,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd2_of_derivatives,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd2_of_contDiff_two,
    ``ModAB.Rational.LocalAnalysis.tendsto_poleValue_of_contDiff_two,
    ``ModAB.Rational.LocalAnalysis.tendsto_error_coefficient_of_contDiff_two,
    ``ModAB.Rational.LocalAnalysis.exists_second_derivative_dd2_arbitrary_of_derivatives,
    ``ModAB.Rational.LocalAnalysis.exists_second_derivative_dd2_arbitrary,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd2_arbitrary_of_contDiff_two,
    ``ModAB.Rational.LocalAnalysis.dd1_comm,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd1_arbitrary_of_contDiff_one,
    ``ModAB.Rational.LocalAnalysis.tendsto_error_coefficient_arbitrary_of_contDiff_two,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd2_arbitrary_of_derivatives,
    ``ModAB.Rational.LocalAnalysis.tendsto_dd1_arbitrary_of_derivatives,
    ``ModAB.Rational.LocalAnalysis.tendsto_divided_differences_of_contDiffAt,
    ``ModAB.Rational.LocalAnalysis.tendsto_poleValue_of_contDiffAt,
    ``ModAB.Rational.LocalAnalysis.tendsto_error_coefficient_of_contDiffAt,
    ``ModAB.Rational.Regularity.regularizedReciprocal_at,
    ``ModAB.Rational.Regularity.regularizedReciprocal_off_root,
    ``ModAB.Rational.Regularity.regularizedReciprocal_continuousAt,
    ``ModAB.Rational.Regularity.regularizedReciprocal_contDiffAt_of_dslope,
    ``ModAB.Rational.Regularity.analyticAt_dslope,
    ``ModAB.Rational.Regularity.regularizedReciprocal_analyticAt,
    ``ModAB.Rational.Regularity.regularizedReciprocal_contDiffAt,
    ``ModAB.Rational.Regularity.reciprocal_second_derivative,
    ``ModAB.Rational.Regularity.iteratedDeriv_shift_at,
    ``ModAB.Rational.Regularity.iteratedDeriv_shift_at_point,
    ``ModAB.Rational.Regularity.hasDerivAt_vanishing_factor,
    ``ModAB.Rational.Regularity.factor_second_derivative,
    ``ModAB.Rational.Regularity.factor_third_derivative_at,
    ``ModAB.Rational.Regularity.dslope_derivative_coefficients_of_contDiff,
    ``ModAB.Rational.Regularity.dslope_derivative_coefficients,
    ``ModAB.Rational.Regularity.regularizedReciprocal_second_coefficient_of_dslope,
    ``ModAB.Rational.Regularity.regularizedReciprocal_second_coefficient,
    ``ModAB.Rational.Regularity.regularizedReciprocal_error_constant_of_dslope,
    ``ModAB.Rational.Regularity.regularizedReciprocal_error_constant,
    ``ModAB.Rational.Regularity.regularizedReciprocal_exists_contDiffOn_interval,
    ``ModAB.Rational.Regularity.regularizedReciprocal_contDiffAt_of_contDiff,
    ``ModAB.Rational.Regularity.regularizedReciprocal_second_coefficient_of_contDiff,
    ``ModAB.Rational.Regularity.regularizedReciprocal_error_constant_of_contDiff,
    ``ModAB.Rational.Regularity.regularizedReciprocal_exists_contDiffOn_interval_of_contDiff,
    ``ModAB.Rational.AnalyticBridge.tendsto_analytic_error_coefficient,
    ``ModAB.Rational.AnalyticBridge.eventually_analytic_denominator_ne_zero,
    ``ModAB.Rational.AnalyticBridge.tendsto_normalized_rational_error,
    ``ModAB.Rational.AnalyticBridge.eventually_raw_denominator_ne_zero,
    ``ModAB.Rational.AnalyticBridge.tendsto_normalized_rawPoint_error,
    ``ModAB.Rational.SmoothBridge.tendsto_smooth_error_coefficient,
    ``ModAB.Rational.SmoothBridge.eventually_smooth_denominator_ne_zero,
    ``ModAB.Rational.SmoothBridge.tendsto_normalized_rational_error,
    ``ModAB.Rational.SmoothBridge.eventually_raw_denominator_ne_zero,
    ``ModAB.Rational.SmoothBridge.tendsto_normalized_rawPoint_error,
    ``ModAB.Rational.OneSided.weight_antitone_of_factors,
    ``ModAB.Rational.OneSided.delta_monotone_of_weight_antitone,
    ``ModAB.Rational.OneSided.exists_delta_weight_limit,
    ``ModAB.Rational.OneSided.exists_positive_weight_limit,
    ``ModAB.Rational.OneSided.exists_positive_weight_limit_of_same_side,
    ``ModAB.Rational.OneSided.tendsto_dd2_fixed_endpoint,
    ``ModAB.Rational.OneSided.tendsto_poleValue_fixed_endpoint,
    ``ModAB.Rational.OneSided.baryRoot_error_eq_mixedCoefficient,
    ``ModAB.Rational.OneSided.tendsto_mixedCoefficient,
    ``ModAB.Rational.OneSided.tendsto_dd1_of_contDiffAt,
    ``ModAB.Rational.OneSided.tendsto_dd2_fixed_endpoint_of_contDiffAt,
    ``ModAB.Rational.OneSided.tendsto_poleValue_fixed_endpoint_of_contDiffAt,
    ``ModAB.Rational.OneSided.tendsto_mixedCoefficient_of_contDiffAt,
    ``ModAB.Rational.OneSidedAssembly.barycentric_one_sided_hasQOrder,
    ``ModAB.Rational.OneSidedAssembly.barycentric_one_sided_hasQOrder_of_bounded_scaling,
    ``ModAB.Rational.OneSidedAssembly.rawPoint_barycentric_bridge,
    ``ModAB.Rational.OneSidedAssembly.rawPoint_one_sided_hasQOrder_of_bounded_scaling,
    ``ModAB.RationalLocalAssembly.bounded_coefficient_recurrence_hasQOrder,
    ``ModAB.RationalLocalAssembly.cycle_cone_hasQOrders,
    ``ModAB.Rational.LocalNeighborhood.nodes_tendsto_of_small_radius,
    ``ModAB.Rational.LocalNeighborhood.uniform_three_of_sequential,
    ``ModAB.Rational.LocalNeighborhood.exists_uniform_raw_coefficient,
    ``ModAB.Rational.LocalNeighborhood.exists_positive_raw_bounds,
    ``ModAB.Rational.LocalNeighborhood.raw_error_eq_coefficient,
    ``ModAB.Rational.LocalState.Valid.admissible,
    ``ModAB.Rational.LocalState.Valid.radius,
    ``ModAB.Rational.LocalState.latest_update,
    ``ModAB.Rational.LocalState.interior_near,
    ``ModAB.Rational.LocalState.fresh_update_right,
    ``ModAB.Rational.LocalState.fresh_update_left,
    ``ModAB.Rational.LocalState.step_preserves,
    ``ModAB.Rational.Alternating.orbit_zero,
    ``ModAB.Rational.Alternating.orbit_succ,
    ``ModAB.Rational.Alternating.valid_orbit,
    ``ModAB.Rational.Alternating.continue_orbit,
    ``ModAB.Rational.Alternating.opposite_shift,
    ``ModAB.Rational.Alternating.history_two_shifts,
    ``ModAB.Rational.Alternating.latest_ne_root,
    ``ModAB.Rational.Alternating.latest_tendsto,
    ``ModAB.Rational.Alternating.opposite_tendsto,
    ``ModAB.Rational.Alternating.history_tendsto,
    ``ModAB.Rational.Alternating.endpoints_min_max,
    ``ModAB.Rational.Alternating.endpoints_tendsto,
    ``ModAB.Rational.Alternating.endpoint_product,
    ``ModAB.Rational.Alternating.raw_recurrence,
    ``ModAB.Rational.Alternating.hasQOrder_of_neighborhood,
    ``ModAB.Rational.Curvature.curvatureNumerator_eq,
    ``ModAB.Rational.Curvature.curvatureNumerator_sign,
    ``ModAB.Rational.Curvature.curvatureChangesSign_eq_false,
    ``ModAB.Rational.Curvature.exists_radius_no_curvature_change,
    ``ModAB.Rational.CycleUniform.tendsto_firstCoefficient,
    ``ModAB.Rational.CycleUniform.tendsto_firstPole,
    ``ModAB.Rational.CycleUniform.tendsto_divided_differences_eventually,
    ``ModAB.Rational.CycleUniform.eventually_ne_of_ratio_zero,
    ``ModAB.Rational.CycleUniform.first_step_limits,
    ``ModAB.Rational.CycleUniform.second_step_limits,
    ``ModAB.Rational.CycleUniform.smallRadius_pos,
    ``ModAB.Rational.CycleUniform.smallRadius_lt_one,
    ``ModAB.Rational.CycleUniform.smallRadius_tendsto,
    ``ModAB.Rational.CycleUniform.exists_uniform_cone_of_sequential,
    ``ModAB.Rational.CycleUniform.Cone.distinct,
    ``ModAB.Rational.CycleUniform.shrinking_cone_limits,
    ``ModAB.Rational.CycleUniform.exists_uniform_negative_coefficients,
    ``ModAB.Rational.CycleUniform.localK_regularizedReciprocal,
    ``ModAB.Rational.CycleUniform.smooth_exists_uniform_negative_coefficients,
    ``ModAB.Rational.CycleState.cone_enlarge_theta,
    ``ModAB.Rational.CycleState.Bounds.restrict,
    ``ModAB.Rational.CycleState.exists_bounds,
    ``ModAB.Rational.CycleState.orbit_cone,
    ``ModAB.Rational.CycleState.coefficient_limits,
    ``ModAB.Rational.CycleState.orbit_hasQOrders,
    ``ModAB.Rational.CycleUniformAlgebra.delta_one_sub,
    ``ModAB.Rational.CycleUniformAlgebra.residual_ratio_from_step,
    ``ModAB.Rational.CycleUniformAlgebra.relativeDenominator_eq,
    ``ModAB.Rational.CycleUniformAlgebra.relativeDenominator_ne_zero,
    ``ModAB.Rational.CycleUniformAlgebra.mixed_denominator_identity,
    ``ModAB.Rational.CycleUniformAlgebra.bary_mixed_pole_identity,
    ``ModAB.Rational.CycleUniformAlgebra.bary_mixed_denominator_ne_zero,
    ``ModAB.Rational.CycleUniformAlgebra.source_denominator_ne_zero_of_mixed_pole,
    ``ModAB.Rational.CycleUniformAlgebra.normalized_second_identity,
    ``ModAB.Rational.CycleUniformAlgebra.second_error_eq,
    ``ModAB.Rational.CycleUniformAlgebra.baryRoot_second_error,
    ``ModAB.Rational.CycleUniformAlgebra.rawPoint_second_error,
    ``ModAB.Rational.CycleUniformAlgebra.source_rawPoint_second_error,
    ``ModAB.Rational.CycleUniformAlgebra.denominator_swap_ne_zero,
    ``ModAB.Rational.CycleUniformAlgebra.rawPoint_swap,
    ``ModAB.Rational.CycleUniformBridge.reciprocal_value_at_error,
    ``ModAB.Rational.CycleUniformBridge.ordinate_ne_zero_of_reciprocal,
    ``ModAB.Rational.CycleUniformBridge.first_source_denominator_ne_zero,
    ``ModAB.Rational.CycleUniformBridge.first_rawPoint_error,
    ``ModAB.Rational.CycleUniformBridge.residual_ratio_eq_ordinateRatio,
    ``ModAB.Rational.CycleUniformBridge.second_source_denominator_ne_zero,
    ``ModAB.Rational.CycleUniformBridge.second_rawPoint_error,
    ``ModAB.Rational.CycleUniformBridge.smooth_exists_uniform_raw_pair,
    ``ModAB.Rational.SecantLocal.normalizedSlope_at,
    ``ModAB.Rational.SecantLocal.normalizedSlope_continuousAt,
    ``ModAB.Rational.SecantLocal.function_eq_normalizedSlope,
    ``ModAB.Rational.SecantLocal.exists_radius_normalizedSlope,
    ``ModAB.Rational.SecantLocal.residual_bounds_of_normalizedSlope,
    ``ModAB.Rational.SecantLocal.sign_product_of_normalizedSlope,
    ``ModAB.Rational.SecantLocal.exists_radius_simple_root_bounds,
    ``ModAB.Rational.SecantLocal.secant_error_of_normalized_values,
    ``ModAB.Rational.SecantLocal.exists_radius_secant_error,
    ``ModAB.Rational.SecantLocal.normalized_product_bound,
    ``ModAB.Rational.SecantLocal.exists_radius_scaled_secant_error,
    ``ModAB.Rational.CycleSafeguards.residual_ratio_bound,
    ``ModAB.Rational.CycleSafeguards.abFactor_near_one,
    ``ModAB.Rational.CycleSafeguards.exists_scaled_neighborhood,
    ``ModAB.Rational.CycleSafeguards.ScaledNeighborhood.shrink,
    ``ModAB.Rational.CycleSafeguards.ScaledNeighborhood.guards,
    ``ModAB.Rational.CycleSafeguards.ScaledNeighborhood.ab_factor,
    ``ModAB.Rational.CycleNeighborhood.exists_pair_neighborhood,
    ``ModAB.Rational.CycleIntermediate.same_side_update_left,
    ``ModAB.Rational.CycleIntermediate.same_side_update_right,
    ``ModAB.Rational.CycleIntermediate.same_side_update,
    ``ModAB.Rational.CycleIntermediate.abFactor_eq_one_sub,
    ``ModAB.Rational.CycleIntermediate.same_side_working_ordinates,
    ``ModAB.Rational.CycleIntermediate.same_side_geometry,
    ``ModAB.Rational.LocalComplete.exists_operational_neighborhood,
    ``ModAB.Rational.LocalComplete.local_alternating_complete,
    ``ModAB.Rational.CycleStateGuards.Geometry.of_valid,
    ``ModAB.Rational.CycleStateGuards.Geometry.of_same_side,
    ``ModAB.Rational.CycleStateGuards.function_ne_zero,
    ``ModAB.Rational.CycleStateGuards.nextMode_of_trial_residual,
    ``ModAB.Rational.CycleStateGuards.guards_latest,
    ``ModAB.Rational.CycleStateUpdate.pair_update_left,
    ``ModAB.Rational.CycleStateUpdate.pair_update_right,
    ``ModAB.Rational.CycleStateUpdate.pair_update_preserves,
    ``ModAB.Rational.CycleTransport.orbit_even_succ,
    ``ModAB.Rational.CycleTransport.even_orbit,
    ``ModAB.Rational.CycleTransport.odd_latest,
    ``ModAB.Rational.CycleTransport.continue_pair,
    ``ModAB.Rational.CycleTransport.source_hasQOrders,
    ``ModAB.Rational.CycleComplete.endpoint_min,
    ``ModAB.Rational.CycleComplete.valid_latest_near,
    ``ModAB.Rational.CycleComplete.valid_opposite_near,
    ``ModAB.Rational.CycleComplete.valid_history_near,
    ``ModAB.Rational.CycleComplete.valid_working,
    ``ModAB.Rational.CycleComplete.pair_law_of_neighborhood,
    ``ModAB.Rational.CycleComplete.local_pair_complete,
    ``ModAB.Rational.OneSidedComplete.monotone_bound_of_eventually,
    ``ModAB.Rational.OneSidedComplete.tendsto_no_cross_cap,
    ``ModAB.Rational.OneSidedComplete.delta_bound_of_same_side_same_sign,
    ``ModAB.Rational.OneSidedComplete.delta_bound_of_barycentric_no_cross,
    ``ModAB.Rational.OneSidedComplete.exists_scaling_bound_of_rawPoint_one_sided,
    ``ModAB.Rational.OneSidedComplete.rawPoint_one_sided_hasQOrder,
    ``ModAB.Rational.OneSidedComplete.abFactor_le_one,
    ``ModAB.Rational.OneSidedComplete.weights_of_AB_updates,
    ``ModAB.Rational.OneSidedComplete.exists_regularized_sign_neighbourhood,
    ``ModAB.Rational.OneSidedComplete.common_sign_product,
    ``ModAB.Rational.OneSidedComplete.opposite_true_residuals,
    ``ModAB.Rational.OneSidedComplete.smooth_rawPoint_one_sided_hasQOrder
  ]
  let allowed : Array Name := #[``propext, ``Classical.choice, ``Quot.sound]
  for name in declarations do
    let dependencies ← Lean.collectAxioms name
    for dependency in dependencies do
      unless allowed.contains dependency do
        throwError "Unexpected axiom {dependency} in {name}"
    logInfo m!"{name}: {dependencies}"
  logInfo m!"PASS: {declarations.size} declarations; only standard logical axioms."
