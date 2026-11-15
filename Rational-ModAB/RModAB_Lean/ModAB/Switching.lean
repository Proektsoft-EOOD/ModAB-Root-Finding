import Mathlib.Basic.Real.Basic
import Mathlib.Tactic
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Inv

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
