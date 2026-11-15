import ModAB.Convergence
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Analysis.Calculus.ContDiff.Deriv
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Topology.Order.OrderClosed

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
