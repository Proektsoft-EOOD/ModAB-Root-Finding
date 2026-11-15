import ModAB.RationalLocalAlgebra

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
