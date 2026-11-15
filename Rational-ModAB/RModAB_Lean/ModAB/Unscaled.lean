import ModAB.Binary64
import ModAB.Cubic
import FloatLib.Floats.Formats.Flocq.Theory.Error.Relative

/-! Full operation-by-operation proof of `prop:switch_raw_error`.
The assumptions concern normal exact intermediate results and absence of overflow,
not already propagated errors of the two sides. The final subtraction is exact.
-/
namespace ModAB.Unscaled
open Binary64
open FloatLib.Floats.Formats.Flocq FloatLib.Numerics
noncomputable section

abbrev eps : ℝ := binary64Epsilon
def minNormal : ℝ := (2:ℝ)^(-1022:ℤ)
def NormalOrZero (x : ℝ) : Prop := x=0 ∨ minNormal≤|x|
def Relative (computed exact : ℝ) : Prop :=
  ∃ δ : ℝ, |δ|≤eps ∧ computed=exact*(1+δ)

theorem eps_bounds : 0<eps ∧ eps<1 ∧ 0<gamma 2 ∧ 0<gamma 3 ∧ gamma 3<1 := by
  norm_num [eps,binary64Epsilon,gamma]

theorem rounded_relative {x : ℝ} (hx : NormalOrZero x) : Relative (Binary64.round x) x := by
  rcases hx with hz | hn
  · subst x; exact ⟨0,by simpa using le_of_lt binary64Epsilon_pos,by simp⟩
  have hp : 0<minNormal := by unfold minNormal; positivity
  have hx0 : x≠0 := by intro h; simp [h] at hn; linarith
  have hr := relative_error_round_FLT_normal (β := binaryRadix) (-1074) 53
    (by norm_num) nearestEven x hx0 (by simpa [minNormal,pow_value] using hn)
  refine ⟨(Binary64.round x-x)/x,?_,?_⟩
  · have hc : bpow binaryRadix (1-53)/2 = eps := by norm_num [pow_value,eps,binary64Epsilon]
    rw [hc] at hr
    simpa [ErrorBounds.relativeError,Binary64.round,exponent,abs_div] using hr
  · field_simp <;> ring

theorem relative_error {a x : ℝ} (h : Relative a x) : Within a x (eps*|x|) := by
  obtain ⟨δ,hδ,rfl⟩ := h
  unfold Within
  have he : x*(1+δ)-x=x*δ := by ring
  rw [he,abs_mul]
  nlinarith [mul_le_mul_of_nonneg_left hδ (abs_nonneg x)]

theorem relative_abs {a x : ℝ} (h : Relative a x) : Relative |a| |x| := by
  obtain ⟨δ,hδ,rfl⟩ := h
  have hd := (abs_le.mp hδ).1
  refine ⟨δ,hδ,?_⟩
  rw [abs_mul,abs_of_pos (by linarith [eps_bounds.2.1] : 0<1+δ)]

theorem normal_halving_exact {x : ℝ} (hx : Representable x)
    (hn : NormalOrZero (x/2)) : Binary64.round (x/2)=x/2 := by
  rcases hn with h0 | hn
  · rw [h0,Binary64.round_zero]
  have hmin : (2:ℝ)^(-1021:ℤ) ≤ |x| := by
    have he : (2:ℝ)^(-1021:ℤ)=2*minNormal := by
      unfold minNormal; rw [show (-1021:ℤ)=1+(-1022) by omega, zpow_add₀ (by norm_num : (2:ℝ)≠0)]; norm_num
    rw [he]
    rw [abs_div,abs_of_pos (by norm_num : (0:ℝ)<2)] at hn
    linarith
  have hm := magnitude_mono_abs binaryRadix (by positivity : (2:ℝ)^(-1021:ℤ)≠0)
    (by simpa [abs_of_pos (by positivity : 0<(2:ℝ)^(-1021:ℤ))] using hmin)
  have hmag : magnitude binaryRadix ((2:ℝ)^(-1021:ℤ)) = -1020 := by
    rw [←pow_value]; simp
  rw [hmag] at hm
  have hf := generic_format_FLT_mul_bpow (β := binaryRadix)
    (-1074) 53 (by norm_num) hx (-1) (by omega)
  exact round_exact (by simpa [pow_value,div_eq_mul_inv] using hf)

theorem factor_two_bound {a b : ℝ} (ha : |a|≤eps) (hb : |b|≤eps) :
    |(1+a)*(1+b)-1|≤2*eps+eps^2 := by
  have hmul := product_error (x:=1) (x':=1+a) (y:=1) (y':=1+b)
    (A:=1+eps) (B:=eps) (C:=eps) (D:=1)
    ((abs_add_le 1 a).trans (by simpa using add_le_add_left ha 1))
    (by simpa [Within] using hb) (by simpa [Within] using ha) (by norm_num)
  convert hmul using 1 <;> (try unfold Within) <;> ring

theorem factor_two_gamma {a b : ℝ} (ha : |a|≤eps) (hb : |b|≤eps) :
    |(1+a)*(1+b)-1|≤gamma 2 := by
  exact (factor_two_bound ha hb).trans (by norm_num [eps,binary64Epsilon,gamma])

theorem factor_three_gamma {a b c : ℝ}
    (ha : |a|≤eps) (hb : |b|≤eps) (hc : |c|≤eps) :
    |(1+a)*(1+b)*(1+c)-1|≤gamma 3 := by
  have he : 0≤eps := le_of_lt eps_bounds.1
  have hab := factor_two_bound ha hb
  have habs : |(1+a)*(1+b)|≤(1+eps)^2 := by
    rw [abs_mul]
    have aa : |1+a|≤1+eps := (abs_add_le _ _).trans (by simpa using add_le_add_left ha 1)
    have bb : |1+b|≤1+eps := (abs_add_le _ _).trans (by simpa using add_le_add_left hb 1)
    nlinarith [mul_le_mul aa bb (abs_nonneg _) (by linarith : 0≤1+eps)]
  have hm := product_error (x:=1) (x':=(1+a)*(1+b)) (y:=1) (y':=1+c)
    habs (by simpa [Within] using hc) hab (by norm_num : |(1:ℝ)|≤1)
  apply (show |(1+a)*(1+b)*(1+c)-1|≤_ by simpa [Within] using hm).trans
  norm_num [eps,binary64Epsilon,gamma]

theorem factor_quotient_gamma {a b c : ℝ}
    (ha : |a|≤eps) (hb : |b|≤eps) (hc : |c|≤eps) :
    |(1+a)*(1+c)/(1+b)-1|≤gamma 3 := by
  have hb0 : 0<1+b := by have h := (abs_le.mp hb).1; linarith [eps_bounds.2.1]
  have hden : 1-eps≤1+b := by have h := (abs_le.mp hb).1; linarith
  have hnum := within_sub (factor_two_bound ha hc) (by simpa [Within] using hb : Within (1+b) 1 eps)
  have hid : (1+a)*(1+c)/(1+b)-1=((1+a)*(1+c)-(1+b))/(1+b) := by field_simp
  rw [hid,abs_div,abs_of_pos hb0]
  apply (div_le_iff₀ hb0).mpr
  have hconst : (2*eps+eps^2)+eps ≤ gamma 3*(1-eps) := by norm_num [eps,binary64Epsilon,gamma]
  have hm := mul_le_mul_of_nonneg_left hden (le_of_lt eps_bounds.2.2.2.1)
  simp only [Within,sub_self,sub_zero] at hnum
  linarith

structure Values where
  (sum diff y v r k leftDiff productY productZ right : ℝ)

def evaluate (f1 f2 f3 : ℝ) : Values :=
  let su := Binary64.round (f1+f2)
  let di := Binary64.round (f2-f1)
  let y := Binary64.round (su/2)
  let v := Binary64.round (|y|/|di|)
  let r := Binary64.round (1-v)
  let k := Binary64.round (r*r)
  let leftDiff := Binary64.round (f3-y)
  let py := Binary64.round (k*|y|)
  let pz := Binary64.round (k*|f3|)
  ⟨su,di,y,v,r,k,leftDiff,py,pz,Binary64.round (py+pz)⟩

/-- Exactly the ten arguments to rounded arithmetic operations in the manuscript. -/
def arguments (f1 f2 f3 : ℝ) : List ℝ :=
  let t:=evaluate f1 f2 f3
  [f1+f2,f2-f1,t.sum/2,|t.y|/|t.diff|,1-t.v,t.r*t.r,
   f3-t.y,t.k*|t.y|,t.k*|f3|,t.productY+t.productZ]

structure NormalRun (f1 f2 f3 : ℝ) : Prop where
  sum : NormalOrZero (f1+f2)
  diff : NormalOrZero (f2-f1)
  half : NormalOrZero ((evaluate f1 f2 f3).sum/2)
  quotient : NormalOrZero (|(evaluate f1 f2 f3).y|/|(evaluate f1 f2 f3).diff|)
  subtract : NormalOrZero (1-(evaluate f1 f2 f3).v)
  square : NormalOrZero ((evaluate f1 f2 f3).r*(evaluate f1 f2 f3).r)
  left : NormalOrZero (f3-(evaluate f1 f2 f3).y)
  productY : NormalOrZero ((evaluate f1 f2 f3).k*|(evaluate f1 f2 f3).y|)
  productZ : NormalOrZero ((evaluate f1 f2 f3).k*|f3|)
  right : NormalOrZero ((evaluate f1 f2 f3).productY+(evaluate f1 f2 f3).productZ)

def NoOverflow (f1 f2 f3 : ℝ) : Prop :=
  ∀ x∈arguments f1 f2 f3, Finite (Binary64.round x)

theorem midpoint_relative {f1 f2 f3 : ℝ} (hn : NormalRun f1 f2 f3) :
    Relative (evaluate f1 f2 f3).y ((f1+f2)/2) := by
  obtain ⟨δ,hd,hs⟩ := rounded_relative hn.sum
  have hh := normal_halving_exact (round_representable (f1+f2)) hn.half
  refine ⟨δ,hd,?_⟩
  change Binary64.round (Binary64.round (f1+f2)/2)=_
  rw [hh,hs]; ring

theorem quotient_error {f1 f2 f3 : ℝ} (hs : f1*f2<0) (hn : NormalRun f1 f2 f3) :
    Within (evaluate f1 f2 f3).v (|(f1+f2)/2|/|f2-f1|)
      (gamma 3*(|(f1+f2)/2|/|f2-f1|)) := by
  obtain ⟨a,ha,hy⟩ := relative_abs (midpoint_relative hn)
  obtain ⟨b,hb,hd⟩ := relative_abs (rounded_relative hn.diff)
  obtain ⟨c,hc,hv⟩ := rounded_relative hn.quotient
  have hp : 0 < |f2-f1| := Cubic.endpoint_difference_pos hs
  have hb0 : 0<1+b := by have h := (abs_le.mp hb).1; linarith [eps_bounds.2.1]
  have hformula : (evaluate f1 f2 f3).v =
      (|(f1+f2)/2|/|f2-f1|)*((1+a)*(1+c)/(1+b)) := by
    change Binary64.round (_/_)=_ at hv
    change Binary64.round (|(evaluate f1 f2 f3).y|/|(evaluate f1 f2 f3).diff|)=_
    rw [hv,hy]
    change |Binary64.round (f2-f1)|=_ at hd
    change (|(f1+f2)/2| *(1+a)/|Binary64.round (f2-f1)|)*(1+c)=_
    rw [hd]; field_simp <;> ring
  unfold Within
  rw [hformula]
  have hq : 0≤|(f1+f2)/2|/|f2-f1| := by positivity
  have hid : (|(f1+f2)/2|/|f2-f1|)*((1+a)*(1+c)/(1+b))-
      |(f1+f2)/2|/|f2-f1| =
      (|(f1+f2)/2|/|f2-f1|)*((1+a)*(1+c)/(1+b)-1) := by ring
  rw [hid,abs_mul,abs_of_nonneg hq]
  nlinarith [mul_le_mul_of_nonneg_left (factor_quotient_gamma ha hb hc) hq]

theorem factor_ranges {f1 f2 f3 : ℝ} (hs : f1*f2<0) (hn : NormalRun f1 f2 f3) :
    let t:=evaluate f1 f2 f3
    0≤t.v ∧ t.v<1 ∧ t.r∈Set.Icc (0:ℝ) 1 ∧ t.k∈Set.Icc (0:ℝ) 1 := by
  let t:=evaluate f1 f2 f3
  have hv0 : 0≤t.v := by
    change 0 ≤ Binary64.round (|t.y|/|t.diff|)
    have h := round_monotone (show (0:ℝ)≤|t.y|/|t.diff| by positivity)
    simpa only [Binary64.round_zero] using h
  have hv := quotient_error hs hn
  have vr := Cubic.ratio_range hs
  have hvl := (abs_le.mp hv).2
  have hv1 : t.v<1 := by
    have hgam := eps_bounds.2.2.2
    have hm := mul_le_mul_of_nonneg_left vr.2 hgam.1.le
    change t.v- (|(f1+f2)/2|/|f2-f1|) ≤ _ at hvl
    nlinarith
  have hr : t.r∈Set.Icc (0:ℝ) 1 := by
    change Binary64.round (1-t.v)∈Set.Icc (0:ℝ) 1
    exact round_interval Binary64.round_zero Binary64.round_one
      (show 1-t.v∈Set.Icc (0:ℝ) 1 by constructor <;> linarith)
  have hk : t.k∈Set.Icc (0:ℝ) 1 := by
    change Binary64.round (t.r*t.r)∈Set.Icc (0:ℝ) 1
    exact round_interval Binary64.round_zero Binary64.round_one
      (show t.r*t.r∈Set.Icc (0:ℝ) 1 by constructor <;> nlinarith [hr.1,hr.2])
  exact ⟨hv0,hv1,hr,hk⟩

theorem factor_error {f1 f2 f3 : ℝ} (hs : f1*f2<0) (hn : NormalRun f1 f2 f3) :
    Within (evaluate f1 f2 f3).k ((1-|(f1+f2)/2|/|f2-f1|)^2) (2*gamma 3) := by
  let t:=evaluate f1 f2 f3
  let v:=|(f1+f2)/2|/|f2-f1|
  obtain ⟨d,hd,hr⟩ := rounded_relative hn.subtract
  obtain ⟨c,hc,hk⟩ := rounded_relative hn.square
  have hmodel : t.k=(1-t.v)^2*((1+d)*(1+d)*(1+c)) := by
    change t.r=(1-t.v)*(1+d) at hr
    change t.k=(t.r*t.r)*(1+c) at hk
    rw [hk,hr]; ring
  have tr := factor_ranges hs hn
  have vr := Cubic.ratio_range hs
  have hh : 0≤(1-t.v)^2 ∧ (1-t.v)^2≤1 := by
    constructor
    · positivity
    · nlinarith [tr.1,tr.2.1]
  have hlocal : Within t.k ((1-t.v)^2) (gamma 3) := by
    unfold Within
    rw [hmodel]
    have hid : (1-t.v)^2*((1+d)*(1+d)*(1+c))-(1-t.v)^2 =
        (1-t.v)^2*((1+d)*(1+d)*(1+c)-1) := by ring
    rw [hid,abs_mul,abs_of_nonneg hh.1]
    have hg := factor_three_gamma hd hd hc
    have hm := mul_le_mul hh.2 hg (abs_nonneg _) (by norm_num : (0:ℝ)≤1)
    simpa using hm
  have hq := quotient_error hs hn
  have hrerr := within_sub (within_refl (1:ℝ)) hq
  have hsquare := square_error
    (x:=1-t.v) (y:=1-v)
    (by constructor <;> linarith [tr.1,tr.2.1])
    (by constructor <;> dsimp [v] <;> linarith [vr.1,vr.2]) hrerr
  have h := within_trans hlocal hsquare
  apply within_mono h
  have hm := mul_le_mul_of_nonneg_left vr.2 (le_of_lt eps_bounds.2.2.2.1)
  dsimp [v] at *
  linarith

/-- The two products and the common final sum give two relative factors per term. -/
theorem distributed_error {a b A B R : ℝ} (ha : 0≤a) (hb : 0≤b)
    (hA : Relative A a) (hB : Relative B b) (hR : Relative R (A+B)) :
    Within R (a+b) (gamma 2*(a+b)) := by
  obtain ⟨d,hd,rfl⟩ := hA
  obtain ⟨e,he,rfl⟩ := hB
  obtain ⟨c,hc,rfl⟩ := hR
  unfold Within
  have hid : (a*(1+d)+b*(1+e))*(1+c)-(a+b)=
      a*((1+d)*(1+c)-1)+b*((1+e)*(1+c)-1) := by ring
  rw [hid]
  calc
    |a*((1+d)*(1+c)-1)+b*((1+e)*(1+c)-1)| ≤
      |a*((1+d)*(1+c)-1)|+|b*((1+e)*(1+c)-1)| := abs_add_le _ _
    _ = a*|(1+d)*(1+c)-1|+b*|(1+e)*(1+c)-1| := by
      rw [abs_mul,abs_mul,abs_of_nonneg ha,abs_of_nonneg hb]
    _ ≤ a*gamma 2+b*gamma 2 := add_le_add
      (mul_le_mul_of_nonneg_left (factor_two_gamma hd hc) ha)
      (mul_le_mul_of_nonneg_left (factor_two_gamma he hc) hb)
    _ = gamma 2*(a+b) := by ring

/-- All side bounds are conclusions from individual concrete rounding operations. -/
theorem side_errors {f1 f2 f3 : ℝ} (hs : f1*f2<0) (hn : NormalRun f1 f2 f3) :
    let t:=evaluate f1 f2 f3
    let H:=max |f2-f1| |f3|
    let k:=(1-|(f1+f2)/2|/|f2-f1|)^2
    Within |t.leftDiff| (originalL f1 f2 f3) ((2*eps+eps^2/2)*H) ∧
    Within t.right (originalR k f1 f2 f3)
      (((3+eps)/2*gamma 2+3*gamma 3+eps/2)*H) := by
  let t:=evaluate f1 f2 f3
  let y:=(f1+f2)/2
  let H:=max |f2-f1| |f3|
  let k:=(1-|(f1+f2)/2|/|f2-f1|)^2
  have hH : 0≤H := (abs_nonneg (f2-f1)).trans (le_max_left _ _)
  have hyH : |y|≤H/2 := (Cubic.midpoint_bound hs).trans (by
    exact div_le_div_of_nonneg_right (le_max_left _ _) (by norm_num))
  have hzH : |f3|≤H := le_max_right _ _
  have he : 0≤eps := le_of_lt eps_bounds.1
  have hy : Within t.y y (eps*|y|) := relative_error (midpoint_relative hn)
  have hyabs : |t.y|≤(1+eps)*|y| := by
    have hh := abs_add_le (t.y-y) y
    simp only [sub_add_cancel] at hh
    have hy' : |t.y-y|≤eps*|y| :=hy
    nlinarith [hy']
  have hyscaled : |t.y|≤(1+eps)*H/2 := by
    have hm := mul_le_mul_of_nonneg_left hyH (by linarith : 0≤1+eps)
    nlinarith
  have tr := factor_ranges hs hn
  have hkrange : |t.k|≤1 := by rw [abs_of_nonneg tr.2.2.2.1]; exact tr.2.2.2.2
  have hk : Within t.k k (2*gamma 3) := factor_error hs hn
  have hleft := within_abs (within_trans
    (relative_error (rounded_relative hn.left)) (within_sub (within_refl f3) hy))
  have hleftabs : |f3-t.y|≤H+(1+eps)*H/2 :=
    (abs_sub _ _).trans (add_le_add hzH hyscaled)
  have hleftbound : eps*|f3-t.y|+(0+eps*|y|)≤(2*eps+eps^2/2)*H := by
    have hm1 := mul_le_mul_of_nonneg_left hleftabs he
    have hm2 := mul_le_mul_of_nonneg_left hyH he
    nlinarith only [hm1,hm2]
  have hsum : Within (|t.y|+|f3|) (|y|+|f3|) (eps*|y|) := by
    simpa using within_add (within_abs hy) (within_refl |f3|)
  have hpert := product_error hkrange hsum hk
    (show |(|y|+|f3|)|≤3*H/2 by rw [abs_of_nonneg (by positivity)]; linarith)
  have hdist := distributed_error
    (a:=t.k*|t.y|) (b:=t.k*|f3|)
    (mul_nonneg tr.2.2.2.1 (abs_nonneg _))
    (mul_nonneg tr.2.2.2.1 (abs_nonneg _))
    (rounded_relative hn.productY) (rounded_relative hn.productZ) (rounded_relative hn.right)
  have hdist' : Within t.right (t.k*(|t.y|+|f3|))
      (gamma 2*(t.k*(|t.y|+|f3|))) := by
    change Within t.right (t.k*|t.y|+t.k*|f3|)
      (gamma 2*(t.k*|t.y|+t.k*|f3|)) at hdist
    simpa only [mul_add] using hdist
  have hsumH : t.k*(|t.y|+|f3|)≤(3+eps)*H/2 := by
    have hmul := mul_le_mul_of_nonneg_right tr.2.2.2.2
      (show 0≤|t.y|+|f3| by positivity)
    nlinarith only [hmul,hyscaled,hzH]
  have hright := within_trans hdist' hpert
  have hrbound : gamma 2*(t.k*(|t.y|+|f3|))+
      (1*(eps*|y|)+(2*gamma 3)*(3*H/2)) ≤
      (((3+eps)/2*gamma 2+3*gamma 3+eps/2)*H) := by
    have hm1 := mul_le_mul_of_nonneg_left hsumH (le_of_lt eps_bounds.2.2.1)
    have hm2 := mul_le_mul_of_nonneg_left hyH he
    nlinarith only [hm1,hm2]
  exact ⟨within_mono hleft hleftbound,within_mono hright hrbound⟩

/-- The full `16 ε H` proposition for the stated concrete rounding sequence.
The no-overflow hypothesis supplies finite-format validity; it is not an error bound. -/
theorem unscaled_error {f1 f2 f3 : ℝ} (hs : f1*f2<0)
    (hn : NormalRun f1 f2 f3) (_hf : NoOverflow f1 f2 f3) :
    let t:=evaluate f1 f2 f3
    let H:=max |f2-f1| |f3|
    let k:=(1-|(f1+f2)/2|/|f2-f1|)^2
    Within (t.right-|t.leftDiff|)
      (originalR k f1 f2 f3-originalL f1 f2 f3) (16*eps*H) := by
  have h := side_errors hs hn
  exact unscaled_error_budget
    ((abs_nonneg (f2-f1)).trans (le_max_left _ _)) h.1 h.2

theorem comparison_agreement {G g m : ℝ} (hm : 0≤m)
    (he : Within g G m) (hmargin : m< |G|) : (0<g ↔ 0<G) := by
  rcases abs_le.mp he with ⟨hlo,hhi⟩
  by_cases hG : 0<G
  · rw [abs_of_pos hG] at hmargin
    constructor
    · exact fun _ => hG
    · intro _; linarith
  · have hGn : G≤0 := le_of_not_gt hG
    rw [abs_of_nonpos hGn] at hmargin
    constructor <;> intro h <;> linarith

theorem unscaled_comparison {f1 f2 f3 : ℝ} (hs : f1*f2<0)
    (hn : NormalRun f1 f2 f3) (hf : NoOverflow f1 f2 f3)
    (hmargin : 16*eps*(max |f2-f1| |f3|) <
      |originalR ((1-|(f1+f2)/2|/|f2-f1|)^2) f1 f2 f3-originalL f1 f2 f3|) :
    |(evaluate f1 f2 f3).leftDiff| < (evaluate f1 f2 f3).right ↔
    originalL f1 f2 f3 < originalR ((1-|(f1+f2)/2|/|f2-f1|)^2) f1 f2 f3 := by
  have heps : 0<eps :=eps_bounds.1
  have hmax : 0≤max |f2-f1| |f3| := (abs_nonneg _).trans (le_max_left _ _)
  have h := comparison_agreement
    (by positivity : 0≤16*eps*(max |f2-f1| |f3|)) (unscaled_error hs hn hf) hmargin
  simpa only [sub_pos] using h

end
end ModAB.Unscaled
