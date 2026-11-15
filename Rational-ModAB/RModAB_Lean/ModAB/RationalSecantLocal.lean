import ModAB.RationalRegularity
import ModAB.Rational

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
