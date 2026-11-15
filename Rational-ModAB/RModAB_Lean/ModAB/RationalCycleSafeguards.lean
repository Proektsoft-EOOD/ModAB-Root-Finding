import ModAB.RationalSecantLocal
import ModAB.RationalCurvature
import ModAB.RationalLocalSafeguards

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
