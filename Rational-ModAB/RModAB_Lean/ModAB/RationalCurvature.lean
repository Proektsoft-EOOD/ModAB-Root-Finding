import ModAB.RationalLocalAnalysis

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
