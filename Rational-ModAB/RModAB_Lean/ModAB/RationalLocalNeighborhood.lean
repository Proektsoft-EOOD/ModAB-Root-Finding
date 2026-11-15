import ModAB.RationalSmoothBridge
import Mathlib.Analysis.SpecificLimits.Basic

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
