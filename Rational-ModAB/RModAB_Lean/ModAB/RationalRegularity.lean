import Mathlib.Analysis.Analytic.IsolatedZeros
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas
import Mathlib.Analysis.Calculus.Deriv.Inv
import ModAB.RationalDSlope

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
