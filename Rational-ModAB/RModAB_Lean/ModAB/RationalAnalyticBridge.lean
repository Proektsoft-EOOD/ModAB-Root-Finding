import ModAB.RationalLocalAnalysis
import ModAB.RationalRegularity

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
