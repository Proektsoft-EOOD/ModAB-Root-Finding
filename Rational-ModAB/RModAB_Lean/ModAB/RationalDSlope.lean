import Mathlib.Analysis.Calculus.DSlope
import Mathlib.Analysis.Calculus.LHopital
import Mathlib.Analysis.Calculus.ContDiff.Deriv
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.Deriv.Pow

/-! Finite-smoothness division at a removable zero, by explicit derivative
formulas and L'Hopital limits. No analyticity or assumed smooth factor is used. -/
namespace ModAB.Rational.DSlope
open Set Filter Topology
open scoped ContDiff
noncomputable section

theorem power_quotient_limit {N F : ℝ → ℝ} {r : ℝ} {k : ℕ}
    (hk : 0 < k) (hN : ContinuousAt N r) (hzero : N r = 0)
    (hF : ContinuousAt F r)
    (hder : ∀ᶠ x in 𝓝[≠] r, HasDerivAt N ((x-r)^(k-1)*F x) x) :
    Tendsto (fun x => N x/(x-r)^k) (𝓝[≠] r) (𝓝 (F r/k)) := by
  have hkn : (k:ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hk)
  have hx : ∀ᶠ x : ℝ in 𝓝[≠] r, x-r ≠ 0 := by
    filter_upwards [self_mem_nhdsWithin] with x hx
    exact sub_ne_zero.mpr hx
  have hdg : ∀ᶠ x : ℝ in 𝓝[≠] r,
      HasDerivAt (fun y : ℝ => (y-r)^k) ((k:ℝ)*(x-r)^(k-1)) x := by
    exact Filter.Eventually.of_forall (fun x => by
      simpa only [mul_one, id_eq] using! ((hasDerivAt_id x).sub_const r).fun_pow k)
  apply HasDerivAt.lhopital_zero_nhdsNE hder hdg
  · filter_upwards [hx] with x hx
    simpa using mul_ne_zero hkn (pow_ne_zero _ hx)
  · simpa [hzero] using hN.tendsto.mono_left nhdsWithin_le_nhds
  · have hc : ContinuousAt (fun x : ℝ => (x-r)^k) r := by fun_prop
    simpa [Nat.ne_of_gt hk] using hc.tendsto.mono_left nhdsWithin_le_nhds
  · have hlim : Tendsto (fun x => F x/(k:ℝ)) (𝓝[≠] r) (𝓝 (F r/(k:ℝ))) :=
      (hF.tendsto.div_const (k:ℝ)).mono_left nhdsWithin_le_nhds
    apply hlim.congr'
    filter_upwards [hx] with x hx
    field_simp [hx, hkn]

def firstFormula (f : ℝ → ℝ) (r x : ℝ) : ℝ :=
  ((x-r)*deriv f x-(f x-f r))/(x-r)^2

def secondFormula (f : ℝ → ℝ) (r x : ℝ) : ℝ :=
  ((x-r)^2*deriv (deriv f) x-2*(x-r)*deriv f x+2*(f x-f r))/(x-r)^3

def first (f : ℝ → ℝ) (r : ℝ) : ℝ → ℝ :=
  Function.update (firstFormula f r) r (deriv (deriv f) r/2)

def second (f : ℝ → ℝ) (r : ℝ) : ℝ → ℝ :=
  Function.update (secondFormula f r) r (deriv (deriv (deriv f)) r/3)

theorem dslope_hasDerivAt_of_ne {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (hx : x ≠ r) :
    HasDerivAt (dslope f r) (first f r x) x := by
  have hn := (hf.hasDerivAt.sub_const (f r)).div
    ((hasDerivAt_id x).sub_const r) (sub_ne_zero.mpr hx)
  have hs : HasDerivAt (slope f r) (firstFormula f r x) x := by
    have he : slope f r = fun y => (f y-f r)/(y-r) := by
      funext y
      simp [slope_def_field]
    rw [he]
    simpa [firstFormula, mul_comm] using! hn
  simpa [first, Function.update_of_ne hx] using
    hs.congr_of_eventuallyEq (dslope_eventuallyEq_slope_of_ne f hx)

theorem dslope_hasDerivAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 f r) :
    HasDerivAt (dslope f r) (deriv (deriv f) r/2) r := by
  have h1 : ContDiffAt ℝ 1 (deriv f) r := hf.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => f x-f r-deriv f r*(x-r)) r := by
    exact (hf.continuousAt.sub continuousAt_const).sub
      (continuousAt_const.mul (continuousAt_id.sub continuousAt_const))
  have hF : ContinuousAt (dslope (deriv f) r) r :=
    continuousAt_dslope_same.mpr (h1.differentiableAt (by norm_num))
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun x => f x-f r-deriv f r*(x-r))
        ((x-r)^(2-1)*dslope (deriv f) r x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hd := ((hx.differentiableAt (by norm_num)).hasDerivAt.sub_const (f r)).sub
      (((hasDerivAt_id x).sub_const r).const_mul (deriv f r))
    convert hd using 1
    · rfl
    · simpa [smul_eq_mul] using sub_smul_dslope (deriv f) r x
  have hl := power_quotient_limit (k:=2) (by norm_num) hN (by simp) hF hder
  rw [dslope_same] at hl
  apply hasDerivAt_iff_tendsto_slope.mpr
  apply hl.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxr : x-r ≠ 0 := sub_ne_zero.mpr hx
  simp only [slope_def_field, dslope_same, dslope_of_ne f hx]
  field_simp [hxr]

theorem first_numerator_derivative {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (h1 : DifferentiableAt ℝ (deriv f) x) :
    HasDerivAt (fun y => (y-r)*deriv f y-(f y-f r))
      ((x-r)*deriv (deriv f) x) x := by
  convert (((hasDerivAt_id x).sub_const r).mul h1.hasDerivAt).sub
    (hf.hasDerivAt.sub_const (f r)) using 1
  · rfl
  · simp only [id_eq]
    ring

theorem first_continuousAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 2 f r) : ContinuousAt (first f r) r := by
  have h1 : ContDiffAt ℝ 1 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 0 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => (x-r)*deriv f x-(f x-f r)) r := by
    exact ((continuousAt_id.sub continuousAt_const).mul h1.continuousAt).sub
      (hf.continuousAt.sub continuousAt_const)
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun y => (y-r)*deriv f y-(f y-f r))
        ((x-r)^(2-1)*deriv (deriv f) x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 1 (deriv f) x := hx.derivWithin (by norm_num)
    simpa using first_numerator_derivative
      (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
  exact continuousAt_update_same.mpr
    (power_quotient_limit (k:=2) (by norm_num) hN (by simp) h2.continuousAt hder)

theorem first_hasDerivAt_of_ne {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 2 f x) (hx : x ≠ r) :
    HasDerivAt (first f r) (second f r x) x := by
  have h1 : ContDiffAt ℝ 1 (deriv f) x := hf.derivWithin (by norm_num)
  have hxr := sub_ne_zero.mpr hx
  have hN := first_numerator_derivative
    (r:=r) (hf.differentiableAt (by norm_num)) (h1.differentiableAt (by norm_num))
  have hQ := hN.div (((hasDerivAt_id x).sub_const r).fun_pow 2) (pow_ne_zero 2 hxr)
  have hformula : HasDerivAt (firstFormula f r) (secondFormula f r x) x := by
    convert hQ using 1
    · rfl
    · simp only [secondFormula, id_eq]
      field_simp [hxr]
      ring
  have he : first f r =ᶠ[𝓝 x] firstFormula f r := by
    filter_upwards [eventually_ne_nhds hx] with y hy
    simp [first, hy]
  simpa [second, hx] using hformula.congr_of_eventuallyEq he

theorem first_hasDerivAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) :
    HasDerivAt (first f r) (deriv (deriv (deriv f)) r/3) r := by
  have h1 : ContDiffAt ℝ 2 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have hF : ContinuousAt (dslope (deriv (deriv f)) r) r :=
    continuousAt_dslope_same.mpr (h2.differentiableAt (by norm_num))
  have hN : ContinuousAt (fun x => (x-r)*deriv f x-(f x-f r)-
      (deriv (deriv f) r/2)*(x-r)^2) r := by
    fun_prop
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun x => (x-r)*deriv f x-(f x-f r)-
        (deriv (deriv f) r/2)*(x-r)^2)
      ((x-r)^(3-1)*dslope (deriv (deriv f)) r x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 2 (deriv f) x := hx.derivWithin (by norm_num)
    have hN := first_numerator_derivative
      (r:=r) (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
    have hD := hN.sub ((((hasDerivAt_id x).sub_const r).fun_pow 2).const_mul
      (deriv (deriv f) r/2))
    convert hD using 1
    · rfl
    · have he := sub_smul_dslope (deriv (deriv f)) r x
      simp only [smul_eq_mul] at he
      have he2 := congrArg (fun z : ℝ => (x-r)*z) he
      norm_num [id_eq]
      nlinarith [he2]
  have hl := power_quotient_limit (k:=3) (by norm_num) hN (by simp) hF hder
  rw [dslope_same] at hl
  apply hasDerivAt_iff_tendsto_slope.mpr
  apply hl.congr'
  filter_upwards [self_mem_nhdsWithin] with x hx
  have hxr : x-r ≠ 0 := sub_ne_zero.mpr hx
  simp only [slope_def_field, first, Function.update_self, Function.update_of_ne hx,
    firstFormula]
  field_simp [hxr]

theorem second_numerator_derivative {f : ℝ → ℝ} {r x : ℝ}
    (hf : DifferentiableAt ℝ f x) (h1 : DifferentiableAt ℝ (deriv f) x)
    (h2 : DifferentiableAt ℝ (deriv (deriv f)) x) :
    HasDerivAt (fun y => (y-r)^2*deriv (deriv f) y-
      2*(y-r)*deriv f y+2*(f y-f r))
      ((x-r)^2*deriv (deriv (deriv f)) x) x := by
  convert ((((hasDerivAt_id x).sub_const r).fun_pow 2).mul h2.hasDerivAt).sub
    ((((hasDerivAt_id x).sub_const r).const_mul 2).mul h1.hasDerivAt) |>.add
      ((hf.hasDerivAt.sub_const (f r)).const_mul 2) using 1
  · rfl
  · simp only [id_eq]
    ring

theorem second_continuousAt_same {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) : ContinuousAt (second f r) r := by
  have h1 : ContDiffAt ℝ 2 (deriv f) r := hf.derivWithin (by norm_num)
  have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) r := h1.derivWithin (by norm_num)
  have h3 : ContDiffAt ℝ 0 (deriv (deriv (deriv f))) r := h2.derivWithin (by norm_num)
  have hN : ContinuousAt (fun x => (x-r)^2*deriv (deriv f) x-
      2*(x-r)*deriv f x+2*(f x-f r)) r := by
    fun_prop
  have hder : ∀ᶠ x in 𝓝[≠] r,
      HasDerivAt (fun y => (y-r)^2*deriv (deriv f) y-
        2*(y-r)*deriv f y+2*(f y-f r))
        ((x-r)^(3-1)*deriv (deriv (deriv f)) x) x := by
    filter_upwards [(hf.eventually (by norm_num)).filter_mono nhdsWithin_le_nhds] with x hx
    have hdx : ContDiffAt ℝ 2 (deriv f) x := hx.derivWithin (by norm_num)
    have hddx : ContDiffAt ℝ 1 (deriv (deriv f)) x := hdx.derivWithin (by norm_num)
    simpa using second_numerator_derivative
      (hx.differentiableAt (by norm_num)) (hdx.differentiableAt (by norm_num))
      (hddx.differentiableAt (by norm_num))
  exact continuousAt_update_same.mpr
    (power_quotient_limit (k:=3) (by norm_num) hN (by simp) h3.continuousAt hder)

theorem dslope_hasDerivAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 2 f x) : HasDerivAt (dslope f r) (first f r x) x := by
  by_cases hx : x=r
  · subst x
    simpa [first] using dslope_hasDerivAt_same hf
  · exact dslope_hasDerivAt_of_ne (hf.differentiableAt (by norm_num)) hx

theorem first_hasDerivAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 3 f x) : HasDerivAt (first f r) (second f r x) x := by
  by_cases hx : x=r
  · subst x
    simpa [second] using first_hasDerivAt_same hf
  · exact first_hasDerivAt_of_ne (hf.of_le (by norm_num)) hx

theorem second_continuousAt {f : ℝ → ℝ} {r x : ℝ}
    (hf : ContDiffAt ℝ 3 f x) : ContinuousAt (second f r) x := by
  by_cases hx : x=r
  · subst x
    exact second_continuousAt_same hf
  · have h1 : ContDiffAt ℝ 2 (deriv f) x := hf.derivWithin (by norm_num)
    have h2 : ContDiffAt ℝ 1 (deriv (deriv f)) x := h1.derivWithin (by norm_num)
    have hc : ContinuousAt (secondFormula f r) x := by
      unfold secondFormula
      apply ContinuousAt.div
      · fun_prop
      · fun_prop
      · exact pow_ne_zero 3 (sub_ne_zero.mpr hx)
    apply hc.congr_of_eventuallyEq
    filter_upwards [eventually_ne_nhds hx] with y hy
    simp [second, hy]

/-- Division by a simple coordinate factor loses exactly one finite derivative.
This theorem only assumes local `C³` regularity, not analyticity. -/
theorem contDiffAt_dslope {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) : ContDiffAt ℝ 2 (dslope f r) r := by
  have hev := hf.eventually (by norm_num)
  obtain ⟨u, hu, hopen, hr⟩ := mem_nhds_iff.mp hev
  have hgood : ∀ x ∈ u, ContDiffAt ℝ 3 f x := fun x hx => hu hx
  have hder0 : ∀ x ∈ u, HasDerivAt (dslope f r) (first f r x) x :=
    fun x hx => dslope_hasDerivAt ((hgood x hx).of_le (by norm_num))
  have hder1 : ∀ x ∈ u, HasDerivAt (first f r) (second f r x) x :=
    fun x hx => first_hasDerivAt (hgood x hx)
  have hcont2 : ContinuousOn (second f r) u :=
    fun x hx => (second_continuousAt (hgood x hx)).continuousWithinAt
  have hc1 : ContDiffOn ℝ 1 (first f r) u := by
    rw [show (1:ℕ∞ω)=0+1 by rfl, contDiffOn_succ_iff_deriv_of_isOpen hopen]
    refine ⟨fun x hx => (hder1 x hx).differentiableAt.differentiableWithinAt, by norm_num, ?_⟩
    exact contDiffOn_zero.mpr (hcont2.congr (fun x hx => (hder1 x hx).deriv))
  have hc0 : ContDiffOn ℝ 2 (dslope f r) u := by
    rw [show (2:ℕ∞ω)=1+1 by norm_num, contDiffOn_succ_iff_deriv_of_isOpen hopen]
    refine ⟨fun x hx => (hder0 x hx).differentiableAt.differentiableWithinAt, by norm_num, ?_⟩
    exact hc1.congr (fun x hx => (hder0 x hx).deriv)
  exact hc0.contDiffAt (hopen.mem_nhds hr)

/-- First coefficient of the removable divided difference. -/
theorem deriv_dslope_at {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 2 f r) :
    deriv (dslope f r) r = deriv (deriv f) r/2 :=
  (dslope_hasDerivAt_same hf).deriv

/-- Second coefficient of the removable divided difference. -/
theorem second_deriv_dslope_at {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r) :
    deriv (deriv (dslope f r)) r = deriv (deriv (deriv f)) r/3 := by
  have he : deriv (dslope f r) =ᶠ[𝓝 r] first f r := by
    filter_upwards [hf.eventually (by norm_num)] with x hx
    exact (dslope_hasDerivAt (r:=r) (hx.of_le (by norm_num))).deriv
  rw [he.deriv_eq]
  exact (first_hasDerivAt_same hf).deriv

end
end ModAB.Rational.DSlope
