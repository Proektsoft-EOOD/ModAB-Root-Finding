import ModAB.RationalSmoothBridge
import ModAB.RationalCycleEntry

/-!
Uniform local coefficients for the negative-coefficient two-step cycle.
The coefficients below are the exact divided-difference expressions. Their
limits are proved from local C² regularity of the removable reciprocal, and
uniform small-cone bounds follow by a sequential contradiction. Thus the
sign and magnitude bounds used by `CycleEntry` are derived, not assumed.
-/
namespace ModAB.Rational.CycleUniform
open Filter Set
open scoped Topology
open LocalAlgebra LocalAnalysis
noncomputable section

def localK (g : ℝ → ℝ) (r : ℝ) : ℝ := (deriv (deriv g) r / 2) / g r

def firstPole (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  poleValue (r+a) (r+s) (r+d) r (g (r+a)) (g (r+s)) (g (r+d))

def firstCoefficient (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  dd2 (r+a) (r+s) (r+d) (g (r+a)) (g (r+s)) (g (r+d)) /
    firstPole g r a s d

def firstError (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  firstCoefficient g r a s d * a * s * d

def ordinateRatio (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  firstCoefficient g r a s d * s * d * (g (r+a) / g (r+firstError g r a s d))

def normalizedGap (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  (1-firstError g r a s d / s) * (1-a/s)

def secondPole (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  let b := firstError g r a s d
  let t := ordinateRatio g r a s d
  poleValue (r+b) (r+a) (r+s) r (g (r+b)) (g (r+a)) (g (r+s)) +
    (t/(1-t)) * g (r+s) * (a/s) * (b/s) / normalizedGap g r a s d

def secondCoefficient (g : ℝ → ℝ) (r a s d : ℝ) : ℝ :=
  let b := firstError g r a s d
  let u := firstCoefficient g r a s d
  let t := ordinateRatio g r a s d
  ((s/d) * dd2 (r+b) (r+a) (r+s) (g (r+b)) (g (r+a)) (g (r+s)) +
    u * g (r+s) * (g (r+a)/g (r+b)) /
      (normalizedGap g r a s d * (1-t))) / secondPole g r a s d

theorem tendsto_firstCoefficient
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstCoefficient g r (a n) (s n) (d n)) atTop (𝓝 (localK g r)) := by
  have h := tendsto_error_coefficient_of_contDiffAt (a:=fun n => r+a n)
    (b:=fun n => r+s n) (c:=fun n => r+d n) hg hgr
    (fun n => by simpa using has n) (fun n => by simpa using hsd n)
    (fun n => by simpa using had n)
    (by simpa using ha.const_add r) (by simpa using hs.const_add r)
    (by simpa using hd.const_add r)
  exact h

theorem tendsto_firstPole
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) := by
  exact tendsto_poleValue_of_contDiffAt hg
    (fun n => by simpa using has n) (fun n => by simpa using hsd n)
    (fun n => by simpa using had n)
    (by simpa using ha.const_add r) (by simpa using hs.const_add r)
    (by simpa using hd.const_add r)

/-- Arbitrarily many coincidences in a finite initial prefix do not affect
the local divided-difference limits. -/
theorem tendsto_divided_differences_eventually
    {g : ℝ → ℝ} {r : ℝ} {a b c : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hdistinct : ∀ᶠ n in atTop, a n ≠ b n ∧ b n ≠ c n ∧ a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) ∧
    Tendsto (fun n => poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g r)) := by
  obtain ⟨N,hN⟩ := eventually_atTop.mp hdistinct
  have hn (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hdd := (tendsto_divided_differences_of_contDiffAt hg
    (fun n => (hn n).1) (fun n => (hn n).2.1) (fun n => (hn n).2.2)
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
    (hc.comp (tendsto_add_atTop_nat N))).2
  have hp := tendsto_poleValue_of_contDiffAt hg
    (fun n => (hn n).1) (fun n => (hn n).2.1) (fun n => (hn n).2.2)
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
    (hc.comp (tendsto_add_atTop_nat N))
  exact ⟨(tendsto_add_atTop_iff_nat N).mp hdd, (tendsto_add_atTop_iff_nat N).mp hp⟩

theorem eventually_ne_of_ratio_zero {a b : ℕ → ℝ}
    (hb : ∀ n, b n ≠ 0) (hlim : Tendsto (fun n => a n / b n) atTop (𝓝 0)) :
    ∀ᶠ n in atTop, a n ≠ b n := by
  filter_upwards [hlim.eventually_ne (by norm_num : (0:ℝ) ≠ 1)] with n hn
  intro h
  apply hn
  simp [h,hb n]

/-- The true first step creates a smaller fresh node and nonzero ordinates
without adding these properties as hypotheses. -/
theorem first_step_limits
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha0 : ∀ n, a n ≠ 0) (hs0 : ∀ n, s n ≠ 0) (hd0 : ∀ n, d n ≠ 0)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0)) :
    Tendsto (fun n => firstError g r (a n) (s n) (d n)) atTop (𝓝 0) ∧
    Tendsto (fun n => firstError g r (a n) (s n) (d n) / a n) atTop (𝓝 0) ∧
    Tendsto (fun n => firstError g r (a n) (s n) (d n) / s n) atTop (𝓝 0) ∧
    (∀ᶠ n in atTop, firstError g r (a n) (s n) (d n) ≠ 0 ∧
      firstError g r (a n) (s n) (d n) ≠ a n ∧
      firstError g r (a n) (s n) (d n) ≠ s n) := by
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  have hb : Tendsto (fun n => firstError g r (a n) (s n) (d n)) atTop (𝓝 0) := by
    simpa only [firstError,mul_zero] using ((hu.mul ha).mul hs).mul hd
  have hba : Tendsto (fun n => firstError g r (a n) (s n) (d n)/a n) atTop (𝓝 0) := by
    have he (n : ℕ) : firstError g r (a n) (s n) (d n)/a n =
        firstCoefficient g r (a n) (s n) (d n)*s n*d n := by
      unfold firstError
      field_simp [ha0 n]
    simp_rw [he]
    simpa only [mul_zero] using (hu.mul hs).mul hd
  have hbs : Tendsto (fun n => firstError g r (a n) (s n) (d n)/s n) atTop (𝓝 0) := by
    have he (n : ℕ) : firstError g r (a n) (s n) (d n)/s n =
        firstCoefficient g r (a n) (s n) (d n)*a n*d n := by
      unfold firstError
      field_simp [hs0 n]
    simp_rw [he]
    simpa only [mul_zero] using (hu.mul ha).mul hd
  refine ⟨hb,hba,hbs,?_⟩
  filter_upwards [hu.eventually_ne hK, eventually_ne_of_ratio_zero ha0 hba,
    eventually_ne_of_ratio_zero hs0 hbs] with n hn hna hns
  exact ⟨mul_ne_zero (mul_ne_zero (mul_ne_zero hn (ha0 n)) (hs0 n)) (hd0 n),hna,hns⟩

/-- The scaled second coefficient has limit K when the two separation ratios
and the radius tend to zero. There is no assumed second error recurrence. -/
theorem second_step_limits
    {g : ℝ → ℝ} {r : ℝ} {a s d : ℕ → ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (has : ∀ n, a n ≠ s n) (hsd : ∀ n, s n ≠ d n) (had : ∀ n, a n ≠ d n)
    (ha0 : ∀ n, a n ≠ 0) (hs0 : ∀ n, s n ≠ 0) (hd0 : ∀ n, d n ≠ 0)
    (ha : Tendsto a atTop (𝓝 0)) (hs : Tendsto s atTop (𝓝 0))
    (hd : Tendsto d atTop (𝓝 0))
    (haslim : Tendsto (fun n => a n/s n) atTop (𝓝 0))
    (hsdlim : Tendsto (fun n => s n/d n) atTop (𝓝 0)) :
    Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 0) ∧
    Tendsto (fun n => normalizedGap g r (a n) (s n) (d n)) atTop (𝓝 1) ∧
    Tendsto (fun n => secondPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) ∧
    Tendsto (fun n => secondCoefficient g r (a n) (s n) (d n)) atTop (𝓝 (localK g r)) := by
  let b := fun n => firstError g r (a n) (s n) (d n)
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  obtain ⟨hb,hba,hbs,hbn⟩ := first_step_limits hg hgr hK has hsd had ha0 hs0 hd0 ha hs hd
  have ha' : Tendsto (fun n => r+a n) atTop (𝓝 r) := by simpa using ha.const_add r
  have hs' : Tendsto (fun n => r+s n) atTop (𝓝 r) := by simpa using hs.const_add r
  have hb' : Tendsto (fun n => r+b n) atTop (𝓝 r) := by simpa using hb.const_add r
  have hga := hg.continuousAt.tendsto.comp ha'
  have hgs := hg.continuousAt.tendsto.comp hs'
  have hgb := hg.continuousAt.tendsto.comp hb'
  have hJ : Tendsto (fun n => g (r+a n)/g (r+b n)) atTop (𝓝 1) := by
    have hh := hga.div hgb hgr
    rw [div_self hgr] at hh
    convert hh using 1
    ext n
    rfl
  have ht : Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 0) := by
    simpa only [ordinateRatio,Function.comp_def,b,mul_zero,zero_mul] using ((hu.mul hs).mul hd).mul hJ
  have hq : Tendsto (fun n => normalizedGap g r (a n) (s n) (d n)) atTop (𝓝 1) := by
    simpa only [normalizedGap,sub_zero,mul_one] using
      ((tendsto_const_nhds (x:=(1:ℝ))).sub hbs).mul ((tendsto_const_nhds (x:=(1:ℝ))).sub haslim)
  have hdist : ∀ᶠ n in atTop, r+b n ≠ r+a n ∧ r+a n ≠ r+s n ∧ r+b n ≠ r+s n := by
    filter_upwards [hbn] with n hn
    exact ⟨by simpa using hn.2.1,by simpa using has n,by simpa using hn.2.2⟩
  obtain ⟨hG,hP⟩ := tendsto_divided_differences_eventually hg hdist hb' ha' hs'
  have hμ : Tendsto (fun n => 1-ordinateRatio g r (a n) (s n) (d n)) atTop (𝓝 1) := by
    simpa using tendsto_const_nhds.sub ht
  have hδ : Tendsto (fun n => ordinateRatio g r (a n) (s n) (d n) /
      (1-ordinateRatio g r (a n) (s n) (d n))) atTop (𝓝 0) := by
    have hh := ht.div hμ (by norm_num : (1:ℝ) ≠ 0)
    rw [zero_div] at hh
    convert hh using 1
  have hcorr := (((hδ.mul hgs).mul haslim).mul hbs).div hq (by norm_num : (1:ℝ) ≠ 0)
  have hP2 : Tendsto (fun n => secondPole g r (a n) (s n) (d n)) atTop (𝓝 (g r)) := by
    simpa only [secondPole,b,Function.comp_def,Pi.div_apply,mul_zero,zero_mul,zero_div,add_zero] using hP.add hcorr
  have hdom := ((hu.mul hgs).mul hJ).div (hq.mul hμ) (by norm_num : (1:ℝ)*1 ≠ 0)
  have hV := ((hsdlim.mul hG).add hdom).div hP2 hgr
  refine ⟨ht,hq,hP2,?_⟩
  simp only [mul_one,zero_mul,zero_add,div_one,mul_div_cancel_right₀ _ hgr] at hV
  convert hV using 1
  ext n
  rfl


/-- A concrete positive null sequence, used only to extract uniform cones. -/
def smallRadius (n : ℕ) : ℝ := 1 / ((n : ℝ)+2)

theorem smallRadius_pos (n : ℕ) : 0 < smallRadius n := by
  unfold smallRadius
  positivity

theorem smallRadius_lt_one (n : ℕ) : smallRadius n < 1 := by
  unfold smallRadius
  apply (div_lt_one (by positivity : (0:ℝ) < (n:ℝ)+2)).mpr
  nlinarith [Nat.cast_nonneg (α:=ℝ) n]

theorem smallRadius_tendsto : Tendsto smallRadius atTop (𝓝 0) := by
  have h : Tendsto (fun n : ℕ => (n:ℝ)+2) atTop atTop :=
    tendsto_natCast_atTop_atTop.atTop_add tendsto_const_nhds
  change Tendsto (fun n : ℕ => 1/((n:ℝ)+2)) atTop (𝓝 0)
  simpa only [one_div,Function.comp_def] using tendsto_inv_atTop_zero.comp h

/-- Sequential extraction of a single uniform neighbourhood; the witness
triples are selected from the failure of uniformity, not supplied externally. -/
theorem exists_uniform_cone_of_sequential {P : ℝ → ℝ → ℝ → Prop}
    (hseq : ∀ a s d : ℕ → ℝ,
      (∀ n, CycleEntry.Cone (smallRadius n) (smallRadius n) (a n) (s n) (d n)) →
      ∀ᶠ n in atTop, P (a n) (s n) (d n)) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d → P a s d := by
  by_contra h
  push Not at h
  choose a s d hc hp using fun n => h (smallRadius n) (smallRadius n)
    (smallRadius_pos n) (smallRadius_lt_one n) (smallRadius_pos n)
  obtain ⟨n,hn⟩ := (hseq a s d hc).exists
  exact hp n hn

theorem Cone.distinct {θ ε a s d : ℝ} (h : CycleEntry.Cone θ ε a s d)
    (hθ : θ < 1) : a ≠ s ∧ s ≠ d ∧ a ≠ d := by
  have has : a ≠ s := by
    intro he
    have hh := h.opposite
    rw [he] at hh
    nlinarith [sq_nonneg s]
  have hsd : s ≠ d := by
    intro he
    have hh := h.opposite
    rw [he] at hh
    linarith [h.same]
  have had : a ≠ d := by
    intro he
    have horder := h.ordered_errors hθ.le
    have hsmall := h.opposite_small
    rw [he] at horder
    have hp := abs_pos.mpr h.history_ne_zero
    nlinarith
  exact ⟨has,hsd,had⟩

/-- A shrinking cone forces both independent separation ratios to zero. -/
theorem shrinking_cone_limits {a s d : ℕ → ℝ}
    (hc : ∀ n, CycleEntry.Cone (smallRadius n) (smallRadius n) (a n) (s n) (d n)) :
    Tendsto a atTop (𝓝 0) ∧ Tendsto s atTop (𝓝 0) ∧ Tendsto d atTop (𝓝 0) ∧
    Tendsto (fun n => a n/s n) atTop (𝓝 0) ∧
    Tendsto (fun n => s n/d n) atTop (𝓝 0) := by
  have hord (n : ℕ) := (hc n).ordered_errors (smallRadius_lt_one n).le
  have hdabs := squeeze_zero (fun n => abs_nonneg (d n))
    (fun n => (hc n).history_bound) smallRadius_tendsto
  have hsabs := squeeze_zero (fun n => abs_nonneg (s n)) (fun n => (hord n).2) hdabs
  have haabs := squeeze_zero (fun n => abs_nonneg (a n)) (fun n => (hord n).1) hsabs
  have hratio₁ : ∀ n, |a n/s n| ≤ smallRadius n := by
    intro n
    rw [abs_div]
    exact (div_le_iff₀ (abs_pos.mpr (hc n).opposite_ne_zero)).mpr (hc n).latest_small
  have hratio₂ : ∀ n, |s n/d n| ≤ smallRadius n := by
    intro n
    rw [abs_div]
    exact (div_le_iff₀ (abs_pos.mpr (hc n).history_ne_zero)).mpr (hc n).opposite_small
  exact ⟨(tendsto_zero_iff_abs_tendsto_zero _).mpr haabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr hsabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr hdabs,
    (tendsto_zero_iff_abs_tendsto_zero _).mpr
      (squeeze_zero (fun n => abs_nonneg _) hratio₁ smallRadius_tendsto),
    (tendsto_zero_iff_abs_tendsto_zero _).mpr
      (squeeze_zero (fun n => abs_nonneg _) hratio₂ smallRadius_tendsto)⟩

/-- The uniform signed coefficient bounds needed for the invariant-cycle
proof follow from local regularity. The same neighbourhood excludes every
singular denominator in the normalized pair formulas. -/
theorem exists_uniform_negative_coefficients
    {g : ℝ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r < 0) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        firstCoefficient g r a s d ∈ Ioo (2*localK g r) (localK g r/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*localK g r) (localK g r/2) ∧
        firstPole g r a s d ≠ 0 ∧ secondPole g r a s d ≠ 0 ∧
        g (r+a) ≠ 0 ∧ g (r+s) ≠ 0 ∧ g (r+d) ≠ 0 ∧
        g (r+firstError g r a s d) ≠ 0 ∧
        firstError g r a s d ≠ 0 ∧ firstError g r a s d ≠ a ∧
        firstError g r a s d ≠ s ∧
        0 < normalizedGap g r a s d ∧ ordinateRatio g r a s d < 1 ∧
        0 < g (r+a) / g (r+firstError g r a s d) := by
  apply exists_uniform_cone_of_sequential
  intro a s d hc
  obtain ⟨ha,hs,hd,haslim,hsdlim⟩ := shrinking_cone_limits hc
  have hdist (n : ℕ) := Cone.distinct (hc n) (smallRadius_lt_one n)
  have has := fun n => (hdist n).1
  have hsd := fun n => (hdist n).2.1
  have had := fun n => (hdist n).2.2
  have ha0 := fun n => (hc n).latest_ne_zero
  have hs0 := fun n => (hc n).opposite_ne_zero
  have hd0 := fun n => (hc n).history_ne_zero
  have hu := tendsto_firstCoefficient hg hgr has hsd had ha hs hd
  have hp := tendsto_firstPole hg has hsd had ha hs hd
  obtain ⟨hb,hba,hbs,hbn⟩ := first_step_limits hg hgr hK.ne has hsd had ha0 hs0 hd0 ha hs hd
  obtain ⟨ht,hq,hp₂,hv⟩ := second_step_limits hg hgr hK.ne has hsd had
    ha0 hs0 hd0 ha hs hd haslim hsdlim
  have hga : Tendsto (fun n => g (r+a n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using ha.const_add r)
  have hgs : Tendsto (fun n => g (r+s n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hs.const_add r)
  have hgd : Tendsto (fun n => g (r+d n)) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hd.const_add r)
  have hgb : Tendsto (fun n => g (r+firstError g r (a n) (s n) (d n))) atTop (𝓝 (g r)) :=
    hg.continuousAt.tendsto.comp (by simpa using hb.const_add r)
  have hJ : Tendsto (fun n => g (r+a n)/g (r+firstError g r (a n) (s n) (d n))) atTop (𝓝 1) := by
    have hh := hga.div hgb hgr
    rw [div_self hgr] at hh
    convert hh using 1
  have hI : Ioo (2*localK g r) (localK g r/2) ∈ 𝓝 (localK g r) :=
    Ioo_mem_nhds (by linarith) (by linarith)
  filter_upwards [hu.eventually hI,hv.eventually hI,hp.eventually_ne hgr,
    hp₂.eventually_ne hgr,hga.eventually_ne hgr,hgs.eventually_ne hgr,
    hgd.eventually_ne hgr,hgb.eventually_ne hgr,hbn,
    (tendsto_order.1 hq).1 0 (by norm_num),
    (tendsto_order.1 ht).2 1 (by norm_num),
    (tendsto_order.1 hJ).1 0 (by norm_num)] with n hu hv hp hp₂ hga hgs hgd hgb hb hq ht hJ
  exact ⟨hu,hv,hp,hp₂,hga,hgs,hgd,hgb,hb.1,hb.2.1,hb.2.2,hq,ht,hJ⟩


/-- The coefficient in the uniform cone theorem is exactly the paper's
C³ coefficient, not a separately chosen constant. -/
theorem localK_regularizedReciprocal
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    localK (Regularity.regularizedReciprocal f r) r = AnalyticBridge.errorConstant f r := by
  have h := Regularity.regularizedReciprocal_error_constant_of_contDiff hf hzero hd
  have hder : iteratedDeriv 2 (Regularity.regularizedReciprocal f r) r =
      deriv (deriv (Regularity.regularizedReciprocal f r)) r := by
    rw [show 2 = 1+1 from rfl, iteratedDeriv_succ', iteratedDeriv_one]
  rw [hder] at h
  simpa only [localK,AnalyticBridge.errorConstant] using h

/-- Full finite-C³ uniform-cone theorem for the actual removable reciprocal. -/
theorem smooth_exists_uniform_negative_coefficients
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hd : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r < 0) :
    let g := Regularity.regularizedReciprocal f r
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        firstCoefficient g r a s d ∈ Ioo (2*AnalyticBridge.errorConstant f r) (AnalyticBridge.errorConstant f r/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*AnalyticBridge.errorConstant f r) (AnalyticBridge.errorConstant f r/2) ∧
        firstPole g r a s d ≠ 0 ∧ secondPole g r a s d ≠ 0 ∧
        g (r+a) ≠ 0 ∧ g (r+s) ≠ 0 ∧ g (r+d) ≠ 0 ∧
        g (r+firstError g r a s d) ≠ 0 ∧
        firstError g r a s d ≠ 0 ∧ firstError g r a s d ≠ a ∧
        firstError g r a s d ≠ s ∧
        0 < normalizedGap g r a s d ∧ ordinateRatio g r a s d < 1 ∧
        0 < g (r+a) / g (r+firstError g r a s d) := by
  dsimp only
  have hg := Regularity.regularizedReciprocal_contDiffAt_of_contDiff hf hd
  have hgr : Regularity.regularizedReciprocal f r r ≠ 0 := by
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hk := localK_regularizedReciprocal hf hzero hd
  have result := exists_uniform_negative_coefficients hg hgr (by rw [hk]; exact hK)
  simpa only [hk] using result

end
end ModAB.Rational.CycleUniform
