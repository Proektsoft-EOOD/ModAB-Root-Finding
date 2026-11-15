import ModAB.RationalOneSidedAssembly
import ModAB.RationalRegularity
import ModAB.RationalDSlope
import ModAB.TopEntry

/-!
The bounded-scaling premise of the earlier connected one-sided theorem is
derived here from the actual raw rational update and its no-crossing sign.
The bounding cap has a finite limit by local C¹ regularity. No convergence
rate, summability of residual ratios, or lower bound on the accumulated weight
is assumed. The hypothesis of an unmodified, convergent, nondegenerate
one-sided tail is retained, exactly as in Proposition 12.
-/
namespace ModAB.Rational.OneSidedComplete
open Filter Set
open scoped Topology
open LocalAlgebra OneSided OneSidedAssembly
noncomputable section

/-- An increasing quantity bounded on a tail is bounded everywhere by the
same constant. This removes any finite-prefix issue in the weight bound. -/
theorem monotone_bound_of_eventually {d : ℕ → ℝ} {M : ℝ}
    (hd : Monotone d) (hM : ∀ᶠ n in atTop, d n ≤ M) : ∀ n, d n ≤ M := by
  obtain ⟨N, hN⟩ := eventually_atTop.mp hM
  intro n
  exact (hd (le_max_left n N)).trans (hN (max n N) (le_max_right n N))

/-- The upper bound extracted from the no-crossing sign has a finite limit
for two coalescing fresh nodes and a fixed opposite endpoint. -/
theorem tendsto_no_cross_cap
    {g : ℝ → ℝ} {a b : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n ≠ b n)
    (hs : s ≠ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => -((s-a n)*(s-b n)) *
      dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s) / g s)
      atTop (𝓝 (-((s-r)*(s-r))*confluentDD2 g r s/g s)) := by
  exact (((tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)).neg.mul
    (tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb)).div_const (g s)

/-- Sign-invariant version of the scaling bound: the pole value and fixed
reciprocal value need only have the same nonzero sign. -/
theorem delta_bound_of_same_side_same_sign {E G P g δ π : ℝ}
    (hE : 0 < E) (hPg : 0 < P*g) (hδ : 0 ≤ δ) (hπ : 0 < π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π)) < 0) :
    δ < -π*G/g := by
  rcases mul_pos_iff.mp hPg with hpos | hneg
  · exact delta_bound_of_same_side hE hpos.1 hpos.2 hδ hπ hsame
  · have hsame' : E*(-G+δ*(-g)/π)/(-P+E*(δ*(-g)/π)) < 0 := by
      have hnum : -G+δ*(-g)/π = -(G+δ*g/π) := by ring
      have hden : -P+E*(δ*(-g)/π) = -(P+E*(δ*g/π)) := by ring
      rw [hnum,hden,mul_neg,neg_div_neg_eq]
      exact hsame
    have h := delta_bound_of_same_side hE (neg_pos.mpr hneg.1)
      (neg_pos.mpr hneg.2) hδ hπ hsame'
    convert h using 1
    ring

/-- Exact no-crossing data give the pointwise scaling bound. The sign of the
divided difference is a consequence, not an assumption. -/
theorem delta_bound_of_barycentric_no_cross
    {f g : ℝ → ℝ} {a b s r μ y : ℝ}
    (hba : b-a ≠ 0) (hsb : s-b ≠ 0) (hsa : s-a ≠ 0)
    (har : a-r ≠ 0) (hbr : b-r ≠ 0) (hsr : s-r ≠ 0)
    (hfs : f s ≠ 0) (hμ : 0 < μ) (hμone : μ ≤ 1)
    (hquotA : (a-r)/f a = g a) (hquotB : (b-r)/f b = g b)
    (hquotS : (s-r)/f s = g s)
    (hH : dd2 a b s (1/f a) (1/f b) (1/(μ*f s)) ≠ 0)
    (hstep : y = baryRoot a b s (f a) (f b) (μ*f s))
    (hE : 0 < (a-r)*(b-r))
    (hP : 0 < poleValue a b s r (g a) (g b) (g s) * g s)
    (hπ : 0 < (s-a)*(s-b))
    (hsame : (y-r)*(s-r) < 0) :
    1/μ-1 < -((s-a)*(s-b))*dd2 a b s (g a) (g b) (g s)/g s := by
  have hδ : 0 ≤ 1/μ-1 := by
    have h : 1 ≤ 1/μ := by
      simpa only [one_div_one] using one_div_le_one_div_of_le hμ hμone
    linarith
  let E := (a-r)*(b-r)
  let G := dd2 a b s (g a) (g b) (g s)
  let P := poleValue a b s r (g a) (g b) (g s)
  let π := (s-a)*(s-b)
  let δ := 1/μ-1
  let z := E*(G+δ*g s/π)/(P+E*(δ*g s/π))
  have hid : y-r = z*(s-r) := by
    rw [hstep, baryRoot_error_eq_mixedCoefficient hba hsb hsa har hbr hsr
      hfs (ne_of_gt hμ) hH hquotA hquotB hquotS]
    dsimp [mixedCoefficient,z,E,G,P,π,δ]
    have heq : (1/μ-1)*g s*(a-r)*(b-r)/((s-a)*(s-b)) =
        ((a-r)*(b-r))*((1/μ-1)*g s/((s-a)*(s-b))) := by ring
    rw [heq]
    ring
  have hz : z < 0 := by
    rw [hid] at hsame
    by_contra hn
    have hn0 : 0 ≤ z := le_of_not_gt hn
    have hnonneg := mul_nonneg hn0 (mul_self_nonneg (s-r))
    nlinarith
  exact delta_bound_of_same_side_same_sign hE hP hδ hπ hz

/-- The actual source-oriented rawPoint tail supplies a global finite bound
on accumulated scaling. No bound on scaling is assumed. -/
theorem exists_scaling_bound_of_rawPoint_one_sided
    {f g : ℝ → ℝ} {x μ : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgrgs : 0 < g r * g s) (hs : s ≠ r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n ≠ r)
    (hadj : ∀ n, x (n+1) ≠ x n) (hxs : ∀ n, x n ≠ s)
    (hμ : ∀ n, 0 < μ n) (hμone : ∀ n, μ n ≤ 1) (hanti : Antitone μ)
    (hfs : f s ≠ 0) (hfx : ∀ n, f (x n) ≠ 0)
    (hquot : ∀ n, (x n-r)/f (x n) = g (x n))
    (hquots : (s-r)/f s = g s)
    (hside : ∀ n, (x n-r)*(s-r) < 0)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0)
    (hstep : ∀ᶠ n in atTop, x (n+2) =
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n))) :
    ∃ M : ℝ, ∀ n, 1/μ n-1 ≤ M := by
  have hxnext : Tendsto (fun n => x (n+1)) atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 1).mpr hx
  have hpole := tendsto_poleValue_fixed_endpoint_of_contDiffAt hg hadj hs hxnext hx
  have hP : ∀ᶠ n in atTop,
      0 < poleValue (x (n+1)) (x n) s r (g (x (n+1))) (g (x n)) (g s) * g s :=
    (hpole.mul_const (g s)).eventually (Ioi_mem_nhds hgrgs)
  have hπlim : Tendsto (fun n => (s-x (n+1))*(s-x n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub hxnext).mul (tendsto_const_nhds.sub hx)
  have hsr : s-r ≠ 0 := sub_ne_zero.mpr hs
  have hπ : ∀ᶠ n in atTop, 0 < (s-x (n+1))*(s-x n) :=
    hπlim.eventually (Ioi_mem_nhds (mul_self_pos.mpr hsr))
  let M := -((s-r)*(s-r))*confluentDD2 g r s/g s + 1
  have hcap := tendsto_no_cross_cap hg hadj hs hxnext hx
  have hcapM : ∀ᶠ n in atTop, -((s-x (n+1))*(s-x n))*
      dd2 (x (n+1)) (x n) s (g (x (n+1))) (g (x n)) (g s) / g s < M :=
    hcap.eventually (Iio_mem_nhds (show
      -((s-r)*(s-r))*confluentDD2 g r s/g s < M by dsimp [M]; linarith))
  refine ⟨M, monotone_bound_of_eventually
    (delta_monotone_of_weight_antitone hμ hanti) ?_⟩
  filter_upwards [hP,hπ,hcapM,hD,hstep] with n hPn hπn hcapn hDn hstepn
  have hbridge := rawPoint_barycentric_bridge (sub_ne_zero.mpr (hadj n).symm)
    (sub_ne_zero.mpr (hxs n).symm) (sub_ne_zero.mpr (hxs (n+1)).symm)
    (hfx (n+1)) (hfx n) hfs (ne_of_gt (hμ n)) hDn
  have hE : 0 < (x (n+1)-r)*(x n-r) := by
    have hm := mul_pos_of_neg_of_neg (hside (n+1)) (hside n)
    by_contra hn
    have hn0 : (x (n+1)-r)*(x n-r) ≤ 0 := le_of_not_gt hn
    have hprod := mul_nonpos_of_nonpos_of_nonneg hn0 (mul_self_nonneg (s-r))
    nlinarith
  have hb := delta_bound_of_barycentric_no_cross
    (sub_ne_zero.mpr (hadj n).symm) (sub_ne_zero.mpr (hxs n).symm)
    (sub_ne_zero.mpr (hxs (n+1)).symm)
    (sub_ne_zero.mpr (hxr (n+1))) (sub_ne_zero.mpr (hxr n)) hsr
    hfs (hμ n) (hμone n) (hquot (n+1)) (hquot n) hquots hbridge.1
    (hstepn.trans hbridge.2) hE hPn hπn (hside (n+2))
  exact hb.le.trans hcapn.le

/-- Complete nondegenerate one-sided order theorem for actual raw proposals.
The positive limiting accumulated weight and golden Q-order are derived from
the tail's noncrossing signs and local regularity, without bounded-scaling
premises. The unmodified-tail and nonzero-leading-coefficient hypotheses are
the intended hypotheses of Proposition 12, not a claim of universal entry. -/
theorem rawPoint_one_sided_hasQOrder
    {f g : ℝ → ℝ} {x μ : ℕ → ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgrgs : 0 < g r * g s) (hs : s ≠ r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n ≠ r)
    (hadj : ∀ n, x (n+1) ≠ x n) (hxs : ∀ n, x n ≠ s)
    (hμ : ∀ n, 0 < μ n) (hμone : ∀ n, μ n ≤ 1) (hanti : Antitone μ)
    (hfs : f s ≠ 0) (hfx : ∀ n, f (x n) ≠ 0)
    (hquot : ∀ n, (x n-r)/f (x n) = g (x n))
    (hquots : (s-r)/f s = g s)
    (hside : ∀ n, (x n-r)*(s-r) < 0)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0)
    (hstep : ∀ᶠ n in atTop, x (n+2) =
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1) ≠ 0) :
    (∃ w : ℝ, 0 < w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  obtain ⟨M,hM⟩ := exists_scaling_bound_of_rawPoint_one_sided hg hgrgs hs hx hxr
    hadj hxs hμ hμone hanti hfs hfx hquot hquots hside hD hstep
  have hgr : g r ≠ 0 := by intro h; rw [h,zero_mul] at hgrgs; exact lt_irrefl _ hgrgs
  exact rawPoint_one_sided_hasQOrder_of_bounded_scaling hg hgr hs hx hxr
    hadj hxs hμ hμone hanti hM hfs hfx hquot hquots hD hstep hC

/-- The actual AB correction is at most one whenever the two fresh residuals
have the same sign. Its positive fallback is included. -/
theorem abFactor_le_one {new old : ℝ} (h : 0 < new/old) :
    ModAB.Convergence.abFactor new old ≤ 1 := by
  unfold ModAB.Convergence.abFactor
  split_ifs <;> linarith

/-- The weight hypotheses used above follow from the source's actual AB
factors, including its positive fallback. -/
theorem weights_of_AB_updates {μ new old : ℕ → ℝ}
    (hμ0 : 0 < μ 0) (hμ01 : μ 0 ≤ 1)
    (hratio : ∀ n, 0 < new n / old n)
    (hstep : ∀ n, μ (n+1) = μ n * ModAB.Convergence.abFactor (new n) (old n)) :
    (∀ n, 0 < μ n) ∧ (∀ n, μ n ≤ 1) ∧ Antitone μ := by
  have hμ : ∀ n, 0 < μ n := by
    intro n
    induction n with
    | zero => exact hμ0
    | succ n ih => rw [hstep]; exact mul_pos ih (ModAB.Convergence.abFactor_positive _ _)
  have hanti : Antitone μ := weight_antitone_of_factors (fun n => (hμ n).le)
    (fun n => abFactor_le_one (hratio n)) hstep
  exact ⟨hμ,fun n => (hanti (Nat.zero_le n)).trans hμ01,hanti⟩

/-- At every C³ simple root, the regularized reciprocal has the same nonzero
sign throughout some neighbourhood. No normalization of the function is
needed: negative derivatives are covered as well. -/
theorem exists_regularized_sign_neighbourhood {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ∀ z : ℝ, |z-r| < ε →
      0 < Regularity.regularizedReciprocal f r r *
        Regularity.regularizedReciprocal f r z := by
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 2 g r :=
    Regularity.regularizedReciprocal_contDiffAt_of_dslope (DSlope.contDiffAt_dslope hf) hd
  have hgr : g r ≠ 0 := by
    dsimp [g]
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hcont : ContinuousAt (fun z => g r * g z) r := continuousAt_const.mul hg.continuousAt
  have hevent : ∀ᶠ z in 𝓝 r, 0 < g r * g z :=
    hcont.tendsto.eventually (Ioi_mem_nhds (mul_self_pos.mpr hgr))
  obtain ⟨ε,hε,hball⟩ := Metric.mem_nhds_iff.mp hevent
  refine ⟨ε,hε,?_⟩
  intro z hz
  exact hball (by simpa only [Metric.mem_ball,Real.dist_eq] using hz)

/-- Algebraic conversion from a common nonzero sign relative to a reference
to a common sign between the two factors themselves. -/
theorem common_sign_product {q a b : ℝ}
    (ha : 0 < q*a) (hb : 0 < q*b) : 0 < a*b := by
  rcases mul_pos_iff.mp ha with ha | ha <;>
    rcases mul_pos_iff.mp hb with hb | hb
  · exact mul_pos ha.2 hb.2
  · linarith [ha.1,hb.1]
  · linarith [ha.1,hb.1]
  · exact mul_pos_of_neg_of_neg ha.2 hb.2

/-- Removing the reciprocal factors translates the geometric opposite-side
condition into the corresponding true-residual sign condition. -/
theorem opposite_true_residuals {a s r fa fs ga gs : ℝ}
    (hfa : fa ≠ 0) (hfs : fs ≠ 0)
    (hqa : (a-r)/fa = ga) (hqs : (s-r)/fs = gs)
    (hg : 0 < ga*gs) (he : (a-r)*(s-r) < 0) : fa*fs < 0 := by
  have hea : a-r = ga*fa := (div_eq_iff hfa).mp hqa
  have hes : s-r = gs*fs := (div_eq_iff hfs).mp hqs
  rw [hea,hes] at he
  by_contra hn
  have hn0 : 0 ≤ fa*fs := le_of_not_gt hn
  have hnonneg := mul_nonneg hg.le hn0
  nlinarith

/-- C³ version of Proposition 12. The neighbourhood, signs of the removable
reciprocal, positivity/monotonicity of the actual AB products, and their
uniform scaling bound are all derived. Only the explicitly stated unmodified
convergent tail and nonzero limiting coefficient remain assumptions. -/
theorem smooth_rawPoint_one_sided_hasQOrder
    {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hzero : f r = 0) (hd : deriv f r ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧ ∀ (s : ℝ) (x μ : ℕ → ℝ),
      s ≠ r → |s-r| < ε →
      (∀ n, |x n-r| < ε) → Tendsto x atTop (𝓝 r) →
      (∀ n, x n ≠ r) → (∀ n, x (n+1) ≠ x n) → (∀ n, x n ≠ s) →
      (∀ n, (x n-r)*(s-r) < 0) →
      0 < μ 0 → μ 0 ≤ 1 →
      (∀ n, μ (n+1) = μ n *
        ModAB.Convergence.abFactor (f (x (n+2))) (f (x (n+1)))) →
      (∀ᶠ n in atTop, ModAB.Rational.denominator (x (n+1)) s
        (f (x (n+1))) (μ n*f s) (x n) (f (x n)) ≠ 0) →
      (∀ᶠ n in atTop, x (n+2) = ModAB.Rational.rawPoint (x (n+1)) s
        (f (x (n+1))) (μ n*f s) (x n) (f (x n))) →
      limitingCoefficient (Regularity.regularizedReciprocal f r) r s
        (⨆ n, 1/μ n-1) ≠ 0 →
      (∃ w : ℝ, 0 < w ∧ Tendsto μ atTop (𝓝 w)) ∧
      ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
        (|limitingCoefficient (Regularity.regularizedReciprocal f r) r s
          (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  obtain ⟨ε,hε,hsign⟩ := exists_regularized_sign_neighbourhood hf hd
  refine ⟨ε,hε,?_⟩
  intro s x μ hs hsmallS hsmallX hx hxr hadj hxs hside hμ0 hμ01 hμstep hD hstep hC
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 1 g r :=
    (Regularity.regularizedReciprocal_contDiffAt_of_dslope
      (DSlope.contDiffAt_dslope hf) hd).of_le (by norm_num)
  have hgs : 0 < g r * g s := hsign s hsmallS
  have hgx : ∀ n, 0 < g r * g (x n) := fun n => hsign (x n) (hsmallX n)
  have hquot : ∀ n, (x n-r)/f (x n) = g (x n) :=
    fun n => (Regularity.regularizedReciprocal_off_root hzero (hxr n)).symm
  have hquots : (s-r)/f s = g s :=
    (Regularity.regularizedReciprocal_off_root hzero hs).symm
  have hfs : f s ≠ 0 := by
    intro h
    have hg0 : g s = 0 := by rw [← hquots,h,div_zero]
    rw [hg0,mul_zero] at hgs
    exact lt_irrefl _ hgs
  have hfx : ∀ n, f (x n) ≠ 0 := by
    intro n h
    have hg0 : g (x n) = 0 := by rw [← hquot n,h,div_zero]
    have hgn := hgx n
    rw [hg0,mul_zero] at hgn
    exact lt_irrefl _ hgn
  have hres : ∀ n, f (x n)*f s < 0 := by
    intro n
    exact opposite_true_residuals (hfx n) hfs (hquot n) hquots
      (common_sign_product (hgx n) hgs) (hside n)
  have hratio : ∀ n, 0 < f (x (n+2))/f (x (n+1)) := by
    intro n
    have hm : 0 < f (x (n+2))*f (x (n+1)) :=
      common_sign_product (q := -f s) (by nlinarith [hres (n+2)])
        (by nlinarith [hres (n+1)])
    exact div_pos_iff.mpr (mul_pos_iff.mp hm)
  obtain ⟨hμ,hμone,hanti⟩ := weights_of_AB_updates hμ0 hμ01 hratio hμstep
  exact rawPoint_one_sided_hasQOrder hg hgs hs hx hxr hadj hxs hμ hμone hanti
    hfs hfx hquot hquots hside hD hstep hC

end
end ModAB.Rational.OneSidedComplete
