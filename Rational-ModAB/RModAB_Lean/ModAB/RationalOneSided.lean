import ModAB.RationalLocalAnalysis
import Mathlib.Topology.Order.MonotoneConvergence

/-!
The one-sided local regime: nonvanishing of the accumulated weight is obtained
from a finite bound on δ, not from an assumed superlinear rate.  The coefficient
limit is proved for actual divided differences with a fixed opposite endpoint.
The entry into this regime and exclusion of the exceptional zero coefficient
are separate hypotheses in the paper's local-order result.
-/
namespace ModAB.Rational.OneSided
open Filter Set
open scoped Topology
open ModAB.Rational.LocalAlgebra ModAB.Rational.LocalAnalysis
noncomputable section

/-- The accumulated weight is decreasing when each AB multiplier is in (0,1]. -/
theorem weight_antitone_of_factors {μ α : ℕ→ℝ}
    (hμ : ∀ n, 0≤μ n) (hα : ∀ n, α n≤1)
    (hstep : ∀ n, μ (n+1)=μ n*α n) : Antitone μ := by
  apply antitone_nat_of_succ_le
  intro n
  rw [hstep]
  nlinarith [hμ n,hα n]

theorem delta_monotone_of_weight_antitone {μ : ℕ→ℝ}
    (hμ : ∀ n, 0<μ n) (hanti : Antitone μ) :
    Monotone (fun n => 1/μ n-1) := by
  intro i j hij
  have h := one_div_le_one_div_of_le (hμ j) (hanti hij)
  linarith

/-- A bounded monotone δ sequence has a finite nonnegative limit, giving a
strictly positive limiting weight without any hypothesis about convergence order. -/
theorem exists_delta_weight_limit {δ : ℕ→ℝ} {M : ℝ}
    (hδ : ∀ n, 0≤δ n) (hmono : Monotone δ) (hM : ∀ n, δ n≤M) :
    ∃ L : ℝ, 0≤L ∧ L≤M ∧ Tendsto δ atTop (𝓝 L) ∧
      Tendsto (fun n => 1/(1+δ n)) atTop (𝓝 (1/(1+L))) ∧
      0<1/(1+L) ∧ 1/(1+M)≤1/(1+L) := by
  have hbdd : BddAbove (range δ) := by
    refine ⟨M,?_⟩
    rintro _ ⟨n,rfl⟩
    exact hM n
  let L : ℝ := ⨆ n, δ n
  have hlim : Tendsto δ atTop (𝓝 L) := tendsto_atTop_ciSup hmono hbdd
  have hL0 : 0≤L := ge_of_tendsto hlim (Eventually.of_forall hδ)
  have hLM : L≤M := le_of_tendsto hlim (Eventually.of_forall hM)
  have hL : 0<1+L := by linarith
  have hw : Tendsto (fun n => 1/(1+δ n)) atTop (𝓝 (1/(1+L))) :=
    tendsto_const_nhds.div (tendsto_const_nhds.add hlim) (ne_of_gt hL)
  exact ⟨L,hL0,hLM,hlim,hw,one_div_pos.mpr hL,(weight_lower_bound hL0 hLM).2⟩

theorem exists_positive_weight_limit {μ : ℕ→ℝ} {M : ℝ}
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) :
    ∃ L : ℝ, 0<L ∧ 1/(1+M)≤L ∧ Tendsto μ atTop (𝓝 L) := by
  have hδ : ∀ n, 0≤1/μ n-1 := by
    intro n
    have h : 1≤1/μ n := by
      simpa only [one_div_one] using one_div_le_one_div_of_le (hμ n) (hμone n)
    linarith
  obtain ⟨L,hL0,hLM,hlim,hw,hpos,hlower⟩ :=
    exists_delta_weight_limit hδ (delta_monotone_of_weight_antitone hμ hanti) hbound
  refine ⟨1/(1+L),hpos,hlower,?_⟩
  simpa only [weight_from_delta (ne_of_gt (hμ _))] using hw

/-- The no-crossing sign test supplies the bound needed for a positive weight
limit.  No convergence rate or summability of residual ratios is assumed. -/
theorem exists_positive_weight_limit_of_same_side
    {μ E G P g π : ℕ→ℝ} {M : ℝ}
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hE : ∀ n, 0<E n) (hP : ∀ n, 0<P n)
    (hg : ∀ n, 0<g n) (hπ : ∀ n, 0<π n)
    (hsame : ∀ n, E n*(G n+(1/μ n-1)*g n/π n)/
      (P n+E n*((1/μ n-1)*g n/π n))<0)
    (hcap : ∀ n, -π n*G n/g n≤M) :
    ∃ L : ℝ, 0<L ∧ 1/(1+M)≤L ∧ Tendsto μ atTop (𝓝 L) := by
  apply exists_positive_weight_limit hμ hμone hanti
  intro n
  have hδ : 0≤1/μ n-1 := by
    have h : 1≤1/μ n := by
      simpa only [one_div_one] using one_div_le_one_div_of_le (hμ n) (hμone n)
    linarith
  exact (delta_bound_of_same_side (hE n) (hP n) (hg n) hδ (hπ n) (hsame n)).le.trans (hcap n)

/-- The second divided difference with two coincident limiting nodes and one
fixed node.  It is defined by the derivative and a secant, not by division by zero. -/
def confluentDD2 (g : ℝ→ℝ) (r s : ℝ) : ℝ :=
  ((g s-g r)/(s-r)-deriv g r)/(s-r)

theorem tendsto_dd2_fixed_endpoint {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (confluentDD2 g r s)) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hga := hg.continuous.continuousAt.tendsto.comp ha
  have hgb := hg.continuous.continuousAt.tendsto.comp hb
  have hfirst := tendsto_dd1_arbitrary_of_contDiff_one hg hab ha hb
  have hsec : Tendsto (fun n => (g s-g (b n))/(s-b n)) atTop
      (𝓝 ((g s-g r)/(s-r))) :=
    (tendsto_const_nhds.sub hgb).div (tendsto_const_nhds.sub hb) hsr
  convert (hsec.sub hfirst).div (tendsto_const_nhds.sub ha) hsr using 1
  · funext n
    rfl
  · rfl

theorem tendsto_poleValue_fixed_endpoint {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (g r)) := by
  have hfirst := tendsto_dd1_arbitrary_of_contDiff_one hg hab ha hb
  have hsecond := tendsto_dd2_fixed_endpoint hg hab hs ha hb
  have h := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (hfirst.mul (ha.sub_const r))).add
      ((hsecond.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [poleValue,Function.comp_def,sub_self,mul_zero,sub_zero,add_zero] using h

/-- Coefficient multiplying the two fresh errors in the exact mixed-ordinate
formula when the opposite endpoint s is fixed. -/
def mixedCoefficient (g : ℝ→ℝ) (r s δ a b : ℝ) : ℝ :=
  (dd2 a b s (g a) (g b) (g s)+δ*g s/((s-a)*(s-b))) /
    (poleValue a b s r (g a) (g b) (g s)+
      δ*g s*(a-r)*(b-r)/((s-a)*(s-b)))*(s-r)

/-- The coefficient whose limit is proved below is exactly the coefficient of
the new rational error, provided g(x)=(x-r)/f(x) at the three supplied nodes. -/
theorem baryRoot_error_eq_mixedCoefficient {g : ℝ→ℝ} {a b s r ya yb ys μ : ℝ}
    (hba : b-a≠0) (hsb : s-b≠0) (hsa : s-a≠0)
    (har : a-r≠0) (hbr : b-r≠0) (hsr : s-r≠0)
    (hys : ys≠0) (hμ : μ≠0)
    (hH : dd2 a b s (1/ya) (1/yb) (1/(μ*ys))≠0)
    (hga : (a-r)/ya=g a) (hgb : (b-r)/yb=g b) (hgs : (s-r)/ys=g s) :
    baryRoot a b s ya yb (μ*ys)-r =
      mixedCoefficient g r s (1/μ-1) a b*((a-r)*(b-r)) := by
  rw [baryRoot_scaled_error hba hsb hsa har hbr hsr hys hμ hH,hga,hgb,hgs]
  unfold mixedCoefficient
  ring

theorem tendsto_mixedCoefficient {g : ℝ→ℝ} {a b δ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiff ℝ 1 g) (hgr : g r≠0) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hδ : Tendsto δ atTop (𝓝 L)) :
    Tendsto (fun n => mixedCoefficient g r s (δ n) (a n) (b n)) atTop
      (𝓝 ((confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r))) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hπ : Tendsto (fun n => (s-a n)*(s-b n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)
  have hπne : (s-r)*(s-r)≠0 := mul_ne_zero hsr hsr
  have hsecond := tendsto_dd2_fixed_endpoint hg hab hs ha hb
  have hpole := tendsto_poleValue_fixed_endpoint hg hab hs ha hb
  have hcorr := (hδ.mul_const (g s)).div hπ hπne
  have hsmall := (((hδ.mul_const (g s)).mul (ha.sub_const r)).mul
    (hb.sub_const r)).div hπ hπne
  have hden : Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s)+
      δ n*g s*(a n-r)*(b n-r)/((s-a n)*(s-b n))) atTop (𝓝 (g r)) := by
    simpa only [sub_self,mul_zero,zero_div,add_zero,Pi.div_apply] using hpole.add hsmall
  simpa only [mixedCoefficient,Pi.div_apply] using ((hsecond.add hcorr).div hden hgr).mul_const (s-r)

/-- Only C¹ regularity near the limiting root is needed for the two moving
nodes; the function may have poles elsewhere. -/
theorem tendsto_dd1_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n))) atTop (𝓝 (deriv g r)) := by
  obtain ⟨δ,hδ,hball⟩ := Metric.mem_nhds_iff.mp (hg.eventually (by simp))
  let J := Ioo (r-δ) (r+δ)
  have hrJ : J ∈ 𝓝 r := Ioo_mem_nhds (by linarith) (by linarith)
  have hJ : ∀ z ∈ J, ContDiffAt ℝ 1 g z := by
    intro z hz
    apply hball
    rw [Metric.mem_ball,Real.dist_eq]
    apply abs_lt.mpr
    dsimp [J] at hz
    constructor <;> linarith [hz.1,hz.2]
  obtain ⟨N,hN⟩ := eventually_atTop.mp ((ha.eventually hrJ).and (hb.eventually hrJ))
  have hm (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hreg (n : ℕ) (z : ℝ)
      (hz : z ∈ Icc (min (a (n+N)) (b (n+N))) (max (a (n+N)) (b (n+N)))) :
      ContDiffAt ℝ 1 g z := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (hm n).2.1).trans_le hz.1,
      hz.2.trans_lt (max_lt (hm n).1.2 (hm n).2.2)⟩
  have hg' : ContDiffAt ℝ 0 (deriv g) r := hg.derivWithin (by norm_num)
  apply (tendsto_add_atTop_iff_nat N).mp
  exact tendsto_dd1_arbitrary_of_derivatives
    (fun n z hz => (hreg n z hz).continuousAt.continuousWithinAt)
    (fun n z hz => ((hreg n z ⟨hz.1.le,hz.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
    hg'.continuousAt (fun n => hab (n+N))
    (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))


/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_dd2_fixed_endpoint_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) s (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (confluentDD2 g r s)) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hga := hg.continuousAt.tendsto.comp ha
  have hgb := hg.continuousAt.tendsto.comp hb
  have hfirst := tendsto_dd1_of_contDiffAt hg hab ha hb
  have hsec : Tendsto (fun n => (g s-g (b n))/(s-b n)) atTop
      (𝓝 ((g s-g r)/(s-r))) :=
    (tendsto_const_nhds.sub hgb).div (tendsto_const_nhds.sub hb) hsr
  convert (hsec.sub hfirst).div (tendsto_const_nhds.sub ha) hsr using 1
  · funext n
    rfl
  · rfl

/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_poleValue_fixed_endpoint_of_contDiffAt {g : ℝ→ℝ} {a b : ℕ→ℝ} {r s : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s))
      atTop (𝓝 (g r)) := by
  have hfirst := tendsto_dd1_of_contDiffAt hg hab ha hb
  have hsecond := tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have h := ((hg.continuousAt.tendsto.comp ha).sub
    (hfirst.mul (ha.sub_const r))).add
      ((hsecond.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [poleValue,Function.comp_def,sub_self,mul_zero,sub_zero,add_zero] using h

/-- Local regularity version; no assumptions on g away from the root. -/
theorem tendsto_mixedCoefficient_of_contDiffAt {g : ℝ→ℝ} {a b δ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hab : ∀ n, a n≠b n) (hs : s≠r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hδ : Tendsto δ atTop (𝓝 L)) :
    Tendsto (fun n => mixedCoefficient g r s (δ n) (a n) (b n)) atTop
      (𝓝 ((confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r))) := by
  have hsr : s-r≠0 := sub_ne_zero.mpr hs
  have hπ : Tendsto (fun n => (s-a n)*(s-b n)) atTop (𝓝 ((s-r)*(s-r))) :=
    (tendsto_const_nhds.sub ha).mul (tendsto_const_nhds.sub hb)
  have hπne : (s-r)*(s-r)≠0 := mul_ne_zero hsr hsr
  have hsecond := tendsto_dd2_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have hpole := tendsto_poleValue_fixed_endpoint_of_contDiffAt hg hab hs ha hb
  have hcorr := (hδ.mul_const (g s)).div hπ hπne
  have hsmall := (((hδ.mul_const (g s)).mul (ha.sub_const r)).mul
    (hb.sub_const r)).div hπ hπne
  have hden : Tendsto (fun n => poleValue (a n) (b n) s r (g (a n)) (g (b n)) (g s)+
      δ n*g s*(a n-r)*(b n-r)/((s-a n)*(s-b n))) atTop (𝓝 (g r)) := by
    simpa only [sub_self,mul_zero,zero_div,add_zero,Pi.div_apply] using hpole.add hsmall
  simpa only [mixedCoefficient,Pi.div_apply] using ((hsecond.add hcorr).div hden hgr).mul_const (s-r)

end
end ModAB.Rational.OneSided
