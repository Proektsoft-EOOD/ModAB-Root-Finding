import ModAB.RationalOneSided
import ModAB.RationalGoldenOrder

/-!
Connected one-sided local-order theorems.  The hypotheses describe a convergent
nonvanishing rational tail with one fixed opposite endpoint.  Its recurrence
coefficient and golden-ratio Q-order are conclusions, not assumed rates.
Entry of the safeguarded program into such a tail is not claimed here.
-/
namespace ModAB.Rational.OneSidedAssembly
open Filter Set
open scoped Topology
open ModAB.Rational.LocalAlgebra ModAB.Rational.OneSided
noncomputable section

def limitingCoefficient (g : ℝ→ℝ) (r s L : ℝ) : ℝ :=
  (confluentDD2 g r s+L*g s/((s-r)*(s-r)))/g r*(s-r)

/-- Actual barycentric rational updates with a fixed opposite endpoint have
golden-ratio Q-order when their derived limiting coefficient is nonzero. -/
theorem barycentric_one_sided_hasQOrder
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s L : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, μ n≠0) (hfs : f s≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hδ : Tendsto (fun n => 1/μ n-1) atTop (𝓝 L))
    (hH : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0)
    (hstep : ∀ᶠ n in atTop,
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s))
    (hC : limitingCoefficient g r s L≠0) :
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s L|^(Real.goldenRatio-1)) := by
  let C : ℕ→ℝ := fun n => mixedCoefficient g r s (1/μ n-1) (x (n+1)) (x n)
  have hxnext : Tendsto (fun n => x (n+1)) atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 1).mpr hx
  have hClim : Tendsto C atTop (𝓝 (limitingCoefficient g r s L)) :=
    tendsto_mixedCoefficient_of_contDiffAt hg hgr hadj hs hxnext hx hδ
  have hrec : ∀ᶠ n in atTop,
      x (n+2)-r=C n*(x (n+1)-r)*(x n-r) := by
    filter_upwards [hH,hstep] with n hHn hstepn
    rw [hstepn]
    have h := baryRoot_error_eq_mixedCoefficient
      (sub_ne_zero.mpr (hadj n).symm) (sub_ne_zero.mpr (hxs n).symm)
      (sub_ne_zero.mpr (hxs (n+1)).symm)
      (sub_ne_zero.mpr (hxr (n+1))) (sub_ne_zero.mpr (hxr n)) (sub_ne_zero.mpr hs)
      hfs (hμ n) hHn (hquot (n+1)) (hquot n) hquots
    simpa only [C,mul_assoc] using h
  exact ModAB.RationalGoldenOrder.eventual_two_error_recurrence_hasQOrder
    (fun n => sub_ne_zero.mpr (hxr n)) hC (by simpa using hx.sub_const r) hClim hrec

/-- Bounded decreasing positive weights supply their own nonzero limiting
weight and δ-limit.  The limiting coefficient is evaluated at sup δ. -/
theorem barycentric_one_sided_hasQOrder_of_bounded_scaling
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s M : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) (hfs : f s≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hH : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0)
    (hstep : ∀ᶠ n in atTop,
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1)≠0) :
    (∃ w : ℝ, 0<w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  have hbdd : BddAbove (range (fun n => 1/μ n-1)) := by
    refine ⟨M,?_⟩
    rintro _ ⟨n,rfl⟩
    exact hbound n
  have hδ := tendsto_atTop_ciSup (delta_monotone_of_weight_antitone hμ hanti) hbdd
  obtain ⟨w,hw,_,hwl⟩ := exists_positive_weight_limit hμ hμone hanti hbound
  exact ⟨⟨w,hw,hwl⟩,barycentric_one_sided_hasQOrder hg hgr hs hx hxr hadj hxs
    (fun n => ne_of_gt (hμ n)) hfs hquot hquots hδ hH hstep hC⟩

/-- Reorder the source's endpoint--opposite--history data into the two-fresh-
nodes--scaled-opposite order used by the local identities. -/
theorem rawPoint_barycentric_bridge {a b s ya yb ys μ : ℝ}
    (hba : b-a≠0) (hsb : s-b≠0) (hsa : s-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hys : ys≠0) (hμ : μ≠0)
    (hD : ModAB.Rational.denominator a s ya (μ*ys) b yb≠0) :
    dd2 a b s (1/ya) (1/yb) (1/(μ*ys))≠0 ∧
    ModAB.Rational.rawPoint a s ya (μ*ys) b yb=baryRoot a b s ya yb (μ*ys) := by
  have hbs : b-s≠0 := by intro h; apply hsb; linarith
  have hh := bary_denominator_ne_zero hsa hbs hba hya (mul_ne_zero hμ hys) hyb hD
  have hswap := dd2_swap_last (v₁:=1/ya) (v₂:=1/yb) (v₃:=1/(μ*ys)) hba hsb hsa
  constructor
  · rw [hswap]
    exact hh
  · rw [rawPoint_eq_baryRoot hsa hbs hba hya (mul_ne_zero hμ hys) hyb hD]
    exact (baryRoot_swap_last hba hsb hsa).symm

/-- Connected source-proposal result: the normalized rawPoint iteration,
bounded accumulated scaling, and the analytic hypotheses imply golden Q-order.
The clamp/closing/fallback must already be inactive on the stated tail. -/
theorem rawPoint_one_sided_hasQOrder_of_bounded_scaling
    {f g : ℝ→ℝ} {x μ : ℕ→ℝ} {r s M : ℝ}
    (hg : ContDiffAt ℝ 1 g r) (hgr : g r≠0) (hs : s≠r)
    (hx : Tendsto x atTop (𝓝 r)) (hxr : ∀ n, x n≠r)
    (hadj : ∀ n, x (n+1)≠x n) (hxs : ∀ n, x n≠s)
    (hμ : ∀ n, 0<μ n) (hμone : ∀ n, μ n≤1) (hanti : Antitone μ)
    (hbound : ∀ n, 1/μ n-1≤M) (hfs : f s≠0) (hfx : ∀ n, f (x n)≠0)
    (hquot : ∀ n, (x n-r)/f (x n)=g (x n))
    (hquots : (s-r)/f s=g s)
    (hD : ∀ᶠ n in atTop,
      ModAB.Rational.denominator (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n))≠0)
    (hstep : ∀ᶠ n in atTop, x (n+2)=
      ModAB.Rational.rawPoint (x (n+1)) s (f (x (n+1))) (μ n*f s) (x n) (f (x n)))
    (hC : limitingCoefficient g r s (⨆ n, 1/μ n-1)≠0) :
    (∃ w : ℝ, 0<w ∧ Tendsto μ atTop (𝓝 w)) ∧
    ModAB.RationalOrder.HasQOrder (fun n => x n-r) Real.goldenRatio
      (|limitingCoefficient g r s (⨆ n, 1/μ n-1)|^(Real.goldenRatio-1)) := by
  have hbridge : ∀ᶠ n in atTop,
      dd2 (x (n+1)) (x n) s (1/f (x (n+1))) (1/f (x n)) (1/(μ n*f s))≠0 ∧
      x (n+2)=baryRoot (x (n+1)) (x n) s (f (x (n+1))) (f (x n)) (μ n*f s) := by
    filter_upwards [hD,hstep] with n hDn hstepn
    have h := rawPoint_barycentric_bridge (sub_ne_zero.mpr (hadj n).symm)
      (sub_ne_zero.mpr (hxs n).symm) (sub_ne_zero.mpr (hxs (n+1)).symm)
      (hfx (n+1)) (hfx n) hfs (ne_of_gt (hμ n)) hDn
    exact ⟨h.1,hstepn.trans h.2⟩
  exact barycentric_one_sided_hasQOrder_of_bounded_scaling hg hgr hs hx hxr hadj hxs
    hμ hμone hanti hbound hfs hquot hquots (hbridge.mono fun _ h => h.1)
    (hbridge.mono fun _ h => h.2) hC

end
end ModAB.Rational.OneSidedAssembly
