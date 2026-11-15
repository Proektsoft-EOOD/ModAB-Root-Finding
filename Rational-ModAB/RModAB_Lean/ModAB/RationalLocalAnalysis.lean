import ModAB.RationalLocalAlgebra
import Mathlib.Analysis.Calculus.LocalExtr.Rolle
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.ContDiff.Deriv
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Topology.Order.Basic

/-!
Analytic bridge for second divided differences.  The mean-value theorem is
proved by applying Rolle's theorem twice to the interpolation remainder.
Derivatives are supplied explicitly, so the statements also apply to local
C² hypotheses after the standard derivative/continuity conversions.
-/
namespace ModAB.Rational.LocalAnalysis
open Set Filter
open scoped Topology
open LocalAlgebra
noncomputable section

/-- The Newton quadratic interpolates its third node. -/
theorem quadratic_at_third {a b c ga gb gc : ℝ}
    (hab : a < b) (hbc : b < c) :
    ga + dd1 a b ga gb * (c-a) +
      dd2 a b c ga gb gc * ((c-a)*(c-b)) = gc := by
  have hab0 : b-a ≠ 0 := ne_of_gt (sub_pos.mpr hab)
  have hbc0 : c-b ≠ 0 := ne_of_gt (sub_pos.mpr hbc)
  have hac0 : c-a ≠ 0 := ne_of_gt (sub_pos.mpr (hab.trans hbc))
  unfold dd2 dd1
  field_simp
  <;> ring

/-- A second divided difference is half a second derivative at an intermediate
point.  Every analytic hypothesis is local to the interval of the three nodes. -/
theorem exists_second_derivative_eq_twice_dd2
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℝ}
    (hab : a < b) (hbc : b < c)
    (hg : ContinuousOn g (Icc a c))
    (hg₁ : ContinuousOn g₁ (Icc a c))
    (hd : ∀ x ∈ Ioo a c, HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ x ∈ Ioo a c, HasDerivAt g₁ (g₂ x) x) :
    ∃ ξ ∈ Ioo a c, g₂ ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  let s := dd1 a b (g a) (g b)
  let D := dd2 a b c (g a) (g b) (g c)
  let q : ℝ → ℝ := fun x => g a+s*(x-a)+D*((x-a)*(x-b))
  let q₁ : ℝ → ℝ := fun x => s+D*(2*x-a-b)
  have hq : ∀ x, HasDerivAt q (q₁ x) x := by
    intro x
    have hxa := (hasDerivAt_id x).sub_const a
    have hxb := (hasDerivAt_id x).sub_const b
    convert ((hasDerivAt_const x (g a)).add (hxa.const_mul s)).add
      ((hxa.mul hxb).const_mul D) using 1 <;> (try funext y) <;> (try dsimp [q, q₁]) <;> ring
  have hq₁ : ∀ x, HasDerivAt q₁ (2*D) x := by
    intro x
    convert (hasDerivAt_const x s).add
      ((((hasDerivAt_id x).const_mul 2).sub_const a).sub_const b |>.const_mul D)
      using 1 <;> (try funext y) <;> (try dsimp [q₁]) <;> ring
  have hqc : Continuous q := continuous_iff_continuousAt.mpr
    (fun x => (hq x).continuousAt)
  have hq₁c : Continuous q₁ := continuous_iff_continuousAt.mpr
    (fun x => (hq₁ x).continuousAt)
  have hqa : q a = g a := by simp [q]
  have hqb : q b = g b := by
    dsimp [q, s]
    unfold dd1
    have h : b-a ≠ 0 := ne_of_gt (sub_pos.mpr hab)
    field_simp
    <;> ring
  have hqcval : q c = g c := quadratic_at_third hab hbc
  let h : ℝ → ℝ := fun x => g x-q x
  let h₁ : ℝ → ℝ := fun x => g₁ x-q₁ x
  have hc : ContinuousOn h (Icc a c) := hg.sub hqc.continuousOn
  have h₁c : ContinuousOn h₁ (Icc a c) := hg₁.sub hq₁c.continuousOn
  have hder : ∀ x ∈ Ioo a c, HasDerivAt h (h₁ x) x :=
    fun x hx => (hd x hx).sub (hq x)
  have hder₁ : ∀ x ∈ Ioo a c, HasDerivAt h₁ (g₂ x-2*D) x :=
    fun x hx => (hd₁ x hx).sub (hq₁ x)
  have ha : h a = 0 := by simp [h, hqa]
  have hb : h b = 0 := by simp [h, hqb]
  have hcv : h c = 0 := by simp [h, hqcval]
  obtain ⟨u, hu, hzu⟩ := exists_hasDerivAt_eq_zero hab
    (hc.mono (Icc_subset_Icc_right hbc.le)) (ha.trans hb.symm)
    (fun x hx => hder x ⟨hx.1, hx.2.trans hbc⟩)
  obtain ⟨v, hv, hzv⟩ := exists_hasDerivAt_eq_zero hbc
    (hc.mono (Icc_subset_Icc_left hab.le)) (hb.trans hcv.symm)
    (fun x hx => hder x ⟨hab.trans hx.1, hx.2⟩)
  have huv : u < v := hu.2.trans hv.1
  have huva : Icc u v ⊆ Icc a c :=
    Icc_subset_Icc hu.1.le hv.2.le
  obtain ⟨ξ, hξ, hz⟩ := exists_hasDerivAt_eq_zero huv
    (h₁c.mono huva) (hzu.trans hzv.symm)
    (fun x hx => hder₁ x ⟨hu.1.trans hx.1, hx.2.trans hv.2⟩)
  exact ⟨ξ, ⟨hu.1.trans hξ.1, hξ.2.trans hv.2⟩, sub_eq_zero.mp hz⟩


/-- First divided differences converge when both distinct ordered nodes coalesce. -/
theorem tendsto_dd1_of_derivatives
    {g g₁ : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hab : ∀ n, a n < b n)
    (hg : ∀ n, ContinuousOn g (Icc (a n) (b n)))
    (hd : ∀ n x, x ∈ Ioo (a n) (b n) → HasDerivAt g (g₁ x) x)
    (hg₁ : ContinuousAt g₁ r)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (g₁ r)) := by
  choose ξ hξ heq using fun n =>
    exists_hasDerivAt_eq_slope g g₁ (hab n) (hg n) (hd n)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hb
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq, dd1] using hg₁.tendsto.comp hlim

/-- Second divided differences converge without a node-separation assumption. -/
theorem tendsto_dd2_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (hg : ∀ n, ContinuousOn g (Icc (a n) (c n)))
    (hg₁ : ∀ n, ContinuousOn g₁ (Icc (a n) (c n)))
    (hd : ∀ n x, x ∈ Ioo (a n) (c n) → HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ n x, x ∈ Ioo (a n) (c n) → HasDerivAt g₁ (g₂ x) x)
    (hg₂ : ContinuousAt g₂ r)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g₂ r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_eq_twice_dd2
    (hab n) (hbc n) (hg n) (hg₁ n) (hd n) (hd₁ n)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hc
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hv := (hg₂.tendsto.comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- C² specialization of the divided-difference limit. -/
theorem tendsto_dd2_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  exact tendsto_dd2_of_derivatives hab hbc
    (fun _ => hg.continuous.continuousOn)
    (fun _ => hg'.continuous.continuousOn)
    (fun _ x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (fun _ x _ => (hg'.differentiable (by norm_num) x).hasDerivAt)
    (hg'.continuous_deriv (by norm_num)).continuousAt ha hc

/-- The interpolation polynomial evaluated at the limit root tends to g(r). -/
theorem tendsto_poleValue_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) (c n) r
      (g (a n)) (g (b n)) (g (c n))) atTop (𝓝 (g r)) := by
  have hb : Tendsto b atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le ha hc
      (fun n => (hab n).le) (fun n => (hbc n).le)
  have h₁ := tendsto_dd1_of_derivatives hab
    (fun _ => hg.continuous.continuousOn)
    (fun _ x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (hg.continuous_deriv (by norm_num)).continuousAt ha hb
  have h₂ := tendsto_dd2_of_contDiff_two hg hab hbc ha hc
  have hv := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero,
    poleValue] using hv

/-- The exact true-data error coefficient has its analytic limit. -/
theorem tendsto_error_coefficient_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g) (hgr : g r ≠ 0)
    (hab : ∀ n, a n < b n) (hbc : ∀ n, b n < c n)
    (ha : Tendsto a atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) :=
  (tendsto_dd2_of_contDiff_two hg hab hbc ha hc).div
    (tendsto_poleValue_of_contDiff_two hg hab hbc ha hc) hgr

/-- The second-derivative mean-value statement is independent of node order. -/
theorem exists_second_derivative_dd2_arbitrary_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℝ}
    (hg : ContinuousOn g (Icc (min a (min b c)) (max a (max b c))))
    (hg₁ : ContinuousOn g₁ (Icc (min a (min b c)) (max a (max b c))))
    (hd : ∀ x ∈ Ioo (min a (min b c)) (max a (max b c)), HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ x ∈ Ioo (min a (min b c)) (max a (max b c)), HasDerivAt g₁ (g₂ x) x)
    (hab : a ≠ b) (hbc : b ≠ c) (hac : a ≠ c) :
    ∃ ξ ∈ Ioo (min a (min b c)) (max a (max b c)),
      g₂ ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  let L := min a (min b c)
  let U := max a (max b c)
  have hs : ∀ u v w : ℝ, L ≤ u → w ≤ U → u < v → v < w →
      ∃ ξ ∈ Ioo u w, g₂ ξ = 2 * dd2 u v w (g u) (g v) (g w) := by
    intro u v w hLu hwU huv hvw
    exact exists_second_derivative_eq_twice_dd2 huv hvw
      (hg.mono (Icc_subset_Icc hLu hwU))
      (hg₁.mono (Icc_subset_Icc hLu hwU))
      (fun x hx => hd x ⟨hLu.trans_lt hx.1, hx.2.trans_le hwU⟩)
      (fun x hx => hd₁ x ⟨hLu.trans_lt hx.1, hx.2.trans_le hwU⟩)
  have hLa : L ≤ a := min_le_left _ _
  have hLb : L ≤ b := (min_le_right _ _).trans (min_le_left _ _)
  have hLc : L ≤ c := (min_le_right _ _).trans (min_le_right _ _)
  have haU : a ≤ U := le_max_left _ _
  have hbU : b ≤ U := (le_max_left _ _).trans (le_max_right _ _)
  have hcU : c ≤ U := (le_max_right _ _).trans (le_max_right _ _)
  have hba0 : b-a ≠ 0 := sub_ne_zero.mpr hab.symm
  have hcb0 : c-b ≠ 0 := sub_ne_zero.mpr hbc.symm
  have hca0 : c-a ≠ 0 := sub_ne_zero.mpr hac.symm
  have hab0 : a-b ≠ 0 := sub_ne_zero.mpr hab
  have hbc0 : b-c ≠ 0 := sub_ne_zero.mpr hbc
  have hac0 : a-c ≠ 0 := sub_ne_zero.mpr hac
  have hswap := dd2_swap_last (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hfirst := dd2_swap_first (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hcyc := dd2_cyclic (v₁ := g a) (v₂ := g b) (v₃ := g c) hba0 hcb0 hca0
  have hcyc₂ := dd2_cyclic (v₁ := g b) (v₂ := g c) (v₃ := g a) hcb0 hac0 hab0
  have hrev := dd2_swap_last (v₁ := g c) (v₂ := g a) (v₃ := g b) hac0 hba0 hbc0
  have lift {u w D : ℝ} (hLu : L ≤ u) (hwU : w ≤ U)
      (ht : ∃ ξ ∈ Ioo u w, g₂ ξ = 2*D)
      (hD : dd2 a b c (g a) (g b) (g c) = D) :
      ∃ ξ ∈ Ioo L U, g₂ ξ = 2*dd2 a b c (g a) (g b) (g c) := by
    obtain ⟨ξ,hξ,heq⟩ := ht
    exact ⟨ξ, ⟨hLu.trans_lt hξ.1, hξ.2.trans_le hwU⟩, by rw [hD]; exact heq⟩
  rcases lt_or_gt_of_ne hab with hab' | hba'
  · rcases lt_or_gt_of_ne hbc with hbc' | hcb'
    · exact lift hLa hcU (hs a b c hLa hcU hab' hbc') rfl
    · rcases lt_or_gt_of_ne hac with hac' | hca'
      · exact lift hLa hbU (hs a c b hLa hbU hac' hcb') hswap
      · exact lift hLc hbU (hs c a b hLc hbU hca' hab') (hcyc.trans hcyc₂)
  · rcases lt_or_gt_of_ne hac with hac' | hca'
    · exact lift hLb hcU (hs b a c hLb hcU hba' hac') hfirst
    · rcases lt_or_gt_of_ne hbc with hbc' | hcb'
      · exact lift hLb haU (hs b c a hLb haU hbc' hca') hcyc
      · exact lift hLc haU (hs c b a hLc haU hcb' hba') ((hcyc.trans hcyc₂).trans hrev)

/-- C² specialization of the unordered mean-value theorem. -/
theorem exists_second_derivative_dd2_arbitrary
    {g : ℝ → ℝ} {a b c : ℝ} (hg : ContDiff ℝ 2 g)
    (hab : a ≠ b) (hbc : b ≠ c) (hac : a ≠ c) :
    ∃ ξ ∈ Ioo (min a (min b c)) (max a (max b c)),
      deriv (deriv g) ξ = 2 * dd2 a b c (g a) (g b) (g c) := by
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  exact exists_second_derivative_dd2_arbitrary_of_derivatives
    hg.continuous.continuousOn hg'.continuous.continuousOn
    (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    (fun x _ => (hg'.differentiable (by norm_num) x).hasDerivAt) hab hbc hac

/-- C² gives the second divided-difference limit for every distinct-node ordering. -/
theorem tendsto_dd2_arbitrary_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_dd2_arbitrary
    hg (hab n) (hbc n) (hac n)
  have hl : Tendsto (fun n => min (a n) (min (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.min (hb.min hc)
  have hu : Tendsto (fun n => max (a n) (max (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.max (hb.max hc)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le hl hu
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hg' : ContDiff ℝ 1 (deriv g) := hg.deriv'
  have hv := (((hg'.continuous_deriv (by norm_num)).continuousAt.tendsto).comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- Symmetry of the first divided difference. -/
theorem dd1_comm (a b ga gb : ℝ) : dd1 a b ga gb = dd1 b a gb ga := by
  unfold dd1
  have hn : ga-gb = -(gb-ga) := by ring
  have hd : a-b = -(b-a) := by ring
  rw [hn, hd, neg_div_neg_eq]

/-- First divided differences converge for either ordering of the nodes. -/
theorem tendsto_dd1_arbitrary_of_contDiff_one
    {g : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 1 g) (hab : ∀ n, a n ≠ b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (deriv g r)) := by
  have hmvt : ∀ n, ∃ ξ ∈ Ioo (min (a n) (b n)) (max (a n) (b n)),
      deriv g ξ = dd1 (a n) (b n) (g (a n)) (g (b n)) := by
    intro n
    rcases lt_or_gt_of_ne (hab n) with h | h
    · simpa only [min_eq_left h.le, max_eq_right h.le, dd1] using
        exists_hasDerivAt_eq_slope g (deriv g) h hg.continuous.continuousOn
          (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
    · rw [min_eq_right h.le, max_eq_left h.le, dd1_comm (a n) (b n)]
      exact exists_hasDerivAt_eq_slope g (deriv g) h hg.continuous.continuousOn
        (fun x _ => (hg.differentiable (by norm_num) x).hasDerivAt)
  choose ξ hξ heq using hmvt
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le (by simpa using ha.min hb)
      (by simpa using ha.max hb) (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq] using
    (hg.continuous_deriv_one.continuousAt.tendsto).comp hlim

/-- The exact coefficient limit is independent of node ordering and relative spacing. -/
theorem tendsto_error_coefficient_arbitrary_of_contDiff_two
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiff ℝ 2 g) (hgr : g r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r))
    (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n =>
      dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) := by
  have h₁ := tendsto_dd1_arbitrary_of_contDiff_one (hg.of_le (by norm_num)) hab ha hb
  have h₂ := tendsto_dd2_arbitrary_of_contDiff_two hg hab hbc hac ha hb hc
  have hv := ((hg.continuous.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  have hp : Tendsto (fun n => poleValue (a n) (b n) (c n) r
      (g (a n)) (g (b n)) (g (c n))) atTop (𝓝 (g r)) := by
    simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero, poleValue] using hv
  exact h₂.div hp hgr

/-- Unordered coalescing-node limit from explicit local derivatives. -/
theorem tendsto_dd2_arbitrary_of_derivatives
    {g g₁ g₂ : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ∀ n, ContinuousOn g (Icc (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n)))))
    (hg₁ : ∀ n, ContinuousOn g₁ (Icc (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n)))))
    (hd : ∀ n x, x ∈ Ioo (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n))) → HasDerivAt g (g₁ x) x)
    (hd₁ : ∀ n x, x ∈ Ioo (min (a n) (min (b n) (c n))) (max (a n) (max (b n) (c n))) → HasDerivAt g₁ (g₂ x) x)
    (hg₂ : ContinuousAt g₂ r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g₂ r / 2)) := by
  choose ξ hξ heq using fun n => exists_second_derivative_dd2_arbitrary_of_derivatives
    (hg n) (hg₁ n) (hd n) (hd₁ n) (hab n) (hbc n) (hac n)
  have hl : Tendsto (fun n => min (a n) (min (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.min (hb.min hc)
  have hu : Tendsto (fun n => max (a n) (max (b n) (c n))) atTop (𝓝 r) :=
    by simpa using ha.max (hb.max hc)
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le hl hu
      (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  have hv := (hg₂.tendsto.comp hlim).div_const 2
  convert hv using 1
  ext n
  simp only [Function.comp_def]
  rw [heq n]
  ring

/-- Local derivative form of the first divided-difference limit. -/
theorem tendsto_dd1_arbitrary_of_derivatives
    {g g₁ : ℝ → ℝ} {a b : ℕ → ℝ} {r : ℝ}
    (hg : ∀ n, ContinuousOn g (Icc (min (a n) (b n)) (max (a n) (b n))))
    (hd : ∀ n x, x ∈ Ioo (min (a n) (b n)) (max (a n) (b n)) → HasDerivAt g (g₁ x) x)
    (hg₁ : ContinuousAt g₁ r) (hab : ∀ n, a n ≠ b n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n)))
      atTop (𝓝 (g₁ r)) := by
  have hmvt : ∀ n, ∃ ξ ∈ Ioo (min (a n) (b n)) (max (a n) (b n)),
      g₁ ξ = dd1 (a n) (b n) (g (a n)) (g (b n)) := by
    intro n
    have hgn := hg n
    have hdn := hd n
    rcases lt_or_gt_of_ne (hab n) with h | h
    · simp only [min_eq_left h.le, max_eq_right h.le] at hgn hdn ⊢
      exact exists_hasDerivAt_eq_slope g g₁ h hgn hdn
    · simp only [min_eq_right h.le, max_eq_left h.le] at hgn hdn ⊢
      rw [dd1_comm (a n) (b n)]
      exact exists_hasDerivAt_eq_slope g g₁ h hgn hdn
  choose ξ hξ heq using hmvt
  have hlim : Tendsto ξ atTop (𝓝 r) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le (by simpa using ha.min hb)
      (by simpa using ha.max hb) (fun n => (hξ n).1.le) (fun n => (hξ n).2.le)
  simpa only [Function.comp_def, heq] using hg₁.tendsto.comp hlim

/-- Only C² regularity at the coalescence point is required; arbitrary finite
prefixes and arbitrary node orderings are allowed. -/
theorem tendsto_divided_differences_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd1 (a n) (b n) (g (a n)) (g (b n))) atTop (𝓝 (deriv g r)) ∧
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (deriv (deriv g) r / 2)) := by
  obtain ⟨δ, hδ, hball⟩ := Metric.mem_nhds_iff.mp (hg.eventually (by simp))
  let J := Ioo (r-δ) (r+δ)
  have hrJ : J ∈ 𝓝 r := Ioo_mem_nhds (by linarith) (by linarith)
  have hJ : ∀ x ∈ J, ContDiffAt ℝ 2 g x := by
    intro x hx
    apply hball
    rw [Metric.mem_ball, Real.dist_eq]
    apply abs_lt.mpr
    dsimp [J] at hx
    constructor <;> linarith [hx.1,hx.2]
  obtain ⟨N,hN⟩ := eventually_atTop.mp ((ha.eventually hrJ).and ((hb.eventually hrJ).and (hc.eventually hrJ)))
  have hm (n : ℕ) := hN (n+N) (Nat.le_add_left N n)
  have hreg₃ (n : ℕ) (x : ℝ)
      (hx : x ∈ Icc (min (a (n+N)) (min (b (n+N)) (c (n+N))))
        (max (a (n+N)) (max (b (n+N)) (c (n+N))))) : ContDiffAt ℝ 2 g x := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (lt_min (hm n).2.1.1 (hm n).2.2.1)).trans_le hx.1,
      hx.2.trans_lt (max_lt (hm n).1.2 (max_lt (hm n).2.1.2 (hm n).2.2.2))⟩
  have hreg₂ (n : ℕ) (x : ℝ)
      (hx : x ∈ Icc (min (a (n+N)) (b (n+N))) (max (a (n+N)) (b (n+N)))) :
      ContDiffAt ℝ 2 g x := by
    apply hJ
    exact ⟨(lt_min (hm n).1.1 (hm n).2.1.1).trans_le hx.1,
      hx.2.trans_lt (max_lt (hm n).1.2 (hm n).2.1.2)⟩
  have hg' : ContDiffAt ℝ 1 (deriv g) r := hg.derivWithin (by norm_num)
  have hg'' : ContDiffAt ℝ 0 (deriv (deriv g)) r := hg'.derivWithin (by norm_num)
  constructor
  · apply (tendsto_add_atTop_iff_nat N).mp
    exact tendsto_dd1_arbitrary_of_derivatives
      (fun n x hx => (hreg₂ n x hx).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₂ n x ⟨hx.1.le,hx.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
      hg'.continuousAt (fun n => hab (n+N))
      (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
  · apply (tendsto_add_atTop_iff_nat N).mp
    exact tendsto_dd2_arbitrary_of_derivatives
      (fun n x hx => (hreg₃ n x hx).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₃ n x hx).derivWithin (m:=1) (by norm_num)).continuousAt.continuousWithinAt)
      (fun n x hx => ((hreg₃ n x ⟨hx.1.le,hx.2.le⟩).differentiableAt (by norm_num)).hasDerivAt)
      (fun n x hx => (((hreg₃ n x ⟨hx.1.le,hx.2.le⟩).derivWithin (m:=1) (by norm_num)).differentiableAt (by norm_num)).hasDerivAt)
      hg''.continuousAt (fun n => hab (n+N)) (fun n => hbc (n+N)) (fun n => hac (n+N))
      (ha.comp (tendsto_add_atTop_nat N)) (hb.comp (tendsto_add_atTop_nat N))
      (hc.comp (tendsto_add_atTop_nat N))

/-- The Newton pole value limit under local C² regularity. -/
theorem tendsto_poleValue_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 (g r)) := by
  obtain ⟨h₁,h₂⟩ := tendsto_divided_differences_of_contDiffAt hg hab hbc hac ha hb hc
  have hv := ((hg.continuousAt.tendsto.comp ha).sub
    (h₁.mul (ha.sub_const r))).add ((h₂.mul (ha.sub_const r)).mul (hb.sub_const r))
  simpa only [Function.comp_def, sub_self, mul_zero, sub_zero, add_zero, poleValue] using hv

/-- The exact coefficient limit under local C² regularity. -/
theorem tendsto_error_coefficient_of_contDiffAt
    {g : ℝ → ℝ} {a b c : ℕ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0)
    (hab : ∀ n, a n ≠ b n) (hbc : ∀ n, b n ≠ c n) (hac : ∀ n, a n ≠ c n)
    (ha : Tendsto a atTop (𝓝 r)) (hb : Tendsto b atTop (𝓝 r)) (hc : Tendsto c atTop (𝓝 r)) :
    Tendsto (fun n => dd2 (a n) (b n) (c n) (g (a n)) (g (b n)) (g (c n)) /
      poleValue (a n) (b n) (c n) r (g (a n)) (g (b n)) (g (c n)))
      atTop (𝓝 ((deriv (deriv g) r / 2) / g r)) :=
  (tendsto_divided_differences_of_contDiffAt hg hab hbc hac ha hb hc).2.div
    (tendsto_poleValue_of_contDiffAt hg hab hbc hac ha hb hc) hgr

end
end ModAB.Rational.LocalAnalysis
