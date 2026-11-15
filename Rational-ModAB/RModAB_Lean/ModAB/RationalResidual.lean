import ModAB.Convergence
import Mathlib.Topology.Sequences

/-!
Residual convergence on the original compact interval. These results do not
assume that the bracket width tends to zero, and they do not assume uniqueness
until the final corollary. They therefore apply to the uncapped residual rule
without silently importing the capped width theorem.
-/
namespace ModAB.RationalResidual
open Set Filter Topology
noncomputable section

/-- Every cluster point of a residual-convergent sequence in the interval is a
zero. The membership hypothesis also supplies the domain of continuity. -/
theorem limit_is_zero {f : ℝ → ℝ} {a b z : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hz : z ∈ Icc a b)
    (hx : ∀ n, x n ∈ Icc a b) (hlim : Tendsto x atTop (𝓝 z))
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0)) : f z = 0 := by
  have hfseq : Tendsto (fun n => f (x n)) atTop (𝓝 (f z)) :=
    (hf z hz).tendsto.comp
      (tendsto_nhdsWithin_iff.mpr ⟨hlim, Filter.Eventually.of_forall hx⟩)
  have he : |f z| = 0 := tendsto_nhds_unique hfseq.abs hres
  exact abs_eq_zero.mp he

/-- Residual convergence gives proximity to the zero set, even if the
enclosing intervals have a positive limiting width or the zero is not unique.
The witness is an actual zero in the original interval. -/
theorem eventually_near_zero_set {f : ℝ → ℝ} {a b : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hx : ∀ n, x n ∈ Icc a b)
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0))
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ n in atTop, ∃ z ∈ Icc a b, f z = 0 ∧ |x n - z| < ε := by
  by_contra h
  have hbad : ∃ᶠ n in atTop,
      ¬ (∃ z ∈ Icc a b, f z = 0 ∧ |x n - z| < ε) := not_eventually.mp h
  obtain ⟨ψ, hψ, hψbad⟩ := extraction_of_frequently_atTop hbad
  obtain ⟨z, hz, φ, hφ, ht⟩ :=
    isCompact_Icc.tendsto_subseq (fun n => hx (ψ n))
  have ht' : Tendsto (fun n => x (ψ (φ n))) atTop (𝓝 z) := ht
  have hr' : Tendsto (fun n => |f (x (ψ (φ n)))|) atTop (𝓝 0) :=
    hres.comp (hψ.tendsto_atTop.comp hφ.tendsto_atTop)
  have hz0 : f z = 0 := limit_is_zero hf hz (fun n => hx (ψ (φ n))) ht' hr'
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp ht' ε hε
  have hd : |x (ψ (φ N)) - z| < ε := by
    simpa only [Real.dist_eq] using hN N (le_refl N)
  exact hψbad (φ N) ⟨z, hz, hz0, hd⟩

/-- With a unique zero in the original interval, residual convergence implies
convergence of the evaluated points themselves. No derivative bound is used. -/
theorem tendsto_unique_zero {f : ℝ → ℝ} {a b r : ℝ} {x : ℕ → ℝ}
    (hf : ContinuousOn f (Icc a b)) (hx : ∀ n, x n ∈ Icc a b)
    (hres : Tendsto (fun n => |f (x n)|) atTop (𝓝 0))
    (hunique : ∀ z ∈ Icc a b, f z = 0 → z = r) :
    Tendsto x atTop (𝓝 r) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  have he := eventually_near_zero_set hf hx hres hε
  obtain ⟨N, hN⟩ := eventually_atTop.mp he
  refine ⟨N, fun n hn => ?_⟩
  obtain ⟨z, hz, hz0, hd⟩ := hN n hn
  have hzr := hunique z hz hz0
  simpa only [hzr, Real.dist_eq] using hd

end
end ModAB.RationalResidual
