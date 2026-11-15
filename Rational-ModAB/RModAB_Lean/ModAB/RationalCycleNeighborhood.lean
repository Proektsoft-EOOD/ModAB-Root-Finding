import ModAB.RationalCycleState
import ModAB.RationalCycleUniformBridge
import ModAB.RationalCycleSafeguards

/-! One common neighbourhood for the exact negative-K pair and every source
safeguard. All fields are derived from C³, a simple root and nonzero curvature.
In particular neither future coefficient bounds nor guard decisions are input.
-/
namespace ModAB.Rational.CycleNeighborhood
open Set Filter
open scoped Topology
open Regularity CycleUniform CycleState CycleSafeguards SecantLocal
noncomputable section

structure PairNeighborhood (f : ℝ → ℝ) (r θ ε C δ : ℝ) : Prop where
  bounds : Bounds (regularizedReciprocal f r) r θ ε C
  guards : ScaledNeighborhood f r C ε δ
  same_sign : ∀ x y, |x-r| ≤ ε → |y-r| ≤ ε →
    (0 < f x*f y ↔ 0 < (x-r)*(y-r))
  raw_pair : ∀ a s d, CycleEntry.Cone θ ε a s d →
    let g := regularizedReciprocal f r
    let K := AnalyticBridge.errorConstant f r
    let b := firstError g r a s d
    let μ := 1-f (r+b)/f (r+a)
    firstCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
    secondCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
    denominator (r+a) (r+s) (f (r+a)) (f (r+s)) (r+d) (f (r+d)) ≠ 0 ∧
    rawPoint (r+a) (r+s) (f (r+a)) (f (r+s)) (r+d) (f (r+d))-r = b ∧
    μ ∈ Ioo (0:ℝ) 1 ∧
    denominator (r+b) (r+s) (f (r+b)) (μ*f (r+s)) (r+a) (f (r+a)) ≠ 0 ∧
    rawPoint (r+b) (r+s) (f (r+b)) (μ*f (r+s)) (r+a) (f (r+a))-r =
      secondCoefficient g r a s d*(a*b*d) ∧
    b ≠ 0 ∧ b ≠ a ∧ b ≠ s ∧ 0 < a*b

/-- Simultaneous constants for the interpolation coefficients, local cone,
actual AB update, scaled secant, curvature and true-residual tests. -/
theorem exists_pair_neighborhood {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : AnalyticBridge.errorConstant f r < 0) :
    ∃ θ ε C δ, PairNeighborhood f r θ ε C δ := by
  obtain ⟨θ,ε₀,hθ0,hθ1,hε₀,hraw⟩ :=
    CycleUniformBridge.smooth_exists_uniform_raw_pair hf hz hd hK
  let C := -2*AnalyticBridge.errorConstant f r
  have hC : 0 < C := by dsimp [C]; linarith
  obtain ⟨δ,hδ,ε₁,hε₁,hguards⟩ := exists_scaled_neighborhood
    (hf.of_le (by norm_num)) hz hd h₂ hC.le
  obtain ⟨ε₂,hε₂,hcompare⟩ := exists_radius_simple_root_bounds
    (hf.differentiableAt (by norm_num)) hz hd
  let ε := min ε₀ (min ε₁ (min (ε₂/2) (min 1 (θ/C))))
  have hε : 0 < ε := by dsimp [ε]; positivity
  have e₀ : ε ≤ ε₀ := by simp [ε,min_le_iff]
  have e₁ : ε ≤ ε₁ := by simp [ε,min_le_iff]
  have e₂ : ε < ε₂ := by
    have hh : ε ≤ ε₂/2 := by simp [ε,min_le_iff]
    linarith
  have eone : ε ≤ 1 := by simp [ε,min_le_iff]
  have eC : ε ≤ θ/C := by simp [ε,min_le_iff]
  have hsmall : C*ε^2 ≤ θ := by
    have hCe : ε*C ≤ θ := (le_div_iff₀ hC).mp eC
    have hs := mul_nonneg (mul_nonneg hC.le hε.le) (sub_nonneg.mpr eone)
    nlinarith
  have hrawε := fun a s d (hc : CycleEntry.Cone θ ε a s d) =>
    hraw a s d (hc.enlarge_radius e₀)
  refine ⟨θ,ε,C,δ,{
    bounds := ⟨hθ0,hθ1,hε,hC,hsmall,?_⟩
    guards := hguards.shrink hε e₁
    same_sign := ?_
    raw_pair := hrawε
  }⟩
  · intro a s d hc
    have hu := (hrawε a s d hc).1
    have hv := (hrawε a s d hc).2.1
    have hU : firstCoefficient (regularizedReciprocal f r) r a s d < 0 := by linarith [hu.2]
    have hV : secondCoefficient (regularizedReciprocal f r) r a s d < 0 := by linarith [hv.2]
    refine ⟨hU,hV,?_,?_⟩
    · rw [abs_of_neg hU]
      dsimp [C]
      linarith [hu.1]
    · rw [abs_of_neg hV]
      dsimp [C]
      linarith [hv.1]
  · intro x y hx hy
    exact (sign_product_of_normalizedSlope hz hd
      (hcompare x (hx.trans_lt e₂)).1 (hcompare y (hy.trans_lt e₂)).1).1

end
end ModAB.Rational.CycleNeighborhood
