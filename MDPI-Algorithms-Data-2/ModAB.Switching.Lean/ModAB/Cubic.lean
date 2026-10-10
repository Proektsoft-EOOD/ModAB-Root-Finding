import ModAB.Switching
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.Deriv.Pow

/-! The manuscript's `lem:switch_cubic`, including its derivative argument. -/
namespace ModAB.Cubic
noncomputable section

def Q (v : ℝ) : ℝ := (1-v)^2*(1+v)

theorem Q_expanded (v : ℝ) : Q v = 1-v-v^2+v^3 := by unfold Q; ring

theorem Q_hasDerivAt (v : ℝ) : HasDerivAt Q ((v-1)*(3*v+1)) v := by
  have h := (((hasDerivAt_const v (1:ℝ)).sub (hasDerivAt_id v)).sub
    ((hasDerivAt_id v).pow 2)).add ((hasDerivAt_id v).pow 3)
  convert h using 1 <;> (try funext x) <;> (try dsimp [Q]) <;> ring

theorem Q_derivative_negative {v : ℝ} (hv : v ∈ Set.Icc (0:ℝ) (1/2)) :
    deriv Q v < 0 := by
  rw [(Q_hasDerivAt v).deriv]
  exact mul_neg_of_neg_of_pos (by linarith [hv.2]) (by linarith [hv.1])

theorem Q_strictAntiOn : StrictAntiOn Q (Set.Icc (0:ℝ) (1/2)) := by
  apply strictAntiOn_of_deriv_neg (convex_Icc _ _)
  · exact (by unfold Q; fun_prop : Continuous Q).continuousOn
  · intro x hx
    apply Q_derivative_negative
    exact interior_subset hx

theorem Q_bounds {v : ℝ} (hv : v ∈ Set.Icc (0:ℝ) (1/2)) :
    (3/8:ℝ) ≤ Q v ∧ Q v ≤ 1 := by
  have lo := Q_strictAntiOn.antitoneOn hv (by norm_num : (1/2:ℝ) ∈ Set.Icc 0 (1/2)) hv.2
  have hi := Q_strictAntiOn.antitoneOn (by norm_num : (0:ℝ) ∈ Set.Icc 0 (1/2)) hv hv.1
  norm_num [Q] at lo hi
  exact ⟨lo,hi⟩

theorem endpoint_difference {f1 f2 : ℝ} (hs : f1*f2 < 0) :
    |f2-f1| = |f1|+|f2| := by
  rcases mul_neg_iff.mp hs with h | h
  · rw [abs_of_pos h.1, abs_of_neg h.2, abs_of_neg (by linarith : f2-f1<0)]
    ring
  · rw [abs_of_neg h.1, abs_of_pos h.2, abs_of_pos (by linarith : 0<f2-f1)]
    ring

theorem endpoint_difference_pos {f1 f2 : ℝ} (hs : f1*f2 < 0) : 0 < |f2-f1| := by
  rw [endpoint_difference hs]
  have h1 : f1 ≠ 0 := by intro h; simp [h] at hs
  exact add_pos_of_pos_of_nonneg (abs_pos.mpr h1) (abs_nonneg _)

theorem midpoint_bound {f1 f2 : ℝ} (hs : f1*f2<0) :
    |(f1+f2)/2| ≤ |f2-f1|/2 := by
  rw [abs_div, abs_of_pos (by norm_num : (0:ℝ)<2), endpoint_difference hs]
  exact div_le_div_of_nonneg_right (abs_add_le _ _) (by norm_num)

theorem ratio_range {f1 f2 : ℝ} (hs : f1*f2<0) :
    |(f1+f2)/2| / |f2-f1| ∈ Set.Icc (0:ℝ) (1/2) := by
  have hd := endpoint_difference_pos hs
  constructor
  · positivity
  · apply (div_le_iff₀ hd).mpr
    linarith [midpoint_bound hs]

/-- Exact distributed expression; no claim about intermediate floating-point overflow. -/
theorem switching_bound {f1 f2 : ℝ} (hs : f1*f2<0) (f3 : ℝ) :
    let D := |f2-f1|
    let y := (f1+f2)/2
    let v := |y|/D
    let H := max D |f3|
    (1-v)^2*|y|+(1-v)^2*|f3| ≤ H*Q v ∧ H*Q v ≤ H := by
  dsimp only
  let D := |f2-f1|
  let y := (f1+f2)/2
  let v := |y|/D
  let H := max D |f3|
  change (1-v)^2*|y|+(1-v)^2*|f3| ≤ H*Q v ∧ H*Q v ≤ H
  have hd : 0<D := endpoint_difference_pos hs
  have hv : v ∈ Set.Icc (0:ℝ) (1/2) := ratio_range hs
  have hDH : D≤H := le_max_left _ _
  have hzH : |f3|≤H := le_max_right _ _
  have hH : 0≤H := (le_of_lt hd).trans hDH
  have hyv : |y|=v*D := by dsimp [v]; field_simp
  have hsum : |y|+|f3|≤H*(v+1) := by
    rw [hyv]
    nlinarith [mul_le_mul_of_nonneg_left hDH hv.1]
  constructor
  · have h := mul_le_mul_of_nonneg_left hsum (sq_nonneg (1-v))
    dsimp [Q]
    nlinarith only [h]
  · exact (mul_le_mul_of_nonneg_left (Q_bounds hv).2 hH).trans_eq (mul_one H)

end
end ModAB.Cubic
