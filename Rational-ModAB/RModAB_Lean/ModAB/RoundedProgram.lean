import ModAB.Binary64

namespace ModAB.Binary64
noncomputable section

lemma abs_le_two_of_unit {x : ℝ} (hx : x ∈ Set.Icc (0:ℝ) 1) : |x| ≤ 2 := by
  rw [abs_of_nonneg hx.1]; linarith [hx.2]

lemma abs_le_two_of_half {x : ℝ} (hx : x ∈ Set.Icc (0:ℝ) (1/2)) : |x| ≤ 2 := by
  rw [abs_of_nonneg hx.1]; linarith [hx.2]

lemma factor_ranges {rho : ℝ} (hr : rho ∈ Set.Icc (0:ℝ) 1) :
    let rh := round rho
    let a := round (1-rh)
    let b := round (1+rh)
    let v := round (round (a/2)/b)
    let r := round (1-v)
    let k := round (r*r)
    rh ∈ Set.Icc (0:ℝ) 1 ∧ a ∈ Set.Icc (0:ℝ) 1 ∧
    b ∈ Set.Icc (1:ℝ) 2 ∧ v ∈ Set.Icc (0:ℝ) (1/2) ∧
    r ∈ Set.Icc (1/2:ℝ) 1 ∧ k ∈ Set.Icc (1/4:ℝ) 1 := by
  dsimp only
  have hrh := round_interval round_zero round_one hr
  have ha := round_interval round_zero round_one
    (show 1-round rho ∈ Set.Icc (0:ℝ) 1 by constructor <;> linarith [hrh.1, hrh.2])
  have hb := round_interval round_one round_two
    (show 1+round rho ∈ Set.Icc (1:ℝ) 2 by constructor <;> linarith [hrh.1, hrh.2])
  have hbpos : 0 < round (1+round rho) := by linarith [hb.1]
  have hvarg : round (round (1-round rho)/2)/round (1+round rho) ∈
      Set.Icc (0:ℝ) (1/2) := by
    rw [factor_halving_exact hr]
    constructor
    · exact div_nonneg (by linarith [ha.1]) (le_of_lt hbpos)
    · apply (div_le_iff₀ hbpos).mpr
      linarith [ha.2, hb.1]
  have hv := round_interval round_zero round_half hvarg
  have hrarg : 1-round (round (round (1-round rho)/2)/round (1+round rho)) ∈
      Set.Icc (1/2:ℝ) 1 := by constructor <;> linarith [hv.1, hv.2]
  have hrr := round_interval round_half round_one hrarg
  have hkarg : round (1-round (round (round (1-round rho)/2)/round (1+round rho))) *
      round (1-round (round (round (1-round rho)/2)/round (1+round rho))) ∈
      Set.Icc (1/4:ℝ) 1 := by constructor <;> nlinarith [hrr.1, hrr.2]
  exact ⟨hrh, ha, hb, hv, hrr, round_interval round_quarter round_one hkarg⟩

/-- Every field and every local error contract is derived from concrete rounding. -/
def factorTrace (rho : ℝ) (hr : rho ∈ Set.Icc (0:ℝ) 1) :
    FactorTrace binary64Epsilon rho := by
  let rh := round rho
  let a := round (1-rh)
  let b := round (1+rh)
  let v := round (round (a/2)/b)
  let r := round (1-v)
  let k := round (r*r)
  have ranges := factor_ranges hr
  change rh ∈ Set.Icc (0:ℝ) 1 ∧ a ∈ Set.Icc (0:ℝ) 1 ∧
    b ∈ Set.Icc (1:ℝ) 2 ∧ v ∈ Set.Icc (0:ℝ) (1/2) ∧
    r ∈ Set.Icc (1/2:ℝ) 1 ∧ k ∈ Set.Icc (1/4:ℝ) 1 at ranges
  have hrh := ranges.1
  have ha := ranges.2.1
  have hb := ranges.2.2.1
  have hv := ranges.2.2.2.1
  have hrr := ranges.2.2.2.2.1
  have hk := ranges.2.2.2.2.2
  have half : round (a/2) = a/2 := factor_halving_exact hr
  have hquot : (a/2)/b ∈ Set.Icc (0:ℝ) (1/2) := by
    have bp : 0 < b := by linarith [hb.1]
    constructor
    · exact div_nonneg (by linarith [ha.1]) (le_of_lt bp)
    · apply (div_le_iff₀ bp).mpr; linarith [ha.2, hb.1]
  refine {
    rhoHat := rh, aHat := a, bHat := b, vHat := v, rHat := r, kHat := k
    rho_range := hrh, a_range := ha, b_range := hb, v_range := hv
    r_range := hrr, k_range := hk
    rho_error := round_error (abs_le_two_of_unit hr)
    a_error := ?_, b_error := ?_, v_error := ?_, r_error := ?_, k_error := ?_ }
  · exact round_error (abs_le_two_of_unit (by constructor <;> linarith [hrh.1, hrh.2]))
  · apply round_error; apply abs_le.mpr; constructor <;> linarith [hrh.1, hrh.2]
  · dsimp only [v]; rw [half]; exact round_error (abs_le_two_of_half hquot)
  · apply round_error; apply abs_le.mpr; constructor <;> linarith [hv.1, hv.2]
  · change Within (round (r*r)) (r^2) _
    rw [pow_two]
    apply round_error; apply abs_le.mpr; constructor <;> nlinarith [hrr.1, hrr.2]

lemma opposite_sum_bound {x y : ℝ} (hx : |x| ≤ 1) (hy : |y| ≤ 1)
    (hs : (x ≤ 0 ∧ 0 ≤ y) ∨ (y ≤ 0 ∧ 0 ≤ x)) : |x+y| ≤ 1 := by
  rcases abs_le.mp hx with ⟨hxl,hxu⟩
  rcases abs_le.mp hy with ⟨hyl,hyu⟩
  apply abs_le.mpr
  rcases hs with ⟨ha,hb⟩ | ⟨ha,hb⟩ <;> constructor <;> linarith

lemma round_signs {x y : ℝ}
    (hs : (x ≤ 0 ∧ 0 ≤ y) ∨ (y ≤ 0 ∧ 0 ≤ x)) :
    (round x ≤ 0 ∧ 0 ≤ round y) ∨ (round y ≤ 0 ∧ 0 ≤ round x) := by
  have hz := round_zero
  rcases hs with ⟨ha,hb⟩ | ⟨ha,hb⟩
  · exact Or.inl ⟨by simpa using round_monotone ha, by simpa using round_monotone hb⟩
  · exact Or.inr ⟨by simpa using round_monotone ha, by simpa using round_monotone hb⟩

/-- Concrete evaluation of all 17 binary arithmetic operations in the valid-input path. -/
def scaledTrace (rho p1 p2 p3 : ℝ) (hr : rho ∈ Set.Icc (0:ℝ) 1)
    (h1 : |p1| ≤ 1) (h2 : |p2| ≤ 1) (h3 : |p3| ≤ 1)
    (hs : (p1 ≤ 0 ∧ 0 ≤ p2) ∨ (p2 ≤ 0 ∧ 0 ≤ p1)) :
    ScaledTrace binary64Epsilon rho p1 p2 p3 := by
  let ft := factorTrace rho hr
  let q1 := round p1
  let q2 := round p2
  let q3 := round p3
  let s := round (q1+q2)
  let h := round (s/2)
  let d := round (q3-h)
  let u := round (ft.kHat*|h|)
  let v := round (ft.kHat*|q3|)
  let w := round (u+v)
  let g := round (w-|d|)
  have hq1 : |q1| ≤ 1 := round_abs_le_one h1
  have hq2 : |q2| ≤ 1 := round_abs_le_one h2
  have hq3 : |q3| ≤ 1 := round_abs_le_one h3
  have hsarg : |q1+q2| ≤ 1 := opposite_sum_bound hq1 hq2 (round_signs hs)
  have hsum : |s| ≤ 1 := round_abs_le_one hsarg
  have hharg : |s/2| ≤ 1/2 := by rw [abs_div]; norm_num; linarith
  have hh : |h| ≤ 1/2 := abs_le.mpr
    (round_interval round_neg_half round_half (abs_le.mp hharg))
  have hdarg : |q3-h| ≤ 2 := by
    rcases abs_le.mp hq3 with ⟨hl,hu⟩
    rcases abs_le.mp hh with ⟨hl',hu'⟩
    apply abs_le.mpr; constructor <;> linarith
  have hd : |d| ≤ 2 := round_abs_le_two hdarg
  have hk := ft.k_range
  have huarg : ft.kHat*|h| ∈ Set.Icc (0:ℝ) (1/2) := by
    constructor
    · exact mul_nonneg (by linarith [hk.1]) (abs_nonneg _)
    · nlinarith [hk.1, hk.2, abs_nonneg h]
  have hvarg : ft.kHat*|q3| ∈ Set.Icc (0:ℝ) 1 := by
    constructor
    · exact mul_nonneg (by linarith [hk.1]) (abs_nonneg _)
    · nlinarith [hk.1, hk.2, abs_nonneg q3]
  have hu : u ∈ Set.Icc (0:ℝ) (1/2) := round_interval round_zero round_half huarg
  have hv : v ∈ Set.Icc (0:ℝ) 1 := round_interval round_zero round_one hvarg
  have hwarg : u+v ∈ Set.Icc (0:ℝ) 2 := by constructor <;> linarith [hu.1,hu.2,hv.1,hv.2]
  have hw : w ∈ Set.Icc (0:ℝ) 2 := round_interval round_zero round_two hwarg
  have hgarg : |w - abs d| ≤ 2 := by
    apply abs_le.mpr; constructor <;> linarith [hw.1, hw.2, abs_nonneg d]
  refine {
    factor := ft, p1Hat := q1, p2Hat := q2, p3Hat := q3
    sumHat := s, hHat := h, diffHat := d, productH := u, productP := v
    rightHat := w, gapHat := g
    p1_range := hq1, p2_range := hq2, p3_range := hq3, sum_range := hsum, h_range := hh
    p1_error := round_error (by linarith), p2_error := round_error (by linarith)
    p3_error := round_error (by linarith), sum_error := round_error (by linarith)
    h_error := round_error (by linarith), diff_error := round_error hdarg
    productH_error := round_error (abs_le_two_of_half huarg)
    productP_error := round_error (abs_le_two_of_unit hvarg)
    right_error := round_error (by rw [abs_of_nonneg hwarg.1]; exact hwarg.2)
    gap_error := round_error hgarg }

/-- No error or range assumptions about an execution trace occur in this theorem. -/
theorem concrete_gap_error (rho p1 p2 p3 : ℝ) (hr : rho ∈ Set.Icc (0:ℝ) 1)
    (h1 : |p1| ≤ 1) (h2 : |p2| ≤ 1) (h3 : |p3| ≤ 1)
    (hs : (p1 ≤ 0 ∧ 0 ≤ p2) ∨ (p2 ≤ 0 ∧ 0 ≤ p1)) :
    Within (scaledTrace rho p1 p2 p3 hr h1 h2 h3 hs).gapHat
      (exactG rho p1 p2 p3) (25*binary64Epsilon) :=
  (scaled_errors_algebraic (le_of_lt binary64Epsilon_pos) hr
    (normalized_midpoint_bound h1 h2 hs) h3 (scaledTrace rho p1 p2 p3 hr h1 h2 h3 hs)).2.2.2.2

/-- The earlier scaled-error proof uses the derivative/mean-value route for κ. -/
theorem concrete_gap_error_mean_value (rho p1 p2 p3 : ℝ) (hr : rho ∈ Set.Icc (0:ℝ) 1)
    (h1 : |p1| ≤ 1) (h2 : |p2| ≤ 1) (h3 : |p3| ≤ 1)
    (hs : (p1 ≤ 0 ∧ 0 ≤ p2) ∨ (p2 ≤ 0 ∧ 0 ≤ p1)) :
    Within (scaledTrace rho p1 p2 p3 hr h1 h2 h3 hs).gapHat
      (exactG rho p1 p2 p3) (25*binary64Epsilon) :=
  (scaled_errors (le_of_lt binary64Epsilon_pos) hr
    (normalized_midpoint_bound h1 h2 hs) h3 (scaledTrace rho p1 p2 p3 hr h1 h2 h3 hs)).2.2.2.2

end
end ModAB.Binary64
