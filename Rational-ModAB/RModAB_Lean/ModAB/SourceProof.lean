import ModAB.Source

namespace ModAB.Source
open IEEE64 Binary64
noncomputable section

def scale (env : Nat → ℝ) : ℝ := max (max |env 0| |env 1|) |env 2|
def rho (env : Nat → ℝ) : ℝ := min |env 0| |env 1| / max |env 0| |env 1|
def normalized (env : Nat → ℝ) (i : Nat) : ℝ := env i / scale env

def Domain (env : Nat → ℝ) : Prop :=
  env 0 ≠ 0 ∧ env 1 ≠ 0 ∧
  ((env 0 ≤ 0 ∧ 0 ≤ env 1) ∨ (env 1 ≤ 0 ∧ 0 ≤ env 0))

theorem input_facts (env : Nat → ℝ) (hd : Domain env) :
    0 < scale env ∧ rho env ∈ Set.Icc (0:ℝ) 1 ∧
    |normalized env 0| ≤ 1 ∧ |normalized env 1| ≤ 1 ∧ |normalized env 2| ≤ 1 ∧
    ((normalized env 0 ≤ 0 ∧ 0 ≤ normalized env 1) ∨
      (normalized env 1 ≤ 0 ∧ 0 ≤ normalized env 0)) := by
  have h0 : 0 < |env 0| := abs_pos.mpr hd.1
  have h1 : 0 < |env 1| := abs_pos.mpr hd.2.1
  have hp : 0 < scale env := lt_of_lt_of_le h0
    ((le_max_left _ _).trans (le_max_left _ _))
  have hb (i : Nat) (hi : |env i| ≤ scale env) : |normalized env i| ≤ 1 := by
    dsimp [normalized]
    rw [abs_div, abs_of_pos hp]
    exact (div_le_one hp).mpr hi
  refine ⟨hp, endpoint_rho_range h0 h1,
    hb 0 ((le_max_left _ _).trans (le_max_left _ _)),
    hb 1 ((le_max_right _ _).trans (le_max_left _ _)), hb 2 (le_max_right _ _), ?_⟩
  rcases hd.2.2 with ⟨ha,hb⟩ | ⟨ha,hb⟩
  · exact Or.inl ⟨div_nonpos_of_nonpos_of_nonneg ha hp.le, div_nonneg hb hp.le⟩
  · exact Or.inr ⟨div_nonpos_of_nonpos_of_nonneg ha hp.le, div_nonneg hb hp.le⟩

def trace (env : Nat → ℝ) (hd : Domain env) :
    ScaledTrace binary64Epsilon (rho env) (normalized env 0) (normalized env 1) (normalized env 2) :=
  let h := input_facts env hd
  scaledTrace (rho env) (normalized env 0) (normalized env 1) (normalized env 2)
    h.2.1 h.2.2.1 h.2.2.2.1 h.2.2.2.2.1 h.2.2.2.2.2

/-- Equality checks the actual nesting of all 17 rounding sites, including both halvings. -/
theorem source_eq_trace (env : Nat → ℝ) (hd : Domain env) :
    gapExpr.real env = (trace env hd).gapHat := by rfl

/-- No overflow and no division by zero anywhere in the translated arithmetic expression. -/
theorem source_safe (env : Nat → ℝ) (hd : Domain env) : gapExpr.Safe env := by
  let t := trace env hd
  obtain ⟨hscale, hrho, hp1, hp2, hp3, hsign⟩ := input_facts env hd
  have srho : rhoExpr.Safe env := by
    refine ⟨⟨True.intro,True.intro⟩, ⟨True.intro,True.intro⟩, ?_, ?_⟩
    · exact abs_le_two_of_unit hrho
    · intro _
      change max |env 0| |env 1| ≠ 0
      exact ne_of_gt (lt_of_lt_of_le (abs_pos.mpr hd.1) (le_max_left _ _))
  have sa : aExpr.Safe env := by
    refine ⟨True.intro, srho, ?_, by intro h; cases h⟩
    change |1-t.factor.rhoHat| ≤ 2
    apply abs_le.mpr; constructor <;> linarith [t.factor.rho_range.1,t.factor.rho_range.2]
  have sb : bExpr.Safe env := by
    refine ⟨True.intro, srho, ?_, by intro h; cases h⟩
    change |1+t.factor.rhoHat| ≤ 2
    apply abs_le.mpr; constructor <;> linarith [t.factor.rho_range.1,t.factor.rho_range.2]
  have sah : (Expr.op .div aExpr .two).Safe env := by
    refine ⟨sa,True.intro,?_,by intro _; norm_num [Expr.real]⟩
    change |t.factor.aHat/2| ≤ 2
    apply abs_le.mpr; constructor <;> linarith [t.factor.a_range.1,t.factor.a_range.2]
  have sv : vExpr.Safe env := by
    refine ⟨sah,sb,?_,?_⟩
    · change |Binary64.round (t.factor.aHat/2)/t.factor.bHat| ≤ 2
      have half : Binary64.round (t.factor.aHat/2) = t.factor.aHat/2 := factor_halving_exact hrho
      rw [half]
      have hb : 0 < t.factor.bHat := by linarith [t.factor.b_range.1]
      rw [abs_of_nonneg (div_nonneg (by linarith [t.factor.a_range.1]) hb.le)]
      apply (div_le_iff₀ hb).mpr
      linarith [t.factor.a_range.2,t.factor.b_range.1]
    · intro _; change t.factor.bHat ≠ 0; linarith [t.factor.b_range.1]
  have sr : rExpr.Safe env := by
    refine ⟨True.intro,sv,?_,by intro h; cases h⟩
    change |1-t.factor.vHat| ≤ 2
    apply abs_le.mpr; constructor <;> linarith [t.factor.v_range.1,t.factor.v_range.2]
  have sk : kExpr.Safe env := by
    refine ⟨sr,sr,?_,by intro h; cases h⟩
    change |t.factor.rHat*t.factor.rHat| ≤ 2
    apply abs_le.mpr; constructor <;> nlinarith [t.factor.r_range.1,t.factor.r_range.2]
  have ss : scaleExpr.Safe env := ⟨⟨True.intro,True.intro⟩,True.intro⟩
  have sp1 : p1Expr.Safe env :=
    ⟨True.intro,ss,by change |normalized env 0| ≤ 2; linarith, fun _ => ne_of_gt hscale⟩
  have sp2 : p2Expr.Safe env :=
    ⟨True.intro,ss,by change |normalized env 1| ≤ 2; linarith, fun _ => ne_of_gt hscale⟩
  have sp3 : p3Expr.Safe env :=
    ⟨True.intro,ss,by change |normalized env 2| ≤ 2; linarith, fun _ => ne_of_gt hscale⟩
  have sadd : (Expr.op .add p1Expr p2Expr).Safe env := by
    refine ⟨sp1,sp2,?_,by intro h; cases h⟩
    change |t.p1Hat+t.p2Hat| ≤ 2
    exact (abs_add_le _ _).trans (by linarith [t.p1_range,t.p2_range])
  have sh : hExpr.Safe env := by
    refine ⟨sadd,True.intro,?_,by intro _; norm_num [Expr.real]⟩
    change |t.sumHat/2| ≤ 2
    rw [abs_div]; norm_num; linarith [t.sum_range]
  have sl : leftExpr.Safe env := by
    refine ⟨sp3,sh,?_,by intro h; cases h⟩
    change |t.p3Hat-t.hHat| ≤ 2
    exact (abs_sub _ _).trans (by linarith [t.p3_range,t.h_range])
  have sm1 : term1Expr.Safe env := by
    refine ⟨sk,sh,?_,by intro h; cases h⟩
    change |t.factor.kHat * abs t.hHat| ≤ 2
    rw [abs_mul, abs_abs, abs_of_nonneg (by linarith [t.factor.k_range.1] : 0 ≤ t.factor.kHat)]
    nlinarith [t.factor.k_range.1,t.factor.k_range.2,t.h_range,abs_nonneg t.hHat]
  have sm2 : term2Expr.Safe env := by
    refine ⟨sk,sp3,?_,by intro h; cases h⟩
    change |t.factor.kHat * abs t.p3Hat| ≤ 2
    rw [abs_mul, abs_abs, abs_of_nonneg (by linarith [t.factor.k_range.1] : 0 ≤ t.factor.kHat)]
    nlinarith [t.factor.k_range.1,t.factor.k_range.2,t.p3_range,abs_nonneg t.p3Hat]
  have hm1 : t.productH ∈ Set.Icc (0:ℝ) (1/2) := by
    change Binary64.round (t.factor.kHat * |t.hHat|) ∈ Set.Icc (0:ℝ) (1/2)
    apply round_interval round_zero round_half
    constructor
    · exact mul_nonneg (by linarith [t.factor.k_range.1]) (abs_nonneg _)
    · nlinarith [t.factor.k_range.1,t.factor.k_range.2,t.h_range,abs_nonneg t.hHat]
  have hm2 : t.productP ∈ Set.Icc (0:ℝ) 1 := by
    change Binary64.round (t.factor.kHat * |t.p3Hat|) ∈ Set.Icc (0:ℝ) 1
    apply round_interval round_zero round_one
    constructor
    · exact mul_nonneg (by linarith [t.factor.k_range.1]) (abs_nonneg _)
    · nlinarith [t.factor.k_range.1,t.factor.k_range.2,t.p3_range,abs_nonneg t.p3Hat]
  have sright : rightExpr.Safe env := by
    refine ⟨sm1,sm2,?_,by intro h; cases h⟩
    change |t.productH+t.productP| ≤ 2
    apply abs_le.mpr; constructor <;> linarith [hm1.1,hm1.2,hm2.1,hm2.2]
  have hright : t.rightHat ∈ Set.Icc (0:ℝ) 2 := by
    change Binary64.round (t.productH+t.productP) ∈ Set.Icc (0:ℝ) 2
    apply round_interval round_zero round_two
    constructor <;> linarith [hm1.1,hm1.2,hm2.1,hm2.2]
  have hdiff : |t.diffHat| ≤ 2 := by
    change |Binary64.round (t.p3Hat-t.hHat)| ≤ 2
    apply round_abs_le_two
    exact (abs_sub _ _).trans (by linarith [t.p3_range,t.h_range])
  refine ⟨sright,sl,?_,by intro h; cases h⟩
  change |t.rightHat - abs t.diffHat| ≤ 2
  apply abs_le.mpr; constructor <;> linarith [hright.1,hright.2,hdiff,abs_nonneg t.diffHat]

/-- End-to-end connection from the executable IEEE word expression to the mathematical trace. -/
theorem source_word_refines (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i))) :
    valid (gapExpr.bits env) ∧ value (gapExpr.bits env) =
      (trace (fun i => value (env i)) hd).gapHat := by
  obtain ⟨hf,he⟩ := expression_refines gapExpr env hi (source_safe _ hd)
  exact ⟨hf,he.trans (source_eq_trace _ hd)⟩

theorem source_error (env : Nat → ℝ) (hd : Domain env) :
    Within (gapExpr.real env)
      (exactG (rho env) (normalized env 0) (normalized env 1) (normalized env 2))
      (25*binary64Epsilon) := by
  rw [source_eq_trace env hd]
  obtain ⟨_,hr,h1,h2,h3,hs⟩ := input_facts env hd
  exact concrete_gap_error _ _ _ _ hr h1 h2 h3 hs

/-- Main binary64 error theorem: only finite input data and the bracket sign condition are assumed. -/
theorem binary64_program_error (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i))) :
    valid (gapExpr.bits env) ∧ Within (value (gapExpr.bits env))
      (exactG (rho (fun i => value (env i)))
        (normalized (fun i => value (env i)) 0)
        (normalized (fun i => value (env i)) 1)
        (normalized (fun i => value (env i)) 2)) (25*binary64Epsilon) := by
  obtain ⟨hf,he⟩ := expression_refines gapExpr env hi (source_safe _ hd)
  exact ⟨hf,by rw [he]; exact source_error _ hd⟩

/-- Source-level error bound retaining the main theorem's mean-value proof route. -/
theorem source_error_mean_value (env : Nat → ℝ) (hd : Domain env) :
    Within (gapExpr.real env)
      (exactG (rho env) (normalized env 0) (normalized env 1) (normalized env 2))
      (25*binary64Epsilon) := by
  rw [source_eq_trace env hd]
  obtain ⟨_,hr,h1,h2,h3,hs⟩ := input_facts env hd
  exact concrete_gap_error_mean_value _ _ _ _ hr h1 h2 h3 hs

/-- The concrete binary64 theorem with the derivative/mean-value κ argument. -/
theorem binary64_program_error_mean_value (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i))) :
    valid (gapExpr.bits env) ∧ Within (value (gapExpr.bits env))
      (exactG (rho (fun i => value (env i)))
        (normalized (fun i => value (env i)) 0)
        (normalized (fun i => value (env i)) 1)
        (normalized (fun i => value (env i)) 2)) (25*binary64Epsilon) := by
  obtain ⟨hf,he⟩ := expression_refines gapExpr env hi (source_safe _ hd)
  exact ⟨hf,by rw [he]; exact source_error_mean_value _ hd⟩

end
end ModAB.Source
