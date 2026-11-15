import ModAB.SourceProof

namespace ModAB.Source
open IEEE64 Binary64
noncomputable section

def mathematicalGap (env : Nat → ℝ) : ℝ :=
  exactG (rho env) (normalized env 0) (normalized env 1) (normalized env 2)

/-- The decimal C# literal is exactly the representable number 2^-48. -/
theorem guard_literal : (3.552713678800500929355621337890625e-15 : ℝ) = 32*binary64Epsilon := by
  norm_num [binary64Epsilon]

theorem guard_exact : Binary64.round (32*binary64Epsilon) = 32*binary64Epsilon := by
  rw [binary64_guard_value]
  have h := round_exact (power_representable (-48) (by norm_num))
  norm_num at h ⊢
  exact h

theorem decision_switch_iff (g : ℝ) : decision g = 3 ↔ 32*binary64Epsilon < g := by
  rw [binary64_guard_value]
  simp only [decision]
  split_ifs <;> simp_all

theorem decision_no_switch_iff (g : ℝ) : decision g = 2 ↔ g < -(32*binary64Epsilon) := by
  have he := binary64Epsilon_pos
  rw [binary64_guard_value]
  simp only [decision]
  split_ifs <;> simp_all <;> linarith

/-- A positive certificate from the binary-word program implies the original strict inequality. -/
theorem binary64_certified_switch (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i)))
    (hc : decision (value (gapExpr.bits env)) = 3) :
    originalL (value (env 0)) (value (env 1)) (value (env 2)) <
      originalR (kappa (rho (fun i => value (env i))))
        (value (env 0)) (value (env 1)) (value (env 2)) := by
  let re : Nat → ℝ := fun i => value (env i)
  have herr := (binary64_program_error env hi hd).2
  have hg := guard_positive binary64Epsilon_pos herr ((decision_switch_iff _).mp hc)
  have hp := (input_facts re hd).1
  apply (normalization_sign hp (kappa (rho re)) (re 0) (re 1) (re 2)).mp
  exact hg

/-- A negative certificate proves the strict reverse inequality. -/
theorem binary64_certified_no_switch (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i)))
    (hc : decision (value (gapExpr.bits env)) = 2) :
    originalR (kappa (rho (fun i => value (env i))))
        (value (env 0)) (value (env 1)) (value (env 2)) <
      originalL (value (env 0)) (value (env 1)) (value (env 2)) := by
  let re : Nat → ℝ := fun i => value (env i)
  have herr := (binary64_program_error env hi hd).2
  have hg := guard_negative binary64Epsilon_pos herr ((decision_no_switch_iff _).mp hc)
  change exactG (rho re) (re 0/scale re) (re 1/scale re) (re 2/scale re) < 0 at hg
  dsimp only [exactG, exactR, exactL, exactH] at hg
  rw [normalization_identity (input_facts re hd).1] at hg
  have h := (div_lt_iff₀ (input_facts re hd).1).mp hg
  change originalR (kappa (rho re)) (re 0) (re 1) (re 2) < originalL (re 0) (re 1) (re 2)
  linarith

/-- The 57ε separation guarantees that the appropriate branch is taken. -/
theorem binary64_margin_guarantees (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i))) :
    (57*binary64Epsilon < mathematicalGap (fun i => value (env i)) →
      decision (value (gapExpr.bits env)) = 3) ∧
    (mathematicalGap (fun i => value (env i)) < -(57*binary64Epsilon) →
      decision (value (gapExpr.bits env)) = 2) := by
  have he := (binary64_program_error env hi hd).2
  exact ⟨fun h => (decision_switch_iff _).mpr (complete_positive he h),
    fun h => (decision_no_switch_iff _).mpr (complete_negative he h)⟩

end
end ModAB.Source
