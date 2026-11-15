import ModAB.Correctness

namespace ModAB.Source
open IEEE64 Binary64
open FloatLib.Floats.Formats.BinaryInterchange

/-- Bit fields of the exact C# guard 2^-48. -/
def guardWord : Word := Model.ofFields FloatFormat.binary64 false 975 0

theorem guardWord_valid : valid guardWord :=
  Model.isFinite_ofFields_ieee _ FloatFormat.isIEEE_binary64 false 975 0 (by decide)

theorem guardWord_value : value guardWord = 1/(2:ℝ)^48 := by
  change Model.toReal (Model.ofFields FloatFormat.binary64 false 975 0) = _
  rw [Model.toReal_ofFields_normal _ _ _ _ (by decide) (by decide) (by decide)]
  norm_num [FloatFormat.binary64, Model.pow2, Binary64.pow_value]

/-- Executable IEEE comparisons; all four returned integers agree with the C# enum. -/
def decisionWord (g : Word) : Nat :=
  if Model.compare g guardWord = some .gt then 3
  else if Model.compare g (Model.neg guardWord) = some .lt then 2 else 1

theorem decisionWord_refines (g : Word) (hg : valid g) : decisionWord g = decision (value g) := by
  have hng : valid (Model.neg guardWord) := by simpa [valid] using guardWord_valid
  have hp := Model.compare_eq_some_gt_iff_toReal_gt_of_isFinite g guardWord hg guardWord_valid
  have hn := Model.compare_eq_some_lt_iff_toReal_lt_of_isFinite g (Model.neg guardWord) hg hng
  have gv := guardWord_value
  have ngv : value (Model.neg guardWord) = -(1/(2:ℝ)^48) := by
    change Model.toReal (Model.neg guardWord) = _
    rw [Model.toReal_neg guardWord guardWord_valid]
    exact congrArg Neg.neg gv
  change (Model.compare g guardWord = some .gt ↔ value guardWord < value g) at hp
  change (Model.compare g (Model.neg guardWord) = some .lt ↔ value g < value (Model.neg guardWord)) at hn
  rw [gv] at hp
  rw [ngv] at hn
  simp only [decisionWord, decision, hp, hn]

/-- The certified positive branch uses executable comparisons as well as executable arithmetic. -/
theorem executable_switch_sound (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i)))
    (hc : decisionWord (gapExpr.bits env) = 3) :
    originalL (value (env 0)) (value (env 1)) (value (env 2)) <
      originalR (kappa (rho (fun i => value (env i))))
        (value (env 0)) (value (env 1)) (value (env 2)) := by
  rw [decisionWord_refines _ (binary64_program_error env hi hd).1] at hc
  exact binary64_certified_switch env hi hd hc

theorem executable_no_switch_sound (env : Nat → Word) (hi : ∀ i, valid (env i))
    (hd : Domain (fun i => value (env i)))
    (hc : decisionWord (gapExpr.bits env) = 2) :
    originalR (kappa (rho (fun i => value (env i))))
        (value (env 0)) (value (env 1)) (value (env 2)) <
      originalL (value (env 0)) (value (env 1)) (value (env 2)) := by
  rw [decisionWord_refines _ (binary64_program_error env hi hd).1] at hc
  exact binary64_certified_no_switch env hi hd hc

end ModAB.Source
