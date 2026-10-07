import ModAB.RoundedProgram
import FloatLib.Floats.Formats.BinaryInterchange.Arithmetic.DivisionSemantics
import FloatLib.Floats.Formats.BinaryInterchange.Arithmetic.SignedSemantics.Subtraction
import FloatLib.Floats.Formats.BinaryInterchange.Arithmetic.Constants
import FloatLib.Floats.Formats.BinaryInterchange.Format.Catalog

namespace ModAB.IEEE64
open FloatLib.Floats.Formats.BinaryInterchange
noncomputable section
abbrev Word := Model FloatFormat.binary64
abbrev value (x : Word) : ℝ := Model.toReal x
abbrev valid (x : Word) : Prop := Model.isFinite x = true

/-- The descriptor has 11 exponent bits and 52 fraction bits; this equality
identifies its gradual-underflow rounding with our concrete real-number map. -/
theorem descriptor_round (x : ℝ) : Model.roundAt FloatFormat.binary64 x = Binary64.round x := by
  rfl

theorem descriptor_max : value (Model.posMaxFinite FloatFormat.binary64) = Binary64.maxFinite := by
  change Model.toReal (Model.posMaxFinite FloatFormat.binary64) = _
  rw [Model.toReal_posMaxFinite]
  norm_num [FloatFormat.binary64, FloatFormat.maxFiniteFracField,
    FloatFormat.fracMaskNat, FloatFormat.maxNormalExponent,
    FloatFormat.maxFiniteExpField, FloatFormat.expMaskNat, Model.pow2, FloatFormat.Encoding.maxFiniteExponent,
    Model.bpow, Binary64.maxFinite, Binary64.pow_value]

theorem two_le_max : (2:ℝ) ≤ value (Model.posMaxFinite FloatFormat.binary64) := by
  rw [descriptor_max]
  exact (Binary64.round_finite (x := 2) (by norm_num)).2.trans' (by simp)

theorem add_refines (x y : Word) (hx : valid x) (hy : valid y)
    (h : |value x+value y| ≤ 2) :
    valid (Model.add x y) ∧ value (Model.add x y) = Binary64.round (value x+value y) := by
  have hf := Model.isFinite_add_of_abs_toReal_add_le_posMaxFinite x y
    FloatFormat.isIEEE_binary64 hx hy (h.trans two_le_max)
  exact ⟨hf, (Model.toReal_add_eq_roundAt x y FloatFormat.isIEEE_binary64 hx hy hf).trans
    (descriptor_round _)⟩

theorem sub_refines (x y : Word) (hx : valid x) (hy : valid y)
    (h : |value x-value y| ≤ 2) :
    valid (Model.sub x y) ∧ value (Model.sub x y) = Binary64.round (value x-value y) := by
  have hf := Model.isFinite_sub_of_abs_toReal_sub_le_posMaxFinite x y
    FloatFormat.isIEEE_binary64 hx hy (h.trans two_le_max)
  exact ⟨hf, (Model.toReal_sub_eq_roundAt x y FloatFormat.isIEEE_binary64 hx hy hf).trans
    (descriptor_round _)⟩

theorem mul_refines (x y : Word) (hx : valid x) (hy : valid y)
    (h : |value x*value y| ≤ 2) :
    valid (Model.mul x y) ∧ value (Model.mul x y) = Binary64.round (value x*value y) := by
  have hb : |value x| * |value y| ≤ value (Model.posMaxFinite FloatFormat.binary64) := by
    rw [← abs_mul]; exact h.trans two_le_max
  have hf := Model.isFinite_mul_of_abs_mul_le_posMaxFinite x y
    FloatFormat.isIEEE_binary64 hx hy hb
  exact ⟨hf, (Model.toReal_mul_eq_roundAt x y FloatFormat.isIEEE_binary64 hx hy hf).trans
    (descriptor_round _)⟩

theorem div_refines (x y : Word) (hx : valid x) (hy : valid y)
    (hy0 : value y ≠ 0) (h : |value x/value y| ≤ 2) :
    valid (Model.div x y) ∧ value (Model.div x y) = Binary64.round (value x/value y) := by
  have hz : Model.isZero y = false := by
    cases he : Model.isZero y with
    | false => rfl
    | true => exact False.elim (hy0 (Model.toReal_eq_zero_of_isZero y he))
  have hf := Model.isFinite_div_of_abs_div_le_posMaxFinite x y
    FloatFormat.isIEEE_binary64 hx hy hz (h.trans two_le_max)
  exact ⟨hf, (Model.toReal_div_eq_roundAt x y FloatFormat.isIEEE_binary64 hx hy hz hf).trans
    (descriptor_round _)⟩
end
end ModAB.IEEE64
