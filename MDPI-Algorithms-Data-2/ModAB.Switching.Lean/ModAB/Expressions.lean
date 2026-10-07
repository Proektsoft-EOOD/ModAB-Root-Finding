import ModAB.IEEE64
import FloatLib.Floats.Formats.BinaryInterchange.Operations.Compare.Proof

namespace ModAB.IEEE64
open FloatLib.Floats.Formats.BinaryInterchange
noncomputable section

abbrev one : Word := Model.posOne FloatFormat.binary64
abbrev two : Word := Model.add one one

lemma one_valid : valid one := Model.isFinite_posOne _
lemma one_value : value one = 1 := Model.toReal_posOne _
lemma two_valid : valid two := (add_refines one one one_valid one_valid (by rw [one_value]; norm_num)).1
lemma two_value : value two = 2 := by
  have h := (add_refines one one one_valid one_valid (by rw [one_value]; norm_num)).2
  simpa only [one_value, show (1:ℝ)+1=2 by norm_num, Binary64.round_two] using h

/- The sign-decoding argument is adapted from FloatLib (MIT):
   see THIRD_PARTY_NOTICES/FloatLib-LICENSE.txt. -/
lemma abs_value (x : Word) (hx : valid x) : value (Model.abs x) = |value x| := by
  have hsigned : FloatFormat.binary64.supportsSignedZero = true := by decide
  obtain ⟨d, hd⟩ := Model.exists_toDyadic?_of_isFinite hx
  have hsign := Model.sign_eq_signBit_of_toDyadic?_some hd
  have hreal : value x = d.toReal := by simp [Model.toReal_eq, hd]
  have hp := FloatLib.Floats.Formats.Flocq.bpow.pos FloatLib.Numerics.binaryRadix d.exponent
  cases hs : Model.signBit x with
  | false =>
    have hn : 0 ≤ value x := by
      rw [hreal]
      simp only [FloatLib.Numerics.Dyadic.toReal, FloatLib.Numerics.Dyadic.signedSignificand,
        hsign.trans hs, Bool.false_eq_true, ↓reduceIte]
      exact mul_nonneg (Nat.cast_nonneg _) hp.le
    simpa [Model.abs, Model.copySign, hsigned, hs] using (abs_of_nonneg hn).symm
  | true =>
    have hn : value x ≤ 0 := by
      rw [hreal]
      simp only [FloatLib.Numerics.Dyadic.toReal, FloatLib.Numerics.Dyadic.signedSignificand,
        hsign.trans hs, ↓reduceIte, Int.cast_neg]
      exact mul_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr (Nat.cast_nonneg _)) hp.le
    have habs : Model.abs x = Model.neg x := by
      simp [Model.abs, Model.copySign, Model.neg, hs]
    change Model.toReal (Model.abs x) = |Model.toReal x|
    rw [habs, Model.toReal_neg x hx, abs_of_nonpos hn]

lemma min_valid (x y : Word) (hx : valid x) (hy : valid y) : valid (Model.minimum x y) := by
  have hn := Model.chooseNaN2_none_of_not_isNaN x y
    (Model.isNaN_eq_false_of_isFinite_eq_true x hx) (Model.isNaN_eq_false_of_isFinite_eq_true y hy)
  unfold Model.minimum
  simp only [Model.withNaNSelection_of_none _ _ hn]
  cases Model.compareNonNaN x y (Model.isNaN_eq_false_of_isFinite_eq_true x hx)
    (Model.isNaN_eq_false_of_isFinite_eq_true y hy) <;>
    simp [valid, apply_ite, hx, hy]

lemma max_valid (x y : Word) (hx : valid x) (hy : valid y) : valid (Model.maximum x y) := by
  have hn := Model.chooseNaN2_none_of_not_isNaN x y
    (Model.isNaN_eq_false_of_isFinite_eq_true x hx) (Model.isNaN_eq_false_of_isFinite_eq_true y hy)
  unfold Model.maximum
  simp only [Model.withNaNSelection_of_none _ _ hn]
  cases Model.compareNonNaN x y (Model.isNaN_eq_false_of_isFinite_eq_true x hx)
    (Model.isNaN_eq_false_of_isFinite_eq_true y hy) <;>
    simp [valid, apply_ite, hx, hy]

inductive Op where | add | sub | mul | div
  deriving DecidableEq

def Op.exact : Op → ℝ → ℝ → ℝ
  | .add => (·+·) | .sub => (·-·) | .mul => (·*·) | .div => (·/·)

def Op.bits : Op → Word → Word → Word
  | .add => Model.add | .sub => Model.sub | .mul => Model.mul | .div => Model.div

inductive Expr where
  | input (index : Nat)
  | one | two
  | abs (x : Expr)
  | min (x y : Expr)
  | max (x y : Expr)
  | op (o : Op) (x y : Expr)

def Expr.real (env : Nat → ℝ) : Expr → ℝ
  | .input i => env i
  | .one => 1 | .two => 2
  | .abs x => |x.real env|
  | .min x y => Min.min (x.real env) (y.real env)
  | .max x y => Max.max (x.real env) (y.real env)
  | .op o x y => Binary64.round (o.exact (x.real env) (y.real env))

def Expr.bits (env : Nat → Word) : Expr → Word
  | .input i => env i
  | .one => IEEE64.one | .two => IEEE64.two
  | .abs x => Model.abs (x.bits env)
  | .min x y => Model.minimum (x.bits env) (y.bits env)
  | .max x y => Model.maximum (x.bits env) (y.bits env)
  | .op o x y => o.bits (x.bits env) (y.bits env)

/-- The safety predicate contains exact magnitude and nonzero-divisor conditions,
not floating-point error assumptions. -/
def Expr.Safe (env : Nat → ℝ) : Expr → Prop
  | .input _ | .one | .two => True
  | .abs x => x.Safe env
  | .min x y | .max x y => x.Safe env ∧ y.Safe env
  | .op o x y => x.Safe env ∧ y.Safe env ∧
      |o.exact (x.real env) (y.real env)| ≤ 2 ∧ (o = .div → y.real env ≠ 0)

/-- Structural refinement from executable finite IEEE words to the rounded real semantics.
The theorem proves finiteness of every subtree while establishing its value. -/
theorem expression_refines (e : Expr) (env : Nat → Word)
    (hi : ∀ i, valid (env i)) (hs : e.Safe (fun i => value (env i))) :
    valid (e.bits env) ∧ value (e.bits env) = e.real (fun i => value (env i)) := by
  induction e with
  | input i => exact ⟨hi i, rfl⟩
  | one => exact ⟨one_valid, one_value⟩
  | two => exact ⟨two_valid, two_value⟩
  | abs x ih =>
    obtain ⟨hx, he⟩ := ih hs
    exact ⟨by simpa only [Expr.bits, valid, Model.isFinite_abs] using hx,
      by simpa only [Expr.bits, Expr.real, he] using abs_value (x.bits env) hx⟩
  | min x y ihx ihy =>
    obtain ⟨hx, ex⟩ := ihx hs.1
    obtain ⟨hy, ey⟩ := ihy hs.2
    exact ⟨min_valid _ _ hx hy, by simpa only [Expr.bits, Expr.real, ex, ey] using
      Model.toReal_minimum_eq_min_of_isFinite (x.bits env) (y.bits env) hx hy⟩
  | max x y ihx ihy =>
    obtain ⟨hx, ex⟩ := ihx hs.1
    obtain ⟨hy, ey⟩ := ihy hs.2
    exact ⟨max_valid _ _ hx hy, by simpa only [Expr.bits, Expr.real, ex, ey] using
      Model.toReal_maximum_eq_max_of_isFinite (x.bits env) (y.bits env) hx hy⟩
  | op o x y ihx ihy =>
    obtain ⟨hx, ex⟩ := ihx hs.1
    obtain ⟨hy, ey⟩ := ihy hs.2.1
    have hb : |o.exact (value (x.bits env)) (value (y.bits env))| ≤ 2 := by
      rw [ex, ey]; exact hs.2.2.1
    have hd : o = .div → value (y.bits env) ≠ 0 := by rw [ey]; exact hs.2.2.2
    cases o with
    | add => simpa only [Expr.bits, Expr.real, Op.bits, Op.exact, ex, ey] using add_refines _ _ hx hy hb
    | sub => simpa only [Expr.bits, Expr.real, Op.bits, Op.exact, ex, ey] using sub_refines _ _ hx hy hb
    | mul => simpa only [Expr.bits, Expr.real, Op.bits, Op.exact, ex, ey] using mul_refines _ _ hx hy hb
    | div => simpa only [Expr.bits, Expr.real, Op.bits, Op.exact, ex, ey] using div_refines _ _ hx hy (hd rfl) hb

end
end ModAB.IEEE64
