import ModAB.Switching
import FloatLib.Floats.Formats.Flocq.Theory.Rounding.Order
import FloatLib.Floats.Formats.Flocq.Theory.Error.Multiplication
import FloatLib.Floats.Formats.Flocq.Theory.Analysis.Neighbors

/-!+Concrete nearest-even rounding for precision 53 and minimum subnormal exponent -1074.
The real-valued grid has no upper exponent bound. `Finite` adds the binary64
upper bound; all operation results used below are proved to lie within it.
-/

namespace ModAB.Binary64

open FloatLib.Floats.Formats.Flocq
open FloatLib.Numerics

noncomputable section

abbrev exponent : ℤ → ℤ := fltExp (-1074) 53
instance : ValidExp exponent := fltValidExp (-1074) 53 (by norm_num)
instance : MonotoneExp exponent := fltMonotoneExp (-1074) 53

abbrev Representable (x : ℝ) : Prop := genericFormat binaryRadix exponent x

def round (x : ℝ) : ℝ :=
  FloatLib.Floats.Formats.Flocq.round (β := binaryRadix) (fexp := exponent) nearestEven x

def maxFinite : ℝ := ((2:ℝ)^53-1)*(2:ℝ)^971
def Finite (x : ℝ) : Prop := Representable x ∧ |x| ≤ maxFinite

lemma pow_value (e : ℤ) : bpow binaryRadix e = (2:ℝ)^e := by
  simp [bpow, binaryRadix, Radix.toReal]

theorem round_representable (x : ℝ) : Representable (round x) :=
  generic_format_round nearestEven x

theorem round_exact {x : ℝ} (hx : Representable x) : round x = x :=
  round_preserves_generic nearestEven x hx

theorem round_monotone : Monotone round := by
  intro x y h
  exact round_mono nearestEven h

theorem power_representable (e : ℤ) (he : -1074 ≤ e) :
    Representable ((2:ℝ)^e) := by
  rw [← pow_value]
  apply generic_format_bpow
  dsimp [exponent, fltExp]
  exact max_le (by omega) he

@[simp] theorem round_zero : round 0 = 0 := round_exact generic_format_zero
@[simp] theorem round_one : round 1 = 1 := by
  have h := round_exact (power_representable 0 (by norm_num))
  simpa using h
@[simp] theorem round_two : round 2 = 2 := by
  have h := round_exact (power_representable 1 (by norm_num))
  simpa using h
@[simp] theorem round_half : round (1/2) = (1/2:ℝ) := by
  have h := round_exact (power_representable (-1) (by norm_num))
  norm_num at h ⊢
  exact h
@[simp] theorem round_quarter : round (1/4) = (1/4:ℝ) := by
  have h := round_exact (power_representable (-2) (by norm_num))
  norm_num at h ⊢
  exact h

theorem round_neg_exact {x : ℝ} (hx : Representable x) : round (-x) = -x :=
  round_exact (generic_format_neg x hx)

@[simp] theorem round_neg_one : round (-1) = -1 := by
  have h := round_neg_exact (power_representable 0 (by norm_num))
  simpa using h
@[simp] theorem round_neg_two : round (-2) = -2 := by
  have h := round_neg_exact (power_representable 1 (by norm_num))
  simpa using h
@[simp] theorem round_neg_half : round (-(1/2)) = -(1/2:ℝ) := by
  have h := round_neg_exact (power_representable (-1) (by norm_num))
  norm_num at h ⊢
  exact h

theorem round_interval {a b x : ℝ} (ha : round a = a) (hb : round b = b)
    (hx : x ∈ Set.Icc a b) : round x ∈ Set.Icc a b := by
  constructor
  · have h := round_monotone hx.1
    rwa [ha] at h
  · have h := round_monotone hx.2
    rwa [hb] at h

theorem round_abs_le_one {x : ℝ} (hx : |x| ≤ 1) : |round x| ≤ 1 :=
  abs_le.mpr (round_interval round_neg_one round_one (abs_le.mp hx))

theorem round_abs_le_two {x : ℝ} (hx : |x| ≤ 2) : |round x| ≤ 2 :=
  abs_le.mpr (round_interval round_neg_two round_two (abs_le.mp hx))

/-- Absolute one-step error, including gradual underflow and the endpoints ±2. -/
theorem round_error {x : ℝ} (hx : |x| ≤ 2) : Within (round x) x binary64Epsilon := by
  by_cases h0 : x = 0
  · subst x; simpa [Within] using le_of_lt binary64Epsilon_pos
  by_cases hp : x = 2
  · subst x; simpa [Within] using le_of_lt binary64Epsilon_pos
  by_cases hn : x = -2
  · subst x; simpa [Within] using le_of_lt binary64Epsilon_pos
  have hlt : |x| < 2 := by
    rcases le_total 0 x with h | h
    · rw [abs_of_nonneg h] at hx ⊢; exact lt_of_le_of_ne hx hp
    · rw [abs_of_nonpos h] at hx ⊢; exact lt_of_le_of_ne hx (by intro he; apply hn; linarith)
  have hmag : magnitude binaryRadix x ≤ 1 :=
    magnitude_le_of_abs_lt_bpow binaryRadix x 1 h0 (by simpa [pow_value] using hlt)
  have hexp : cexp binaryRadix exponent x ≤ -52 := by
    dsimp [cexp, exponent, fltExp]
    exact max_le (by omega) (by norm_num)
  have hpw := (bpow_le_bpow_iff binaryRadix _ _).mpr hexp
  have hu : ulp binaryRadix exponent x ≤ (2:ℝ)^(-52:ℤ) := by
    rw [ulp.of_ne_zero _ _ _ h0]
    simpa [pow_value] using hpw
  have he := error_bound_ulp (β := binaryRadix) (fexp := exponent) nearestEven x
  change |round x-x| ≤ _ at he
  apply le_trans he
  have h := div_le_div_of_nonneg_right hu (by norm_num : (0:ℝ) ≤ 2)
  convert h using 1 <;> norm_num [binary64Epsilon]

theorem round_finite {x : ℝ} (hx : |x| ≤ 2) : Finite (round x) := by
  refine ⟨round_representable x, le_trans (round_abs_le_two hx) ?_⟩
  have hp : (1:ℝ) ≤ 2^971 := one_le_pow₀ (by norm_num)
  dsimp [maxFinite]
  norm_num only [show (2:ℝ)^53 = 9007199254740992 by norm_num]
  calc
    (2:ℝ) ≤ 9007199254740991 := by norm_num
    _ ≤ 9007199254740991 * 2^971 := by simpa using mul_le_mul_of_nonneg_left hp (by norm_num : (0:ℝ) ≤ 9007199254740991)

lemma epsilon_power : binary64Epsilon = (2:ℝ)^(-53:ℤ) := by
  norm_num [binary64Epsilon]

@[simp] theorem round_epsilon : round binary64Epsilon = binary64Epsilon := by
  rw [epsilon_power]
  exact round_exact (power_representable (-53) (by norm_num))

/-- There is no representable number strictly between 1-2^-53 and 1. -/
theorem gap_below_one {x : ℝ} (hx : Representable x) (h0 : 0 ≤ x) (h1 : x < 1) :
    binary64Epsilon ≤ 1-x := by
  by_cases hh : x < 1/2
  · have he : binary64Epsilon < (1/2:ℝ) := by norm_num [binary64Epsilon]
    linarith
  have hlo : (1/2:ℝ) ≤ x := le_of_not_gt hh
  have hx0 : x ≠ 0 := by linarith
  have hm : magnitude binaryRadix x = 0 := by
    apply magnitude_eq_of_bpow_bounds binaryRadix x 0 hx0
    · simpa [pow_value, abs_of_nonneg h0] using hlo
    · simpa [pow_value, abs_of_nonneg h0] using h1
  obtain ⟨n, hn⟩ := scaled_mantissa_int_of_generic x hx
  have hrepr := scaled_mantissa_mul_bpow (β := binaryRadix) (fexp := exponent) x
  rw [hn] at hrepr
  have hc : cexp binaryRadix exponent x = -53 := by simp [cexp, exponent, fltExp, hm]
  rw [hc, pow_value] at hrepr
  have hnlt : (n:ℝ) < 9007199254740992 := by
    norm_num at hrepr
    nlinarith
  have hnZ : n < 9007199254740992 := by exact_mod_cast hnlt
  have hnZ' : n ≤ 9007199254740991 := by omega
  have hnR : (n:ℝ) ≤ 9007199254740991 := by exact_mod_cast hnZ'
  norm_num [binary64Epsilon] at hrepr ⊢
  nlinarith

theorem exact_half_of_ge_epsilon {x : ℝ} (hx : Representable x)
    (hlo : binary64Epsilon ≤ x) : round (x/2) = x/2 := by
  have xp : 0 < x := lt_of_lt_of_le binary64Epsilon_pos hlo
  have hmag := magnitude_mono_pos binaryRadix binary64Epsilon_pos hlo
  have hm : magnitude binaryRadix binary64Epsilon = -52 := by
    rw [epsilon_power, ← pow_value]
    simp
  rw [hm] at hmag
  have hf := generic_format_FLT_mul_bpow (β := binaryRadix)
    (-1074) 53 (by norm_num) hx (-1) (by omega)
  have hrepr : Representable (x/2) := by simpa [pow_value, div_eq_mul_inv] using hf
  exact round_exact hrepr

/-- The first halving in the C# evaluation is exact; the midpoint halving need not be. -/
theorem factor_halving_exact {rho : ℝ} (hr : rho ∈ Set.Icc (0:ℝ) 1) :
    round (round (1-round rho)/2) = round (1-round rho)/2 := by
  have hrh := round_interval round_zero round_one hr
  by_cases heq : round rho = 1
  · simp [heq]
  have hlt : round rho < 1 := lt_of_le_of_ne hrh.2 heq
  have hg := gap_below_one (round_representable rho) hrh.1 hlt
  have hga := round_monotone hg
  rw [round_epsilon] at hga
  exact exact_half_of_ge_epsilon (round_representable _) hga

end
end ModAB.Binary64
