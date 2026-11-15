import ModAB.RationalUncapped
import ModAB.Rational

/-!
Quantitative local interfaces for the rational interpolation phase.

This file proves the safeguard consequences of explicit local error and
residual bounds. It does not assert that `C³` regularity supplies those bounds,
that an arbitrary global run enters this neighbourhood, or that fixed positive
stopping tolerances are asymptotically inactive. The clamp, zero-tolerance
closing operation and uncapped continuation decision are the existing project
definitions. The recurrence results are an exact-real local model.
-/
namespace ModAB.RationalLocalSafeguards
open Set Filter Topology
open ModAB.Convergence
noncomputable section

/-- Three local errors give a uniform contraction against either endpoint. -/
theorem triple_error_bound {u v d z B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε)
    (hu : |u| ≤ ε) (hv : |v| ≤ ε) (hd : |d| ≤ ε)
    (hz : |z| ≤ B * |u| * |v| * |d|) :
    |z| ≤ (B * ε ^ 2) * min |u| |v| := by
  have hvd : |v| * |d| ≤ ε ^ 2 := by
    calc |v| * |d| ≤ ε * ε := mul_le_mul hv hd (abs_nonneg _) hε
         _ = ε ^ 2 := by ring
  have hud : |u| * |d| ≤ ε ^ 2 := by
    calc |u| * |d| ≤ ε * ε := mul_le_mul hu hd (abs_nonneg _) hε
         _ = ε ^ 2 := by ring
  have hleft : |z| ≤ (B * ε ^ 2) * |u| := by
    calc |z| ≤ (B * |u|) * (|v| * |d|) := by nlinarith [hz]
         _ ≤ (B * |u|) * ε ^ 2 :=
           mul_le_mul_of_nonneg_left hvd (mul_nonneg hB (abs_nonneg _))
         _ = (B * ε ^ 2) * |u| := by ring
  have hright : |z| ≤ (B * ε ^ 2) * |v| := by
    calc |z| ≤ (B * |v|) * (|u| * |d|) := by nlinarith [hz]
         _ ≤ (B * |v|) * ε ^ 2 :=
           mul_le_mul_of_nonneg_left hud (mul_nonneg hB (abs_nonneg _))
         _ = (B * ε ^ 2) * |v| := by ring
  rw [mul_min_of_nonneg |u| |v| (mul_nonneg hB (sq_nonneg ε))]
  exact le_min hleft hright

/-- The usual quadratic secant bound has the matching minimum-error form. -/
theorem pair_error_bound {u v z S ε : ℝ} (hS : 0 ≤ S)
    (hu : |u| ≤ ε) (hv : |v| ≤ ε)
    (hz : |z| ≤ S * |u| * |v|) :
    |z| ≤ (S*ε) * min |u| |v| := by
  have hε : 0 ≤ ε := (abs_nonneg _).trans hu
  have hl := mul_le_mul_of_nonneg_left hv (mul_nonneg hS (abs_nonneg u))
  have hr := mul_le_mul_of_nonneg_left hu (mul_nonneg hS (abs_nonneg v))
  rw [mul_min_of_nonneg |u| |v| (mul_nonneg hS hε)]
  apply le_min <;> nlinarith

/-- The secant error bound implies the lower bound on its actual clamp radius. -/
theorem radius_lower_of_secant_error {a b r ya yb t : ℝ}
    (ha : a ≤ r) (hb : r ≤ b)
    (hs : |ModAB.Rational.secant a b ya yb - r| ≤
      t * min |a-r| |b-r|) :
    (1-t) * min |a-r| |b-r| ≤ ModAB.Rational.radius a b ya yb := by
  have hleft := min_le_left |a-r| |b-r|
  have hright := min_le_right |a-r| |b-r|
  have habsa : |a-r| = r-a := by rw [abs_of_nonpos (sub_nonpos.mpr ha)]; ring
  have habsb : |b-r| = b-r := abs_of_nonneg (sub_nonneg.mpr hb)
  rcases abs_le.mp hs with ⟨hsl, hsr⟩
  unfold ModAB.Rational.radius
  apply le_min
  · nlinarith
  · nlinarith

/-- This is an inactivity theorem for the existing `0.7` clamp, rather than
an independent replacement implementation of a clamp. -/
theorem limited_eq_raw_of_local_estimates {a b r ya yb raw q t : ℝ}
    (ha : a ≤ r) (hb : r ≤ b)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|)
    (hsec : |ModAB.Rational.secant a b ya yb-r| ≤
      t * min |a-r| |b-r|)
    (hsmall : q+t ≤ (7/10:ℝ) * (1-t)) :
    ModAB.Rational.limited a b ya yb raw = raw := by
  have hm : 0 ≤ min |a-r| |b-r| := le_min (abs_nonneg _) (abs_nonneg _)
  have hrad := radius_lower_of_secant_error ha hb hsec
  have hdev : |raw-ModAB.Rational.secant a b ya yb| ≤
      (q+t) * min |a-r| |b-r| := by
    have htri := abs_sub_le raw r (ModAB.Rational.secant a b ya yb)
    rw [abs_sub_comm r (ModAB.Rational.secant a b ya yb)] at htri
    nlinarith
  have hsmallm := mul_le_mul_of_nonneg_right hsmall hm
  have hbound : |raw-ModAB.Rational.secant a b ya yb| ≤
      (7/10:ℝ) * ModAB.Rational.radius a b ya yb := by nlinarith
  rcases abs_le.mp hbound with ⟨hl, hr⟩
  exact ModAB.Rational.limited_eq_self (by linarith) (by linarith)

/-- The local point bound gives the literal factor-42 residual exception. -/
theorem forty_two_residual_of_local_bounds {f : ℝ → ℝ}
    {a b r raw lower upper q : ℝ}
    (hupper : 0 ≤ upper)
    (hm : 0 < min |a-r| |b-r|)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|)
    (hfraw : |f raw| ≤ upper * |raw-r|)
    (hflower : lower * min |a-r| |b-r| ≤ min |f a| |f b|)
    (hsmall : 42 * upper * q < lower) :
    42 * |f raw| < min |f a| |f b| := by
  have hp := mul_le_mul_of_nonneg_left hraw hupper
  have hs := mul_lt_mul_of_pos_right hsmall hm
  nlinarith

/-- The residual exception prevents fallback regardless of the width clock. -/
theorem nextMode_interpolation_of_forty_two {B S m w R : ℝ} (t : ℕ)
    (hR : 42 * R < m) :
    ModAB.RationalUncapped.nextMode B (1/42) S m w R t =
      .interpolation S (t+1) := by
  have hr : R < (1/42:ℝ) * m :=
    (ModAB.RationalUncapped.forty_two_test R m).mpr hR
  unfold ModAB.RationalUncapped.nextMode
  split_ifs <;> simp_all

/-- In the zero-tolerance idealization, closing does not move an interior trial. -/
theorem closeCandidate_zero {a b x : ℝ} (hx : x ∈ Icc a b)
    (side : EndpointSide) : ModAB.Rational.closeCandidate a b 0 x side = x := by
  cases side <;> simp only [ModAB.Rational.closeCandidate] <;> split_ifs with h
  · exfalso; linarith [hx.1, h.1]
  · rfl
  · exfalso; linarith [hx.2, h.1]
  · rfl

/-- A sufficiently small error relative to both endpoints is strictly interior. -/
theorem interior_of_local_error {a b r raw q : ℝ}
    (ha : a < r) (hb : r < b) (hq : q < 1)
    (hraw : |raw-r| ≤ q * min |a-r| |b-r|) : raw ∈ Ioo a b := by
  have hm : 0 < min |a-r| |b-r| := by
    apply lt_min <;> apply abs_pos.mpr <;> linarith
  have hsmall : |raw-r| < min |a-r| |b-r| :=
    hraw.trans_lt (by simpa using mul_lt_mul_of_pos_right hq hm)
  rcases abs_lt.mp hsmall with ⟨hl, hr⟩
  have hleft := min_le_left |a-r| |b-r|
  have hright := min_le_right |a-r| |b-r|
  have habsa : |a-r| = r-a := by rw [abs_of_neg (sub_neg.mpr ha)]; ring
  have habsb : |b-r| = b-r := abs_of_pos (sub_pos.mpr hb)
  constructor <;> linarith

/-- Operational bridge with every analytic input displayed: the existing
guarded candidate is the raw rational point, zero-tolerance closing leaves it
unchanged, and the existing factor-42 policy continues interpolation. The
curvature decision and nonzero rational denominator are explicit premises. -/
theorem local_operational_bridge {f : ℝ → ℝ}
    {a b r xp q t lower upper B S w : ℝ} (clock : ℕ) (side : EndpointSide)
    (ha : a < r) (hb : r < b) (hq : q < 1)
    (hden : ModAB.Rational.denominator a b (f a) (f b) xp (f xp) ≠ 0)
    (hraw : |ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)-r| ≤
      q * min |a-r| |b-r|)
    (hsec : |ModAB.Rational.secant a b (f a) (f b)-r| ≤
      t * min |a-r| |b-r|)
    (hclamp : q+t ≤ (7/10:ℝ)*(1-t)) (hupper : 0 ≤ upper)
    (hfraw : |f (ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp))| ≤
      upper * |ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)-r|)
    (hflower : lower * min |a-r| |b-r| ≤ min |f a| |f b|)
    (hres : 42*upper*q < lower) :
    let raw := ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp)
    let trial := ModAB.Rational.candidate a b (f a) (f b) xp (f xp) true
    trial = raw ∧ ModAB.Rational.closeCandidate a b 0 trial side = raw ∧
      ModAB.RationalUncapped.nextMode B (1/42) S
        (min |f a| |f b|) w |f trial| clock = .interpolation S (clock+1) := by
  dsimp only
  have hc : ModAB.Rational.candidate a b (f a) (f b) xp (f xp) true =
      ModAB.Rational.rawPoint a b (f a) (f b) xp (f xp) := by
    rw [ModAB.Rational.candidate, ite_eq_left ⟨rfl, hden⟩]
    exact limited_eq_raw_of_local_estimates ha.le hb.le hraw hsec hclamp
  rw [hc]
  have hi := interior_of_local_error ha hb hq hraw
  have hm : 0 < min |a-r| |b-r| := by
    apply lt_min <;> apply abs_pos.mpr <;> linarith
  exact ⟨rfl, closeCandidate_zero ⟨hi.1.le, hi.2.le⟩ side,
    nextMode_interpolation_of_forty_two clock
      (forty_two_residual_of_local_bounds hupper hm hraw hfraw hflower hres)⟩

/-- Positive local coefficient means the new point crosses from the most
recent side to the opposite side, even before any asymptotic limit is used. -/
theorem positive_coefficient_crossing {en es ed ez H : ℝ}
    (hH : 0 < H) (hnd : 0 < en*ed) (hns : en*es < 0)
    (he : ez = H*en*es*ed) : 0 < ez*es ∧ ez*en < 0 := by
  have hs : es ≠ 0 := by intro hs; simp [hs] at hns
  have hh : 0 < H*(en*ed) := mul_pos hH hnd
  rw [he]
  constructor
  · have h := mul_pos hh (sq_pos_of_ne_zero hs)
    convert h using 1
    ring
  · have h := mul_neg_of_pos_of_neg hh hns
    convert h using 1
    ring

/-- Crossing from a fresh left endpoint restores two true ordinates. -/
theorem crossing_left_to_right_restores_true (s : OrdinateState) (y : ℝ)
    (hp : s.previous = some .left) (hf : ModAB.Rational.Fresh s) :
    (recordUpdate s .right y).leftWork = (recordUpdate s .right y).leftTrue ∧
    (recordUpdate s .right y).rightWork = (recordUpdate s .right y).rightTrue := by
  have h : s.leftWork = s.leftTrue := by simpa [ModAB.Rational.Fresh, hp] using hf
  simp [recordUpdate, hp, h]

/-- The symmetric crossing restores two true ordinates as well. -/
theorem crossing_right_to_left_restores_true (s : OrdinateState) (y : ℝ)
    (hp : s.previous = some .right) (hf : ModAB.Rational.Fresh s) :
    (recordUpdate s .left y).leftWork = (recordUpdate s .left y).leftTrue ∧
    (recordUpdate s .left y).rightWork = (recordUpdate s .left y).rightTrue := by
  have h : s.rightWork = s.rightTrue := by simpa [ModAB.Rational.Fresh, hp] using hf
  simp [recordUpdate, hp, h]

/-- A bounded triple-product recurrence preserves a sufficiently small initial
neighbourhood. This is an explicit recurrence interface, not a source-code
refinement or an assumed differentiability theorem. -/
theorem recurrence_local_invariant {e : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 ≤ 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hstep : ∀ n, |e (n+3)| ≤ B * |e (n+2)| * |e (n+1)| * |e n|) :
    ∀ n, |e n| ≤ ε := by
  have hq0 : 0 ≤ B*ε^2 := mul_nonneg hB (sq_nonneg ε)
  have htrip : ∀ n, |e n| ≤ ε ∧ |e (n+1)| ≤ ε ∧ |e (n+2)| ≤ ε := by
    intro n
    induction n with
    | zero => exact ⟨h0, h1, h2⟩
    | succ n ih =>
      have hc := triple_error_bound hB hε ih.2.2 ih.2.1 ih.1 (hstep n)
      have hm : min |e (n+2)| |e (n+1)| ≤ ε := (min_le_left _ _).trans ih.2.2
      have hb := mul_le_mul_of_nonneg_left hm hq0
      have hqe : (B*ε^2)*ε ≤ ε := by
        simpa using mul_le_mul_of_nonneg_right hq hε
      have he : |e (n+3)| ≤ ε := hc.trans (hb.trans hqe)
      exact ⟨ih.2.1, ih.2.2, by simpa [Nat.add_assoc] using he⟩
  exact fun n => (htrip n).1

/-- The same local recurrence entails convergence, not merely boundedness. -/
theorem recurrence_tendsto_zero {e : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hstep : ∀ n, |e (n+3)| ≤ B * |e (n+2)| * |e (n+1)| * |e n|) :
    Tendsto e atTop (𝓝 0) := by
  have hi := recurrence_local_invariant hB hε hq.le h0 h1 h2 hstep
  have hq0 : 0 ≤ B*ε^2 := mul_nonneg hB (sq_nonneg ε)
  have ht : Tendsto (fun n => |e (n+2)|) atTop (𝓝 0) := by
    apply ModAB.RationalUncapped.geometric_residual hq0 hq
    · exact fun _ => abs_nonneg _
    · intro n
      have hc := triple_error_bound hB hε (hi (n+2)) (hi (n+1)) (hi n) (hstep n)
      have hm := mul_le_mul_of_nonneg_left (min_le_left |e (n+2)| |e (n+1)|) hq0
      simpa [Nat.add_assoc] using hc.trans hm
  have habs : Tendsto (fun n => |e n|) atTop (𝓝 0) :=
    (tendsto_add_atTop_iff_nat 2).mp ht
  apply Metric.tendsto_atTop.mpr
  intro δ hδ
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp habs δ hδ
  exact ⟨N, fun n hn => by simpa [Real.dist_eq] using hN n hn⟩

/-- A bounded error coefficient is enough for the local convergence premise
needed by an order theorem; no convergence assumption on the errors is used. -/
theorem bounded_coefficient_recurrence_tendsto {e H : ℕ → ℝ} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B*ε^2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hH : ∀ n, |H n| ≤ B)
    (he : ∀ n, e (n+3) = H n * e (n+2) * e (n+1) * e n) :
    Tendsto e atTop (𝓝 0) := by
  apply recurrence_tendsto_zero hB hε hq h0 h1 h2
  intro n
  rw [he n, abs_mul, abs_mul, abs_mul]
  gcongr
  exact hH n

end
end ModAB.RationalLocalSafeguards
