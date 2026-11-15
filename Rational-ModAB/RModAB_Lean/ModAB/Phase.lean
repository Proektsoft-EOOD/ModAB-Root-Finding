import ModAB.Switching

/-! The width-policy projection of a continuing AB phase. The exception constructor
records the cap; requiring residual progress can only remove permitted executions.
Neither the geometric envelope nor a convergence conclusion is a field of the model. -/
namespace ModAB.Convergence
noncomputable section

def threshold (B S : ℝ) (t : ℕ) : ℝ := 2*B*S/(2:ℝ)^t
def envelopeConstant (B : ℝ) (M : ℕ) : ℝ := (2:ℝ)^(M+1)*B

inductive ActivePhase (B : ℝ) (M : ℕ) (S : ℝ) : ℕ → ℕ → ℝ → Prop
  | start : ActivePhase B M S 0 0 S
  | width {t c : ℕ} {w w' : ℝ} :
      ActivePhase B M S t c w → 0≤w' → w'≤w → w'≤threshold B S (t+1) →
      ActivePhase B M S (t+1) 0 w'
  | exception {t c : ℕ} {w w' : ℝ} :
      ActivePhase B M S t c w → 0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      ActivePhase B M S (t+1) (c+1) w'

theorem active_nonnegative {B S w : ℝ} {M t c : ℕ} (hS : 0≤S)
    (h : ActivePhase B M S t c w) : 0≤w := by
  induction h with
  | start => exact hS
  | width _ hw _ _ _ => exact hw
  | exception _ hw _ _ _ _ => exact hw

theorem active_nonincrease {B S w : ℝ} {M t c : ℕ}
    (h : ActivePhase B M S t c w) : w≤S := by
  induction h with
  | start => exact le_refl _
  | width _ _ hw _ ih => exact hw.trans ih
  | exception _ _ hw _ _ ih => exact hw.trans ih

/-- The manuscript's last successful width test j=t-c, including j=0. -/
theorem last_width_success {B S w : ℝ} {M t c : ℕ} (hB : (1/2:ℝ)≤B) (hS : 0≤S)
    (h : ActivePhase B M S t c w) :
    ∃ j : ℕ, j+c=t ∧ c≤M ∧ w≤threshold B S j := by
  induction h with
  | start =>
      refine ⟨0,rfl,Nat.zero_le _,?_⟩
      simp only [threshold,pow_zero,div_one]
      nlinarith [mul_nonneg hS (show 0≤2*B-1 by linarith)]
  | @width t c w w' old hw hmono htest ih =>
      exact ⟨t+1,by omega,Nat.zero_le _,htest⟩
  | @exception t c w w' old hw hmono htest hcap ih =>
      obtain ⟨j,hj,hc,hbound⟩:=ih
      exact ⟨j,by omega,by omega,hmono.trans hbound⟩

theorem threshold_shift (B S : ℝ) (j c : ℕ) :
    threshold B S j=(2:ℝ)^c*threshold B S (j+c) := by
  unfold threshold
  rw [pow_add]
  field_simp

/-- `lem:adaptive_envelope`, following its last-successful-test argument. -/
theorem active_envelope {B S w : ℝ} {M t c : ℕ} (hB : (1/2:ℝ)≤B) (hS : 0≤S)
    (h : ActivePhase B M S t c w) : w≤envelopeConstant B M*S/(2:ℝ)^t := by
  obtain ⟨j,hj,hc,hw⟩:=last_width_success hB hS h
  have hp : (2:ℝ)^c≤(2:ℝ)^M := pow_le_pow_right₀ (by norm_num) hc
  have ht : 0≤threshold B S t := by unfold threshold; positivity
  calc
    w≤threshold B S j := hw
    _=(2:ℝ)^c*threshold B S t := by rw [threshold_shift B S j c,hj]
    _≤(2:ℝ)^M*threshold B S t := mul_le_mul_of_nonneg_right hp ht
    _=envelopeConstant B M*S/(2:ℝ)^t := by unfold threshold envelopeConstant; rw [pow_succ]; ring

theorem C_at_least_one {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) : 1≤envelopeConstant B M := by
  have hp : (1:ℝ)≤2^M := one_le_pow₀ (by norm_num)
  unfold envelopeConstant
  rw [pow_succ]
  nlinarith

theorem C32 : envelopeConstant 2 3=32 := by norm_num [envelopeConstant]

theorem envelope_to_dyadic {B S w : ℝ} {M t c d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : ActivePhase B M S t c w) (ht : d≤t) : w≤S/(2:ℝ)^(t-d) := by
  have he := active_envelope hB hS h
  have hm := div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_right hC hS)
    (by positivity : 0≤(2:ℝ)^t)
  have ht' : t=d+(t-d) := by omega
  calc
    w≤envelopeConstant B M*S/(2:ℝ)^t := he
    _≤(2:ℝ)^d*S/(2:ℝ)^t := hm
    _=S/(2:ℝ)^(t-d) := by conv_lhs => rw [ht',pow_add]; field_simp

def phaseCredit (d t : ℕ) : ℕ := if t≤d+1 then 1 else t-d

/-- Completed phase: t−1 continuing steps, one last AB update, then bisection. -/
theorem completed_phase {B S u v : ℝ} {M t c d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (ht : 1≤t) (hphase : ActivePhase B M S (t-1) c u) (huv : v≤u) :
    v/2≤S/(2:ℝ)^(phaseCredit d t) ∧ t+1≤(d+2)*phaseCredit d t := by
  by_cases hshort : t≤d+1
  · simp only [phaseCredit,if_pos hshort,pow_one]
    exact ⟨div_le_div_of_nonneg_right (huv.trans (active_nonincrease hphase)) (by norm_num),by omega⟩
  · have hlong : d≤t-1 := by omega
    have hw := envelope_to_dyadic hB hS hC hphase hlong
    have he : t-d=(t-1-d)+1 := by omega
    have hratio : (S/(2:ℝ)^(t-1-d))/2=S/(2:ℝ)^(t-d) := by rw [he,pow_succ]; ring
    simp only [phaseCredit,if_neg hshort]
    constructor
    · exact (div_le_div_of_nonneg_right (huv.trans hw) (by norm_num)).trans_eq hratio
    · have hh : 2≤t-d := by omega
      have hm := Nat.mul_le_mul_left (d+1) (show 1≤t-d by omega)
      have htde : t=d+(t-d) := by omega
      nlinarith

inductive Tail (B : ℝ) (M : ℕ) (S : ℝ) : ℕ → ℝ → Prop
  | active {t c w} : ActivePhase B M S t c w → Tail B M S t w
  | stopped {t c u w} : ActivePhase B M S t c u → 0≤w → w≤u → Tail B M S (t+1) w

theorem tail_nonnegative {B S w : ℝ} {M t : ℕ} (hS : 0≤S)
    (h : Tail B M S t w) : 0≤w := by
  cases h with
  | active h => exact active_nonnegative hS h
  | stopped _ hw _ => exact hw

/-- `lem:tail`: the same four cases (active/stopped, short/long). -/
theorem tail_bound {B S w : ℝ} {M t d : ℕ}
    (hB : (1/2:ℝ)≤B) (hS : 0≤S) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : Tail B M S t w) :
    ∃ k : ℕ, w≤S/(2:ℝ)^k ∧ t≤(d+2)*k+(d+1) := by
  cases h with
  | @active t c w hp =>
      by_cases ht : t≤d
      · exact ⟨0,by simpa using active_nonincrease hp,by omega⟩
      · refine ⟨t-d,envelope_to_dyadic hB hS hC hp (by omega),?_⟩
        have hm := Nat.mul_le_mul_left (d+1) (show 1≤t-d by omega)
        have htde : t=d+(t-d) := by omega
        nlinarith
  | @stopped t c u w hp hw hwu =>
      by_cases ht : t≤d
      · exact ⟨0,by simpa using hwu.trans (active_nonincrease hp),by omega⟩
      · refine ⟨t-d,hwu.trans (envelope_to_dyadic hB hS hC hp (by omega)),?_⟩
        have hm := Nat.zero_le ((d+1)*(t-d))
        have htde : t=d+(t-d) := by omega
        nlinarith

end
end ModAB.Convergence
