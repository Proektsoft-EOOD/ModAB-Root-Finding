import ModAB.Phase
import Mathlib.Analysis.SpecialFunctions.Log.Base
import Mathlib.Analysis.SpecificLimits.Basic

/-! Exact-arithmetic width histories, their decomposition, and the manuscript's
amortized global estimate. The model permits earlier fallback and an arbitrary
switching test; the residual test only restricts its exception transitions.
-/
namespace ModAB.Convergence
open Filter Topology
noncomputable section

inductive Mode where
  | bisection
  | ab (S : ℝ) (t c : ℕ)

/-- One constructor per permitted width-policy transition. -/
inductive History (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Mode → Prop
  | initial : History B M W0 0 W0 .bisection
  | bisectStay {n w} : History B M W0 n w .bisection →
      History B M W0 (n+1) (w/2) .bisection
  | bisectStart {n w} : History B M W0 n w .bisection →
      History B M W0 (n+1) (w/2) (.ab (w/2) 0 0)
  | abWidth {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → w'≤threshold B S (t+1) →
      History B M W0 (n+1) w' (.ab S (t+1) 0)
  | abException {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      History B M W0 (n+1) w' (.ab S (t+1) (c+1))
  | abStop {n w w' S t c} : History B M W0 n w (.ab S t c) →
      0≤w' → w'≤w → History B M W0 (n+1) w' .bisection

/-- Completed parts: isolated bisections and AB blocks ending in bisection. -/
inductive CompletedPrefix (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Prop
  | initial : CompletedPrefix B M W0 0 W0
  | bisect {n w} : CompletedPrefix B M W0 n w → CompletedPrefix B M W0 (n+1) (w/2)
  | phase {n w t c u v} : CompletedPrefix B M W0 n w →
      ActivePhase B M w t c u → 0≤v → v≤u →
      CompletedPrefix B M W0 (n+t+2) (v/2)

theorem prefix_nonnegative {B W0 w : ℝ} {M n : ℕ} (hW0 : 0≤W0)
    (h : CompletedPrefix B M W0 n w) : 0≤w := by
  induction h with
  | initial => exact hW0
  | bisect _ ih => positivity
  | phase _ _ hv _ _ => positivity

inductive Decomposed (B : ℝ) (M : ℕ) (W0 : ℝ) : ℕ → ℝ → Mode → Prop
  | ready {n w} : CompletedPrefix B M W0 n w → Decomposed B M W0 n w .bisection
  | active {n S t c w} : CompletedPrefix B M W0 n S → ActivePhase B M S t c w →
      Decomposed B M W0 (n+t) w (.ab S t c)
  | stopped {n S t c u w} : CompletedPrefix B M W0 n S → ActivePhase B M S t c u →
      0≤w → w≤u → Decomposed B M W0 (n+t+1) w .bisection

/-- The partition in the manuscript is derived from transitions, not assumed. -/
theorem history_decomposes {B W0 w : ℝ} {M n : ℕ} {mode : Mode}
    (h : History B M W0 n w mode) : Decomposed B M W0 n w mode := by
  induction h with
  | initial => exact .ready .initial
  | bisectStay old ih =>
      cases ih with
      | ready hp => exact .ready (.bisect hp)
      | stopped hp ha hw hmono =>
          simpa only [Nat.add_assoc] using Decomposed.ready (CompletedPrefix.phase hp ha hw hmono)
  | bisectStart old ih =>
      cases ih with
      | ready hp => simpa using Decomposed.active (CompletedPrefix.bisect hp) ActivePhase.start
      | stopped hp ha hw hmono =>
          simpa only [Nat.add_zero,Nat.add_assoc] using
            Decomposed.active (CompletedPrefix.phase hp ha hw hmono) ActivePhase.start
  | abWidth old hw hmono ht ih =>
      cases ih with
      | active hp ha =>
          simpa only [Nat.add_assoc] using Decomposed.active hp (ActivePhase.width ha hw hmono ht)
  | abException old hw hmono ht hc ih =>
      cases ih with
      | active hp ha =>
          simpa only [Nat.add_assoc] using Decomposed.active hp (ActivePhase.exception ha hw hmono ht hc)
  | abStop old hw hmono ih =>
      cases ih with
      | active hp ha => exact .stopped hp ha hw hmono

theorem dyadic_compose (W0 : ℝ) (h k : ℕ) :
    (W0/(2:ℝ)^h)/(2:ℝ)^k=W0/(2:ℝ)^(h+k) := by rw [pow_add]; ring

/-- Sum of the credits of completed blocks; same induction as summing their bounds. -/
theorem prefix_credits {B W0 w : ℝ} {M n d : ℕ}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : CompletedPrefix B M W0 n w) :
    ∃ H : ℕ, w≤W0/(2:ℝ)^H ∧ n≤(d+2)*H := by
  induction h with
  | initial => exact ⟨0,by simp,by simp⟩
  | @bisect n w old ih =>
      obtain ⟨H,hw,hn⟩:=ih
      refine ⟨H+1,?_,?_⟩
      · simpa only [pow_one] using
          (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by simpa using dyadic_compose W0 H 1)
      · nlinarith
  | @phase n w t c u v old ha hv hmono ih =>
      obtain ⟨H,hw,hn⟩:=ih
      have hphase : ActivePhase B M w ((t+1)-1) c u := by simpa using ha
      obtain ⟨hb,hcost⟩:=completed_phase hB (prefix_nonnegative hW0 old) hC
        (by omega : 1≤t+1) hphase hmono
      let k:=phaseCredit d (t+1)
      refine ⟨H+k,?_,?_⟩
      · exact hb.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · change n+t+2≤(d+2)*(H+phaseCredit d (t+1))
        nlinarith

theorem history_credits {B W0 w : ℝ} {M n d : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : History B M W0 n w mode) :
    ∃ H : ℕ, w≤W0/(2:ℝ)^H ∧ n≤(d+2)*H+(d+1) := by
  have hd:=history_decomposes h
  cases hd with
  | ready hp =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      exact ⟨H,hw,by omega⟩
  | @active n S t c w hp ha =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      obtain ⟨k,hk,ht⟩:=tail_bound hB (prefix_nonnegative hW0 hp) hC (Tail.active ha)
      refine ⟨H+k,?_,?_⟩
      · exact hk.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · nlinarith
  | @stopped n S t c u w hp ha hw0 hmono =>
      obtain ⟨H,hw,hn⟩:=prefix_credits hB hW0 hC hp
      obtain ⟨k,hk,ht⟩:=tail_bound hB (prefix_nonnegative hW0 hp) hC (Tail.stopped ha hw0 hmono)
      refine ⟨H+k,?_,?_⟩
      · exact hk.trans ((div_le_div_of_nonneg_right hw (by positivity : 0≤(2:ℝ)^k)).trans_eq
          (dyadic_compose W0 H k))
      · nlinarith

/-- Integer step in `thm:quantitative`: n≤qH+q−1 implies floor(n/q)≤H. -/
theorem credit_floor {n q H : ℕ} (hq : 0<q) (hn : n≤q*H+(q-1)) : n/q≤H := by
  apply Nat.le_of_lt_succ
  apply (Nat.div_lt_iff_lt_mul hq).mpr
  have hq' : q-1+1=q := by omega
  nlinarith

theorem dyadic_antitone {W0 : ℝ} (hW0 : 0≤W0) {h k : ℕ} (hh : h≤k) :
    W0/(2:ℝ)^k≤W0/(2:ℝ)^h := by
  exact div_le_div_of_nonneg_left hW0 (by positivity)
    (pow_le_pow_right₀ (by norm_num) hh)

/-- `thm:quantitative`; natural-number division is the displayed floor. -/
theorem global_width_bound {B W0 w : ℝ} {M n d : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (h : History B M W0 n w mode) : w≤W0/(2:ℝ)^(n/(d+2)) := by
  obtain ⟨H,hw,hn⟩:=history_credits hB hW0 hC h
  have hc : n/(d+2)≤H := credit_floor (by omega) (by simpa using hn)
  exact hw.trans (dyadic_antitone hW0 hc)

def dyadicExponent (C : ℝ) : ℕ := Nat.ceil (Real.logb 2 C)

theorem ceil_log_bound {C : ℝ} (hC : 0<C) : C≤(2:ℝ)^(dyadicExponent C) := by
  have h := (Real.logb_le_iff_le_rpow (by norm_num : (1:ℝ)<2) hC).mp
    (Nat.le_ceil (Real.logb 2 C))
  simpa only [Real.rpow_natCast, dyadicExponent] using h

theorem global_width_bound_log {B W0 w : ℝ} {M n : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (h : History B M W0 n w mode) :
    w≤W0/(2:ℝ)^(n/(dyadicExponent (envelopeConstant B M)+2)) :=
  global_width_bound hB hW0 (ceil_log_bound (lt_of_lt_of_le (by norm_num) (C_at_least_one hB))) h

/-- `cor:C32`, without a numerical approximation to logarithms. -/
theorem C32_global_width {W0 w : ℝ} {n : ℕ} {mode : Mode} (hW0 : 0≤W0)
    (h : History 2 3 W0 n w mode) : w≤W0/(2:ℝ)^(n/7) := by
  exact global_width_bound (d:=5) (by norm_num) hW0 (by norm_num [envelopeConstant]) h

theorem finite_termination {B W0 w τ : ℝ} {M d K : ℕ} {mode : Mode}
    (hB : (1/2:ℝ)≤B) (hW0 : 0≤W0) (hC : envelopeConstant B M≤(2:ℝ)^d)
    (hK : W0/(2:ℝ)^K≤τ) (h : History B M W0 ((d+2)*K) w mode) : w≤τ := by
  have hb:=global_width_bound hB hW0 hC h
  have hdiv : ((d+2)*K)/(d+2)=K := Nat.mul_div_right K (by omega)
  rw [hdiv] at hb
  exact hb.trans hK

def toleranceExponent (W0 τ : ℝ) : ℕ := Nat.ceil (Real.logb 2 (W0/τ))

theorem tolerance_exponent_sufficient {W0 τ : ℝ} (hW0 : 0≤W0) (hτ : 0<τ) :
    W0/(2:ℝ)^(toleranceExponent W0 τ)≤τ := by
  by_cases hz : W0=0
  · simp [hz,hτ.le]
  · have hW : 0<W0 := lt_of_le_of_ne hW0 (Ne.symm hz)
    have hp := ceil_log_bound (div_pos hW hτ)
    change W0/τ≤(2:ℝ)^(toleranceExponent W0 τ) at hp
    apply (div_le_iff₀ (by positivity : 0<(2:ℝ)^(toleranceExponent W0 τ))).mpr
    exact (div_le_iff₀ hτ).mp hp |>.trans_eq (mul_comm _ _)

end
end ModAB.Convergence
