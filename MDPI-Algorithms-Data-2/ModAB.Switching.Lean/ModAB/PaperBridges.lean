import ModAB.Decision
import ModAB.Cubic

/-! Explicit bridges from the executable scaled criterion to the manuscript's
original coefficient, for both endpoint signs and exactly three finite inputs. -/
namespace ModAB.Source
open IEEE64
noncomputable section

theorem domain_iff_product_negative (env : ℕ→ℝ) : Domain env ↔ env 0*env 1<0 := by
  constructor
  · rintro ⟨h0,h1,h⟩
    rcases h with ⟨ha,hb⟩|⟨hb,ha⟩
    · exact mul_neg_of_neg_of_pos (lt_of_le_of_ne ha h0) (lt_of_le_of_ne hb (Ne.symm h1))
    · exact mul_neg_of_pos_of_neg (lt_of_le_of_ne ha (Ne.symm h0)) (lt_of_le_of_ne hb h1)
  · intro h
    have h0 : env 0≠0 := by intro hz; simp [hz] at h
    have h1 : env 1≠0 := by intro hz; simp [hz] at h
    refine ⟨h0,h1,?_⟩
    rcases mul_neg_iff.mp h with h|h
    · exact Or.inr ⟨h.2.le,h.1.le⟩
    · exact Or.inl ⟨h.1.le,h.2.le⟩

theorem original_factor_identity {f1 f2 : ℝ} (hs : f1*f2<0) :
    (1-|(f1+f2)/2|/|f2-f1|)^2=kappa (min |f1| |f2|/max |f1| |f2|) := by
  rw [Cubic.endpoint_difference hs]
  rcases mul_neg_iff.mp hs with h|h
  · simpa only [abs_of_pos h.1,abs_of_neg h.2,neg_neg,add_comm,min_comm,max_comm] using
      endpoint_factor (neg_pos.mpr h.2) h.1
  · simpa only [abs_of_neg h.1,abs_of_pos h.2,neg_neg] using
      endpoint_factor (neg_pos.mpr h.1) h.2

def inputs (x y z : Word) : ℕ→Word
  | 0=>x
  | 1=>y
  | _=>z

theorem inputs_valid {x y z : Word} (hx : valid x) (hy : valid y) (hz : valid z) :
    ∀i,valid (inputs x y z i) := by
  intro i
  cases i with
  | zero => exact hx
  | succ i => cases i with
    | zero => exact hy
    | succ i => exact hz

/-- Positive certificate with the original k and no unused environment assumptions. -/
theorem three_input_switch {x y z : Word} (hx : valid x) (hy : valid y) (hz : valid z)
    (hs : value x*value y<0) (hc : decisionWord (gapExpr.bits (inputs x y z))=3) :
    originalL (value x) (value y) (value z) <
      originalR ((1-|(value x+value y)/2|/|value y-value x|)^2) (value x) (value y) (value z) := by
  have hd : Domain (fun i=>value (inputs x y z i)) := (domain_iff_product_negative _).mpr hs
  have h:=executable_switch_sound (inputs x y z) (inputs_valid hx hy hz) hd hc
  simpa only [rho,inputs,original_factor_identity hs] using h

theorem three_input_no_switch {x y z : Word} (hx : valid x) (hy : valid y) (hz : valid z)
    (hs : value x*value y<0) (hc : decisionWord (gapExpr.bits (inputs x y z))=2) :
    originalR ((1-|(value x+value y)/2|/|value y-value x|)^2) (value x) (value y) (value z) <
      originalL (value x) (value y) (value z) := by
  have hd : Domain (fun i=>value (inputs x y z i)) := (domain_iff_product_negative _).mpr hs
  have h:=executable_no_switch_sound (inputs x y z) (inputs_valid hx hy hz) hd hc
  simpa only [rho,inputs,original_factor_identity hs] using h

/-- The 57ε completeness statement now also ends at executable comparisons. -/
theorem executable_margin_guarantees (env : ℕ→Word) (hi : ∀i,valid (env i))
    (hd : Domain (fun i=>value (env i))) :
    (57*binary64Epsilon<mathematicalGap (fun i=>value (env i)) → decisionWord (gapExpr.bits env)=3) ∧
    (mathematicalGap (fun i=>value (env i))< -(57*binary64Epsilon) → decisionWord (gapExpr.bits env)=2) := by
  rw [decisionWord_refines _ (binary64_program_error env hi hd).1]
  exact binary64_margin_guarantees env hi hd

end
end ModAB.Source
