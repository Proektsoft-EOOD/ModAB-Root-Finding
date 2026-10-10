import ModAB.Bracketing
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-! The derivative-based residual-to-distance estimate `eq:residual_bound`.
The mean-value point is exhibited before the derivative lower bound is used. -/
namespace ModAB.Residual
open Set
noncomputable section

theorem mean_value_lower_bound {f : ℝ→ℝ} {a b m : ℝ}
    (hab : a<b) (hc : ContinuousOn f (Icc a b))
    (hd : ∀x∈Ioo a b,DifferentiableAt ℝ f x)
    (hl : ∀x∈Ioo a b,m≤|deriv f x|) :
    m*(b-a)≤|f b-f a| := by
  obtain ⟨c,hc',hs⟩:=exists_hasDerivAt_eq_slope f (deriv f) hab hc
    (fun x hx=>(hd x hx).hasDerivAt)
  have hm:=hl c hc'
  rw [hs,abs_div,abs_of_pos (sub_pos.mpr hab)] at hm
  exact (le_div_iff₀ (sub_pos.mpr hab)).mp hm

theorem residual_distance {f : ℝ→ℝ} {a b x z m : ℝ}
    (hx : x∈Icc a b) (hz : z∈Icc a b) (hzero : f z=0)
    (hm : 0<m) (hc : ContinuousOn f (Icc a b))
    (hd : ∀t∈Ioo a b,DifferentiableAt ℝ f t)
    (hl : ∀t∈Ioo a b,m≤|deriv f t|) :
    |x-z|≤|f x|/m := by
  apply (le_div_iff₀ hm).mpr
  rcases lt_trichotomy x z with h|h|h
  · have hi : Icc x z⊆Icc a b:=fun t ht=>⟨hx.1.trans ht.1,ht.2.trans hz.2⟩
    have hj : Ioo x z⊆Ioo a b:=fun t ht=>⟨lt_of_le_of_lt hx.1 ht.1,lt_of_lt_of_le ht.2 hz.2⟩
    have hv:=mean_value_lower_bound h (hc.mono hi)
      (fun t ht=>hd t (hj ht)) (fun t ht=>hl t (hj ht))
    rw [hzero,zero_sub,abs_neg] at hv
    rw [abs_of_neg (sub_neg.mpr h)]
    nlinarith only [hv]
  · subst z; simp
  · have hi : Icc z x⊆Icc a b:=fun t ht=>⟨hz.1.trans ht.1,ht.2.trans hx.2⟩
    have hj : Ioo z x⊆Ioo a b:=fun t ht=>⟨lt_of_le_of_lt hz.1 ht.1,lt_of_lt_of_le ht.2 hx.2⟩
    have hv:=mean_value_lower_bound h (hc.mono hi)
      (fun t ht=>hd t (hj ht)) (fun t ht=>hl t (hj ht))
    rw [hzero,sub_zero] at hv
    rw [abs_of_pos (sub_pos.mpr h)]
    nlinarith only [hv]

end
end ModAB.Residual
