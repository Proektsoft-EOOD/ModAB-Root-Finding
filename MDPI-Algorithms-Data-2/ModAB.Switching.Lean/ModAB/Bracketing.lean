import ModAB.Schedule
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Instances.Real.Lemmas

/-! Exact-arithmetic interval invariants and the final nested-interval argument. -/
namespace ModAB.Convergence
open Set Filter Topology
noncomputable section

/-- The convex-combination identity used in `lem:interiority`. -/
theorem secant_convex {a b u v : ℝ} (hu : u<0) (hv : 0<v) :
    let alpha:=v/(v+|u|)
    0<alpha ∧ alpha<1 ∧ (a*v-b*u)/(v-u)=alpha*a+(1-alpha)*b := by
  dsimp
  rw [abs_of_neg hu]
  have hd : 0<v-u := by linarith
  refine ⟨div_pos hv hd,(div_lt_one hd).mpr (by linarith),?_⟩
  simp only [← sub_eq_add_neg]
  field_simp [ne_of_gt hd]
  <;> ring

theorem secant_interior {a b u v : ℝ} (hab : a<b) (huv : u*v<0) :
    (a*v-b*u)/(v-u) ∈ Ioo a b := by
  rcases mul_neg_iff.mp huv with h|h
  · obtain ⟨halpha0,halpha1,he⟩:=secant_convex (a:=b) (b:=a) h.2 h.1
    have hs : (a*v-b*u)/(v-u)=(b*u-a*v)/(u-v) := by
      have hd : v-u≠0 := by linarith [h.1,h.2]
      have hd' : u-v≠0 := by linarith [h.1,h.2]
      field_simp
      <;> ring
    rw [hs,he]
    constructor <;> nlinarith [mul_pos halpha0 (sub_pos.mpr hab),mul_pos (sub_pos.mpr halpha1) (sub_pos.mpr hab)]
  · obtain ⟨halpha0,halpha1,he⟩:=secant_convex (a:=a) (b:=b) h.1 h.2
    rw [he]
    constructor <;> nlinarith [mul_pos halpha0 (sub_pos.mpr hab),mul_pos (sub_pos.mpr halpha1) (sub_pos.mpr hab)]

/-- A modified ordinate has the sign of its nonzero true residual. -/
def SameSign (u f : ℝ) : Prop := 0<u*f

theorem positive_correction {u f θ : ℝ} (h : SameSign u f) (hθ : 0<θ) :
    SameSign (θ*u) f := by unfold SameSign at *; nlinarith [mul_pos hθ h]

theorem inserted_sign {f : ℝ} (hf : f≠0) : SameSign f f := by
  exact mul_self_pos.mpr hf

theorem working_opposite {u v fa fb : ℝ}
    (hu : SameSign u fa) (hv : SameSign v fb) (hf : fa*fb<0) : u*v<0 := by
  unfold SameSign at hu hv
  have hp : 0<(u*v)*(fa*fb) := by nlinarith [mul_pos hu hv]
  rcases mul_pos_iff.mp hp with h|h
  · linarith [h.2]
  · exact h.1

/-- Nonzero interior evaluation leaves exactly one sign-changing subinterval. -/
theorem sign_update {fa fb fx : ℝ} (hs : fa*fb<0) (hx : fx≠0) :
    (fa*fx<0 ∧ ¬fx*fb<0) ∨ (fx*fb<0 ∧ ¬fa*fx<0) := by
  rcases mul_neg_iff.mp hs with h|h
  · rcases lt_or_gt_of_ne hx with hx|hx
    · exact Or.inl ⟨mul_neg_of_pos_of_neg h.1 hx,not_lt_of_gt (mul_pos_of_neg_of_neg hx h.2)⟩
    · exact Or.inr ⟨mul_neg_of_pos_of_neg hx h.2,not_lt_of_gt (mul_pos h.1 hx)⟩
  · rcases lt_or_gt_of_ne hx with hx|hx
    · exact Or.inr ⟨mul_neg_of_neg_of_pos hx h.2,not_lt_of_gt (mul_pos_of_neg_of_neg h.1 hx)⟩
    · exact Or.inl ⟨mul_neg_of_neg_of_pos h.1 hx,not_lt_of_gt (mul_pos hx h.2)⟩

/-- Both possible updates are nested and do not increase the width. -/
theorem interval_update {a b x : ℝ} (hx : x∈Icc a b) :
    Icc a x ⊆ Icc a b ∧ Icc x b ⊆ Icc a b ∧
    0≤x-a ∧ x-a≤b-a ∧ 0≤b-x ∧ b-x≤b-a := by
  refine ⟨fun y hy=>⟨hy.1,hy.2.trans hx.2⟩,
    fun y hy=>⟨hx.1.trans hy.1,hy.2⟩,?_,?_,?_,?_⟩ <;> linarith [hx.1,hx.2]

theorem bisection_width (a b : ℝ) :
    (a+b)/2-a=(b-a)/2 ∧ b-(a+b)/2=(b-a)/2 := by constructor <;> ring

theorem bracket_contains_zero {f : ℝ→ℝ} {a b : ℝ} (hab : a≤b)
    (hf : ContinuousOn f (Icc a b)) (hs : f a*f b<0) :
    ∃ z∈Icc a b, f z=0 := by
  rcases mul_neg_iff.mp hs with h|h
  · exact intermediate_value_Icc' hab hf ⟨h.2.le,h.1.le⟩
  · exact intermediate_value_Icc hab hf ⟨h.1.le,h.2.le⟩

theorem residual_progress {r old retained η : ℝ} (hold : 0<old)
    (hret : old≤retained) (hη : η<1) (hnew : r<η*old) :
    min retained r=r ∧ r<η*old := by
  have hr : r<retained := by nlinarith
  exact ⟨min_eq_right hr.le,hnew⟩

/-- Repeated substitution; the strict statement requires at least one step. -/
theorem residual_progress_iterated {R : ℕ→ℝ} {η : ℝ} (hη : 0<η)
    (h : ∀ n, R (n+1)<η*R n) {k : ℕ} (hk : 0<k) : R k<η^k*R 0 := by
  cases k with
  | zero => omega
  | succ k =>
      induction k with
      | zero => simpa using h 0
      | succ k ih =>
          have hmul:=mul_lt_mul_of_pos_left (ih (by omega)) hη
          calc R (k+1+1)<η*R (k+1):=h (k+1)
               _<η*(η^(k+1)*R 0):=hmul
               _=η^(k+1+1)*R 0:=by rw [pow_succ]; ring

theorem interval_distance {a b x y : ℝ} (hx : x∈Icc a b) (hy : y∈Icc a b) :
    |x-y|≤b-a := by rw [abs_le]; constructor <;> linarith [hx.1,hx.2,hy.1,hy.2]

/-- Every choice inside shrinking intervals converges to their common point. -/
theorem chosen_points_tendsto {a b : ℕ→ℝ} {x : ℝ} {z : ℕ→ℝ}
    (hx : ∀n,x∈Icc (a n) (b n)) (hz : ∀n,z n∈Icc (a n) (b n))
    (hw : Tendsto (fun n=>b n-a n) atTop (𝓝 0)) : Tendsto z atTop (𝓝 x) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw ε hε
  refine ⟨N,fun n hn=>?_⟩
  rw [Real.dist_eq]
  have hpos : 0≤b n-a n := by linarith [(hx n).1,(hx n).2]
  have hsmall:=hN n hn
  rw [Real.dist_eq,sub_zero,abs_of_nonneg hpos] at hsmall
  exact (interval_distance (hz n) (hx n)).trans_lt hsmall

/-- Last paragraph of `thm:global_simple`, including the chosen-zero sequence. -/
theorem nested_intervals_root {f : ℝ→ℝ} {a b : ℕ→ℝ}
    (hnested : ∀n,Icc (a (n+1)) (b (n+1))⊆Icc (a n) (b n))
    (hzero : ∀n,∃z∈Icc (a n) (b n),f z=0)
    (hf : ContinuousOn f (Icc (a 0) (b 0)))
    (hw : Tendsto (fun n=>b n-a n) atTop (𝓝 0)) :
    ∃x, (⋂n,Icc (a n) (b n))={x} ∧ f x=0 := by
  have hne : ∀n,(Icc (a n) (b n)).Nonempty := fun n=>by
    obtain ⟨z,hz,_⟩:=hzero n; exact ⟨z,hz⟩
  obtain ⟨x,hx⟩:=IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed
    (fun n=>Icc (a n) (b n)) hnested hne isCompact_Icc (fun _=>isClosed_Icc)
  have hxall : ∀n,x∈Icc (a n) (b n) := mem_iInter.mp hx
  have hinter : (⋂n,Icc (a n) (b n))={x} := by
    apply Set.eq_singleton_iff_unique_mem.mpr
    refine ⟨hx,?_⟩
    intro y hy
    have hyall:=mem_iInter.mp hy
    have hlim:=chosen_points_tendsto hxall hyall hw
    exact tendsto_nhds_unique tendsto_const_nhds hlim
  choose z hz hval using hzero
  have hzlim:=chosen_points_tendsto hxall hz hw
  have hmono : Antitone (fun n=>Icc (a n) (b n)):=antitone_nat_of_succ_le hnested
  have hz0 : ∀n,z n∈Icc (a 0) (b 0):=fun n=>hmono (Nat.zero_le n) (hz n)
  have hcont : Tendsto (fun n=>f (z n)) atTop (𝓝 (f x)) :=
    (hf x (hxall 0)).tendsto.comp (tendsto_nhdsWithin_iff.mpr ⟨hzlim,Filter.Eventually.of_forall hz0⟩)
  have heq : f x=0 := by
    have he : (fun n=>f (z n))=(fun _ : ℕ=>0) := funext hval
    rw [he] at hcont
    exact (tendsto_nhds_unique tendsto_const_nhds hcont).symm
  exact ⟨x,hinter,heq⟩

end
end ModAB.Convergence
