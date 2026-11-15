import ModAB.Global
import ModAB.Bracketing

/-! Assembly of the interval invariants, the width-policy execution, and the
qualitative and quantitative convergence theorems. All arithmetic here is real. -/
namespace ModAB.Convergence
open Set Filter Topology
noncomputable section

/-- The actual adaptive decision after a nonterminating AB update. -/
def adaptiveMode (B S η residual previousMinimum w : ℝ) (M t c : ℕ) : Mode :=
  if w≤threshold B S (t+1) then .ab S (t+1) 0
  else if residual<η*previousMinimum ∧ c<M then .ab S (t+1) (c+1)
  else .bisection

/-- The residual condition restricts the permitted exceptions; no estimate is assumed. -/
theorem adaptive_step {B S η residual previousMinimum w w' : ℝ} {M t c : ℕ}
    (hw' : 0≤w') (hmono : w'≤w) :
    Step B M w (.ab S t c) w' (adaptiveMode B S η residual previousMinimum w' M t c) := by
  unfold adaptiveMode
  split_ifs with ht hr
  · exact .abWidth hw' hmono ht
  · exact .abException hw' hmono (lt_of_not_ge ht) hr.2
  · exact .abStop hw' hmono

/-- Both sign orientations of the secant proposal are admitted. -/
theorem proposal_interior {a b u v fa fb : ℝ} (hab : a<b)
    (hu : SameSign u fa) (hv : SameSign v fb) (hf : fa*fb<0) :
    (a*v-b*u)/(v-u)∈Ioo a b := secant_interior hab (working_opposite hu hv hf)

/-- One nonterminating exact-arithmetic sign update. -/
def Updated (f : ℝ→ℝ) (a b x a' b' : ℝ) : Prop :=
  (a',b') = if f a*f x<0 then (a,x) else (x,b)

theorem updated_invariants {f : ℝ→ℝ} {a b x a' b' : ℝ}
    (hx : x∈Ioo a b) (hs : f a*f b<0) (hz : f x≠0)
    (hu : Updated f a b x a' b') :
    a'<b' ∧ f a'*f b'<0 ∧ Icc a' b'⊆Icc a b := by
  unfold Updated at hu
  have hi:=interval_update ⟨hx.1.le,hx.2.le⟩
  by_cases hax : f a*f x<0
  · rw [if_pos hax,Prod.mk.injEq] at hu
    rcases hu with ⟨rfl,rfl⟩
    exact ⟨hx.1,hax,hi.1⟩
  · rw [if_neg hax,Prod.mk.injEq] at hu
    have hxb : f x*f b<0 := by
      rcases sign_update hs hz with h|h
      · exact False.elim (hax h.1)
      · exact h.1
    rcases hu with ⟨rfl,rfl⟩
    exact ⟨hx.2,hxb,hi.2.1⟩

/-- An infinite nonterminating execution, with the sign update made explicit.
The proposal is any interior point and the modes obey the width policy. This
is the generality used in Remark `rem:proof_independence`; `proposal_interior`
verifies the AB proposal, and a midpoint is interior as well. -/
structure ExactRun (B : ℝ) (M : ℕ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→Mode
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  initial_mode : mode 0=.bisection
  interior : ∀n,point n∈Ioo (a n) (b n)
  nonzero : ∀n,f (point n)≠0
  update : ∀n,Updated f (a n) (b n) (point n) (a (n+1)) (b (n+1))
  width_step : ∀n,Step B M (b n-a n) (mode n) (b (n+1)-a (n+1)) (mode (n+1))

theorem exact_invariants {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) :
    ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 ∧
      Icc (r.a (n+1)) (r.b (n+1))⊆Icc (r.a n) (r.b n) := by
  have hs : ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 := by
    intro n
    induction n with
    | zero => exact ⟨r.initial_order,r.initial_sign⟩
    | succ n ih => exact (updated_invariants (r.interior n) ih.2 (r.nonzero n) (r.update n)).imp_right And.left
  intro n
  exact ⟨(hs n).1,(hs n).2,(updated_invariants (r.interior n) (hs n).2 (r.nonzero n) (r.update n)).2.2⟩

def ExactRun.widthRun {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) : Run B M where
  width n:=r.b n-r.a n
  mode:=r.mode
  nonneg_initial:=by linarith [r.initial_order]
  initial:=r.initial_mode
  next:=r.width_step

theorem exact_zeros {B : ℝ} {M : ℕ} {f : ℝ→ℝ} (r : ExactRun B M f) :
    ∀n,∃z∈Icc (r.a n) (r.b n),f z=0 := by
  have hi:=exact_invariants r
  have hmono : Antitone (fun n=>Icc (r.a n) (r.b n)):=antitone_nat_of_succ_le (fun n=>(hi n).2.2)
  intro n
  exact bracket_contains_zero (hi n).1.le (r.continuous.mono (hmono (Nat.zero_le n))) (hi n).2.1

/-- Full nonterminating alternative of `thm:global_simple`. Exact-zero termination
is the other alternative and needs no limit argument. -/
theorem global_convergence {B : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (r : ExactRun B M f) :
    Tendsto (fun n=>r.b n-r.a n) atTop (𝓝 0) ∧
    ∃x,(⋂n,Icc (r.a n) (r.b n))={x} ∧ f x=0 := by
  have hw:=width_tendsto_zero hB r.widthRun
  exact ⟨hw,nested_intervals_root (fun n=>(exact_invariants r n).2.2) (exact_zeros r) r.continuous hw⟩

theorem exact_global_width {B : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (r : ExactRun B M f) (n : ℕ) :
    r.b n-r.a n≤(r.b 0-r.a 0)/(2:ℝ)^(n/(dyadicExponent (envelopeConstant B M)+2)) :=
  global_width_bound_log hB r.widthRun.nonneg_initial (history_of_run r.widthRun n)

/-- `cor:finite_termination`: a continued run reaches the width tolerance by qK. -/
theorem exact_finite_tolerance {B τ : ℝ} {M : ℕ} {f : ℝ→ℝ}
    (hB : (1/2:ℝ)≤B) (hτ : 0<τ) (r : ExactRun B M f) :
    let q:=dyadicExponent (envelopeConstant B M)+2
    let K:=toleranceExponent (r.b 0-r.a 0) τ
    r.b (q*K)-r.a (q*K)≤τ := by
  dsimp
  exact finite_termination hB r.widthRun.nonneg_initial
    (ceil_log_bound (lt_of_lt_of_le (by norm_num) (C_at_least_one hB)))
    (tolerance_exponent_sufficient r.widthRun.nonneg_initial hτ)
    (history_of_run r.widthRun _)

/-- `eq:distance_zero_set`, expressed as an explicit nearby zero. -/
theorem distance_to_zero {f : ℝ→ℝ} {a b x τ : ℝ} (hx : x∈Icc a b)
    (hz : ∃z∈Icc a b,f z=0) (hw : b-a≤τ) : ∃z,f z=0 ∧ |x-z|≤τ := by
  obtain ⟨z,hz,hfz⟩:=hz
  exact ⟨z,hfz,(interval_distance hx hz).trans hw⟩

end
end ModAB.Convergence
