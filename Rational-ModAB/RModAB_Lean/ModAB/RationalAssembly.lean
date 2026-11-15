import ModAB.RationalUncapped
import ModAB.Rational

/-!
Operational assembly of the exact-real uncapped method.  Candidate membership
and the scheduling transition are proved from the proposal, closing operation,
true-sign update, and mode decision.  They are not assumed as independent
invariants.  The source proposal relation includes the right-endpoint fallback.
This is an exact-real model of those branches, not a full IEEE refinement.
-/
namespace ModAB.RationalAssembly
open Set Filter Topology
open ModAB.Convergence ModAB.Rational
noncomputable section

def RecentEndpoint (a b latest : ℝ) : EndpointSide→Prop
  | .left => latest=a
  | .right => latest=b

def stepPoint (a b tolerance pre : ℝ) (recent : EndpointSide) :
    ModAB.RationalUncapped.Mode→ℝ
  | .bisection => (a+b)/2
  | .interpolation _ _ => closeCandidate a b tolerance pre recent

/-- The source proposal relation, including its endpoint fallback, supplies
closed membership after the guarded closing operation. -/
theorem source_closed_update {f : ℝ→ℝ} {a b ya yb pre tolerance x a' b' : ℝ}
    (recent : EndpointSide) (hab : a<b) (hs : f a*f b<0)
    (hya : SameSign ya (f a)) (hyb : SameSign yb (f b))
    (hp : SourceProposal a b ya yb pre) (ht : 0<tolerance)
    (hx : x=closeCandidate a b tolerance pre recent)
    (hz : f x≠0) (hu : Updated f a b x a' b') :
    a'<b' ∧ f a'*f b'<0 ∧ Icc a' b'⊆Icc a b := by
  have hpre:=sourceProposal_mem_closed hab (working_opposite hya hyb hs) hp
  have hc : x∈Icc a b := by
    rw [hx]
    exact closeCandidate_mem_closed recent ht hpre
  exact ModAB.RationalBridge.updated_invariants_closed hab hc hs hz hu

/-- Operational data.  The phase clock advances through `nextMode`; neither
bracket membership nor width contraction is postulated.  The working-sign
premises are supplied by positive AB corrections in the ordinate-state layer. -/
structure SourceRun (B η : ℝ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→ModAB.RationalUncapped.Mode
  ya : ℕ→ℝ
  yb : ℕ→ℝ
  preproposal : ℕ→ℝ
  tolerance : ℕ→ℝ
  recent : ℕ→EndpointSide
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  initial_point : point 0=a 0 ∨ point 0=b 0
  initial_mode : mode 0=.bisection
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  working_left : ∀n,SameSign (ya n) (f (a n))
  working_right : ∀n,SameSign (yb n) (f (b n))
  proposal : ∀n,SourceProposal (a n) (b n) (ya n) (yb n) (preproposal n)
  tolerance_pos : ∀n,0<tolerance n
  recent_correct : ∀n S t,mode n=.interpolation S t→
    RecentEndpoint (a n) (b n) (point n) (recent n)
  point_rule : ∀n,point (n+1)=stepPoint (a n) (b n) (tolerance n)
    (preproposal n) (recent n) (mode n)
  nonzero : ∀n,f (point (n+1))≠0
  update : ∀n,Updated f (a n) (b n) (point (n+1)) (a (n+1)) (b (n+1))
  bisect_mode : ∀n,mode n=.bisection→
    mode (n+1)=.bisection ∨ mode (n+1)=.interpolation ((b n-a n)/2) 0
  continue_mode : ∀n S t,mode n=.interpolation S t→
    mode (n+1)=ModAB.RationalUncapped.nextMode B η S
      (min |f (a n)| |f (b n)|) (b (n+1)-a (n+1)) |f (point (n+1))| t

theorem proposal_mem_of_invariants {B η : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (n : ℕ) (hab : s.a n<s.b n)
    (hs : f (s.a n)*f (s.b n)<0) : s.point (n+1)∈Icc (s.a n) (s.b n) := by
  rw [s.point_rule n]
  cases hm : s.mode n with
  | bisection =>
    simp only [stepPoint]
    constructor <;> linarith
  | interpolation S t =>
    simp only [stepPoint]
    exact closeCandidate_mem_closed (s.recent n) (s.tolerance_pos n)
      (sourceProposal_mem_closed hab
        (working_opposite (s.working_left n) (s.working_right n) hs) (s.proposal n))

theorem source_invariants {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) :
    ∀n,s.a n<s.b n ∧ f (s.a n)*f (s.b n)<0 := by
  intro n
  induction n with
  | zero => exact ⟨s.initial_order,s.initial_sign⟩
  | succ n ih =>
    have h:=ModAB.RationalBridge.updated_invariants_closed ih.1
      (proposal_mem_of_invariants s n ih.1 ih.2) ih.2 (s.nonzero n) (s.update n)
    exact ⟨h.1,h.2.1⟩

theorem source_candidate_mem {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    s.point (n+1)∈Icc (s.a n) (s.b n) :=
  proposal_mem_of_invariants s n (source_invariants s n).1 (source_invariants s n).2

theorem source_width_step {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    0≤s.b (n+1)-s.a (n+1) ∧ s.b (n+1)-s.a (n+1)≤s.b n-s.a n := by
  have hx:=source_candidate_mem s n
  have hu:=s.update n
  unfold Updated at hu
  split_ifs at hu with hs
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2]
    constructor <;> linarith [hx.1,hx.2]
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2]
    constructor <;> linarith [hx.1,hx.2]

theorem source_bisection_width {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f)
    (n : ℕ) (hm : s.mode n=.bisection) :
    s.b (n+1)-s.a (n+1)=(s.b n-s.a n)/2 := by
  have hp : s.point (n+1)=(s.a n+s.b n)/2 := by
    simpa [stepPoint,hm] using s.point_rule n
  have hu:=s.update n
  unfold Updated at hu
  split_ifs at hu with hs
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2,hp]
    ring
  · rw [Prod.mk.injEq] at hu
    rw [hu.1,hu.2,hp]
    ring

theorem source_policy {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) (n : ℕ) :
    ModAB.RationalUncapped.Step B η (min |f (s.a n)| |f (s.b n)|)
      (s.b n-s.a n) |f (s.point n)| (s.mode n)
      (s.b (n+1)-s.a (n+1)) |f (s.point (n+1))| (s.mode (n+1)) := by
  cases hm : s.mode n with
  | bisection =>
    have hw:=source_bisection_width s n hm
    rcases s.bisect_mode n hm with hn|hn
    · rw [hw,hn]
      exact .bisectStay _ _ _
    · rw [hw,hn]
      exact .bisectStart _ _ _
  | interpolation S t =>
    rw [s.continue_mode n S t hm]
    exact ModAB.RationalUncapped.decision_step
      (source_width_step s n).1 (source_width_step s n).2

/-- The abstract uncapped run is built from the operational proposal and tests.
Its candidate-membership and policy fields are derived theorems above. -/
def SourceRun.exactRun {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f) :
    ModAB.RationalUncapped.ExactRun B η f where
  a:=s.a
  b:=s.b
  point:=s.point
  mode:=s.mode
  initial_order:=s.initial_order
  initial_sign:=s.initial_sign
  initial_point:=s.initial_point
  continuous:=s.continuous
  candidate_mem:=source_candidate_mem s
  nonzero:=s.nonzero
  update:=s.update
  policy:=source_policy s

theorem source_global_residual {B η : ℝ} {f : ℝ→ℝ} (s : SourceRun B η f)
    (hη0 : 0≤η) (hη1 : η<1) :
    Tendsto (fun n=>|f (s.point n)|) atTop (𝓝 0) :=
  ModAB.RationalUncapped.global_residual_convergence hη0 hη1 s.exactRun

theorem source_unique_root_convergence {B η z : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (hη0 : 0≤η) (hη1 : η<1)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z) :
    Tendsto s.point atTop (𝓝 z) :=
  ModAB.RationalUncapped.unique_root_convergence hη0 hη1 s.exactRun hunique

/-- Actual closing-operation equations give the accepted-point spacing. -/
theorem source_closing_separation {B η : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (n : ℕ) {S : ℝ} {t : ℕ}
    (hm : s.mode n=.interpolation S t) (hw : s.tolerance n<s.b n-s.a n) :
    s.tolerance n/2≤|s.point (n+1)-s.point n| := by
  have hp : s.point (n+1)=closeCandidate (s.a n) (s.b n) (s.tolerance n)
      (s.preproposal n) (s.recent n) := by
    simpa [stepPoint,hm] using s.point_rule n
  have hrecent:=s.recent_correct n S t hm
  cases he : s.recent n with
  | left =>
    have hold : s.point n=s.a n := by simpa [RecentEndpoint,he] using hrecent
    have hsep:=closeCandidate_left_separation (x:=s.preproposal n) (s.tolerance_pos n) hw
    rw [hp,he,hold]
    exact (show s.tolerance n/2≤
      closeCandidate (s.a n) (s.b n) (s.tolerance n) (s.preproposal n) .left-s.a n by linarith).trans
      (le_abs_self _)
  | right =>
    have hold : s.point n=s.b n := by simpa [RecentEndpoint,he] using hrecent
    have hsep:=closeCandidate_right_separation (x:=s.preproposal n) (s.tolerance_pos n) hw
    rw [hp,he,hold,abs_sub_comm]
    exact (show s.tolerance n/2≤
      s.b n-closeCandidate (s.a n) (s.b n) (s.tolerance n) (s.preproposal n) .right by linarith).trans
      (le_abs_self _)

/-- Dynamic-tolerance width termination for the actual closing operation in
the exact-real model, under uniqueness and a positive absolute tolerance floor.
There is no assumed spacing inequality in this theorem. -/
theorem source_finite_width_tolerance {B η z τ0 : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B η f) (hη0 : 0≤η) (hη1 : η<1)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z)
    (hτ0 : 0<τ0) (hfloor : ∀n,τ0≤s.tolerance n) :
    ∃n,s.b n-s.a n≤s.tolerance n := by
  by_cases hb : ∀N,∃n,N≤n ∧ ModAB.RationalUncapped.isBisection (s.mode n)=true
  · have hw:=ModAB.RationalUncapped.infinite_bisections s.exactRun.policyRun hb
    obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw τ0 hτ0
    have hsmall:=hN N (le_refl N)
    change |s.b N-s.a N-0|<τ0 at hsmall
    have hn : s.b N-s.a N<τ0 := by
      exact (le_abs_self _).trans_lt (by simpa using hsmall)
    exact ⟨N,hn.le.trans (hfloor N)⟩
  push Not at hb
  obtain ⟨N,hN⟩:=hb
  have hp:=source_unique_root_convergence s hη0 hη1 hunique
  by_contra hn
  push Not at hn
  apply ModAB.RationalUncapped.no_convergent_separated_tail hp
    (show 0<τ0/2 by linarith)
  apply eventually_atTop.mpr
  refine ⟨N,fun n hNn=>?_⟩
  cases hm : s.mode n with
  | bisection =>
    have hfalse:=hN n hNn
    simp [hm,ModAB.RationalUncapped.isBisection] at hfalse
  | interpolation S t =>
    exact (show τ0/2≤s.tolerance n/2 by linarith [hfloor n]).trans
      (source_closing_separation s n hm (hn n))

theorem eta42_source_finite_width {B z τ0 : ℝ} {f : ℝ→ℝ}
    (s : SourceRun B (1/42) f)
    (hunique : ∀x∈Icc (s.a 0) (s.b 0),f x=0→x=z)
    (hτ0 : 0<τ0) (hfloor : ∀n,τ0≤s.tolerance n) :
    ∃n,s.b n-s.a n≤s.tolerance n :=
  source_finite_width_tolerance s (by norm_num) (by norm_num) hunique hτ0 hfloor

end
end ModAB.RationalAssembly
