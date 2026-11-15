import ModAB.RationalBridge
import ModAB.RationalResidual
import ModAB.RationalSpacing

/-!
Uncapped continuation for a bracket-preserving trial method.

This module does not impose the old cap on residual exceptions.  A continuing
phase advances its clock after every interpolation update and may continue only
by satisfying the geometric width test or the true-residual decrease test.
Switching, curvature checks, and working ordinates are deliberately unrestricted.
The conclusions concern exact arithmetic; no binary64 claim is made here.
-/
namespace ModAB.RationalUncapped
open Set Filter Topology
open ModAB.Convergence
noncomputable section

inductive Mode where
  | bisection
  | interpolation (scale : ℝ) (clock : ℕ)

/-- `m` is the smaller true endpoint residual before the update. -/
inductive Step (B η m : ℝ) : ℝ → ℝ → Mode → ℝ → ℝ → Mode → Prop
  | bisectStay (w R R') : Step B η m w R .bisection (w/2) R' .bisection
  | bisectStart (w R R') :
      Step B η m w R .bisection (w/2) R' (.interpolation (w/2) 0)
  | width {w w' R R' S t} : 0≤w' → w'≤w →
      w'≤threshold B S (t+1) →
      Step B η m w R (.interpolation S t) w' R' (.interpolation S (t+1))
  | residual {w w' R R' S t} : 0≤w' → w'≤w →
      threshold B S (t+1)<w' → R'<η*m →
      Step B η m w R (.interpolation S t) w' R' (.interpolation S (t+1))
  | stop {w w' R R' S t} : 0≤w' → w'≤w →
      Step B η m w R (.interpolation S t) w' R' .bisection

/-- The actual uncapped continuation decision.  Neither successful width tests
nor residual exceptions reset the phase clock. -/
def nextMode (B η S m w' R' : ℝ) (t : ℕ) : Mode :=
  if w'≤threshold B S (t+1) then .interpolation S (t+1)
  else if R'<η*m then .interpolation S (t+1)
  else .bisection

theorem decision_step {B η S m w R w' R' : ℝ} {t : ℕ}
    (h0 : 0≤w') (hmono : w'≤w) :
    Step B η m w R (.interpolation S t) w' R' (nextMode B η S m w' R' t) := by
  unfold nextMode
  split_ifs with hw hr
  · exact .width h0 hmono hw
  · exact .residual h0 hmono (lt_of_not_ge hw) hr
  · exact .stop h0 hmono

theorem forty_two_test (R m : ℝ) : R<(1/42:ℝ)*m ↔ 42*R<m := by
  constructor <;> intro h <;> linarith

/-- Residual `R n` belongs to the point evaluated most recently.  That point is
an endpoint, so the smaller endpoint residual is at most `R n`. -/
structure Run (B η : ℝ) where
  width : ℕ→ℝ
  residual : ℕ→ℝ
  minimum : ℕ→ℝ
  mode : ℕ→Mode
  nonneg_initial : 0≤width 0
  nonneg_residual : ∀n,0≤residual n
  minimum_le : ∀n,minimum n≤residual n
  next : ∀n,Step B η (minimum n) (width n) (residual n) (mode n)
    (width (n+1)) (residual (n+1)) (mode (n+1))

theorem step_nonnegative {B η m w R w' R' : ℝ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : 0≤w) : 0≤w' := by
  cases hs with
  | bisectStay => positivity
  | bisectStart => positivity
  | width h _ _ => exact h
  | residual h _ _ _ => exact h
  | stop h _ => exact h

theorem run_nonnegative {B η : ℝ} (r : Run B η) (n : ℕ) : 0≤r.width n := by
  induction n with
  | zero => exact r.nonneg_initial
  | succ n ih => exact step_nonnegative (r.next n) ih

theorem step_nonincrease {B η m w R w' R' : ℝ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : 0≤w) : w'≤w := by
  cases hs with
  | bisectStay => linarith
  | bisectStart => linarith
  | width _ h _ => exact h
  | residual _ h _ _ => exact h
  | stop _ h => exact h

theorem width_antitone {B η : ℝ} (r : Run B η) : Antitone r.width :=
  antitone_nat_of_succ_le (fun n=>step_nonincrease (r.next n) (run_nonnegative r n))

def isBisection : Mode→Bool
  | .bisection=>true
  | .interpolation _ _=>false

def bisectionCount (m : ℕ→Mode) : ℕ→ℕ
  | 0=>0
  | n+1=>bisectionCount m n + if isBisection (m n) then 1 else 0

theorem count_monotone (m : ℕ→Mode) : Monotone (bisectionCount m) := by
  apply monotone_nat_of_le_succ
  intro n
  simp only [bisectionCount]
  omega

theorem step_count_bound {B η m W0 w R w' R' : ℝ} {k : ℕ} {s s' : Mode}
    (hs : Step B η m w R s w' R' s') (hw : w≤W0/(2:ℝ)^k) :
    w'≤W0/(2:ℝ)^(k+if isBisection s then 1 else 0) := by
  cases hs with
  | bisectStay =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | bisectStart =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | width _ hm _ => simpa [isBisection] using hm.trans hw
  | residual _ hm _ _ => simpa [isBisection] using hm.trans hw
  | stop _ hm => simpa [isBisection] using hm.trans hw

theorem bisection_count_bound {B η : ℝ} (r : Run B η) (n : ℕ) :
    r.width n≤r.width 0/(2:ℝ)^(bisectionCount r.mode n) := by
  induction n with
  | zero => simp [bisectionCount]
  | succ n ih => exact step_count_bound (r.next n) ih

theorem count_unbounded {m : ℕ→Mode}
    (h : ∀N,∃n,N≤n ∧ isBisection (m n)=true) : ∀k,∃n,k≤bisectionCount m n := by
  intro k
  induction k with
  | zero => exact ⟨0,by simp [bisectionCount]⟩
  | succ k ih =>
      obtain ⟨n,hn⟩:=ih
      obtain ⟨j,hnj,hj⟩:=h n
      refine ⟨j+1,?_⟩
      have hm:=count_monotone m hnj
      simp only [bisectionCount,hj,ite_true]
      omega

theorem infinite_bisections {B η : ℝ} (r : Run B η)
    (h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  have hc : Tendsto (bisectionCount r.mode) atTop atTop := by
    apply tendsto_atTop_atTop.mpr
    intro k
    obtain ⟨n,hn⟩:=count_unbounded h k
    exact ⟨n,fun j hj=>hn.trans (count_monotone r.mode hj)⟩
  exact squeeze_zero (run_nonnegative r) (bisection_count_bound r)
    ((dyadic_tendsto _).comp hc)

/-- There is no reset of a continuing phase's scale or clock. -/
theorem continuing_step {B η m w R w' R' S : ℝ} {t : ℕ} {s' : Mode}
    (hs : Step B η m w R (.interpolation S t) w' R' s')
    (hno : isBisection s'=false) :
    s'=.interpolation S (t+1) ∧
      (w'≤threshold B S (t+1) ∨ R'<η*m) := by
  cases hs with
  | width _ _ ht => exact ⟨rfl,Or.inl ht⟩
  | residual _ _ _ hr => exact ⟨rfl,Or.inr hr⟩
  | stop => simp [isBisection] at hno

theorem continuing_suffix {B η S : ℝ} {N t : ℕ} (r : Run B η)
    (hm : r.mode N=.interpolation S t)
    (hno : ∀n,N≤n→isBisection (r.mode n)=false) :
    ∀k,r.mode (N+k)=.interpolation S (t+k) := by
  intro k
  induction k with
  | zero => simpa using hm
  | succ k ih =>
      have hs:=r.next (N+k)
      rw [ih] at hs
      have he:=(continuing_step hs (hno (N+k+1) (by omega))).1
      simpa [Nat.add_assoc] using he

theorem small_widths_tendsto {B η : ℝ} (r : Run B η)
    (h : ∀ε:ℝ,0<ε→∃n,r.width n<ε) : Tendsto r.width atTop (𝓝 0) := by
  apply Metric.tendsto_atTop.mpr
  intro ε hε
  obtain ⟨N,hN⟩:=h ε hε
  refine ⟨N,fun n hn=>?_⟩
  rw [Real.dist_eq,sub_zero,abs_of_nonneg (run_nonnegative r n)]
  exact (width_antitone r hn).trans_lt hN

theorem geometric_residual {R : ℕ→ℝ} {η : ℝ} (hη0 : 0≤η) (hη1 : η<1)
    (hR : ∀n,0≤R n) (hstep : ∀n,R (n+1)≤η*R n) :
    Tendsto R atTop (𝓝 0) := by
  have hb : ∀n,R n≤η^n*R 0 := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
      calc R (n+1)≤η*R n:=hstep n
           _≤η*(η^n*R 0):=mul_le_mul_of_nonneg_left ih hη0
           _=η^(n+1)*R 0:=by rw [pow_succ]; ring
  exact squeeze_zero hR hb
    (by simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one hη0 hη1).mul_const (R 0))

/-- The uncapped policy guarantees width convergence or residual convergence.
This is derived from the transition system, not postulated as a run property. -/
theorem width_or_residual {B η : ℝ} (hη0 : 0≤η) (hη1 : η<1) (r : Run B η) :
    Tendsto r.width atTop (𝓝 0) ∨ Tendsto r.residual atTop (𝓝 0) := by
  by_cases hsmall : ∀ε:ℝ,0<ε→∃n,r.width n<ε
  · exact Or.inl (small_widths_tendsto r hsmall)
  push Not at hsmall
  obtain ⟨ε,hε,hfloor⟩:=hsmall
  by_cases hb : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true
  · exact Or.inl (infinite_bisections r hb)
  push Not at hb
  obtain ⟨N,hN⟩:=hb
  have hno : ∀n,N≤n→isBisection (r.mode n)=false := by
    intro n hn
    exact Bool.eq_false_iff.mpr (hN n hn)
  obtain ⟨S,t,hm⟩ : ∃S t,r.mode N=.interpolation S t := by
    cases he : r.mode N with
    | bisection => simpa [he,isBisection] using hno N (le_refl N)
    | interpolation S t => exact ⟨S,t,rfl⟩
  have hmode:=continuing_suffix r hm hno
  have hlim : Tendsto (fun k:ℕ=>threshold B S (t+k+1)) atTop (𝓝 0) := by
    convert dyadic_tendsto (2*B*S/(2:ℝ)^(t+1)) using 1
    funext k
    simp only [threshold,show t+k+1=(t+1)+k by omega,pow_add,div_mul_eq_div_div]
  obtain ⟨K,hK⟩:=Metric.tendsto_atTop.mp hlim ε hε
  have hcontract : ∀k,r.residual (N+K+k+1)≤η*r.residual (N+K+k) := by
    intro k
    have hs:=r.next (N+(K+k))
    rw [hmode (K+k)] at hs
    have hor:=(continuing_step hs (hno (N+(K+k)+1) (by omega))).2
    have ht : threshold B S (t+(K+k)+1)<ε := by
      have hk:=hK (K+k) (by omega)
      rw [Real.dist_eq,sub_zero] at hk
      exact (le_abs_self _).trans_lt hk
    rcases hor with hw|hr
    · have hf:=hfloor (N+(K+k)+1)
      linarith
    · have hh:=hr.le.trans (mul_le_mul_of_nonneg_left (r.minimum_le (N+(K+k))) hη0)
      simpa [Nat.add_assoc] using hh
  have hr : Tendsto (fun k=>r.residual (N+K+k)) atTop (𝓝 0) :=
    geometric_residual hη0 hη1 (fun k=>r.nonneg_residual _) (by simpa [Nat.add_assoc] using hcontract)
  exact Or.inr ((tendsto_add_atTop_iff_nat (N+K)).mp (by simpa [Nat.add_comm] using hr))

/-- Any positive width-or-residual stopping tolerances are attained.  Continuity
is unnecessary for this scheduling statement. -/
theorem finite_combined_tolerance {B η τw τf : ℝ} (hη0 : 0≤η) (hη1 : η<1)
    (hτw : 0<τw) (hτf : 0<τf) (r : Run B η) :
    ∃n,r.width n<τw ∨ r.residual n<τf := by
  rcases width_or_residual hη0 hη1 r with hw|hr
  · obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hw τw hτw
    refine ⟨N,Or.inl ?_⟩
    have h:=hN N (le_refl N)
    simpa [Real.dist_eq,abs_of_nonneg (run_nonnegative r N)] using h
  · obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp hr τf hτf
    refine ⟨N,Or.inr ?_⟩
    have h:=hN N (le_refl N)
    simpa [Real.dist_eq,abs_of_nonneg (r.nonneg_residual N)] using h

/-- A nonterminating exact run.  `point 0` is an initial endpoint and `point
(n+1)` is the new trial inserted by update `n`. -/
structure ExactRun (B η : ℝ) (f : ℝ→ℝ) where
  a : ℕ→ℝ
  b : ℕ→ℝ
  point : ℕ→ℝ
  mode : ℕ→Mode
  initial_order : a 0<b 0
  initial_sign : f (a 0)*f (b 0)<0
  initial_point : point 0=a 0 ∨ point 0=b 0
  continuous : ContinuousOn f (Icc (a 0) (b 0))
  candidate_mem : ∀n,point (n+1)∈Icc (a n) (b n)
  nonzero : ∀n,f (point (n+1))≠0
  update : ∀n,Updated f (a n) (b n) (point (n+1)) (a (n+1)) (b (n+1))
  policy : ∀n,Step B η (min |f (a n)| |f (b n)|)
    (b n-a n) |f (point n)| (mode n)
    (b (n+1)-a (n+1)) |f (point (n+1))| (mode (n+1))

theorem exact_invariants {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) :
    ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 ∧
      Icc (r.a (n+1)) (r.b (n+1))⊆Icc (r.a n) (r.b n) := by
  have hs : ∀n,r.a n<r.b n ∧ f (r.a n)*f (r.b n)<0 := by
    intro n
    induction n with
    | zero => exact ⟨r.initial_order,r.initial_sign⟩
    | succ n ih =>
      exact (ModAB.RationalBridge.updated_invariants_closed ih.1 (r.candidate_mem n)
        ih.2 (r.nonzero n) (r.update n)).imp_right And.left
  intro n
  exact ⟨(hs n).1,(hs n).2,
    (ModAB.RationalBridge.updated_invariants_closed (hs n).1 (r.candidate_mem n)
      (hs n).2 (r.nonzero n) (r.update n)).2.2⟩

theorem point_endpoint {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n=r.a n ∨ r.point n=r.b n := by
  cases n with
  | zero => exact r.initial_point
  | succ n =>
    have hu:=r.update n
    unfold Updated at hu
    split_ifs at hu with hs
    · exact Or.inr (congrArg Prod.snd hu).symm
    · exact Or.inl (congrArg Prod.fst hu).symm

theorem point_mem {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n∈Icc (r.a n) (r.b n) := by
  rcases point_endpoint r n with h|h
  · rw [h]; exact ⟨le_refl _,(exact_invariants r n).1.le⟩
  · rw [h]; exact ⟨(exact_invariants r n).1.le,le_refl _⟩

theorem point_mem_initial {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) (n : ℕ) :
    r.point n∈Icc (r.a 0) (r.b 0) := by
  have hm : Antitone (fun n=>Icc (r.a n) (r.b n)) :=
    antitone_nat_of_succ_le (fun n=>(exact_invariants r n).2.2)
  exact hm (Nat.zero_le n) (point_mem r n)

def ExactRun.policyRun {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) : Run B η where
  width n:=r.b n-r.a n
  residual n:=|f (r.point n)|
  minimum n:=min |f (r.a n)| |f (r.b n)|
  mode:=r.mode
  nonneg_initial:=by linarith [r.initial_order]
  nonneg_residual n:=abs_nonneg _
  minimum_le n:=by
    rcases point_endpoint r n with h|h
    · rw [h]; exact min_le_left _ _
    · rw [h]; exact min_le_right _ _
  next:=r.policy

theorem exact_zeros {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f) :
    ∀n,∃z∈Icc (r.a n) (r.b n),f z=0 := by
  have hi:=exact_invariants r
  have hm : Antitone (fun n=>Icc (r.a n) (r.b n)):=
    antitone_nat_of_succ_le (fun n=>(hi n).2.2)
  intro n
  exact bracket_contains_zero (hi n).1.le
    (r.continuous.mono (hm (Nat.zero_le n))) (hi n).2.1

theorem residual_of_width {B η : ℝ} {f : ℝ→ℝ} (r : ExactRun B η f)
    (hw : Tendsto (fun n=>r.b n-r.a n) atTop (𝓝 0)) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) := by
  obtain ⟨z,hz,hfz⟩:=nested_intervals_root
    (fun n=>(exact_invariants r n).2.2) (exact_zeros r) r.continuous hw
  have hzall : ∀n,z∈Icc (r.a n) (r.b n) := by
    apply mem_iInter.mp
    rw [hz]
    exact mem_singleton _
  have hp:=chosen_points_tendsto hzall (point_mem r) hw
  have hf : Tendsto (fun n=>f (r.point n)) atTop (𝓝 (f z)) :=
    (r.continuous z (hzall 0)).tendsto.comp
      (tendsto_nhdsWithin_iff.mpr ⟨hp,Filter.Eventually.of_forall (point_mem_initial r)⟩)
  simpa [hfz] using hf.abs

/-- Full residual convergence for continuous sign-bracketed functions and the
actual uncapped width-or-residual continuation schedule. -/
theorem global_residual_convergence {B η : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) := by
  rcases width_or_residual hη0 hη1 r.policyRun with hw|hr
  · exact residual_of_width r hw
  · exact hr

theorem finite_residual_tolerance {B η τ : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (hτ : 0<τ) (r : ExactRun B η f) :
    ∃N,∀n,N≤n→|f (r.point n)|<τ := by
  obtain ⟨N,hN⟩:=Metric.tendsto_atTop.mp (global_residual_convergence hη0 hη1 r) τ hτ
  refine ⟨N,fun n hn=>?_⟩
  simpa [Real.dist_eq] using hN n hn

/-- The concrete factor used by the uncapped implementation. -/
theorem eta42_global_residual {B : ℝ} {f : ℝ→ℝ} (r : ExactRun B (1/42) f) :
    Tendsto (fun n=>|f (r.point n)|) atTop (𝓝 0) :=
  global_residual_convergence (by norm_num) (by norm_num) r

theorem eventually_near_zero_set {B η ε : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (hε : 0<ε) (r : ExactRun B η f) :
    ∀ᶠn in atTop,∃z∈Icc (r.a 0) (r.b 0),f z=0 ∧ |r.point n-z|<ε :=
  ModAB.RationalResidual.eventually_near_zero_set r.continuous (point_mem_initial r)
    (global_residual_convergence hη0 hη1 r) hε

theorem unique_root_convergence {B η z : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f)
    (hunique : ∀x∈Icc (r.a 0) (r.b 0),f x=0→x=z) :
    Tendsto r.point atTop (𝓝 z) :=
  ModAB.RationalResidual.tendsto_unique_zero r.continuous (point_mem_initial r)
    (global_residual_convergence hη0 hη1 r) hunique

/-- Conditional width termination from a positive closing distance.  The
spacing premise refers to accepted/evaluated points; it must be checked for
the implementation's closing-step and rounding guards. -/
theorem width_stop_with_closing {B η z δ τ : ℝ} {f : ℝ→ℝ}
    (hη0 : 0≤η) (hη1 : η<1) (r : ExactRun B η f)
    (hunique : ∀x∈Icc (r.a 0) (r.b 0),f x=0→x=z) (hδ : 0<δ)
    (hclose : ∀ᶠn in atTop,τ<r.b n-r.a n→δ≤|r.point (n+1)-r.point n|) :
    ∃n,r.b n-r.a n≤τ :=
  exists_width_stop_of_converges (unique_root_convergence hη0 hη1 r hunique) hδ hclose

end
end ModAB.RationalUncapped
