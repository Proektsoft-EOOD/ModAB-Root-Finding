import ModAB.RationalLocalState

/-! Inductive local trajectories of the actual guarded interpolation update.
No convergence or error recurrence is assumed. The two discarded-side records
are refreshed by the true-sign transition and the coefficient limit follows
from the finite-C³ raw-step theorem. -/
namespace ModAB.Rational.Alternating
open Set Filter
open scoped Topology
open LocalState LocalNeighborhood
noncomputable section

def orbit (f : ℝ → ℝ) (allow : State → Bool) (initial : State) : ℕ → State
  | 0 => initial
  | n+1 => step f allow (orbit f allow initial n)

@[simp] theorem orbit_zero (f : ℝ → ℝ) (allow : State → Bool) (s : State) :
    orbit f allow s 0 = s := rfl

@[simp] theorem orbit_succ (f : ℝ → ℝ) (allow : State → Bool) (s : State) (n : ℕ) :
    orbit f allow s (n+1) = step f allow (orbit f allow s n) := rfl

theorem valid_orbit {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    ∀ n, Valid f r ε (orbit f allow s n) := by
  intro n
  induction n with
  | zero => exact h
  | succ n ih => exact (step_preserves N ih 2).1

theorem continue_orbit {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (B : ℝ) (n : ℕ) :
    nextMode B f allow (orbit f allow s n) =
      .interpolation (orbit f allow s n).scale ((orbit f allow s n).clock+1) :=
  (step_preserves N (valid_orbit N h n) B).2.2.2.1

theorem opposite_shift {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    opposite (orbit f allow s (n+1)) = latest (orbit f allow s n) :=
  (step_preserves N (valid_orbit N h n) 2).2.2.2.2.1

theorem history_two_shifts {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    selectedHistory (orbit f allow s (n+2)) = latest (orbit f allow s n) := by
  have hh := (step_preserves N (valid_orbit N h (n+1)) 2).2.2.2.2.2.1
  simpa only [← orbit_succ, Nat.add_assoc, Nat.reduceAdd, opposite_shift N h] using hh

theorem latest_ne_root {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : latest s ≠ r := by
  cases hs : s.recent <;> simp only [latest,hs]
  · exact ne_of_lt h.left
  · exact ne_of_gt h.right

/-- The source-state iteration is geometrically bounded and converges. -/
theorem latest_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => latest (orbit f allow s n)) atTop (𝓝 r) := by
  have habs : Tendsto (fun n => |latest (orbit f allow s n)-r|) atTop (𝓝 0) :=
    ModAB.RationalUncapped.geometric_residual
      (mul_nonneg N.bound_nonneg (sq_nonneg ε)) N.contraction
      (fun _ => abs_nonneg _) (fun n => (step_preserves N (valid_orbit N h n) 2).2.2.2.2.2.2)
  apply Metric.tendsto_atTop.mpr
  intro δ hδ
  obtain ⟨n,hn⟩ := Metric.tendsto_atTop.mp habs δ hδ
  exact ⟨n,fun k hk => by simpa [Real.dist_eq] using hn k hk⟩

theorem opposite_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => opposite (orbit f allow s n)) atTop (𝓝 r) := by
  apply (tendsto_add_atTop_iff_nat 1).mp
  simpa only [opposite_shift N h] using latest_tendsto N h

theorem history_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => selectedHistory (orbit f allow s n)) atTop (𝓝 r) := by
  apply (tendsto_add_atTop_iff_nat 2).mp
  simpa only [history_two_shifts N h] using latest_tendsto N h

theorem endpoints_min_max {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) :
    min (latest s) (opposite s) = s.a ∧ max (latest s) (opposite s) = s.b := by
  have hab : s.a ≤ s.b := le_of_lt (h.left.trans h.right)
  cases hs : s.recent <;> simp [latest,opposite,hs,min_eq_left hab,
    min_eq_right hab,max_eq_left hab,max_eq_right hab]

theorem endpoints_tendsto {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) :
    Tendsto (fun n => (orbit f allow s n).a) atTop (𝓝 r) ∧
    Tendsto (fun n => (orbit f allow s n).b) atTop (𝓝 r) := by
  have hl := latest_tendsto N h
  have ho := opposite_tendsto N h
  constructor
  · have heq : (fun n => (orbit f allow s n).a) =
        fun n => min (latest (orbit f allow s n)) (opposite (orbit f allow s n)) :=
      funext (fun n => (endpoints_min_max (valid_orbit N h n)).1.symm)
    rw [heq]
    simpa only [min_self] using hl.min ho
  · have heq : (fun n => (orbit f allow s n).b) =
        fun n => max (latest (orbit f allow s n)) (opposite (orbit f allow s n)) :=
      funext (fun n => (endpoints_min_max (valid_orbit N h n)).2.symm)
    rw [heq]
    simpa only [max_self] using hl.max ho

theorem endpoint_product (s : State) (r : ℝ) :
    (s.a-r)*(s.b-r) = (latest s-r)*(opposite s-r) := by
  cases hs : s.recent <;> simp [latest,opposite,hs,mul_comm]

/-- The three-iterate recurrence is a consequence of history updates, not an
orbit hypothesis. It starts once both discarded-side records have refreshed. -/
theorem raw_recurrence {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s) (n : ℕ) :
    latest (orbit f allow s (n+3))-r =
      rawCoefficient f r (orbit f allow s (n+2)).a (orbit f allow s (n+2)).b
        (selectedHistory (orbit f allow s (n+2))) *
      (latest (orbit f allow s (n+2))-r) *
      (latest (orbit f allow s (n+1))-r) *
      (latest (orbit f allow s n)-r) := by
  have hv := valid_orbit N h (n+2)
  have hp := step_preserves N hv 2
  have hx : latest (orbit f allow s (n+3)) =
      rawPoint (orbit f allow s (n+2)).a (orbit f allow s (n+2)).b
        (f (orbit f allow s (n+2)).a) (f (orbit f allow s (n+2)).b)
        (selectedHistory (orbit f allow s (n+2)))
        (f (selectedHistory (orbit f allow s (n+2)))) := by
    rw [show n+3 = (n+2)+1 by omega, orbit_succ,step,latest_update]
    exact hp.2.2.1.trans hp.2.1
  rw [hx,raw_error_eq_coefficient hv.admissible.2.2.2.1 hv.admissible.2.2.2.2.1
    hv.admissible.2.2.2.2.2]
  have he := endpoint_product (orbit f allow s (n+2)) r
  have ho := opposite_shift N h (n+1)
  simp only [Nat.add_assoc,Nat.reduceAdd] at ho
  rw [ho] at he
  rw [history_two_shifts N h]
  linear_combination rawCoefficient f r (orbit f allow s (n+2)).a
    (orbit f allow s (n+2)).b (latest (orbit f allow s n)) *
    (latest (orbit f allow s n)-r) * he

/-- Complete trajectory-level Q-order theorem. Local estimates ensure all
operational tests, state invariance and convergence; finite C³ regularity
supplies the coefficient limit for the actual rational proposal. -/
theorem hasQOrder_of_neighborhood {f : ℝ → ℝ} {r ε : ℝ} {allow : State → Bool}
    (N : Neighborhood f r ε allow) {s : State} (h : Valid f r ε s)
    (hf : ContDiffAt ℝ 3 f r) (hz : f r=0) (hd : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r ≠ 0) :
    ModAB.RationalOrder.HasQOrder (fun n => latest (orbit f allow s n)-r)
      ModAB.RationalOrder.tau
      (|AnalyticBridge.errorConstant f r| ^ ((ModAB.RationalOrder.tau-1)/2)) := by
  let aa := fun n => (orbit f allow s (n+2)).a
  let bb := fun n => (orbit f allow s (n+2)).b
  let dd := fun n => selectedHistory (orbit f allow s (n+2))
  have had : ∀ n, Admissible r (aa n) (bb n) (dd n) :=
    fun n => (valid_orbit N h (n+2)).admissible
  have ha : Tendsto aa atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (endpoints_tendsto N h).1
  have hb : Tendsto bb atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (endpoints_tendsto N h).2
  have hh : Tendsto dd atTop (𝓝 r) :=
    (tendsto_add_atTop_iff_nat 2).mpr (history_tendsto N h)
  have hcoeff := SmoothBridge.tendsto_normalized_rawPoint_error hf hz hd
    (fun n => (had n).1) (fun n => (had n).2.1) (fun n => (had n).2.2.1)
    (fun n => (had n).2.2.2.1) (fun n => (had n).2.2.2.2.1)
    (fun n => (had n).2.2.2.2.2) ha hb hh
  apply ModAB.RationalOrder.three_error_recurrence_hasQOrder
    (fun n => sub_ne_zero.mpr (latest_ne_root (valid_orbit N h n))) hK
    (by simpa using (latest_tendsto N h).sub_const r) hcoeff
  exact raw_recurrence N h

end
end ModAB.Rational.Alternating
