import ModAB.RationalCycleStateUpdate
import ModAB.RationalAlternating

/-!
Transport from a pointwise two-step source law to the entire safeguarded
trajectory. The law is about one local input state; invariance, iteration
identities, convergence and Q-order of the source trajectory are derived.
-/
namespace ModAB.Rational.CycleTransport
open Set Filter
open scoped Topology
open LocalState CycleState CycleUniform
noncomputable section

def view (r : ℝ) (s : State) : PairState :=
  ⟨latest s-r, opposite s-r, selectedHistory s-r⟩

def pairStep (f : ℝ → ℝ) (allow : State → Bool) (s : State) : State :=
  step f allow (step f allow s)

structure PairLaw (f g : ℝ → ℝ) (r θ ε : ℝ) (allow : State → Bool) : Prop where
  advance : ∀ s, Valid f r ε s →
    CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d →
    Valid f r ε (pairStep f allow s) ∧
    view r (pairStep f allow s) = nextPair g r (view r s) ∧
    latest (step f allow s)-r = (nextPair g r (view r s)).s ∧
    ∀ B, nextMode B f allow s = .interpolation s.scale (s.clock+1) ∧
      nextMode B f allow (step f allow s) =
        .interpolation (step f allow s).scale ((step f allow s).clock+1)

theorem orbit_even_succ (f : ℝ → ℝ) (allow : State → Bool) (s : State) (n : ℕ) :
    Alternating.orbit f allow s (2*(n+1)) =
      pairStep f allow (Alternating.orbit f allow s (2*n)) := by
  rw [show 2*(n+1) = (2*n+1)+1 by omega]
  rfl

theorem even_orbit {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d) :
    ∀ n, Valid f r ε (Alternating.orbit f allow s (2*n)) ∧
      view r (Alternating.orbit f allow s (2*n)) = CycleState.orbit g r (view r s) n := by
  intro n
  induction n with
  | zero => exact ⟨h,rfl⟩
  | succ n ih =>
    have hc : CycleEntry.Cone θ ε
        (view r (Alternating.orbit f allow s (2*n))).a
        (view r (Alternating.orbit f allow s (2*n))).s
        (view r (Alternating.orbit f allow s (2*n))).d := by
      rw [ih.2]
      exact orbit_cone N hp n
    have hh := L.advance _ ih.1 hc
    rw [orbit_even_succ]
    refine ⟨hh.1,?_⟩
    rw [hh.2.1,ih.2]
    rfl

theorem odd_latest {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d)
    (n : ℕ) :
    latest (Alternating.orbit f allow s (2*n+1))-r =
      (CycleState.orbit g r (view r s) (n+1)).s := by
  have he := even_orbit N L h hp n
  have hc : CycleEntry.Cone θ ε
      (view r (Alternating.orbit f allow s (2*n))).a
      (view r (Alternating.orbit f allow s (2*n))).s
      (view r (Alternating.orbit f allow s (2*n))).d := by
    rw [he.2]
    exact orbit_cone N hp n
  have hh := (L.advance _ he.1 hc).2.2.1
  rw [Alternating.orbit_succ,hh,he.2]
  rfl

/-- Both actual continuation decisions pass at every pair, even though the
odd state has a modified retained ordinate and need not be fresh. -/
theorem continue_pair {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d)
    (B : ℝ) (n : ℕ) :
    nextMode B f allow (Alternating.orbit f allow s (2*n)) =
      .interpolation (Alternating.orbit f allow s (2*n)).scale
        ((Alternating.orbit f allow s (2*n)).clock+1) ∧
    nextMode B f allow (Alternating.orbit f allow s (2*n+1)) =
      .interpolation (Alternating.orbit f allow s (2*n+1)).scale
        ((Alternating.orbit f allow s (2*n+1)).clock+1) := by
  have he := even_orbit N L h hp n
  have hc : CycleEntry.Cone θ ε
      (view r (Alternating.orbit f allow s (2*n))).a
      (view r (Alternating.orbit f allow s (2*n))).s
      (view r (Alternating.orbit f allow s (2*n))).d := by
    rw [he.2]
    exact orbit_cone N hp n
  exact (L.advance _ he.1 hc).2.2.2 B

/-- The actual evaluated subsequences have the proved orders. There are no
error-convergence or coefficient-limit assumptions in this transport. -/
theorem source_hasQOrders {f g : ℝ → ℝ} {r θ ε C : ℝ} {allow : State → Bool} {s : State}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (L : PairLaw f g r θ ε allow)
    (h : Valid f r ε s) (hp : CycleEntry.Cone θ ε (view r s).a (view r s).s (view r s).d) :
    ModAB.RationalOrder.HasQOrder
      (fun n => latest (Alternating.orbit f allow s (2*n))-r) 3 |localK g r| ∧
    ModAB.RationalOrder.HasQOrder
      (fun n => latest (Alternating.orbit f allow s (2*n+1))-r) 3 |localK g r| ∧
    Tendsto (fun n => |latest (Alternating.orbit f allow s (2*n+1))-r| /
      |latest (Alternating.orbit f allow s (2*n))-r| ^ (9/5:ℝ)) atTop
      (𝓝 (|localK g r| ^ (2/5:ℝ))) ∧
    Tendsto (fun n => |latest (Alternating.orbit f allow s (2*(n+1)))-r| /
      |latest (Alternating.orbit f allow s (2*n+1))-r| ^ (5/3:ℝ)) atTop
      (𝓝 (|localK g r| ^ (1/3:ℝ))) := by
  have hev (n : ℕ) : latest (Alternating.orbit f allow s (2*n))-r =
      (CycleState.orbit g r (view r s) n).a := by
    exact congrArg PairState.a (even_orbit N L h hp n).2
  have hod := odd_latest N L h hp
  simp_rw [hev,hod]
  exact orbit_hasQOrders hg hgr hK N hp

end
end ModAB.Rational.CycleTransport
