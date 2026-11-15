import ModAB.RationalCycleUniform
import ModAB.RationalLocalAssembly

/-!
The negative-coefficient local pair orbit is constructed from the interpolation
coefficients themselves. Cone invariance, convergence and coefficient limits
are derived, rather than supplied as fields of a hypothetical orbit.
-/
namespace ModAB.Rational.CycleState
open Filter Set
open scoped Topology
open CycleUniform
noncomputable section

structure PairState where
  a : ℝ
  s : ℝ
  d : ℝ

def nextPair (g : ℝ → ℝ) (r : ℝ) (p : PairState) : PairState :=
  let b := firstError g r p.a p.s p.d
  ⟨secondCoefficient g r p.a p.s p.d * p.a * b * p.d, b, p.s⟩

def orbit (g : ℝ → ℝ) (r : ℝ) (p : PairState) : ℕ → PairState
  | 0 => p
  | n+1 => nextPair g r (orbit g r p n)

structure Bounds (g : ℝ → ℝ) (r θ ε C : ℝ) : Prop where
  theta_pos : 0 < θ
  theta_lt : θ < 1
  radius_pos : 0 < ε
  bound_pos : 0 < C
  small : C * ε^2 ≤ θ
  coefficients : ∀ a s d, CycleEntry.Cone θ ε a s d →
    firstCoefficient g r a s d < 0 ∧ secondCoefficient g r a s d < 0 ∧
    |firstCoefficient g r a s d| ≤ C ∧ |secondCoefficient g r a s d| ≤ C

theorem cone_enlarge_theta {θ θ' ε a s d : ℝ}
    (h : CycleEntry.Cone θ ε a s d) (hθ : θ ≤ θ') :
    CycleEntry.Cone θ' ε a s d := by
  exact ⟨h.same,h.opposite,
    h.latest_small.trans (mul_le_mul_of_nonneg_right hθ (abs_nonneg s)),
    h.opposite_small.trans (mul_le_mul_of_nonneg_right hθ (abs_nonneg d)),h.history_bound⟩

/-- Analytic guard neighbourhoods can be intersected with the coefficient
cone without losing the derived coefficient bounds. -/
theorem Bounds.restrict {g : ℝ → ℝ} {r θ ε C θ' ε' : ℝ}
    (N : Bounds g r θ ε C) (hθ0 : 0 < θ') (hθ : θ' ≤ θ)
    (hε0 : 0 < ε') (hε : ε' ≤ ε) (hsmall : C*ε'^2 ≤ θ') :
    Bounds g r θ' ε' C := by
  refine ⟨hθ0,hθ.trans_lt N.theta_lt,hε0,N.bound_pos,hsmall,?_⟩
  intro a s d hc
  exact N.coefficients _ _ _ ((cone_enlarge_theta hc hθ).enlarge_radius hε)

theorem exists_bounds {g : ℝ → ℝ} {r : ℝ}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r < 0) :
    ∃ θ ε C, Bounds g r θ ε C := by
  obtain ⟨θ,δ,hθ0,hθ1,hδ,hall⟩ := exists_uniform_negative_coefficients hg hgr hK
  let C := -2*localK g r
  have hC : 0 < C := by dsimp [C]; linarith
  let ε := min δ (min 1 (θ/C))
  have hε : 0 < ε := lt_min hδ (lt_min (by norm_num) (div_pos hθ0 hC))
  have hεδ : ε ≤ δ := min_le_left _ _
  have hε1 : ε ≤ 1 := (min_le_right _ _).trans (min_le_left _ _)
  have hεC : ε ≤ θ/C := (min_le_right _ _).trans (min_le_right _ _)
  have hCε : C*ε ≤ θ := by
    have hh := (le_div_iff₀ hC).mp hεC
    nlinarith
  refine ⟨θ,ε,C,⟨hθ0,hθ1,hε,hC,?_,?_⟩⟩
  · have hsq : ε^2 ≤ ε := by nlinarith
    exact (mul_le_mul_of_nonneg_left hsq hC.le).trans hCε
  · intro a s d hc
    have hh := hall a s d (hc.enlarge_radius hεδ)
    have hu : firstCoefficient g r a s d < 0 := by linarith [hh.1.2]
    have hv : secondCoefficient g r a s d < 0 := by linarith [hh.2.1.2]
    refine ⟨hu,hv,?_,?_⟩
    · rw [abs_of_neg hu]
      dsimp [C]
      linarith [hh.1.1]
    · rw [abs_of_neg hv]
      dsimp [C]
      linarith [hh.2.1.1]

theorem orbit_cone {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    ∀ n, CycleEntry.Cone θ ε (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d := by
  intro n
  induction n with
  | zero => exact hp
  | succ n ih =>
    obtain ⟨hu,hv,huC,hvC⟩ := N.coefficients _ _ _ ih
    have hc := CycleEntry.pair_preserves_cone N.theta_pos N.theta_lt
      N.bound_pos.le N.radius_pos.le N.small ih hu hv huC hvC
      (show firstError g r (orbit g r p n).a (orbit g r p n).s
        (orbit g r p n).d = _ by rfl) (show _ = _ by rfl)
    have he : θ*ε ≤ ε := by nlinarith [N.theta_lt,N.radius_pos]
    exact hc.enlarge_radius he

/-- The established abstract recurrence interface is now populated by a
concrete recursively defined interpolation orbit. -/
def toRun {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    CycleEntry.Run θ ε C where
  a n := (orbit g r p n).a
  s n := (orbit g r p n).s
  d n := (orbit g r p n).d
  u n := firstCoefficient g r (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d
  v n := secondCoefficient g r (orbit g r p n).a (orbit g r p n).s (orbit g r p n).d
  initial := hp
  u_neg n := (N.coefficients _ _ _ (orbit_cone N hp n)).1
  v_neg n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.1
  u_bound n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.2.1
  v_bound n := (N.coefficients _ _ _ (orbit_cone N hp n)).2.2.2
  first_step _n := rfl
  second_step _n := rfl
  drop_step _n := rfl

/-- Local regularity and the derived cone imply both coefficient limits;
neither convergence nor asymptotic separation is assumed. -/
theorem coefficient_limits {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    Tendsto (toRun N hp).u atTop (𝓝 (localK g r)) ∧
    Tendsto (toRun N hp).v atTop (𝓝 (localK g r)) := by
  let R := toRun N hp
  have hc := R.all_cones N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have he := R.errors_tendsto N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have hs := R.separation_tendsto N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small
  have hd (n : ℕ) := CycleUniform.Cone.distinct (hc n) N.theta_lt
  have hu := tendsto_firstCoefficient hg hgr (fun n => (hd n).1)
    (fun n => (hd n).2.1) (fun n => (hd n).2.2) he.1 he.2.1 he.2.2
  have hv := second_step_limits hg hgr hK (fun n => (hd n).1)
    (fun n => (hd n).2.1) (fun n => (hd n).2.2)
    (fun n => (hc n).latest_ne_zero) (fun n => (hc n).opposite_ne_zero)
    (fun n => (hc n).history_ne_zero) he.1 he.2.1 he.2.2 hs.1 hs.2
  exact ⟨hu,hv.2.2.2⟩

theorem orbit_hasQOrders {g : ℝ → ℝ} {r θ ε C : ℝ} {p : PairState}
    (hg : ContDiffAt ℝ 2 g r) (hgr : g r ≠ 0) (hK : localK g r ≠ 0)
    (N : Bounds g r θ ε C) (hp : CycleEntry.Cone θ ε p.a p.s p.d) :
    ModAB.RationalOrder.HasQOrder (fun n => (orbit g r p n).a) 3 |localK g r| ∧
    ModAB.RationalOrder.HasQOrder (fun n => (orbit g r p (n+1)).s) 3 |localK g r| ∧
    Tendsto (fun n => |(orbit g r p (n+1)).s| /
      |(orbit g r p n).a| ^ (9/5:ℝ)) atTop (𝓝 (|localK g r| ^ (2/5:ℝ))) ∧
    Tendsto (fun n => |(orbit g r p (n+1)).a| /
      |(orbit g r p (n+1)).s| ^ (5/3:ℝ)) atTop (𝓝 (|localK g r| ^ (1/3:ℝ))) := by
  have hh := coefficient_limits hg hgr hK N hp
  exact ModAB.RationalLocalAssembly.cycle_cone_hasQOrders (toRun N hp)
    N.theta_pos N.theta_lt N.bound_pos.le N.radius_pos.le N.small hK hh.1 hh.2

end
end ModAB.Rational.CycleState
