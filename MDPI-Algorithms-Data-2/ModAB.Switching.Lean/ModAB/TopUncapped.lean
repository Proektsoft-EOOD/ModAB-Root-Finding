import ModAB.Uncapped
import ModAB.TopEntry

/-! The f101 counterexample with the switching bisection recorded (top convention).
The original bottom-convention proof remains in Uncapped.lean. -/
namespace ModAB.Examples
open Set Filter Topology
noncomputable section

def topCounterPair : ℕ→ℝ×ℝ
  | 0=>(2664111/8512544,24001/36692)
  | n+1=>let s:=topCounterPair n; let r:=phi s.1 s.2; (s.1*r,r)

def topCounterError (n : ℕ) : ℝ := (topCounterPair n).1
def topCounterRatio (n : ℕ) : ℝ := (topCounterPair n).2

theorem topCounter_recurrence (n : ℕ) :
    topCounterRatio (n+1)=phi (topCounterError n) (topCounterRatio n) ∧
    topCounterError (n+1)=topCounterError n*topCounterRatio (n+1) := by
  constructor <;> rfl

theorem topCounter_invariant (n : ℕ) :
    0<topCounterError n ∧ topCounterError n≤2/5 ∧ topCounterRatio n∈Icc (1/2:ℝ) (7/10) := by
  induction n with
  | zero=>norm_num [topCounterError,topCounterRatio,topCounterPair]
  | succ n ih=>
    have hp:=phi_region ⟨ih.1.le,ih.2.1⟩ ih.2.2
    have hr : topCounterRatio (n+1)∈Icc (1/2:ℝ) (7/10):=by
      rw [(topCounter_recurrence n).1]
      constructor <;> linarith [hp.1,hp.2]
    have hpos : 0<topCounterRatio (n+1):=by linarith [hr.1]
    rw [(topCounter_recurrence n).2]
    exact ⟨mul_pos ih.1 hpos,by nlinarith [hr.2],hr⟩

theorem topCounter_error_shrinks (n : ℕ) :
    topCounterError (n+1)≤(7/10)*topCounterError n := by
  rw [(topCounter_recurrence n).2]
  nlinarith [(topCounter_invariant n).1,(topCounter_invariant (n+1)).2.2.2]

theorem topCounterError_tendsto : Tendsto topCounterError atTop (𝓝 0) :=
  geometric_tendsto (fun n=>(topCounter_invariant n).1.le) (by norm_num)
    (by norm_num) topCounter_error_shrinks

theorem topCounter_width_limit : Tendsto (fun n=>1+topCounterError n) atTop (𝓝 1) := by
  simpa using topCounterError_tendsto.const_add 1

def topCounterPrevious : ℕ→ℝ
  | 0=>111/232
  | n+1=>topCounterError n

theorem topCounter_ratio_identity (n : ℕ) :
    topCounterError n=topCounterRatio n*topCounterPrevious n := by
  cases n with
  | zero=>norm_num [topCounterError,topCounterRatio,topCounterPair,topCounterPrevious]
  | succ n=>simpa only [topCounterPrevious,mul_comm] using (topCounter_recurrence n).2

theorem topCounterPrevious_positive (n : ℕ) : 0<topCounterPrevious n := by
  cases n with
  | zero=>norm_num [topCounterPrevious]
  | succ n=>exact (topCounter_invariant n).1

def topCounterState (n : ℕ) : LeftState :=
  ⟨topCounterError n,(1+topCounterError n)*(topCounterPrevious n+topCounterError n),true⟩

theorem topCounter_state_step (n : ℕ) : leftStep 2 (topCounterState n)=topCounterState (n+1) := by
  have hi:=topCounter_invariant n
  have hp:=topCounterPrevious_positive n
  have hr : 0<topCounterRatio n:=by linarith [hi.2.2.1]
  have hv : 0<(1+topCounterError n)*(topCounterPrevious n+topCounterError n):=
    mul_pos (by linarith [hi.1]) (by linarith [hi.1,hp])
  have ht:=counter_trial_identity hi.1 hp hr (topCounter_ratio_identity n)
  have ht' : leftTrial 2 (topCounterError n)
      ((1+topCounterError n)*(topCounterPrevious n+topCounterError n))=topCounterError (n+1) := by
    rw [ht,(topCounter_recurrence n).2,(topCounter_recurrence n).1]
  have hc:=counter_correction_identity hi.1 hv ht'.symm
  simp only [leftStep,topCounterState,ite_true]
  rw [ht',hc]
  rfl

theorem topCounter_first_steps :
    (leftOrbit 2 ⟨1,7,true⟩ 1).e=3/4 ∧
    (leftOrbit 2 ⟨1,7,true⟩ 2).e=111/232 ∧
    leftOrbit 2 ⟨1,7,true⟩ 3=topCounterState 0 := by
  norm_num [leftOrbit,leftStep,leftTrial,topCounterState,topCounterError,topCounterRatio,topCounterPair,topCounterPrevious]

theorem topCounter_orbit (n : ℕ) : leftOrbit 2 ⟨1,7,true⟩ (n+3)=topCounterState n := by
  induction n with
  | zero=>exact topCounter_first_steps.2.2
  | succ n ih=>
    change leftStep 2 (leftOrbit 2 ⟨1,7,true⟩ (n+3))=topCounterState (n+1)
    rw [ih]
    exact topCounter_state_step n

theorem topCounter_actual_update (n : ℕ) :
    ModAB.Convergence.Updated counterFunction (-topCounterError n) 1
      (-topCounterError (n+1)) (-topCounterError (n+1)) 1 := by
  exact left_updated (p:=2) (topCounter_invariant n).1 (topCounter_invariant (n+1)).1
    (counterFunction_left (topCounter_invariant n).1.le)
    (counterFunction_left (topCounter_invariant (n+1)).1.le)

theorem topCounter_residual_ratio (n : ℕ) :
    |counterFunction (-topCounterError n)| /
      min |counterFunction (-topCounterPrevious n)| |counterFunction 1|≤49/100 := by
  have hp:=topCounterPrevious_positive n
  have hi:=topCounter_invariant n
  have hprev : topCounterPrevious n≤1:=by
    cases n with
    | zero=>norm_num [topCounterPrevious]
    | succ k=>exact le_trans (topCounter_invariant k).2.1 (by norm_num)
  rw [counterFunction_left hi.1.le,counterFunction_left hp.le]
  have hf : counterFunction 1=7:=by norm_num [counterFunction,gluedQuadratic]
  rw [hf,abs_neg,abs_neg,abs_of_nonneg (sq_nonneg _),abs_of_nonneg (sq_nonneg _)]
  norm_num only [abs_of_pos (by norm_num : (0:ℝ)<7)]
  have hmin : min ((topCounterPrevious n)^2) (7:ℝ)=(topCounterPrevious n)^2:=by
    apply min_eq_left
    nlinarith [mul_nonneg hp.le (sub_nonneg.mpr hprev)]
  rw [hmin,topCounter_ratio_identity n]
  have hid : (topCounterRatio n*topCounterPrevious n)^2/(topCounterPrevious n)^2=(topCounterRatio n)^2:=by
    field_simp
  rw [hid]
  nlinarith [hi.2.2.1,hi.2.2.2]

theorem topCounter_hybrid_entry :
    counterFunction (-3)=-9 ∧ counterFunction 1=7 ∧ counterFunction (-1)=-1 ∧
    ModAB.originalL (-9) 7 (-1)<ModAB.originalR ((1-(1:ℝ)/16)^2) (-9) 7 (-1) := by
  norm_num [counterFunction,gluedQuadratic,ModAB.originalL,ModAB.originalR]

theorem topCounter_admissible_step (n : ℕ) :
    -topCounterError (n+1)∈Ioo (-topCounterError n) 1 ∧
    counterFunction (-topCounterError (n+1))≠0 ∧
    0<1-(topCounterError (n+1))^2/(topCounterError n)^2 := by
  have he:=(topCounter_invariant n).1
  have hn:=(topCounter_invariant (n+1)).1
  have hs:=topCounter_error_shrinks n
  have hlt : topCounterError (n+1)<topCounterError n:=by linarith
  constructor
  · constructor <;> linarith
  · constructor
    · intro hz
      have heq:=(counterFunction_zero_iff _).mp hz
      linarith
    · have hsq : (topCounterError (n+1))^2<(topCounterError n)^2:=by nlinarith
      have hratio : (topCounterError (n+1))^2/(topCounterError n)^2<1:=
        (div_lt_one (sq_pos_of_pos he)).mpr hsq
      linarith

theorem topCounter_residual_exception (n : ℕ) :
    |counterFunction (-topCounterError n)| < (1/2)*
      min |counterFunction (-topCounterPrevious n)| |counterFunction 1| := by
  have hp:=topCounterPrevious_positive n
  have hm : 0<min |counterFunction (-topCounterPrevious n)| |counterFunction 1| := by
    apply lt_min
    · rw [counterFunction_left hp.le,abs_neg,abs_of_nonneg (sq_nonneg _)]
      exact sq_pos_of_pos hp
    · norm_num [counterFunction,gluedQuadratic]
  have hb:=topCounter_residual_ratio n
  have hm' := (div_le_iff₀ hm).mp hb
  nlinarith only [hm,hm']

theorem topCounter_width_test_fails (n : ℕ) :
    ModAB.Convergence.threshold 2 2 (n+3)<1+topCounterError n := by
  have hp : (2:ℝ)^3≤2^(n+3):=pow_le_pow_right₀ (by norm_num) (by omega)
  have hpow : 0<(2:ℝ)^(n+3):=by positivity
  have ht : ModAB.Convergence.threshold 2 2 (n+3)≤1:=by
    unfold ModAB.Convergence.threshold
    apply (div_le_iff₀ hpow).mpr
    norm_num at hp ⊢
    exact hp
  linarith [(topCounter_invariant n).1]

def topCounterFullError (j : ℕ) : ℝ := (leftOrbit 2 ⟨1,7,true⟩ j).e

theorem topCounter_full_tail (n : ℕ) : topCounterFullError (n+3)=topCounterError n := by
  simp only [topCounterFullError,topCounter_orbit,topCounterState]

theorem topCounterFullError_positive (j : ℕ) : 0<topCounterFullError j := by
  rcases j with _|_|_|n
  · norm_num [topCounterFullError,leftOrbit]
  · norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial]
  · exact topCounter_full_tail n ▸ (topCounter_invariant n).1

theorem topCounter_full_width_limit : Tendsto (fun j=>1+topCounterFullError j) atTop (𝓝 1) := by
  apply (tendsto_add_atTop_iff_nat 3).mp
  simpa only [topCounter_full_tail] using topCounter_width_limit

theorem topCounterFullError_tendsto : Tendsto topCounterFullError atTop (𝓝 0) := by
  apply (tendsto_add_atTop_iff_nat 3).mp
  simpa only [topCounter_full_tail] using topCounterError_tendsto

theorem topCounter_full_strict_decrease (j : ℕ) : topCounterFullError (j+1)<topCounterFullError j := by
  rcases j with _|_|_|n
  · norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial]
  · rw [show n+3+1=(n+1)+3 by omega,topCounter_full_tail,topCounter_full_tail]
    linarith [topCounter_error_shrinks n,(topCounter_invariant n).1]

theorem topCounter_full_ordinate_positive (j : ℕ) : 0<(leftOrbit 2 ⟨1,7,true⟩ j).v := by
  rcases j with _|_|_|n
  · norm_num [leftOrbit]
  · norm_num [leftOrbit,leftStep,leftTrial]
  · norm_num [leftOrbit,leftStep,leftTrial]
  · rw [topCounter_orbit]
    change 0<(1+topCounterError n)*(topCounterPrevious n+topCounterError n)
    exact mul_pos (by linarith [(topCounter_invariant n).1])
      (by linarith [(topCounter_invariant n).1,topCounterPrevious_positive n])

theorem topCounter_full_actual_update (j : ℕ) :
    ModAB.Convergence.Updated counterFunction (-topCounterFullError j) 1
      (-topCounterFullError (j+1)) (-topCounterFullError (j+1)) 1 ∧
    -topCounterFullError (j+1)∈Ioo (-topCounterFullError j) 1 ∧
    counterFunction (-topCounterFullError (j+1))≠0 := by
  have he:=topCounterFullError_positive j
  have hn:=topCounterFullError_positive (j+1)
  have hd:=topCounter_full_strict_decrease j
  refine ⟨left_updated (p:=2) he hn (counterFunction_left he.le)
    (counterFunction_left hn.le),⟨by linarith,by linarith⟩,?_⟩
  intro hz
  have heq:=(counterFunction_zero_iff _).mp hz
  linarith

/-- Exactly the two continuation tests, with the consecutive-exception cap removed. -/
def topUncappedContinues (j : ℕ) : Prop :=
  1+topCounterFullError j≤ModAB.Convergence.threshold 2 2 j ∨
  |counterFunction (-topCounterFullError j)|<(1/2)*
    min |counterFunction (-topCounterFullError (j-1))| |counterFunction 1|

theorem topCounter_uncapped_continues (j : ℕ) (hj : 1≤j) : topUncappedContinues j := by
  rcases j with _|_|_|n
  · omega
  · left; norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · left; norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · right
    have hp : topCounterFullError (n+3-1)=topCounterPrevious n := by
      cases n with
      | zero=>norm_num [topCounterFullError,leftOrbit,leftStep,leftTrial,topCounterPrevious]
      | succ n=>
        rw [show n+1+3-1=n+3 by omega,topCounter_full_tail]
        rfl
    rw [topCounter_full_tail,hp]
    exact topCounter_residual_exception n

theorem topCounter_no_width_termination {τ : ℝ} (hτ : τ≤1) (j : ℕ) :
    τ<1+topCounterFullError j := by linarith [topCounterFullError_positive j]

theorem topCounter_mixed_tolerance_fails {ea er : ℝ} (ha : ea<1) :
    ∀ᶠj in atTop,ea+er*|(-topCounterFullError j)|<1+topCounterFullError j := by
  have ht := ((topCounterFullError_tendsto.neg).abs.const_mul er).const_add ea
  have ht' : Tendsto (fun j=>ea+er*|(-topCounterFullError j)|) atTop (𝓝 ea):=by
    simpa using ht
  have hevent := (tendsto_order.mp ht').2 1 ha
  filter_upwards [hevent] with j hj
  linarith [topCounterFullError_positive j]

theorem topCounter_global_bound_violation :
    (4:ℝ)/(2:ℝ)^(14/7)<1+topCounterFullError 13 := by
  have h:=topCounterFullError_positive 13
  norm_num
  linarith

/-- All conclusions are derived from the concrete secants, not assumed in a run record. -/
theorem topUncapped_counterexample :
    ContDiff ℝ 1 counterFunction ∧ (∀x,counterFunction x=0↔x=0) ∧
    (∀j,1≤j→topUncappedContinues j) ∧
    Tendsto (fun j=>1+topCounterFullError j) atTop (𝓝 1) ∧
    (∀n,|counterFunction (-topCounterError n)| /
      min |counterFunction (-topCounterPrevious n)| |counterFunction 1|≤49/100) :=
  ⟨counterFunction_C1,counterFunction_zero_iff,topCounter_uncapped_continues,
    topCounter_full_width_limit,topCounter_residual_ratio⟩

/-- The concrete switching bisection and its recorded initial orbit state. -/
theorem topCounter_top_entry :
    ModAB.Convergence.Updated counterFunction (-3) 1 (-1) (-1) 1 ∧
    topLeftSeed (-9) (7) = ⟨1,7,true⟩ := by
  constructor
  · norm_num [ModAB.Convergence.Updated,counterFunction_formula]
  · exact topLeftSeed_eq _ _

/-- The positive correction holds from the first chord, not only in the tail. -/
theorem topCounter_full_correction (j : ℕ) :
    0 < 1-(topCounterFullError (j+1))^2/(topCounterFullError j)^2 := by
  have he := topCounterFullError_positive j
  have hn := topCounterFullError_positive (j+1)
  have hd := topCounter_full_strict_decrease j
  have hs : (topCounterFullError (j+1))^2 < (topCounterFullError j)^2 := by nlinarith
  linarith [(div_lt_one (sq_pos_of_pos he)).mpr hs]

end
end ModAB.Examples
