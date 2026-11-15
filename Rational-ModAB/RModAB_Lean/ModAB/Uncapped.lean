import ModAB.ExamplesBase
import Mathlib.Analysis.Calculus.Deriv.Inv

/-! `prop:uncapped_counterexample`: the actual AB trajectory and the paper's
two partial derivatives, invariant region, and residual/width conclusions. -/
namespace ModAB.Examples
open Set Filter Topology
noncomputable section

def phi (e r : ℝ) : ℝ := (1+e*(1+r))/(1+r+e*(1+2*r))

theorem phi_den_positive {e r : ℝ} (he : 0≤e) (hr : 0<r) :
    0<1+r+e*(1+2*r) := by positivity

theorem phi_hasDerivAt_ratio {e r : ℝ} (he : 0≤e) (hr : 0<r) :
    HasDerivAt (phi e) (-(1+e)^2/(1+r+e*(1+2*r))^2) r := by
  have hn:= (hasDerivAt_const r (1:ℝ)).add
    (((hasDerivAt_const r (1:ℝ)).add (hasDerivAt_id r)).const_mul e)
  have hd:= ((hasDerivAt_const r (1:ℝ)).add (hasDerivAt_id r)).add
    (((hasDerivAt_const r (1:ℝ)).add ((hasDerivAt_id r).const_mul 2)).const_mul e)
  have h:=hn.div hd (ne_of_gt (phi_den_positive he hr))
  convert h using 1 <;> (try funext x) <;> (try dsimp [phi]) <;> field_simp <;> ring

theorem phi_hasDerivAt_error {e r : ℝ} (he : 0≤e) (hr : 0<r) :
    HasDerivAt (fun t=>phi t r) (r^2/(1+r+e*(1+2*r))^2) e := by
  have hn:= (hasDerivAt_const e (1:ℝ)).add ((hasDerivAt_id e).mul_const (1+r))
  have hd:= (hasDerivAt_const e (1+r)).add ((hasDerivAt_id e).mul_const (1+2*r))
  have h:=hn.div hd (ne_of_gt (phi_den_positive he hr))
  convert h using 1 <;> (try funext x) <;> (try dsimp [phi]) <;> field_simp <;> ring

theorem phi_ratio_strictAntiOn {e : ℝ} (he : 0≤e) :
    StrictAntiOn (phi e) (Icc (1/2:ℝ) (7/10)) := by
  apply strictAntiOn_of_deriv_neg (convex_Icc _ _)
  · intro r hr
    exact (phi_hasDerivAt_ratio he (by linarith [hr.1])).continuousAt.continuousWithinAt
  · intro r hr
    have hi:=interior_subset hr
    rw [(phi_hasDerivAt_ratio he (by linarith [hi.1])).deriv]
    apply div_neg_of_neg_of_pos
    · have hp : 0<(1+e)^2:=by positivity
      linarith
    · exact sq_pos_of_pos (phi_den_positive he (by linarith [hi.1]))

theorem phi_error_strictMonoOn {r : ℝ} (hr : 0<r) :
    StrictMonoOn (fun e=>phi e r) (Icc (0:ℝ) (2/5)) := by
  apply strictMonoOn_of_deriv_pos (convex_Icc _ _)
  · intro e he
    exact (phi_hasDerivAt_error he.1 hr).continuousAt.continuousWithinAt
  · intro e he
    have hi:=interior_subset he
    rw [(phi_hasDerivAt_error hi.1 hr).deriv]
    exact div_pos (sq_pos_of_pos hr) (sq_pos_of_pos (phi_den_positive hi.1 hr))

theorem phi_region {e r : ℝ} (he : e∈Icc (0:ℝ) (2/5))
    (hr : r∈Icc (1/2:ℝ) (7/10)) :
    (10/17:ℝ)≤phi e r ∧ phi e r≤16/23 := by
  have hl:= (phi_error_strictMonoOn (r:=7/10) (by norm_num)).monotoneOn
    (by norm_num : (0:ℝ)∈Icc 0 (2/5)) he he.1
  have hl':= (phi_ratio_strictAntiOn he.1).antitoneOn hr
    (by norm_num : (7/10:ℝ)∈Icc (1/2) (7/10)) hr.2
  have hu:= (phi_ratio_strictAntiOn he.1).antitoneOn
    (by norm_num : (1/2:ℝ)∈Icc (1/2) (7/10)) hr hr.1
  have hu':= (phi_error_strictMonoOn (r:=1/2) (by norm_num)).monotoneOn he
    (by norm_num : (2/5:ℝ)∈Icc 0 (2/5)) he.2
  rw [show phi 0 (7/10)=10/17 by norm_num [phi]] at hl
  rw [show phi (2/5) (1/2)=16/23 by norm_num [phi]] at hu'
  exact ⟨hl.trans hl',hu.trans hu'⟩

def counterPair : ℕ→ℝ×ℝ
  | 0=>(146325/384296,1951/3176)
  | n+1=>let s:=counterPair n; let r:=phi s.1 s.2; (s.1*r,r)

def counterError (n : ℕ) : ℝ := (counterPair n).1
def counterRatio (n : ℕ) : ℝ := (counterPair n).2

theorem counter_recurrence (n : ℕ) :
    counterRatio (n+1)=phi (counterError n) (counterRatio n) ∧
    counterError (n+1)=counterError n*counterRatio (n+1) := by
  constructor <;> rfl

theorem counter_invariant (n : ℕ) :
    0<counterError n ∧ counterError n≤2/5 ∧ counterRatio n∈Icc (1/2:ℝ) (7/10) := by
  induction n with
  | zero=>norm_num [counterError,counterRatio,counterPair]
  | succ n ih=>
    have hp:=phi_region ⟨ih.1.le,ih.2.1⟩ ih.2.2
    have hr : counterRatio (n+1)∈Icc (1/2:ℝ) (7/10):=by
      rw [(counter_recurrence n).1]
      constructor <;> linarith [hp.1,hp.2]
    have hpos : 0<counterRatio (n+1):=by linarith [hr.1]
    rw [(counter_recurrence n).2]
    exact ⟨mul_pos ih.1 hpos,by nlinarith [hr.2],hr⟩

theorem counter_error_shrinks (n : ℕ) :
    counterError (n+1)≤(7/10)*counterError n := by
  rw [(counter_recurrence n).2]
  nlinarith [(counter_invariant n).1,(counter_invariant (n+1)).2.2.2]

theorem counterError_tendsto : Tendsto counterError atTop (𝓝 0) :=
  geometric_tendsto (fun n=>(counter_invariant n).1.le) (by norm_num)
    (by norm_num) counter_error_shrinks

theorem counter_width_limit : Tendsto (fun n=>1+counterError n) atTop (𝓝 1) := by
  simpa using counterError_tendsto.const_add 1

def counterPrevious : ℕ→ℝ
  | 0=>75/121
  | n+1=>counterError n

theorem counter_ratio_identity (n : ℕ) :
    counterError n=counterRatio n*counterPrevious n := by
  cases n with
  | zero=>norm_num [counterError,counterRatio,counterPair,counterPrevious]
  | succ n=>simpa only [counterPrevious,mul_comm] using (counter_recurrence n).2

theorem counterPrevious_positive (n : ℕ) : 0<counterPrevious n := by
  cases n with
  | zero=>norm_num [counterPrevious]
  | succ n=>exact (counter_invariant n).1

def counterState (n : ℕ) : LeftState :=
  ⟨counterError n,(1+counterError n)*(counterPrevious n+counterError n),true⟩

theorem counter_trial_identity {e p r : ℝ} (he : 0<e) (hp : 0<p) (hr : 0<r)
    (hrel : e=r*p) :
    leftTrial 2 e ((1+e)*(p+e))=e*phi e r := by
  have hd:=ne_of_gt (phi_den_positive he.le hr)
  have hv : (1+e)*(p+e)+e^2≠0:=by positivity
  simp only [leftTrial,phi]
  field_simp
  rw [hrel]
  ring

theorem counter_correction_identity {e t v : ℝ} (he : 0<e)
    (hv : 0<v) (ht : t=leftTrial 2 e v) :
    v*(1-t^2/e^2)=(1+t)*(e+t) := by
  have hd : v+e^2≠0:=by positivity
  have hbalance : (e-t)*v=e^2*(1+t) := by
    rw [ht]
    dsimp [leftTrial]
    field_simp
    ring
  have hmul:=congrArg (fun x : ℝ=>x*(e+t)) hbalance
  field_simp
  nlinarith only [hmul]

theorem counter_state_step (n : ℕ) : leftStep 2 (counterState n)=counterState (n+1) := by
  have hi:=counter_invariant n
  have hp:=counterPrevious_positive n
  have hr : 0<counterRatio n:=by linarith [hi.2.2.1]
  have hv : 0<(1+counterError n)*(counterPrevious n+counterError n):=
    mul_pos (by linarith [hi.1]) (by linarith [hi.1,hp])
  have ht:=counter_trial_identity hi.1 hp hr (counter_ratio_identity n)
  have ht' : leftTrial 2 (counterError n)
      ((1+counterError n)*(counterPrevious n+counterError n))=counterError (n+1) := by
    rw [ht,(counter_recurrence n).2,(counter_recurrence n).1]
  have hc:=counter_correction_identity hi.1 hv ht'.symm
  simp only [leftStep,counterState,ite_true]
  rw [ht',hc]
  rfl

theorem counter_first_steps :
    (leftOrbit 2 ⟨1,7,false⟩ 1).e=3/4 ∧
    (leftOrbit 2 ⟨1,7,false⟩ 2).e=75/121 ∧
    leftOrbit 2 ⟨1,7,false⟩ 3=counterState 0 := by
  norm_num [leftOrbit,leftStep,leftTrial,counterState,counterError,counterRatio,counterPair,counterPrevious]

theorem counter_orbit (n : ℕ) : leftOrbit 2 ⟨1,7,false⟩ (n+3)=counterState n := by
  induction n with
  | zero=>exact counter_first_steps.2.2
  | succ n ih=>
    change leftStep 2 (leftOrbit 2 ⟨1,7,false⟩ (n+3))=counterState (n+1)
    rw [ih]
    exact counter_state_step n

theorem counter_actual_update (n : ℕ) :
    ModAB.Convergence.Updated counterFunction (-counterError n) 1
      (-counterError (n+1)) (-counterError (n+1)) 1 := by
  exact left_updated (p:=2) (counter_invariant n).1 (counter_invariant (n+1)).1
    (counterFunction_left (counter_invariant n).1.le)
    (counterFunction_left (counter_invariant (n+1)).1.le)

theorem counter_residual_ratio (n : ℕ) :
    |counterFunction (-counterError n)| /
      min |counterFunction (-counterPrevious n)| |counterFunction 1|≤49/100 := by
  have hp:=counterPrevious_positive n
  have hi:=counter_invariant n
  have hprev : counterPrevious n≤1:=by
    cases n with
    | zero=>norm_num [counterPrevious]
    | succ k=>exact le_trans (counter_invariant k).2.1 (by norm_num)
  rw [counterFunction_left hi.1.le,counterFunction_left hp.le]
  have hf : counterFunction 1=7:=by norm_num [counterFunction,gluedQuadratic]
  rw [hf,abs_neg,abs_neg,abs_of_nonneg (sq_nonneg _),abs_of_nonneg (sq_nonneg _)]
  norm_num only [abs_of_pos (by norm_num : (0:ℝ)<7)]
  have hmin : min ((counterPrevious n)^2) (7:ℝ)=(counterPrevious n)^2:=by
    apply min_eq_left
    nlinarith [mul_nonneg hp.le (sub_nonneg.mpr hprev)]
  rw [hmin,counter_ratio_identity n]
  have hid : (counterRatio n*counterPrevious n)^2/(counterPrevious n)^2=(counterRatio n)^2:=by
    field_simp
  rw [hid]
  nlinarith [hi.2.2.1,hi.2.2.2]

theorem counter_hybrid_entry :
    counterFunction (-3)=-9 ∧ counterFunction 1=7 ∧ counterFunction (-1)=-1 ∧
    ModAB.originalL (-9) 7 (-1)<ModAB.originalR ((1-(1:ℝ)/16)^2) (-9) 7 (-1) := by
  norm_num [counterFunction,gluedQuadratic,ModAB.originalL,ModAB.originalR]

theorem counter_admissible_step (n : ℕ) :
    -counterError (n+1)∈Ioo (-counterError n) 1 ∧
    counterFunction (-counterError (n+1))≠0 ∧
    0<1-(counterError (n+1))^2/(counterError n)^2 := by
  have he:=(counter_invariant n).1
  have hn:=(counter_invariant (n+1)).1
  have hs:=counter_error_shrinks n
  have hlt : counterError (n+1)<counterError n:=by linarith
  constructor
  · constructor <;> linarith
  · constructor
    · intro hz
      have heq:=(counterFunction_zero_iff _).mp hz
      linarith
    · have hsq : (counterError (n+1))^2<(counterError n)^2:=by nlinarith
      have hratio : (counterError (n+1))^2/(counterError n)^2<1:=
        (div_lt_one (sq_pos_of_pos he)).mpr hsq
      linarith

theorem counter_residual_exception (n : ℕ) :
    |counterFunction (-counterError n)| < (1/2)*
      min |counterFunction (-counterPrevious n)| |counterFunction 1| := by
  have hp:=counterPrevious_positive n
  have hm : 0<min |counterFunction (-counterPrevious n)| |counterFunction 1| := by
    apply lt_min
    · rw [counterFunction_left hp.le,abs_neg,abs_of_nonneg (sq_nonneg _)]
      exact sq_pos_of_pos hp
    · norm_num [counterFunction,gluedQuadratic]
  have hb:=counter_residual_ratio n
  have hm' := (div_le_iff₀ hm).mp hb
  nlinarith only [hm,hm']

theorem counter_width_test_fails (n : ℕ) :
    ModAB.Convergence.threshold 2 2 (n+3)<1+counterError n := by
  have hp : (2:ℝ)^3≤2^(n+3):=pow_le_pow_right₀ (by norm_num) (by omega)
  have hpow : 0<(2:ℝ)^(n+3):=by positivity
  have ht : ModAB.Convergence.threshold 2 2 (n+3)≤1:=by
    unfold ModAB.Convergence.threshold
    apply (div_le_iff₀ hpow).mpr
    norm_num at hp ⊢
    exact hp
  linarith [(counter_invariant n).1]

def counterFullError (j : ℕ) : ℝ := (leftOrbit 2 ⟨1,7,false⟩ j).e

theorem counter_full_tail (n : ℕ) : counterFullError (n+3)=counterError n := by
  simp only [counterFullError,counter_orbit,counterState]

theorem counterFullError_positive (j : ℕ) : 0<counterFullError j := by
  rcases j with _|_|_|n
  · norm_num [counterFullError,leftOrbit]
  · norm_num [counterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [counterFullError,leftOrbit,leftStep,leftTrial]
  · exact counter_full_tail n ▸ (counter_invariant n).1

theorem counter_full_width_limit : Tendsto (fun j=>1+counterFullError j) atTop (𝓝 1) := by
  apply (tendsto_add_atTop_iff_nat 3).mp
  simpa only [counter_full_tail] using counter_width_limit

theorem counterFullError_tendsto : Tendsto counterFullError atTop (𝓝 0) := by
  apply (tendsto_add_atTop_iff_nat 3).mp
  simpa only [counter_full_tail] using counterError_tendsto

theorem counter_full_strict_decrease (j : ℕ) : counterFullError (j+1)<counterFullError j := by
  rcases j with _|_|_|n
  · norm_num [counterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [counterFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [counterFullError,leftOrbit,leftStep,leftTrial]
  · rw [show n+3+1=(n+1)+3 by omega,counter_full_tail,counter_full_tail]
    linarith [counter_error_shrinks n,(counter_invariant n).1]

theorem counter_full_ordinate_positive (j : ℕ) : 0<(leftOrbit 2 ⟨1,7,false⟩ j).v := by
  rcases j with _|_|_|n
  · norm_num [leftOrbit]
  · norm_num [leftOrbit,leftStep,leftTrial]
  · norm_num [leftOrbit,leftStep,leftTrial]
  · rw [counter_orbit]
    change 0<(1+counterError n)*(counterPrevious n+counterError n)
    exact mul_pos (by linarith [(counter_invariant n).1])
      (by linarith [(counter_invariant n).1,counterPrevious_positive n])

theorem counter_full_actual_update (j : ℕ) :
    ModAB.Convergence.Updated counterFunction (-counterFullError j) 1
      (-counterFullError (j+1)) (-counterFullError (j+1)) 1 ∧
    -counterFullError (j+1)∈Ioo (-counterFullError j) 1 ∧
    counterFunction (-counterFullError (j+1))≠0 := by
  have he:=counterFullError_positive j
  have hn:=counterFullError_positive (j+1)
  have hd:=counter_full_strict_decrease j
  refine ⟨left_updated (p:=2) he hn (counterFunction_left he.le)
    (counterFunction_left hn.le),⟨by linarith,by linarith⟩,?_⟩
  intro hz
  have heq:=(counterFunction_zero_iff _).mp hz
  linarith

/-- Exactly the two continuation tests, with the consecutive-exception cap removed. -/
def uncappedContinues (j : ℕ) : Prop :=
  1+counterFullError j≤ModAB.Convergence.threshold 2 2 j ∨
  |counterFunction (-counterFullError j)|<(1/2)*
    min |counterFunction (-counterFullError (j-1))| |counterFunction 1|

theorem counter_uncapped_continues (j : ℕ) (hj : 1≤j) : uncappedContinues j := by
  rcases j with _|_|_|n
  · omega
  · left; norm_num [counterFullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · left; norm_num [counterFullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · right
    have hp : counterFullError (n+3-1)=counterPrevious n := by
      cases n with
      | zero=>norm_num [counterFullError,leftOrbit,leftStep,leftTrial,counterPrevious]
      | succ n=>
        rw [show n+1+3-1=n+3 by omega,counter_full_tail]
        rfl
    rw [counter_full_tail,hp]
    exact counter_residual_exception n

theorem counter_no_width_termination {τ : ℝ} (hτ : τ≤1) (j : ℕ) :
    τ<1+counterFullError j := by linarith [counterFullError_positive j]

theorem counter_mixed_tolerance_fails {ea er : ℝ} (ha : ea<1) :
    ∀ᶠj in atTop,ea+er*|(-counterFullError j)|<1+counterFullError j := by
  have ht := ((counterFullError_tendsto.neg).abs.const_mul er).const_add ea
  have ht' : Tendsto (fun j=>ea+er*|(-counterFullError j)|) atTop (𝓝 ea):=by
    simpa using ht
  have hevent := (tendsto_order.mp ht').2 1 ha
  filter_upwards [hevent] with j hj
  linarith [counterFullError_positive j]

theorem counter_global_bound_violation :
    (4:ℝ)/(2:ℝ)^(14/7)<1+counterFullError 13 := by
  have h:=counterFullError_positive 13
  norm_num
  linarith

/-- All conclusions are derived from the concrete secants, not assumed in a run record. -/
theorem uncapped_counterexample :
    ContDiff ℝ 1 counterFunction ∧ (∀x,counterFunction x=0↔x=0) ∧
    (∀j,1≤j→uncappedContinues j) ∧
    Tendsto (fun j=>1+counterFullError j) atTop (𝓝 1) ∧
    (∀n,|counterFunction (-counterError n)| /
      min |counterFunction (-counterPrevious n)| |counterFunction 1|≤49/100) :=
  ⟨counterFunction_C1,counterFunction_zero_iff,counter_uncapped_continues,
    counter_full_width_limit,counter_residual_ratio⟩

end
end ModAB.Examples
