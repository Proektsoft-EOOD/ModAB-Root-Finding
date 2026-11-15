import ModAB.OneSided
import ModAB.TopEntry

/-! A simple-root counterexample to uncapped width termination at eta=1/42.
Every iterate is an actual top-convention AB secant. Fast point convergence
coexists with width tending to 1. The cap M=3 forces a return after secant 6. -/
namespace ModAB.Examples
open Set Filter Topology
noncomputable section

def eta42Function (x : ℝ) : ℝ := x+gluedQuadratic 0 (1/10) x

theorem eta42_formula (x : ℝ) :
    eta42Function x = if x≤0 then x else x+x^2/10 := by
  unfold eta42Function gluedQuadratic
  split_ifs <;> ring

theorem eta42_C1 : ContDiff ℝ 1 eta42Function :=
  contDiff_id.add (gluedQuadratic_C1 0 (1/10))

theorem eta42_left {e : ℝ} (he : 0≤e) : eta42Function (-e) = -e := by
  rw [eta42_formula,if_pos (by linarith)]

theorem eta42_right : eta42Function 1 = 11/10 := by norm_num [eta42_formula]

theorem eta42_zero_iff (x : ℝ) : eta42Function x=0 ↔ x=0 := by
  rw [eta42_formula]
  split_ifs with h
  · rfl
  · constructor
    · intro hz; nlinarith [sq_nonneg x]
    · intro hz; simp [hz]

theorem eta42_simple_root : HasDerivAt eta42Function 1 0 := by
  convert (hasDerivAt_id (0:ℝ)).add
    (gluedQuadratic_hasDerivAt_zero 0 (1/10)) using 1 <;> try norm_num
  funext x
  simp [eta42Function]

theorem eta42_strictMono : StrictMono eta42Function := by
  intro a b hab
  rw [eta42_formula,eta42_formula]
  split_ifs with ha hb
  · exact hab
  · nlinarith [sq_nonneg b]
  · linarith
  · nlinarith [mul_pos (sub_pos.mpr hab) (show 0<b+a by linarith)]

/-- The strict switching test is satisfied at the initial midpoint -1. -/
theorem eta42_hybrid_entry :
    eta42Function (-3)=-3 ∧ eta42Function 1=11/10 ∧ eta42Function (-1)=-1 ∧
    ModAB.originalL (-3) (11/10) (-1) <
      ModAB.originalR ((1-(19/20)/(41/10))^2) (-3) (11/10) (-1) := by
  norm_num [eta42_formula,ModAB.originalL,ModAB.originalR]

/-- n=0 denotes the first true secant, after the switching bisection. -/
def eta42Error : ℕ→ℝ
  | 0 => 1/21
  | n+1 => quadraticMap (eta42Error n)

def eta42State (n : ℕ) : LeftState := ⟨eta42Error n,1+eta42Error n,true⟩

theorem eta42_positive (n : ℕ) : 0<eta42Error n := by
  induction n with
  | zero => norm_num [eta42Error]
  | succ n ih => exact quadraticMap_positive ih

theorem eta42_step (n : ℕ) : leftStep 1 (eta42State n) = eta42State (n+1) :=
  quadratic_step (eta42_positive n)

theorem eta42_first_steps :
    leftOrbit 1 ⟨1,11/10,true⟩ 1 = eta42State 0 ∧
    eta42Error 1=1/483 ∧ eta42Error 2=1/234255 := by
  norm_num [leftOrbit,leftStep,leftTrial,eta42State,eta42Error,quadraticMap]

theorem eta42_orbit (n : ℕ) : leftOrbit 1 ⟨1,11/10,true⟩ (n+1)=eta42State n := by
  induction n with
  | zero => exact eta42_first_steps.1
  | succ n ih =>
    change leftStep 1 (leftOrbit 1 ⟨1,11/10,true⟩ (n+1)) = eta42State (n+1)
    rw [ih,eta42_step]

theorem eta42_shrinks (n : ℕ) : eta42Error (n+1) ≤ (1/2)*eta42Error n :=
  quadraticMap_half (eta42_positive n).le

theorem eta42_strict_decrease (n : ℕ) : eta42Error (n+1)<eta42Error n := by
  linarith [eta42_shrinks n,eta42_positive n]

theorem eta42_antitone : Antitone eta42Error :=
  antitone_nat_of_succ_le (fun n=>(eta42_strict_decrease n).le)

theorem eta42_bound (n : ℕ) : eta42Error n ≤ 1/21 :=
  eta42_antitone (Nat.zero_le n)

theorem eta42_next_bound (n : ℕ) : eta42Error (n+1) ≤ 1/483 := by
  rw [← eta42_first_steps.2.1]
  exact eta42_antitone (by omega)

theorem eta42_tendsto : Tendsto eta42Error atTop (𝓝 0) :=
  geometric_tendsto (fun n=>(eta42_positive n).le) (by norm_num)
    (by norm_num) eta42_shrinks

/-- Reciprocal invariant yielding the closed form without astronomical numerals. -/
theorem eta42_reciprocal (n : ℕ) : 1+1/eta42Error n = (22:ℝ)^(2^n) := by
  induction n with
  | zero => norm_num [eta42Error]
  | succ n ih =>
    have he := ne_of_gt (eta42_positive n)
    calc
      1+1/eta42Error (n+1) = (1+1/eta42Error n)^2 := by
        simp only [eta42Error,quadraticMap]
        field_simp
        <;> ring
      _ = ((22:ℝ)^(2^n))^2 := by rw [ih]
      _ = (22:ℝ)^(2^(n+1)) := by rw [← pow_mul, pow_succ]

theorem eta42_closed_form (n : ℕ) : eta42Error n = 1/((22:ℝ)^(2^n)-1) := by
  rw [← eta42_reciprocal]
  simp

theorem eta42_quadratic_order :
    Tendsto (fun n=>eta42Error (n+1)/(eta42Error n)^2) atTop (𝓝 1) := by
  have hq (n : ℕ) : eta42Error (n+1)/(eta42Error n)^2=1/(1+2*eta42Error n) := by
    have he := ne_of_gt (eta42_positive n)
    simp only [eta42Error,quadraticMap]
    field_simp
  simp_rw [hq]
  have ht := (tendsto_const_nhds : Tendsto (fun _ : ℕ=>(1:ℝ)) atTop (𝓝 1)).div
    ((eta42_tendsto.const_mul 2).const_add 1) (by norm_num : (1:ℝ)+2*0≠0)
  change Tendsto (fun n=>(1:ℝ)/(1+2*eta42Error n)) atTop (𝓝 (1/(1+2*(0:ℝ)))) at ht
  simpa only [mul_zero,add_zero,div_one,one_div] using ht

def eta42FullError (j : ℕ) : ℝ := (leftOrbit 1 ⟨1,11/10,true⟩ j).e

theorem eta42_full_tail (n : ℕ) : eta42FullError (n+1)=eta42Error n := by
  simp only [eta42FullError,eta42_orbit,eta42State]

theorem eta42_full_positive (j : ℕ) : 0<eta42FullError j := by
  cases j with
  | zero => norm_num [eta42FullError,leftOrbit]
  | succ n => rw [eta42_full_tail]; exact eta42_positive n

theorem eta42_full_decrease (j : ℕ) : eta42FullError (j+1)<eta42FullError j := by
  cases j with
  | zero => norm_num [eta42FullError,leftOrbit,leftStep,leftTrial]
  | succ n => rw [eta42_full_tail,eta42_full_tail]; exact eta42_strict_decrease n

theorem eta42_full_bound (j : ℕ) : eta42FullError j ≤ 1 := by
  cases j with
  | zero => norm_num [eta42FullError,leftOrbit]
  | succ n => rw [eta42_full_tail]; linarith [eta42_bound n]

theorem eta42_full_tendsto : Tendsto eta42FullError atTop (𝓝 0) := by
  apply (tendsto_add_atTop_iff_nat 1).mp
  simpa only [eta42_full_tail] using eta42_tendsto

theorem eta42_width_limit : Tendsto (fun j=>1+eta42FullError j) atTop (𝓝 1) := by
  simpa using eta42_full_tendsto.const_add 1

theorem eta42_full_ordinate_positive (j : ℕ) :
    0 < (leftOrbit 1 ⟨1,11/10,true⟩ j).v := by
  cases j with
  | zero => norm_num [leftOrbit]
  | succ n => rw [eta42_orbit]; change 0<1+eta42Error n; linarith [eta42_positive n]

/-- The un-clipped correction in leftStep coincides with GetABFactor on this orbit. -/
theorem eta42_positive_correction (j : ℕ) :
    0 < 1-eta42FullError (j+1)/eta42FullError j ∧
    ModAB.Convergence.abFactor (-eta42FullError (j+1)) (-eta42FullError j) =
      1-eta42FullError (j+1)/eta42FullError j := by
  have h := (div_lt_one (eta42_full_positive j)).mpr (eta42_full_decrease j)
  have hp : 0<1-eta42FullError (j+1)/eta42FullError j := by linarith
  exact ⟨hp,by simp [ModAB.Convergence.abFactor,hp]⟩

theorem eta42_actual_update (j : ℕ) :
    ModAB.Convergence.Updated eta42Function (-eta42FullError j) 1
      (-eta42FullError (j+1)) (-eta42FullError (j+1)) 1 ∧
    -eta42FullError (j+1)∈Ioo (-eta42FullError j) 1 ∧
    eta42Function (-eta42FullError (j+1))≠0 := by
  have he:=eta42_full_positive j
  have hn:=eta42_full_positive (j+1)
  have hd:=eta42_full_decrease j
  refine ⟨left_updated (p:=1) he hn (by simpa using eta42_left he.le)
    (by simpa using eta42_left hn.le),⟨by linarith,by linarith⟩,?_⟩
  rw [eta42_left hn.le]; linarith

/-- At and after the third secant, the true-residual ratio is at most 1/485. -/
theorem eta42_residual_ratio (n : ℕ) :
    eta42FullError (n+3)/eta42FullError (n+2) ≤ 1/485 := by
  rw [show n+3=(n+2)+1 by omega,show n+2=(n+1)+1 by omega,
    eta42_full_tail,eta42_full_tail]
  change quadraticMap (eta42Error (n+1))/eta42Error (n+1) ≤ 1/485
  have he:=eta42_positive (n+1)
  have hd : 0<1+2*eta42Error (n+1) := by positivity
  have hq : quadraticMap (eta42Error (n+1))/eta42Error (n+1) =
      eta42Error (n+1)/(1+2*eta42Error (n+1)) := by
    unfold quadraticMap; field_simp
  rw [hq]
  apply (div_le_iff₀ hd).mpr
  linarith [eta42_next_bound n]

theorem eta42_residual_exception (n : ℕ) :
    |eta42Function (-eta42FullError (n+3))| < (1/42)*
      min |eta42Function (-eta42FullError (n+2))| |eta42Function 1| := by
  have he:=eta42_full_positive (n+2)
  rw [eta42_left (eta42_full_positive (n+3)).le,eta42_left he.le,eta42_right,
    abs_neg,abs_neg,abs_of_pos (eta42_full_positive (n+3)),abs_of_pos he]
  norm_num only [abs_of_pos (by norm_num : (0:ℝ)<11/10)]
  rw [min_eq_left (show eta42FullError (n+2)≤11/10 by linarith [eta42_full_bound (n+2)])]
  have h:=(div_le_iff₀ he).mp (eta42_residual_ratio n)
  nlinarith only [h,he]

theorem eta42_width_test_fails (n : ℕ) :
    ModAB.Convergence.threshold 2 2 (n+3)<1+eta42FullError (n+3) := by
  have hp : (2:ℝ)^3≤2^(n+3):=pow_le_pow_right₀ (by norm_num) (by omega)
  have ht : ModAB.Convergence.threshold 2 2 (n+3)≤1:=by
    unfold ModAB.Convergence.threshold
    apply (div_le_iff₀ (by positivity : 0<(2:ℝ)^(n+3))).mpr
    norm_num at hp ⊢
    exact hp
  linarith [eta42_full_positive (n+3)]

def eta42UncappedContinues (j : ℕ) : Prop :=
  1+eta42FullError j≤ModAB.Convergence.threshold 2 2 j ∨
  |eta42Function (-eta42FullError j)| < (1/42)*
    min |eta42Function (-eta42FullError (j-1))| |eta42Function 1|

theorem eta42_uncapped_continues (j : ℕ) (hj : 1≤j) : eta42UncappedContinues j := by
  rcases j with _|_|_|n
  · omega
  · left; norm_num [eta42FullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · left; norm_num [eta42FullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]
  · right
    rw [show n+3-1=n+2 by omega]
    exact eta42_residual_exception n

theorem eta42_no_width_termination {τ : ℝ} (hτ : τ≤1) (j : ℕ) :
    τ<1+eta42FullError j := by linarith [eta42_full_positive j]

theorem eta42_mixed_tolerance_fails {ea er : ℝ} (ha : ea<1) :
    ∀ᶠj in atTop,ea+er*|(-eta42FullError j)|<1+eta42FullError j := by
  have ht := ((eta42_full_tendsto.neg).abs.const_mul er).const_add ea
  have ht' : Tendsto (fun j=>ea+er*|(-eta42FullError j)|) atTop (𝓝 ea):=by simpa using ht
  have hevent := (tendsto_order.mp ht').2 1 ha
  filter_upwards [hevent] with j hj
  linarith [eta42_full_positive j]

/-- The exact implemented continuation decision after true secant t+1. -/
def eta42NextMode (t c : ℕ) : ModAB.Convergence.Mode :=
  ModAB.Convergence.adaptiveMode 2 2 (1/42)
    |eta42Function (-eta42FullError (t+1))|
    (min |eta42Function (-eta42FullError t)| |eta42Function 1|)
    (1+eta42FullError (t+1)) 3 t c

theorem eta42_capped_exception (n c : ℕ) (hc : c<3) :
    eta42NextMode (n+2) c = .ab 2 (n+3) (c+1) := by
  have hw := not_le.mpr (eta42_width_test_fails n)
  have hr := eta42_residual_exception n
  simp only [eta42NextMode,ModAB.Convergence.adaptiveMode]
  rw [if_neg hw,if_pos ⟨hr,hc⟩]

theorem eta42_cap_forces_fallback (n : ℕ) : eta42NextMode (n+2) 3 = .bisection := by
  have hw := not_le.mpr (eta42_width_test_fails n)
  simp only [eta42NextMode,ModAB.Convergence.adaptiveMode]
  rw [if_neg hw]
  norm_num

/-- Exactly two width successes, then three exceptions, then fallback. -/
theorem eta42_capped_schedule :
    eta42NextMode 0 0 = .ab 2 1 0 ∧ eta42NextMode 1 0 = .ab 2 2 0 ∧
    eta42NextMode 2 0 = .ab 2 3 1 ∧ eta42NextMode 3 1 = .ab 2 4 2 ∧
    eta42NextMode 4 2 = .ab 2 5 3 ∧ eta42NextMode 5 3 = .bisection := by
  refine ⟨?_,?_,eta42_capped_exception 0 0 (by norm_num),
    eta42_capped_exception 1 1 (by norm_num),eta42_capped_exception 2 2 (by norm_num),
    eta42_cap_forces_fallback 3⟩
  all_goals norm_num [eta42NextMode,ModAB.Convergence.adaptiveMode,
    eta42FullError,leftOrbit,leftStep,leftTrial,ModAB.Convergence.threshold]

/-- Summary: this refutes uncapped width termination, not point convergence. -/
theorem eta42_uncapped_counterexample :
    ContDiff ℝ 1 eta42Function ∧ StrictMono eta42Function ∧
    (∀x,eta42Function x=0↔x=0) ∧ HasDerivAt eta42Function 1 0 ∧
    (∀j,1≤j→eta42UncappedContinues j) ∧
    Tendsto eta42FullError atTop (𝓝 0) ∧
    Tendsto (fun j=>1+eta42FullError j) atTop (𝓝 1) :=
  ⟨eta42_C1,eta42_strictMono,eta42_zero_iff,eta42_simple_root,
    eta42_uncapped_continues,eta42_full_tendsto,eta42_width_limit⟩

/-- The concrete switching bisection and its recorded initial orbit state. -/
theorem eta42_top_entry :
    ModAB.Convergence.Updated eta42Function (-3) 1 (-1) (-1) 1 ∧
    topLeftSeed (-3) (11/10) = ⟨1,11/10,true⟩ := by
  constructor
  · norm_num [ModAB.Convergence.Updated,eta42_formula]
  · exact topLeftSeed_eq _ _

end
end ModAB.Examples
