import ModAB.ExamplesBase

/-! `prop:one_sided`: the first two unmodified/corrected secants, the exact
quadratic recurrence, its error quotient, and the nonshrinking right endpoint. -/
namespace ModAB.Examples
open Filter Topology
noncomputable section

def quadraticMap (e : ℝ) : ℝ := e^2/(1+2*e)

def quadraticError : ℕ→ℝ
  | 0=>1/7
  | n+1=>quadraticMap (quadraticError n)

def quadraticState (n : ℕ) : LeftState := ⟨quadraticError n,1+quadraticError n,true⟩

theorem quadraticMap_positive {e : ℝ} (he : 0<e) : 0<quadraticMap e := by
  unfold quadraticMap
  positivity

theorem quadraticMap_half {e : ℝ} (he : 0≤e) : quadraticMap e≤(1/2)*e := by
  unfold quadraticMap
  apply (div_le_iff₀ (by positivity : 0<1+2*e)).mpr
  nlinarith

theorem quadraticError_positive (n : ℕ) : 0<quadraticError n := by
  induction n with
  | zero=>norm_num [quadraticError]
  | succ n ih=>exact quadraticMap_positive ih

theorem one_trial_formula (e v : ℝ) : leftTrial 1 e v=e*(v-1)/(v+e) := by
  simp only [leftTrial,pow_one]
  congr 1
  ring

theorem one_correction_identity {e v : ℝ} (he : 0<e) (hv : 1<v) :
    v*(1-leftTrial 1 e v/e)=v*(1+e)/(v+e) ∧
    v*(1+e)/(v+e)=1+leftTrial 1 e v := by
  have hd : v+e≠0:=by linarith
  rw [one_trial_formula]
  constructor <;> field_simp <;> ring

theorem quadratic_step {e : ℝ} (he : 0<e) :
    leftStep 1 ⟨e,1+e,true⟩=⟨quadraticMap e,1+quadraticMap e,true⟩ := by
  have hd : 1+2*e≠0:=ne_of_gt (by positivity)
  have hd' : 1+e+e≠0:=by linarith
  have ht : leftTrial 1 e (1+e)=quadraticMap e:=by
    rw [one_trial_formula]
    dsimp [quadraticMap]
    congr 1 <;> ring
  have hc:=one_correction_identity he (show 1<1+e by linarith)
  simp only [leftStep,pow_one,ite_true]
  rw [hc.1,hc.2,ht]

theorem one_first_steps :
    leftOrbit 1 ⟨1,2,false⟩ 1=⟨1/3,2,true⟩ ∧
    leftOrbit 1 ⟨1,2,false⟩ 2=quadraticState 0 := by
  norm_num [leftOrbit,leftStep,leftTrial,quadraticState,quadraticError]

theorem one_orbit (n : ℕ) : leftOrbit 1 ⟨1,2,false⟩ (n+2)=quadraticState n := by
  induction n with
  | zero=>exact one_first_steps.2
  | succ n ih=>
    change leftStep 1 (leftOrbit 1 ⟨1,2,false⟩ (n+2))=quadraticState (n+1)
    rw [ih]
    exact quadratic_step (quadraticError_positive n)

theorem one_actual_update (n : ℕ) :
    ModAB.Convergence.Updated oneFunction (-quadraticError n) 1
      (-quadraticError (n+1)) (-quadraticError (n+1)) 1 := by
  exact left_updated (p:=1) (quadraticError_positive n) (quadraticError_positive (n+1))
    (by simpa using oneFunction_left (quadraticError_positive n).le)
    (by simpa using oneFunction_left (quadraticError_positive (n+1)).le)

theorem quadraticError_bound (n : ℕ) :
    quadraticError n≤(1/2:ℝ)^n*(1/7) := by
  simpa only [quadraticError] using geometric_bound
    (e:=quadraticError) (by norm_num : (0:ℝ)≤1/2)
    (fun k=>quadraticMap_half (quadraticError_positive k).le) n

theorem quadraticError_tendsto : Tendsto quadraticError atTop (𝓝 0) :=
  geometric_tendsto (fun n=>(quadraticError_positive n).le) (by norm_num)
    (by norm_num) (fun n=>quadraticMap_half (quadraticError_positive n).le)

theorem quadratic_error_quotient (n : ℕ) :
    quadraticError (n+1)/(quadraticError n)^2=1/(1+2*quadraticError n) := by
  have he:=ne_of_gt (quadraticError_positive n)
  simp only [quadraticError,quadraticMap]
  field_simp

theorem quadratic_order :
    Tendsto (fun n=>quadraticError (n+1)/(quadraticError n)^2) atTop (𝓝 1) := by
  simp_rw [quadratic_error_quotient]
  have ht := (tendsto_const_nhds : Tendsto (fun _ : ℕ=>(1:ℝ)) atTop (𝓝 1)).div
    ((quadraticError_tendsto.const_mul 2).const_add 1) (by norm_num : (1:ℝ)+2*0≠0)
  change Tendsto (fun n=>(1:ℝ)/(1+2*quadraticError n)) atTop (𝓝 (1/(1+2*(0:ℝ)))) at ht
  simpa only [mul_zero,add_zero,div_one,one_div] using ht

theorem one_width_limit : Tendsto (fun n=>1+quadraticError n) atTop (𝓝 1) := by
  simpa using quadraticError_tendsto.const_add 1

theorem one_subsequent_values :
    quadraticError 1=1/63 ∧ quadraticError 2=1/4095 ∧ quadraticError 3=1/16777215 := by
  norm_num [quadraticError,quadraticMap]

theorem one_hybrid_entry :
    oneFunction (-3)=-3 ∧ oneFunction 1=2 ∧ oneFunction (-1)=-1 ∧
    ModAB.originalL (-3) 2 (-1)<ModAB.originalR ((1-(1/2)/5)^2) (-3) 2 (-1) := by
  norm_num [oneFunction,gluedQuadratic,ModAB.originalL,ModAB.originalR]

theorem one_admissible_step (n : ℕ) :
    -quadraticError (n+1)∈Set.Ioo (-quadraticError n) 1 ∧
    oneFunction (-quadraticError (n+1))≠0 ∧
    0<1-quadraticError (n+1)/quadraticError n := by
  have he:=quadraticError_positive n
  have hn:=quadraticError_positive (n+1)
  have hs:=quadraticMap_half he.le
  have hlt : quadraticError (n+1)<quadraticError n:=by
    change quadraticMap (quadraticError n)<quadraticError n
    linarith
  constructor
  · constructor <;> linarith
  · constructor
    · rw [oneFunction_left hn.le]; linarith
    · have hr := (div_lt_one he).mpr hlt
      linarith

theorem one_full_width_limit :
    Tendsto (fun j=>1+(leftOrbit 1 ⟨1,2,false⟩ j).e) atTop (𝓝 1) := by
  apply (tendsto_add_atTop_iff_nat 2).mp
  simpa only [one_orbit,quadraticState] using one_width_limit

def oneFullError (j : ℕ) : ℝ := (leftOrbit 1 ⟨1,2,false⟩ j).e

theorem oneFullError_positive (j : ℕ) : 0<oneFullError j := by
  rcases j with _|_|n
  · norm_num [oneFullError,leftOrbit]
  · norm_num [oneFullError,leftOrbit,leftStep,leftTrial]
  · simp only [oneFullError,one_orbit,quadraticState]
    exact quadraticError_positive n

theorem one_full_strict_decrease (j : ℕ) : oneFullError (j+1)<oneFullError j := by
  rcases j with _|_|n
  · norm_num [oneFullError,leftOrbit,leftStep,leftTrial]
  · norm_num [oneFullError,leftOrbit,leftStep,leftTrial]
  · unfold oneFullError
    rw [show n+2+1=(n+1)+2 by omega,one_orbit,one_orbit]
    change quadraticMap (quadraticError n)<quadraticError n
    linarith [quadraticMap_half (quadraticError_positive n).le,quadraticError_positive n]

theorem one_full_ordinate_positive (j : ℕ) : 0<(leftOrbit 1 ⟨1,2,false⟩ j).v := by
  rcases j with _|_|n
  · norm_num [leftOrbit]
  · norm_num [leftOrbit,leftStep,leftTrial]
  · rw [one_orbit]
    change 0<1+quadraticError n
    linarith [quadraticError_positive n]

theorem one_full_actual_update (j : ℕ) :
    ModAB.Convergence.Updated oneFunction (-oneFullError j) 1
      (-oneFullError (j+1)) (-oneFullError (j+1)) 1 ∧
    -oneFullError (j+1)∈Set.Ioo (-oneFullError j) 1 ∧
    oneFunction (-oneFullError (j+1))≠0 := by
  have he:=oneFullError_positive j
  have hn:=oneFullError_positive (j+1)
  have hd:=one_full_strict_decrease j
  refine ⟨left_updated (p:=1) he hn (by simpa using oneFunction_left he.le)
    (by simpa using oneFunction_left hn.le),⟨by linarith,by linarith⟩,?_⟩
  rw [oneFunction_left hn.le]
  linarith

theorem one_sided_quadratic_phase :
    ContDiff ℝ 1 oneFunction ∧
    (∀n,leftOrbit 1 ⟨1,2,false⟩ (n+2)=quadraticState n) ∧
    Tendsto quadraticError atTop (𝓝 0) ∧
    Tendsto (fun n=>quadraticError (n+1)/(quadraticError n)^2) atTop (𝓝 1) ∧
    Tendsto (fun j=>1+(leftOrbit 1 ⟨1,2,false⟩ j).e) atTop (𝓝 1) :=
  ⟨oneFunction_C1,one_orbit,quadraticError_tendsto,quadratic_order,one_full_width_limit⟩

end
end ModAB.Examples
