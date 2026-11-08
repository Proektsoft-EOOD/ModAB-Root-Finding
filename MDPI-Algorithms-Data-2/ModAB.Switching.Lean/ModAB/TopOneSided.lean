import ModAB.OneSided
import ModAB.TopEntry

/-! The top-convention hybrid entry in the discussion after Proposition 2. -/
namespace ModAB.Examples
open Set Filter Topology
noncomputable section

def topOneError : ℕ→ℝ
  | 0 => 1/3
  | n+1 => quadraticMap (topOneError n)

def topOneState (n : ℕ) : LeftState := ⟨topOneError n,1+topOneError n,true⟩

theorem topOne_positive (n : ℕ) : 0<topOneError n := by
  induction n with
  | zero => norm_num [topOneError]
  | succ n ih => exact quadraticMap_positive ih

theorem topOne_step (n : ℕ) : leftStep 1 (topOneState n)=topOneState (n+1) :=
  quadratic_step (topOne_positive n)

theorem topOne_first_steps :
    leftOrbit 1 ⟨1,2,true⟩ 1=topOneState 0 ∧
    topOneError 1=1/15 ∧ topOneError 2=1/255 ∧ topOneError 3=1/65535 := by
  norm_num [leftOrbit,leftStep,leftTrial,topOneState,topOneError,quadraticMap]

theorem topOne_orbit (n : ℕ) : leftOrbit 1 ⟨1,2,true⟩ (n+1)=topOneState n := by
  induction n with
  | zero => exact topOne_first_steps.1
  | succ n ih =>
    change leftStep 1 (leftOrbit 1 ⟨1,2,true⟩ (n+1))=topOneState (n+1)
    rw [ih,topOne_step]

theorem topOne_tendsto : Tendsto topOneError atTop (𝓝 0) :=
  geometric_tendsto (fun n=>(topOne_positive n).le) (by norm_num)
    (by norm_num) (fun n=>quadraticMap_half (topOne_positive n).le)

theorem topOne_quadratic_order :
    Tendsto (fun n=>topOneError (n+1)/(topOneError n)^2) atTop (𝓝 1) := by
  have hq (n : ℕ) : topOneError (n+1)/(topOneError n)^2=1/(1+2*topOneError n) := by
    have he:=ne_of_gt (topOne_positive n)
    simp only [topOneError,quadraticMap]
    field_simp
  simp_rw [hq]
  have ht := (tendsto_const_nhds : Tendsto (fun _ : ℕ=>(1:ℝ)) atTop (𝓝 1)).div
    ((topOne_tendsto.const_mul 2).const_add 1) (by norm_num : (1:ℝ)+2*0≠0)
  change Tendsto (fun n=>(1:ℝ)/(1+2*topOneError n)) atTop (𝓝 (1/(1+2*(0:ℝ)))) at ht
  simpa only [mul_zero,add_zero,div_one,one_div] using ht

def topOneFullError (j : ℕ) : ℝ := (leftOrbit 1 ⟨1,2,true⟩ j).e

theorem topOne_full_tail (n : ℕ) : topOneFullError (n+1)=topOneError n := by
  simp only [topOneFullError,topOne_orbit,topOneState]

theorem topOne_full_positive (j : ℕ) : 0<topOneFullError j := by
  cases j with
  | zero => norm_num [topOneFullError,leftOrbit]
  | succ n => rw [topOne_full_tail]; exact topOne_positive n

theorem topOne_full_decrease (j : ℕ) : topOneFullError (j+1)<topOneFullError j := by
  cases j with
  | zero => norm_num [topOneFullError,leftOrbit,leftStep,leftTrial]
  | succ n =>
    rw [topOne_full_tail,topOne_full_tail]
    change quadraticMap (topOneError n)<topOneError n
    linarith [quadraticMap_half (topOne_positive n).le,topOne_positive n]

theorem topOne_width_limit : Tendsto (fun j=>1+topOneFullError j) atTop (𝓝 1) := by
  apply (tendsto_add_atTop_iff_nat 1).mp
  simpa only [topOne_full_tail,add_zero] using topOne_tendsto.const_add 1

theorem topOne_actual_update (j : ℕ) :
    ModAB.Convergence.Updated oneFunction (-topOneFullError j) 1
      (-topOneFullError (j+1)) (-topOneFullError (j+1)) 1 ∧
    -topOneFullError (j+1)∈Ioo (-topOneFullError j) 1 ∧
    oneFunction (-topOneFullError (j+1))≠0 ∧
    0<1-topOneFullError (j+1)/topOneFullError j := by
  have he:=topOne_full_positive j
  have hn:=topOne_full_positive (j+1)
  have hd:=topOne_full_decrease j
  refine ⟨left_updated (p:=1) he hn (by simpa using oneFunction_left he.le)
    (by simpa using oneFunction_left hn.le),⟨by linarith,by linarith⟩,?_,?_⟩
  · rw [oneFunction_left hn.le]; linarith
  · linarith [(div_lt_one he).mpr hd]

/-- The concrete switching bisection and its recorded initial orbit state. -/
theorem topOne_top_entry :
    ModAB.Convergence.Updated oneFunction (-3) 1 (-1) (-1) 1 ∧
    topLeftSeed (-3) (2) = ⟨1,2,true⟩ := by
  constructor
  · norm_num [ModAB.Convergence.Updated,oneFunction_formula]
  · exact topLeftSeed_eq _ _

end
end ModAB.Examples
