import ModAB.Schedule

/-! Infinite exact-arithmetic runs. The qualitative proof deliberately follows
both bisection cases in `thm:global_simple`, independently of amortized counting. -/
namespace ModAB.Convergence
open Filter Topology
noncomputable section

inductive Step (B : ℝ) (M : ℕ) : ℝ → Mode → ℝ → Mode → Prop
  | bisectStay (w) : Step B M w .bisection (w/2) .bisection
  | bisectStart (w) : Step B M w .bisection (w/2) (.ab (w/2) 0 0)
  | abWidth {w w' S t c} : 0≤w' → w'≤w → w'≤threshold B S (t+1) →
      Step B M w (.ab S t c) w' (.ab S (t+1) 0)
  | abException {w w' S t c} : 0≤w' → w'≤w → threshold B S (t+1)<w' → c<M →
      Step B M w (.ab S t c) w' (.ab S (t+1) (c+1))
  | abStop {w w' S t c} : 0≤w' → w'≤w → Step B M w (.ab S t c) w' .bisection

structure Run (B : ℝ) (M : ℕ) where
  width : ℕ→ℝ
  mode : ℕ→Mode
  nonneg_initial : 0≤width 0
  initial : mode 0=.bisection
  next : ∀n,Step B M (width n) (mode n) (width (n+1)) (mode (n+1))

theorem step_history {B W0 w w' : ℝ} {M n : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (h : History B M W0 n w m) :
    History B M W0 (n+1) w' m' := by
  cases hs with
  | bisectStay => exact .bisectStay h
  | bisectStart => exact .bisectStart h
  | abWidth h0 hm ht => exact .abWidth h h0 hm ht
  | abException h0 hm ht hc => exact .abException h h0 hm ht hc
  | abStop h0 hm => exact .abStop h h0 hm

theorem step_nonnegative {B w w' : ℝ} {M : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (hw : 0≤w) : 0≤w' := by
  cases hs with
  | bisectStay => positivity
  | bisectStart => positivity
  | abWidth h0 _ _ => exact h0
  | abException h0 _ _ _ => exact h0
  | abStop h0 _ => exact h0

theorem history_of_run {B : ℝ} {M : ℕ} (r : Run B M) (n : ℕ) :
    History B M (r.width 0) n (r.width n) (r.mode n) := by
  induction n with
  | zero => rw [r.initial]; exact .initial
  | succ n ih => exact step_history (r.next n) ih

theorem run_nonnegative {B : ℝ} {M : ℕ} (r : Run B M) : ∀n,0≤r.width n := by
  intro n
  induction n with
  | zero => exact r.nonneg_initial
  | succ n ih => exact step_nonnegative (r.next n) ih

def isBisection : Mode→Bool
  | .bisection=>true
  | .ab _ _ _=>false

def bisectionCount (m : ℕ→Mode) : ℕ→ℕ
  | 0=>0
  | n+1=>bisectionCount m n + if isBisection (m n) then 1 else 0

theorem count_monotone (m : ℕ→Mode) : Monotone (bisectionCount m) := by
  apply monotone_nat_of_le_succ
  intro n
  simp only [bisectionCount]
  omega

theorem step_count_bound {B W0 w w' : ℝ} {M k : ℕ} {m m' : Mode}
    (hs : Step B M w m w' m') (hw : w≤W0/(2:ℝ)^k) :
    w'≤W0/(2:ℝ)^(k+if isBisection m then 1 else 0) := by
  cases hs with
  | bisectStay =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | bisectStart =>
      simp only [isBisection,ite_true,pow_succ]
      exact (div_le_div_of_nonneg_right hw (by norm_num : (0:ℝ)≤2)).trans_eq (by ring)
  | abWidth _ hm _ => simpa [isBisection] using hm.trans hw
  | abException _ hm _ _ => simpa [isBisection] using hm.trans hw
  | abStop _ hm => simpa [isBisection] using hm.trans hw

/-- Each bisection halves the width; all other steps only decrease it. -/
theorem bisection_count_bound {B : ℝ} {M : ℕ} (r : Run B M) (n : ℕ) :
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

theorem dyadic_tendsto (S : ℝ) : Tendsto (fun n:ℕ=>S/(2:ℝ)^n) atTop (𝓝 0) := by
  have h:= (tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num : (0:ℝ)≤1/2)
    (by norm_num : (1/2:ℝ)<1)).const_mul S
  simpa only [mul_zero,div_pow,one_pow,div_eq_mul_inv,one_mul,inv_pow] using h

/-- Case 1 of the manuscript's qualitative proof. -/
theorem infinite_bisections {B : ℝ} {M : ℕ} (r : Run B M)
    (h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  have hc : Tendsto (bisectionCount r.mode) atTop atTop := by
    apply tendsto_atTop_atTop.mpr
    intro k
    obtain ⟨n,hn⟩:=count_unbounded h k
    exact ⟨n,fun j hj=>hn.trans (count_monotone r.mode hj)⟩
  exact squeeze_zero (run_nonnegative r) (bisection_count_bound r) ((dyadic_tendsto _).comp hc)

theorem active_step {B S w w' : ℝ} {M t c : ℕ} {m' : Mode}
    (hs : Step B M w (.ab S t c) w' m') (ha : ActivePhase B M S t c w)
    (hno : isBisection m'=false) :
    ∃c',m'=.ab S (t+1) c' ∧ ActivePhase B M S (t+1) c' w' := by
  cases hs with
  | abWidth h0 hm ht => exact ⟨0,rfl,.width ha h0 hm ht⟩
  | abException h0 hm ht hc => exact ⟨c+1,rfl,.exception ha h0 hm ht hc⟩
  | abStop => simp [isBisection] at hno

theorem active_history {B W0 w : ℝ} {M n : ℕ} {m : Mode}
    (hW0 : 0≤W0) (h : History B M W0 n w m) (hno : isBisection m=false) :
    ∃S t c,m=.ab S t c ∧ 0≤S ∧ ActivePhase B M S t c w := by
  cases history_decomposes h with
  | ready => simp [isBisection] at hno
  | stopped => simp [isBisection] at hno
  | active hp ha => exact ⟨_,_,_,rfl,prefix_nonnegative hW0 hp,ha⟩

/-- In an infinite suffix without bisection, the AB phase cannot end. -/
theorem continuing_suffix {B S : ℝ} {M N t c : ℕ} (r : Run B M)
    (hm : r.mode N=.ab S t c) (ha : ActivePhase B M S t c (r.width N))
    (hno : ∀n,N≤n→isBisection (r.mode n)=false) :
    ∀k,∃c',r.mode (N+k)=.ab S (t+k) c' ∧
      ActivePhase B M S (t+k) c' (r.width (N+k)) := by
  intro k
  induction k with
  | zero => exact ⟨c,by simpa using hm,by simpa using ha⟩
  | succ k ih =>
      obtain ⟨c',hm',ha'⟩:=ih
      have hs:=r.next (N+k)
      rw [hm'] at hs
      obtain ⟨c'',he,hactive⟩:=active_step hs ha' (hno (N+k+1) (by omega))
      exact ⟨c'',by simpa [Nat.add_assoc] using he,by simpa [Nat.add_assoc] using hactive⟩

/-- Case 2: extract the remaining phase from the history and apply its envelope. -/
theorem finite_bisections {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M)
    (h : ¬∀N,∃n,N≤n ∧ isBisection (r.mode n)=true) :
    Tendsto r.width atTop (𝓝 0) := by
  push_neg at h
  obtain ⟨N,hN⟩:=h
  have hno : ∀n,N≤n→isBisection (r.mode n)=false := by
    intro n hn
    exact Bool.eq_false_iff.mpr (hN n hn)
  obtain ⟨S,t,c,hm,hS,ha⟩:=active_history r.nonneg_initial (history_of_run r N) (hno N (le_refl N))
  have hs:=continuing_suffix r hm ha hno
  have hlim : Tendsto (fun k=>r.width (N+k)) atTop (𝓝 0) := by
    apply squeeze_zero (fun k=>run_nonnegative r (N+k))
      (fun k=>?_) (dyadic_tendsto (envelopeConstant B M*S/(2:ℝ)^t))
    obtain ⟨c',hm',ha'⟩:=hs k
    have hb:=active_envelope hB hS ha'
    simpa only [pow_add,div_mul_eq_div_div] using hb
  exact (tendsto_add_atTop_iff_nat N).mp (by simpa [Nat.add_comm] using hlim)

/-- `thm:global_simple`, width part, with precisely the paper's two cases. -/
theorem width_tendsto_zero {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M) :
    Tendsto r.width atTop (𝓝 0) := by
  by_cases h : ∀N,∃n,N≤n ∧ isBisection (r.mode n)=true
  · exact infinite_bisections r h
  · exact finite_bisections hB r h

/-- Independent consistency check using the quantitative theorem. -/
theorem width_tendsto_zero_quantitative {B : ℝ} {M : ℕ} (hB : (1/2:ℝ)≤B) (r : Run B M) :
    Tendsto r.width atTop (𝓝 0) := by
  let d:=dyadicExponent (envelopeConstant B M)
  have he : Tendsto (fun n:ℕ=>n/(d+2)) atTop atTop:=Nat.tendsto_div_const_atTop (by omega)
  exact squeeze_zero (run_nonnegative r)
    (fun n=>global_width_bound_log hB r.nonneg_initial (history_of_run r n))
    ((dyadic_tendsto _).comp he)

end
end ModAB.Convergence
