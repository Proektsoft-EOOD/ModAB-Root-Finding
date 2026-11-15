import ModAB.RationalCycleNeighborhood
import ModAB.RationalCycleStateGuards
import ModAB.RationalCycleTransport

/-! Complete negative-coefficient local theorem for the guarded source-state
iteration. The raw coefficient, AB factor, both accepted trials and both
continuation decisions are derived from local C³ data. -/
namespace ModAB.Rational.CycleComplete
open Set Filter
open scoped Topology
open LocalState LocalComplete CycleUniform CycleNeighborhood CycleStateGuards
open CycleIntermediate CycleTransport CycleStateUpdate
noncomputable section

theorem endpoint_min (s : State) (r : ℝ) :
    min |latest s-r| |opposite s-r| = min |s.a-r| |s.b-r| := by
  cases hs : s.recent <;> simp [latest,opposite,hs,min_comm]

theorem valid_latest_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |latest s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [latest,hs]
  · exact h.near_left
  · exact h.near_right

theorem valid_opposite_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |opposite s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [opposite,hs]
  · exact h.near_right
  · exact h.near_left

theorem valid_history_near {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) : |selectedHistory s-r| ≤ ε := by
  cases hs : s.recent <;> simp only [selectedHistory,hs]
  · exact h.near_left_history
  · exact h.near_right_history

theorem valid_working {f : ℝ → ℝ} {r ε : ℝ} {s : State}
    (h : Valid f r ε s) :
    latestWork s = f (latest s) ∧ oppositeWork s = f (opposite s) := by
  cases hs : s.recent <;> simp [latestWork,oppositeWork,latest,opposite,hs,h.fresh,freshOrdinates]

/-- One fresh local cone state executes a genuine same-side update followed
by a genuine crossing, including the two actual safeguarded proposals. -/
theorem pair_law_of_neighborhood {f : ℝ → ℝ} {r θ ε C δ : ℝ}
    (N : PairNeighborhood f r θ ε C δ) :
    PairLaw f (Regularity.regularizedReciprocal f r) r θ ε (allowRational f) := by
  refine ⟨?_⟩
  intro st hs hc
  let g := Regularity.regularizedReciprocal f r
  let a := latest st-r
  let s := opposite st-r
  let d := selectedHistory st-r
  let b := firstError g r a s d
  let v := secondCoefficient g r a s d
  let c := v*a*b*d
  let x := r+b
  let y := r+c
  let μ := 1-f x/f (latest st)
  have hcone : CycleEntry.Cone θ ε a s d := hc
  have hh := N.raw_pair a s d hcone
  dsimp only at hh
  have hea : r+a = latest st := by dsimp [a]; ring
  have hes : r+s = opposite st := by dsimp [s]; ring
  have hed : r+d = selectedHistory st := by dsimp [d]; ring
  obtain ⟨hU,hV,hD₁,hR₁,hμ,hD₂,hR₂,hb0,hba,hbs,hab⟩ := hh
  have hD₁' : denominator (latest st) (opposite st) (f (latest st))
      (f (opposite st)) (selectedHistory st) (f (selectedHistory st)) ≠ 0 := by
    simpa only [hea,hes,hed] using hD₁
  have hR₁' : rawPoint (latest st) (opposite st) (f (latest st))
      (f (opposite st)) (selectedHistory st) (f (selectedHistory st))-r = b := by
    simpa only [hea,hes,hed] using hR₁
  have hμ' : μ ∈ Ioo (0:ℝ) 1 := by
    simpa only [hea] using hμ
  have hD₂' : denominator x (opposite st) (f x) (μ*f (opposite st))
      (latest st) (f (latest st)) ≠ 0 := by
    simpa only [hea,hes] using hD₂
  have hR₂' : rawPoint x (opposite st) (f x) (μ*f (opposite st))
      (latest st) (f (latest st))-r = c := by
    simpa only [x,μ,c,v,b,g,hea,hes,mul_assoc] using hR₂
  obtain ⟨hu,hv,huC,hvC⟩ := N.bounds.coefficients a s d hcone
  have hq0 : 0 ≤ C*ε^2 := mul_nonneg N.bounds.bound_pos.le (sq_nonneg ε)
  have hq1 : C*ε^2 < 1 := N.guards.small.trans (by norm_num)
  have haε : |a| ≤ ε := valid_latest_near hs
  have hsε : |s| ≤ ε := valid_opposite_near hs
  have hdε : |d| ≤ ε := valid_history_near hs
  have hbprod : |b| ≤ C*|a| *|s| *|d| := by
    dsimp [b,firstError]
    simp only [abs_mul]
    gcongr
  have hbmin : |b| ≤ (C*ε^2)*min |a| |s| :=
    ModAB.RationalLocalSafeguards.triple_error_bound N.bounds.bound_pos.le
      N.bounds.radius_pos.le haε hsε hdε hbprod
  have hraw₁ : |rawPoint (latest st) (opposite st) (f (latest st))
      ((1:ℝ)*f (opposite st)) (selectedHistory st) (f (selectedHistory st))-r| ≤
      (C*ε^2)*min |latest st-r| |opposite st-r| := by
    simpa only [one_mul,hR₁'] using hbmin
  have hw := valid_working hs
  have guard₁ := guards_latest N.guards (Geometry.of_valid hs) hw.1
    (by simpa using hw.2) (by norm_num : (1:ℝ) ≠ 0)
    (by simpa using N.guards.factor_pos.le)
    (by simpa using hD₁') hraw₁
  have htrial₁ : trial f (allowRational f) st = x := by
    have hh := guard₁.1
    rw [one_mul] at hh
    dsimp only [x]
    linarith only [hR₁',hh]
  have hbminOrd : |x-r| ≤ (C*ε^2)*min |st.a-r| |st.b-r| := by
    simpa only [x,add_sub_cancel_left,←endpoint_min st r] using hbmin
  have hxi := ModAB.RationalLocalSafeguards.interior_of_local_error hs.left hs.right hq1 hbminOrd
  have hxn := interior_near hs.left hs.right hs.near_left hs.near_right hxi
  have hsame : 0 < (latest st-r)*(x-r) := by
    simpa only [x,add_sub_cancel_left] using hab
  have hi := same_side_update hs hxi hsame N.same_sign
  let mid := update f st x
  have hstep₁ : step f (allowRational f) st = mid := by
    simp only [step,htrial₁]; rfl
  have hgeom := Geometry.of_same_side hs hxi hsame N.same_sign
  have hmidlatest : latest mid = x := hi.2.1
  have hmidopp : opposite mid = opposite st := hi.2.2.1
  have hmidhist : selectedHistory mid = latest st := hi.2.2.2.1
  have hbLatest : |x-r| ≤ (C*ε^2)*|latest st-r| := by
    exact (show |x-r| ≤ (C*ε^2)*min |latest st-r| |opposite st-r| by
      simpa only [x,add_sub_cancel_left] using hbmin).trans
        (mul_le_mul_of_nonneg_left (min_le_left _ _) hq0)
  have habfactor := N.guards.ab_factor (Alternating.latest_ne_root hs)
    (valid_latest_near hs) hxn hbLatest
  have hmidwork : latestWork mid = f (latest mid) := by rw [hmidlatest]; exact hi.2.2.2.2.1
  have hmidoppwork : oppositeWork mid = μ*f (opposite mid) := by
    rw [hmidopp,hi.2.2.2.2.2.1,habfactor.1]
  have hμnear : |μ-1| ≤ δ := by simpa only [habfactor.1] using habfactor.2
  have hb_le_min : |b| ≤ min |a| |s| := hbmin.trans
    (mul_le_of_le_one_left (le_min (abs_nonneg a) (abs_nonneg s)) hq1.le)
  have hb_le_s : |b| ≤ |s| := hb_le_min.trans (min_le_right _ _)
  have hcBound : |c| ≤ (C*ε^2)*|b| := by
    exact CycleEntry.product_bound N.bounds.bound_pos.le N.bounds.radius_pos.le hvC haε hdε
  have hraw₂ : |rawPoint (latest mid) (opposite mid) (f (latest mid))
      (μ*f (opposite mid)) (selectedHistory mid) (f (selectedHistory mid))-r| ≤
      (C*ε^2)*min |latest mid-r| |opposite mid-r| := by
    rw [hmidlatest,hmidopp,hmidhist,hR₂']
    have hh : |b| ≤ |opposite st-r| := hb_le_s
    simpa only [x,add_sub_cancel_left,min_eq_left hh] using hcBound
  have hden₂ : denominator (latest mid) (opposite mid) (f (latest mid))
      (μ*f (opposite mid)) (selectedHistory mid) (f (selectedHistory mid)) ≠ 0 := by
    simpa only [hmidlatest,hmidopp,hmidhist] using hD₂'
  have guard₂ := guards_latest N.guards hgeom hmidwork hmidoppwork
    hμ'.1.ne' hμnear hden₂ hraw₂
  have htrial₂ : trial f (allowRational f) mid = y := by
    have hh := guard₂.1
    rw [hmidlatest,hmidopp,hmidhist] at hh
    dsimp only [y]
    linarith only [hR₂',hh]
  have hpair : pairStep f (allowRational f) st = update f mid y := by
    unfold pairStep
    rw [hstep₁]
    change update f mid (trial f (allowRational f) mid) = _
    rw [htrial₂]
  have hcyOrd : |y-r| ≤ (C*ε^2)*min |st.a-r| |st.b-r| := by
    have hcmin := hcBound.trans (mul_le_mul_of_nonneg_left hb_le_min hq0)
    simpa only [y,add_sub_cancel_left,←endpoint_min st r] using hcmin
  have hyi := ModAB.RationalLocalSafeguards.interior_of_local_error hs.left hs.right hq1 hcyOrd
  have hcross : (latest st-r)*(y-r) < 0 := by
    have hneg := mul_neg_of_neg_of_pos (mul_neg_of_neg_of_pos hv hcone.same) hab
    change a*((r+c)-r) < 0
    dsimp only [c]
    convert hneg using 1; ring
  have hout := pair_update_preserves hs hxi hyi hsame hcross N.same_sign
  rw [hpair]
  refine ⟨hout.1,?_,?_,?_⟩
  · dsimp only [view,CycleState.nextPair]
    rw [hout.2.1,hout.2.2.1,hout.2.2.2]
    simp only [y,x,add_sub_cancel_left]
    rfl
  · rw [hstep₁,hmidlatest]
    dsimp [x,CycleState.nextPair,view,a,s,d,b,g]
    ring
  · intro B
    refine ⟨guard₁.2 B,?_⟩
    rw [hstep₁]
    exact guard₂.2 B


/-- Propositions 10--11 in operational form. Finite C³ smoothness and the
initial geometric cone choose one common neighbourhood. Every later pair
executes both actual guarded source steps, stays in the cone and has the
proved two distinct consecutive-step orders. No future-state, guard,
recurrence, coefficient-limit or convergence premise is present. -/
theorem local_pair_complete {f : ℝ → ℝ} {r : ℝ}
    (hf : ContDiffAt ℝ 3 f r) (hz : f r = 0) (hd : deriv f r ≠ 0)
    (h₂ : deriv (deriv f) r ≠ 0) (hK : AnalyticBridge.errorConstant f r < 0) :
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ st : State, Valid f r ε st →
        CycleEntry.Cone θ ε (view r st).a (view r st).s (view r st).d →
        (∀ n, Valid f r ε (Alternating.orbit f (allowRational f) st (2*n)) ∧
          CycleEntry.Cone θ ε
            (view r (Alternating.orbit f (allowRational f) st (2*n))).a
            (view r (Alternating.orbit f (allowRational f) st (2*n))).s
            (view r (Alternating.orbit f (allowRational f) st (2*n))).d) ∧
        (∀ B n,
          nextMode B f (allowRational f) (Alternating.orbit f (allowRational f) st (2*n)) =
            .interpolation (Alternating.orbit f (allowRational f) st (2*n)).scale
              ((Alternating.orbit f (allowRational f) st (2*n)).clock+1) ∧
          nextMode B f (allowRational f) (Alternating.orbit f (allowRational f) st (2*n+1)) =
            .interpolation (Alternating.orbit f (allowRational f) st (2*n+1)).scale
              ((Alternating.orbit f (allowRational f) st (2*n+1)).clock+1)) ∧
        ModAB.RationalOrder.HasQOrder
          (fun n => latest (Alternating.orbit f (allowRational f) st (2*n))-r)
          3 |AnalyticBridge.errorConstant f r| ∧
        ModAB.RationalOrder.HasQOrder
          (fun n => latest (Alternating.orbit f (allowRational f) st (2*n+1))-r)
          3 |AnalyticBridge.errorConstant f r| ∧
        Tendsto (fun n => |latest (Alternating.orbit f (allowRational f) st (2*n+1))-r| /
          |latest (Alternating.orbit f (allowRational f) st (2*n))-r| ^ (9/5:ℝ)) atTop
          (𝓝 (|AnalyticBridge.errorConstant f r| ^ (2/5:ℝ))) ∧
        Tendsto (fun n => |latest (Alternating.orbit f (allowRational f) st (2*(n+1)))-r| /
          |latest (Alternating.orbit f (allowRational f) st (2*n+1))-r| ^ (5/3:ℝ)) atTop
          (𝓝 (|AnalyticBridge.errorConstant f r| ^ (1/3:ℝ))) := by
  obtain ⟨θ,ε,C,δ,N⟩ := exists_pair_neighborhood hf hz hd h₂ hK
  let g := Regularity.regularizedReciprocal f r
  have hg : ContDiffAt ℝ 2 g r :=
    Regularity.regularizedReciprocal_contDiffAt_of_contDiff hf hd
  have hgr : g r ≠ 0 := by
    dsimp only [g]
    rw [Regularity.regularizedReciprocal_at]
    exact one_div_ne_zero hd
  have hk : localK g r = AnalyticBridge.errorConstant f r :=
    localK_regularizedReciprocal hf hz hd
  have L := pair_law_of_neighborhood N
  refine ⟨θ,ε,N.bounds.theta_pos,N.bounds.theta_lt,N.bounds.radius_pos,?_⟩
  intro st hs hc
  refine ⟨?_,?_,?_⟩
  · intro n
    have he := even_orbit N.bounds L hs hc n
    refine ⟨he.1,?_⟩
    rw [he.2]
    exact CycleState.orbit_cone N.bounds hc n
  · exact fun B n => continue_pair N.bounds L hs hc B n
  · have hh := source_hasQOrders hg hgr (by rw [hk]; exact hK.ne)
      N.bounds L hs hc
    simpa only [hk] using hh

end
end ModAB.Rational.CycleComplete
