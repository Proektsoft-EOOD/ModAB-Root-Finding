import ModAB.TopEntry

/-! Independent exact-arithmetic specification of the rational proposal.
The only algorithmic input taken from the new description is its normalized
candidate and secant-relative limit. These results establish algebra and
admissibility; they assert neither an asymptotic order nor an uncapped width
schedule. All bracket decisions must still use true residuals. -/
namespace ModAB.Rational
open Set
noncomputable section

/-- Secant through working ordinates, which preserve the true endpoint signs. -/
def secant (a b ya yb : ℝ) : ℝ := (a*yb-b*ya)/(yb-ya)

/-- Distance of the secant to its nearer endpoint. -/
def radius (a b ya yb : ℝ) : ℝ :=
  min (secant a b ya yb-a) (b-secant a b ya yb)

/-- The description's normalized denominator, q=ya/yb, qp=fp/yb. -/
def normalizedDenominator (w u q qp : ℝ) : ℝ :=
  w*(qp-q)-u*qp*(1-q)

def normalizedNumerator (w u q qp : ℝ) : ℝ := w*u*q*(qp-1)

def normalizedOffset (w u q qp : ℝ) : ℝ :=
  normalizedNumerator w u q qp / normalizedDenominator w u q qp

def denominator (a b ya yb xp fp : ℝ) : ℝ :=
  normalizedDenominator (b-a) (xp-a) (ya/yb) (fp/yb)

def rawPoint (a b ya yb xp fp : ℝ) : ℝ :=
  a+normalizedOffset (b-a) (xp-a) (ya/yb) (fp/yb)

/-- Real-valued projection; machine NaN/infinity handling is a separate obligation. -/
def limited (a b ya yb raw : ℝ) : ℝ :=
  let s := secant a b ya yb
  let d := radius a b ya yb
  max (s-(7/10:ℝ)*d) (min raw (s+(7/10:ℝ)*d))

/-- Guarded exact-real completion, not the literal exceptional path of C#.
Curvature may disable the proposal. A zero rational denominator falls back
to the secant here; C# instead uses IEEE infinity/NaN and its right-endpoint
fallback. Every admitted finite rational value is projected. -/
def candidate (a b ya yb xp fp : ℝ) (allowRational : Bool) : ℝ :=
  if allowRational = true ∧ denominator a b ya yb xp fp ≠ 0
  then limited a b ya yb (rawPoint a b ya yb xp fp)
  else secant a b ya yb

theorem secant_interior {a b ya yb : ℝ} (hab : a<b) (hy : ya*yb<0) :
    secant a b ya yb ∈ Ioo a b :=
  ModAB.Convergence.secant_interior hab hy

theorem working_ratio_negative {ya yb : ℝ} (hy : ya*yb<0) : ya/yb<0 := by
  rcases mul_neg_iff.mp hy with h|h
  · exact div_neg_of_pos_of_neg h.1 h.2
  · exact div_neg_of_neg_of_pos h.1 h.2

theorem working_second_ne_zero {ya yb : ℝ} (hy : ya*yb<0) : yb≠0 := by
  intro h
  simp [h] at hy

theorem secant_normalized {a b ya yb : ℝ} (hy : ya*yb<0) :
    secant a b ya yb = a-(b-a)*(ya/yb)/(1-ya/yb) := by
  have hv := working_second_ne_zero hy
  have hq := working_ratio_negative hy
  have hc : 1-ya/yb≠0 := by linarith
  have hd : yb-ya≠0 := by
    intro h
    have he : yb=ya := sub_eq_zero.mp h
    rw [he] at hy
    nlinarith [sq_nonneg ya]
  unfold secant
  field_simp [hv,hc,hd]
  ring

theorem radius_pos {a b ya yb : ℝ} (hab : a<b) (hy : ya*yb<0) :
    0<radius a b ya yb := by
  have hs := secant_interior hab hy
  unfold radius
  exact lt_min (sub_pos.mpr hs.1) (sub_pos.mpr hs.2)

theorem radius_le_left (a b ya yb : ℝ) :
    radius a b ya yb ≤ secant a b ya yb-a := min_le_left _ _

theorem radius_le_right (a b ya yb : ℝ) :
    radius a b ya yb ≤ b-secant a b ya yb := min_le_right _ _

theorem limited_bounds {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    secant a b ya yb-(7/10:ℝ)*radius a b ya yb ≤ limited a b ya yb raw ∧
    limited a b ya yb raw ≤ secant a b ya yb+(7/10:ℝ)*radius a b ya yb := by
  have hd := radius_pos hab hy
  unfold limited
  dsimp
  exact ⟨le_max_left _ _,max_le (by nlinarith) (min_le_right _ _)⟩

theorem limited_margins {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    a+(3/10:ℝ)*radius a b ya yb ≤ limited a b ya yb raw ∧
    limited a b ya yb raw ≤ b-(3/10:ℝ)*radius a b ya yb := by
  have hb := limited_bounds (raw:=raw) hab hy
  have hl := radius_le_left a b ya yb
  have hr := radius_le_right a b ya yb
  constructor <;> nlinarith [hb.1,hb.2]

theorem limited_interior {a b ya yb raw : ℝ} (hab : a<b) (hy : ya*yb<0) :
    limited a b ya yb raw ∈ Ioo a b := by
  have hd := radius_pos hab hy
  have hm := limited_margins (raw:=raw) hab hy
  constructor <;> nlinarith [hm.1,hm.2]

theorem limited_eq_self {a b ya yb raw : ℝ}
    (hl : secant a b ya yb-(7/10:ℝ)*radius a b ya yb ≤ raw)
    (hr : raw ≤ secant a b ya yb+(7/10:ℝ)*radius a b ya yb) :
    limited a b ya yb raw=raw := by
  unfold limited
  dsimp
  rw [min_eq_left hr,max_eq_right hl]

theorem candidate_interior {a b ya yb xp fp : ℝ} (allowRational : Bool)
    (hab : a<b) (hy : ya*yb<0) :
    candidate a b ya yb xp fp allowRational ∈ Ioo a b := by
  unfold candidate
  split_ifs
  · exact limited_interior hab hy
  · exact secant_interior hab hy

theorem candidate_fallback_denominator {a b ya yb xp fp : ℝ}
    (allowRational : Bool) (h : denominator a b ya yb xp fp=0) :
    candidate a b ya yb xp fp allowRational=secant a b ya yb := by
  simp [candidate,h]

theorem candidate_fallback_curvature (a b ya yb xp fp : ℝ) :
    candidate a b ya yb xp fp false=secant a b ya yb := by
  simp [candidate]

/-- A conservative exact-real proposal relation for the source branches.
The projected cases cover finite proposals and the two limiting projections of
infinities; the right endpoint admits the source's NaN/rounding fallback.
This is not a theorem about all IEEE operations in the implementation. -/
def SourceProposal (a b ya yb x : ℝ) : Prop :=
  x=secant a b ya yb ∨ (∃raw,x=limited a b ya yb raw) ∨ x=b

theorem sourceProposal_mem_closed {a b ya yb x : ℝ}
    (hab : a<b) (hy : ya*yb<0) (h : SourceProposal a b ya yb x) :
    x ∈ Icc a b := by
  rcases h with rfl|⟨raw,rfl⟩|rfl
  · have hs := secant_interior hab hy
    exact ⟨hs.1.le,hs.2.le⟩
  · have hs := limited_interior (raw:=raw) hab hy
    exact ⟨hs.1.le,hs.2.le⟩
  · exact ⟨hab.le,le_rfl⟩

theorem sourceProposal_interior_of_ne_right {a b ya yb x : ℝ}
    (hab : a<b) (hy : ya*yb<0) (h : SourceProposal a b ya yb x)
    (hne : x≠b) : x ∈ Ioo a b := by
  rcases h with rfl|⟨raw,rfl⟩|rfl
  · exact secant_interior hab hy
  · exact limited_interior hab hy
  · exact False.elim (hne rfl)

theorem candidate_is_sourceProposal (a b ya yb xp fp : ℝ) (allowRational : Bool) :
    SourceProposal a b ya yb (candidate a b ya yb xp fp allowRational) := by
  unfold candidate
  split_ifs
  · exact Or.inr (Or.inl ⟨_,rfl⟩)
  · exact Or.inl rfl

/-- Source-order closing operation, applied AFTER the secant-relative clamp.
The supplied tolerance is computed from the pre-closing point. The source's
second guard prevents a displacement past the opposite endpoint. -/
def closeCandidate (a b tolerance x : ℝ) (mostRecent : ModAB.Convergence.EndpointSide) : ℝ :=
  match mostRecent with
  | .left => if x-a<tolerance/2 ∧ a+tolerance/2<b then a+tolerance/2 else x
  | .right => if b-x<tolerance/2 ∧ a<b-tolerance/2 then b-tolerance/2 else x

theorem closeCandidate_interior {a b tolerance x : ℝ}
    (mostRecent : ModAB.Convergence.EndpointSide)
    (ht : 0<tolerance) (hx : x∈Ioo a b) :
    closeCandidate a b tolerance x mostRecent ∈ Ioo a b := by
  cases mostRecent <;> dsimp only [closeCandidate] <;> split_ifs with h
  · exact ⟨by linarith,h.2⟩
  · exact hx
  · exact ⟨h.2,by linarith⟩
  · exact hx

theorem closeCandidate_mem_closed {a b tolerance x : ℝ}
    (mostRecent : ModAB.Convergence.EndpointSide)
    (ht : 0<tolerance) (hx : x∈Icc a b) :
    closeCandidate a b tolerance x mostRecent ∈ Icc a b := by
  cases mostRecent <;> dsimp only [closeCandidate] <;> split_ifs with h
  · exact ⟨by linarith,h.2.le⟩
  · exact hx
  · exact ⟨h.2.le,by linarith⟩
  · exact hx

theorem closeCandidate_left_separation {a b tolerance x : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    a+tolerance/2 ≤ closeCandidate a b tolerance x .left := by
  have hg : a+tolerance/2<b := by linarith
  simp only [closeCandidate,hg,and_true]
  split_ifs with h <;> linarith

theorem closeCandidate_right_separation {a b tolerance x : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    closeCandidate a b tolerance x .right ≤ b-tolerance/2 := by
  have hg : a<b-tolerance/2 := by linarith
  simp only [closeCandidate,hg,and_true]
  split_ifs with h <;> linarith

theorem closeCandidate_right_fallback_interior {a b tolerance : ℝ}
    (ht : 0<tolerance) (hw : tolerance<b-a) :
    closeCandidate a b tolerance b .right ∈ Ioo a b := by
  have hg : a<b-tolerance/2 := by linarith
  have hz : b-b<tolerance/2 := by linarith
  change (if b-b<tolerance/2 ∧ a<b-tolerance/2 then b-tolerance/2 else b) ∈ Ioo a b
  rw [ite_eq_left ⟨hz,hg⟩]
  exact ⟨hg,by linarith⟩

/-- With an outside history node and opposite endpoint ordinates, the normalized
numerator and denominator cannot vanish together in exact arithmetic. This says
nothing about NaN caused by missing history or floating-point overflow. -/
theorem numerator_ne_zero_of_denominator_zero {w u q qp : ℝ}
    (hw : 0<w) (hu : u<0 ∨ w<u) (hq : q<0)
    (hD : normalizedDenominator w u q qp=0) :
    normalizedNumerator w u q qp≠0 := by
  have hwn : w≠0 := ne_of_gt hw
  have hun : u≠0 := by rcases hu with hu|hu <;> linarith
  have hqn : q≠0 := ne_of_lt hq
  intro hn
  unfold normalizedNumerator at hn
  have hqp : qp-1=0 := (mul_eq_zero.mp hn).resolve_left
    (mul_ne_zero (mul_ne_zero hwn hun) hqn)
  have he : qp=1 := by linarith
  have hwu : w-u≠0 := by rcases hu with hu|hu <;> linarith
  have hc : 1-q≠0 := by linarith
  have hid : normalizedDenominator w u q qp=(w-u)*(1-q) := by
    rw [he]
    unfold normalizedDenominator
    ring
  rw [hid] at hD
  exact (mul_ne_zero hwu hc) hD

/-- Independent coefficient construction: the normalized proposal satisfies
all three cross-multiplied interpolation equations. Genuine function-value
interpolation additionally requires nonzero model denominators at the nodes. -/
theorem normalized_interpolation_equations {w u q qp : ℝ}
    (hw : w≠0) (hq : q≠0) (hD : normalizedDenominator w u q qp≠0) :
    let t := normalizedOffset w u q qp
    let p := 1+t*(1-q)/(q*w)
    let Q := -t/q
    (-t=q*Q ∧ w-t=p*w+Q ∧ u-t=qp*(p*u+Q)) := by
  dsimp [normalizedOffset,normalizedNumerator,normalizedDenominator]
  unfold normalizedDenominator at hD
  constructor
  · field_simp [hw,hq,hD]
  constructor
  · field_simp [hw,hq,hD]
    ring
  · field_simp [hw,hq,hD]
    ring

/-- A history record deliberately keeps discarded TRUE residuals. -/
structure HistoryState where
  left : ℝ
  right : ℝ
  ordinates : ModAB.Convergence.OrdinateState
  dropped : Option (ℝ×ℝ)

/-- Genuine endpoint replacement: this helper records the removed endpoint.
For source executions that reuse the same coordinate, use recordSourceHistory. -/
def recordHistory (s : HistoryState) (moved : ModAB.Convergence.EndpointSide)
    (x y : ℝ) : HistoryState :=
  match moved with
  | .left => ⟨x,s.right,ModAB.Convergence.recordUpdate s.ordinates moved y,
      some (s.left,s.ordinates.leftTrue)⟩
  | .right => ⟨s.left,x,ModAB.Convergence.recordUpdate s.ordinates moved y,
      some (s.right,s.ordinates.rightTrue)⟩

theorem recordHistory_true_values {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    (recordHistory s moved x (f x)).ordinates.leftTrue =
      f (recordHistory s moved x (f x)).left ∧
    (recordHistory s moved x (f x)).ordinates.rightTrue =
      f (recordHistory s moved x (f x)).right := by
  cases moved <;> simp [recordHistory,ModAB.Convergence.recordUpdate,hl,hr]

theorem recordHistory_dropped_true {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    ∃xp,(recordHistory s moved x (f x)).dropped=some (xp,f xp) := by
  cases moved with
  | left => exact ⟨s.left,by simp [recordHistory,hl]⟩
  | right => exact ⟨s.right,by simp [recordHistory,hr]⟩

/-- Truth of an optional history record: an absent history is allowed. -/
def TrueHistory (f : ℝ→ℝ) (history : Option (ℝ×ℝ)) : Prop :=
  ∀xp yp,history=some (xp,yp) → yp=f xp

/-- Source-style history guard. Reusing an endpoint preserves the previous
history, but it does not skip the AB ordinate update. In the source a repeated
same-side reuse can still apply an AB correction before this guard. This local
helper uses the old AB-side update and does not encode the source's independent
lastMoved/sideMoved plateau/reset controller. -/
def recordSourceHistory (s : HistoryState) (moved : ModAB.Convergence.EndpointSide)
    (x y : ℝ) : HistoryState :=
  match moved with
  | .left => ⟨x,s.right,ModAB.Convergence.recordUpdate s.ordinates moved y,
      if x=s.left then s.dropped else some (s.left,s.ordinates.leftTrue)⟩
  | .right => ⟨s.left,x,ModAB.Convergence.recordUpdate s.ordinates moved y,
      if x=s.right then s.dropped else some (s.right,s.ordinates.rightTrue)⟩

theorem recordSourceHistory_true_values {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right) :
    (recordSourceHistory s moved x (f x)).ordinates.leftTrue =
      f (recordSourceHistory s moved x (f x)).left ∧
    (recordSourceHistory s moved x (f x)).ordinates.rightTrue =
      f (recordSourceHistory s moved x (f x)).right := by
  cases moved <;> simp [recordSourceHistory,ModAB.Convergence.recordUpdate,hl,hr]

theorem recordSourceHistory_history_true {f : ℝ→ℝ} (s : HistoryState)
    (moved : ModAB.Convergence.EndpointSide) (x : ℝ)
    (hl : s.ordinates.leftTrue=f s.left) (hr : s.ordinates.rightTrue=f s.right)
    (hd : TrueHistory f s.dropped) :
    TrueHistory f (recordSourceHistory s moved x (f x)).dropped := by
  intro xp yp hh
  cases moved with
  | left =>
      by_cases hx : x=s.left
      · exact hd xp yp (by simpa [recordSourceHistory,hx] using hh)
      · have he : s.left=xp ∧ s.ordinates.leftTrue=yp := by
          simpa [recordSourceHistory,hx] using hh
        calc yp = s.ordinates.leftTrue := he.2.symm
             _ = f s.left := hl
             _ = f xp := congrArg f he.1
  | right =>
      by_cases hx : x=s.right
      · exact hd xp yp (by simpa [recordSourceHistory,hx] using hh)
      · have he : s.right=xp ∧ s.ordinates.rightTrue=yp := by
          simpa [recordSourceHistory,hx] using hh
        calc yp = s.ordinates.rightTrue := he.2.symm
             _ = f s.right := hr
             _ = f xp := congrArg f he.1

theorem recordSourceHistory_reuse_left (s : HistoryState) (y : ℝ) :
    (recordSourceHistory s .left s.left y).dropped=s.dropped := by
  simp [recordSourceHistory]

theorem recordSourceHistory_reuse_right (s : HistoryState) (y : ℝ) :
    (recordSourceHistory s .right s.right y).dropped=s.dropped := by
  simp [recordSourceHistory]

theorem recordSourceHistory_left_eq_recordHistory {s : HistoryState} {x y : ℝ}
    (hx : x≠s.left) :
    recordSourceHistory s .left x y=recordHistory s .left x y := by
  simp [recordSourceHistory,recordHistory,hx]

theorem recordSourceHistory_right_eq_recordHistory {s : HistoryState} {x y : ℝ}
    (hx : x≠s.right) :
    recordSourceHistory s .right x y=recordHistory s .right x y := by
  simp [recordSourceHistory,recordHistory,hx]

/-- The endpoint recorded as most recently evaluated has its true ordinate. -/
def Fresh (s : ModAB.Convergence.OrdinateState) : Prop :=
  match s.previous with
  | some .left => s.leftWork=s.leftTrue
  | some .right => s.rightWork=s.rightTrue
  | none => s.leftWork=s.leftTrue ∧ s.rightWork=s.rightTrue

theorem recordUpdate_fresh (s : ModAB.Convergence.OrdinateState)
    (moved : ModAB.Convergence.EndpointSide) (y : ℝ) :
    Fresh (ModAB.Convergence.recordUpdate s moved y) := by
  cases moved <;> simp [Fresh,ModAB.Convergence.recordUpdate]

theorem recordUpdate_at_most_one_modified (s : ModAB.Convergence.OrdinateState)
    (moved : ModAB.Convergence.EndpointSide) (y : ℝ) :
    (ModAB.Convergence.recordUpdate s moved y).leftWork=
      (ModAB.Convergence.recordUpdate s moved y).leftTrue ∨
    (ModAB.Convergence.recordUpdate s moved y).rightWork=
      (ModAB.Convergence.recordUpdate s moved y).rightTrue := by
  cases moved <;> simp [ModAB.Convergence.recordUpdate]

/-- When a correction is applied on the same side, the old implementation's
working denominator equals the true replaced residual by the freshness invariant. -/
theorem correction_true_denominator (s : ModAB.Convergence.OrdinateState)
    (h : Fresh s) (y : ℝ) :
    (s.previous=some .left → ModAB.Convergence.abFactor y s.leftWork=
      ModAB.Convergence.abFactor y s.leftTrue) ∧
    (s.previous=some .right → ModAB.Convergence.abFactor y s.rightWork=
      ModAB.Convergence.abFactor y s.rightTrue) := by
  constructor
  · intro hp
    have he : s.leftWork=s.leftTrue := by simpa [Fresh,hp] using h
    rw [he]
  · intro hp
    have he : s.rightWork=s.rightTrue := by simpa [Fresh,hp] using h
    rw [he]

/-- Same-sign replacement residuals make the AB multiplier a strict contraction;
together with the old positivity theorem this yields a factor in (0,1). -/
theorem abFactor_lt_one {new old : ℝ} (h : 0<new*old) :
    ModAB.Convergence.abFactor new old<1 := by
  have hr : 0<new/old := by
    rcases mul_pos_iff.mp h with h|h
    · exact div_pos h.1 h.2
    · exact div_pos_of_neg_of_neg h.1 h.2
  unfold ModAB.Convergence.abFactor
  split_ifs <;> linarith

theorem abFactor_mem_Ioo {new old : ℝ} (h : 0<new*old) :
    ModAB.Convergence.abFactor new old ∈ Ioo (0:ℝ) 1 :=
  ⟨ModAB.Convergence.abFactor_positive new old,abFactor_lt_one h⟩

/-- Mixed residuals do not inherit exactness even for f(x)=x: the left ordinate
is multiplied by 1/2, the right ordinate and the dropped value remain true.
The computed point is not the root 0. This is a data-level example, not a claim
that a full run must encounter this state. -/
theorem scaled_linear_candidate_example :
    candidate (-1) 1 (-(1/2:ℝ)) 1 2 2 true=-(1/4:ℝ) := by
  norm_num [candidate,denominator,rawPoint,normalizedOffset,
    normalizedNumerator,normalizedDenominator,limited,radius,secant]

end
end ModAB.Rational
