import ModAB.RationalCycleUniform
import ModAB.RationalCycleUniformAlgebra

/-!
Connect the uniformly controlled coefficients to the actual function values
and the source-normalized raw proposal. Nonvanishing guards follow from the
normalized pole values and regularized reciprocal values.
-/
namespace ModAB.Rational.CycleUniformBridge
open Set
open LocalAlgebra Regularity CycleUniform CycleUniformAlgebra
noncomputable section

theorem reciprocal_value_at_error {f : ℝ → ℝ} {r e : ℝ}
    (hzero : f r = 0) (he : e ≠ 0) :
    e / f (r+e) = regularizedReciprocal f r (r+e) := by
  symm
  simpa only [add_sub_cancel_left] using
    regularizedReciprocal_off_root hzero (show r+e ≠ r by simpa using he)

theorem ordinate_ne_zero_of_reciprocal {f : ℝ → ℝ} {r e : ℝ}
    (hzero : f r = 0) (he : e ≠ 0)
    (hg : regularizedReciprocal f r (r+e) ≠ 0) : f (r+e) ≠ 0 := by
  intro hz
  apply hg
  rw [←reciprocal_value_at_error hzero he,hz,div_zero]

/-- Uniform control of the first pole and reciprocal values supplies the
actual raw-step denominator guard. -/
theorem first_source_denominator_ne_zero {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (has : a ≠ s) (hsd : s ≠ d) (had : a ≠ d)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgd : regularizedReciprocal f r (r+d) ≠ 0)
    (hP : firstPole (regularizedReciprocal f r) r a s d ≠ 0) :
    ModAB.Rational.denominator (r+a) (r+s) (f (r+a)) (f (r+s))
      (r+d) (f (r+d)) ≠ 0 := by
  have hsa : (r+s)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr has.symm
  have hds : (r+d)-(r+s) ≠ 0 := by simpa using sub_ne_zero.mpr hsd.symm
  have hda : (r+d)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr had.symm
  have hp := pole_identity_from_values (r:=r) (x₁:=r+a) (x₂:=r+s) (x₃:=r+d)
    (y₁:=f (r+a)) (y₂:=f (r+s)) (y₃:=f (r+d)) hsa hds hda
    (by simpa using ha) (by simpa using hs) (by simpa using hd)
  simp only [add_sub_cancel_left] at hp
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hs,
    reciprocal_value_at_error hzero hd] at hp
  have hH : dd2 (r+a) (r+s) (r+d) (1/f (r+a)) (1/f (r+s)) (1/f (r+d)) ≠ 0 := by
    intro hz
    apply hP
    change poleValue (r+a) (r+s) (r+d) r
      (regularizedReciprocal f r (r+a)) (regularizedReciprocal f r (r+s))
      (regularizedReciprocal f r (r+d)) = 0
    rw [←hp,hz,zero_mul]
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfd := ordinate_ne_zero_of_reciprocal hzero hd hgd
  intro hz
  apply hH
  rw [bary_denominator_eq hsa hds hda hfa hfs hfd,hz]
  simp

/-- The first coefficient is the coefficient of the implemented raw proposal,
not merely a separate divided-difference model. -/
theorem first_rawPoint_error {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (has : a ≠ s) (hsd : s ≠ d) (had : a ≠ d)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgd : regularizedReciprocal f r (r+d) ≠ 0)
    (hP : firstPole (regularizedReciprocal f r) r a s d ≠ 0) :
    ModAB.Rational.rawPoint (r+a) (r+s) (f (r+a)) (f (r+s))
      (r+d) (f (r+d))-r = firstError (regularizedReciprocal f r) r a s d := by
  have hD := first_source_denominator_ne_zero hzero ha hs hd has hsd had hga hgs hgd hP
  have hsa : (r+s)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr has.symm
  have hds : (r+d)-(r+s) ≠ 0 := by simpa using sub_ne_zero.mpr hsd.symm
  have hda : (r+d)-(r+a) ≠ 0 := by simpa using sub_ne_zero.mpr had.symm
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfd := ordinate_ne_zero_of_reciprocal hzero hd hgd
  rw [rawPoint_eq_baryRoot hsa hds hda hfa hfs hfd hD]
  rw [baryRoot_true_error hsa hds hda (by simpa using ha) (by simpa using hs)
    (by simpa using hd) (bary_denominator_ne_zero hsa hds hda hfa hfs hfd hD)]
  simp only [add_sub_cancel_left]
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hs,
    reciprocal_value_at_error hzero hd]
  unfold firstError firstCoefficient firstPole
  ring

/-- The exact residual ratio used by AB equals the normalized local expression. -/
theorem residual_ratio_eq_ordinateRatio {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0) :
    f (r+firstError (regularizedReciprocal f r) r a s d) / f (r+a) =
      ordinateRatio (regularizedReciprocal f r) r a s d := by
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  have h := residual_ratio_from_step ha hb hfa hfb
    (show firstError (regularizedReciprocal f r) r a s d =
      firstCoefficient (regularizedReciprocal f r) r a s d*a*s*d from rfl)
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb] at h
  exact h

/-- The second normalized pole is exactly the pole in the implementation's
AB-corrected data. -/
theorem second_source_denominator_ne_zero {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hba : firstError (regularizedReciprocal f r) r a s d ≠ a)
    (hbs : firstError (regularizedReciprocal f r) r a s d ≠ s) (has : a ≠ s)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0)
    (ht : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0)
    (hP : secondPole (regularizedReciprocal f r) r a s d ≠ 0) :
    let b := firstError (regularizedReciprocal f r) r a s d
    ModAB.Rational.denominator (r+b) (r+s) (f (r+b))
      ((1-f (r+b)/f (r+a))*f (r+s)) (r+a) (f (r+a)) ≠ 0 := by
  dsimp only
  let g := regularizedReciprocal f r
  let b := firstError g r a s d
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  apply source_denominator_ne_zero_of_mixed_pole ha hb hs
    (sub_ne_zero.mpr hba.symm) (sub_ne_zero.mpr has.symm) (sub_ne_zero.mpr hbs.symm)
    hfa hfb hfs (by rw [hratio]; exact ht)
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb,
    reciprocal_value_at_error hzero hs,hratio]
  exact hP

/-- Actual second raw proposal equals the uniformly controlled coefficient
times a*b*d. The AB factor is formed from the true residuals. -/
theorem second_rawPoint_error {f : ℝ → ℝ} {r a s d : ℝ}
    (hzero : f r = 0) (ha : a ≠ 0) (hs : s ≠ 0) (hd : d ≠ 0)
    (hb : firstError (regularizedReciprocal f r) r a s d ≠ 0)
    (hba : firstError (regularizedReciprocal f r) r a s d ≠ a)
    (hbs : firstError (regularizedReciprocal f r) r a s d ≠ s) (has : a ≠ s)
    (hga : regularizedReciprocal f r (r+a) ≠ 0)
    (hgs : regularizedReciprocal f r (r+s) ≠ 0)
    (hgb : regularizedReciprocal f r
      (r+firstError (regularizedReciprocal f r) r a s d) ≠ 0)
    (ht : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0)
    (hP : secondPole (regularizedReciprocal f r) r a s d ≠ 0) :
    let g := regularizedReciprocal f r
    let b := firstError g r a s d
    ModAB.Rational.rawPoint (r+b) (r+s) (f (r+b))
      ((1-f (r+b)/f (r+a))*f (r+s)) (r+a) (f (r+a))-r =
      secondCoefficient g r a s d*(a*b*d) := by
  dsimp only
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hfa := ordinate_ne_zero_of_reciprocal hzero ha hga
  have hfs := ordinate_ne_zero_of_reciprocal hzero hs hgs
  have hfb := ordinate_ne_zero_of_reciprocal hzero hb hgb
  have hD := second_source_denominator_ne_zero hzero ha hs hb hba hbs has hga hgs hgb ht hP
  have h := source_rawPoint_second_error ha hb hs hd (sub_ne_zero.mpr hba.symm)
    (sub_ne_zero.mpr has.symm) (sub_ne_zero.mpr hbs.symm) hfa hfb hfs
    (by rw [hratio]; exact ht)
    (show firstError (regularizedReciprocal f r) r a s d =
      firstCoefficient (regularizedReciprocal f r) r a s d*a*s*d from rfl) hD
  rw [reciprocal_value_at_error hzero ha,reciprocal_value_at_error hzero hb,
    reciprocal_value_at_error hzero hs] at h
  exact h

/-- Finite-C³ hypotheses alone give a cone on which both actual rational
proposals have bounded negative coefficients, the AB multiplier lies strictly
between zero and one, and both implementation denominator guards hold.
No divided-difference estimate or asymptotic recurrence is a premise. -/
theorem smooth_exists_uniform_raw_pair
    {f : ℝ → ℝ} {r : ℝ} (hf : ContDiffAt ℝ 3 f r)
    (hzero : f r = 0) (hder : deriv f r ≠ 0)
    (hK : AnalyticBridge.errorConstant f r < 0) :
    let g := regularizedReciprocal f r
    let K := AnalyticBridge.errorConstant f r
    ∃ θ ε : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < ε ∧
      ∀ a s d, CycleEntry.Cone θ ε a s d →
        let b := firstError g r a s d
        let μ := 1-f (r+b)/f (r+a)
        firstCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
        secondCoefficient g r a s d ∈ Ioo (2*K) (K/2) ∧
        ModAB.Rational.denominator (r+a) (r+s) (f (r+a)) (f (r+s))
          (r+d) (f (r+d)) ≠ 0 ∧
        ModAB.Rational.rawPoint (r+a) (r+s) (f (r+a)) (f (r+s))
          (r+d) (f (r+d))-r = b ∧
        μ ∈ Ioo (0:ℝ) 1 ∧
        ModAB.Rational.denominator (r+b) (r+s) (f (r+b))
          (μ*f (r+s)) (r+a) (f (r+a)) ≠ 0 ∧
        ModAB.Rational.rawPoint (r+b) (r+s) (f (r+b))
          (μ*f (r+s)) (r+a) (f (r+a))-r =
          secondCoefficient g r a s d*(a*b*d) ∧
        b ≠ 0 ∧ b ≠ a ∧ b ≠ s ∧ 0 < a*b := by
  dsimp only
  obtain ⟨θ,ε,hθ0,hθ1,hε,huniform⟩ :=
    smooth_exists_uniform_negative_coefficients hf hzero hder hK
  refine ⟨θ,ε,hθ0,hθ1,hε,?_⟩
  intro a s d hc
  obtain ⟨hU,hV,hP₁,hP₂,hga,hgs,hgd,hgb,hb,hba,hbs,hq,ht,hJ⟩ :=
    huniform a s d hc
  have ha := hc.latest_ne_zero
  have hs := hc.opposite_ne_zero
  have hd := hc.history_ne_zero
  obtain ⟨has,hsd,had⟩ := CycleUniform.Cone.distinct hc hθ1
  have hratio := residual_ratio_eq_ordinateRatio hzero ha hb hga hgb
  have hmu : 1-ordinateRatio (regularizedReciprocal f r) r a s d ≠ 0 := by
    linarith
  have hu : firstCoefficient (regularizedReciprocal f r) r a s d < 0 := by
    linarith [hU.2]
  have hsdneg : s*d < 0 := by
    have hprod : a^2*(s*d) < 0 := by
      convert mul_neg_of_neg_of_pos hc.opposite hc.same using 1; ring
    by_contra hh
    have hn := mul_nonneg (sq_nonneg a) (le_of_not_gt hh)
    linarith
  have hupos : 0 < firstCoefficient (regularizedReciprocal f r) r a s d*s*d := by
    convert mul_pos_of_neg_of_neg hu hsdneg using 1; ring
  have htpos : 0 < ordinateRatio (regularizedReciprocal f r) r a s d :=
    mul_pos hupos hJ
  have hmubounds : 1-f (r+firstError (regularizedReciprocal f r) r a s d)/f (r+a)
      ∈ Ioo (0:ℝ) 1 := by
    rw [hratio]
    exact ⟨by linarith,by linarith⟩
  have habpos : 0 < a*firstError (regularizedReciprocal f r) r a s d := by
    unfold firstError
    convert mul_pos (sq_pos_of_ne_zero ha) hupos using 1; ring
  exact ⟨hU,hV,
    first_source_denominator_ne_zero hzero ha hs hd has hsd had hga hgs hgd hP₁,
    first_rawPoint_error hzero ha hs hd has hsd had hga hgs hgd hP₁,
    hmubounds,
    second_source_denominator_ne_zero hzero ha hs hb hba hbs has hga hgs hgb hmu hP₂,
    second_rawPoint_error hzero ha hs hd hb hba hbs has hga hgs hgb hmu hP₂,
    hb,hba,hbs,habpos⟩

end
end ModAB.Rational.CycleUniformBridge
