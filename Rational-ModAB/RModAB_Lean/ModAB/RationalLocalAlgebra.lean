import ModAB.Rational

/-!
Exact identities behind the local rational-step analysis.  The nodes and values
are arbitrary real data.  The theorems state every required nonvanishing
condition; they do not assume differentiability or an asymptotic recurrence.
The positivity result at the end is the finite-data crossing mechanism used in
the paper's local K > 0 argument.  Passing from C³ hypotheses to the signs of
divided differences is a separate analytic statement.
-/
namespace ModAB.Rational.LocalAlgebra
noncomputable section

/-- First divided difference, with the function values supplied explicitly. -/
def dd1 (x₁ x₂ v₁ v₂ : ℝ) : ℝ := (v₂-v₁)/(x₂-x₁)

/-- Second divided difference in Newton form. -/
def dd2 (x₁ x₂ x₃ v₁ v₂ v₃ : ℝ) : ℝ :=
  (dd1 x₂ x₃ v₂ v₃-dd1 x₁ x₂ v₁ v₂)/(x₃-x₁)

/-- Barycentric quotient specifying the root of the linear-fractional model.
Its denominator is required nonzero in all theorems that interpret it as a root. -/
def baryRoot (x₁ x₂ x₃ y₁ y₂ y₃ : ℝ) : ℝ :=
  dd2 x₁ x₂ x₃ (x₁/y₁) (x₂/y₂) (x₃/y₃) /
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)

/-- The Newton quadratic through g evaluated at the proposed root r. -/
def poleValue (x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ) : ℝ :=
  g₁-dd1 x₁ x₂ g₁ g₂*(x₁-r)+
    dd2 x₁ x₂ x₃ g₁ g₂ g₃*(x₁-r)*(x₂-r)

theorem dd2_add (x₁ x₂ x₃ a₁ a₂ a₃ b₁ b₂ b₃ : ℝ) :
    dd2 x₁ x₂ x₃ (a₁+b₁) (a₂+b₂) (a₃+b₃) =
      dd2 x₁ x₂ x₃ a₁ a₂ a₃+dd2 x₁ x₂ x₃ b₁ b₂ b₃ := by
  unfold dd2 dd1
  ring

theorem dd2_smul (x₁ x₂ x₃ c a₁ a₂ a₃ : ℝ) :
    dd2 x₁ x₂ x₃ (c*a₁) (c*a₂) (c*a₃) =
      c*dd2 x₁ x₂ x₃ a₁ a₂ a₃ := by
  unfold dd2 dd1
  ring

theorem dd2_barycentric {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ =
      v₁/((x₁-x₂)*(x₁-x₃))+
      v₂/((x₂-x₁)*(x₂-x₃))+v₃/((x₃-x₁)*(x₃-x₂)) := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  have h₁₃ : x₁-x₃≠0 := by intro h; apply h₃₁; linarith
  have h₂₃ : x₂-x₃≠0 := by intro h; apply h₃₂; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂,h₁₃,h₂₃]
  ring

theorem dd2_affine {x₁ x₂ x₃ p q : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) :
    dd2 x₁ x₂ x₃ (p*x₁+q) (p*x₂+q) (p*x₃+q)=0 := by
  have hleft : dd1 x₁ x₂ (p*x₁+q) (p*x₂+q)=p := by
    unfold dd1
    field_simp [h₂₁]
    ring
  have hright : dd1 x₂ x₃ (p*x₂+q) (p*x₃+q)=p := by
    unfold dd1
    field_simp [h₃₂]
    ring
  simp [dd2,hleft,hright]

theorem dd2_swap_last {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₁ x₃ x₂ v₁ v₃ v₂ := by
  have h₂₃ : x₂-x₃≠0 := by intro h; apply h₃₂; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₂₃]
  ring

theorem dd2_swap_first {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₂ x₁ x₃ v₂ v₁ v₃ := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂]
  ring

theorem dd2_cyclic {x₁ x₂ x₃ v₁ v₂ v₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ v₃ = dd2 x₂ x₃ x₁ v₂ v₃ v₁ := by
  have h₁₂ : x₁-x₂≠0 := by intro h; apply h₂₁; linarith
  have h₁₃ : x₁-x₃≠0 := by intro h; apply h₃₁; linarith
  unfold dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁₂,h₁₃]
  ring

theorem baryRoot_swap_last {x₁ x₂ x₃ y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ y₃ = baryRoot x₁ x₃ x₂ y₁ y₃ y₂ := by
  unfold baryRoot
  rw [dd2_swap_last h₂₁ h₃₂ h₃₁,dd2_swap_last h₂₁ h₃₂ h₃₁]

theorem dd2_change_third {x₁ x₂ x₃ v₁ v₂ v₃ Δ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) :
    dd2 x₁ x₂ x₃ v₁ v₂ (v₃+Δ) =
      dd2 x₁ x₂ x₃ v₁ v₂ v₃+Δ/((x₃-x₁)*(x₃-x₂)) := by
  unfold dd2 dd1
  field_simp [h₃₂,h₃₁]
  ring

/-- This linear identity does not require nonzero ordinate or node denominators. -/
theorem dd2_error_decomposition (x₁ x₂ x₃ y₁ y₂ y₃ r : ℝ) :
    dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) =
      dd2 x₁ x₂ x₃ (x₁/y₁) (x₂/y₂) (x₃/y₃)-
        r*dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃) := by
  unfold dd2 dd1
  ring

theorem baryRoot_error_identity {x₁ x₂ x₃ y₁ y₂ y₃ r : ℝ}
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)≠0) :
    (baryRoot x₁ x₂ x₃ y₁ y₂ y₃-r)*
      dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃) =
        dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) := by
  rw [dd2_error_decomposition]
  unfold baryRoot
  field_simp [hH]

/-- The pole-part identity is exact, not a Taylor approximation. -/
theorem pole_identity {x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0) :
    dd2 x₁ x₂ x₃ (g₁/(x₁-r)) (g₂/(x₂-r)) (g₃/(x₃-r))*
      ((x₁-r)*(x₂-r)*(x₃-r)) = poleValue x₁ x₂ x₃ r g₁ g₂ g₃ := by
  unfold poleValue dd2 dd1
  field_simp [h₂₁,h₃₂,h₃₁,h₁,h₂,h₃]
  ring

theorem pole_identity_from_values {x₁ x₂ x₃ r y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0) :
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)*((x₁-r)*(x₂-r)*(x₃-r)) =
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) := by
  have h := pole_identity (g₁:=(x₁-r)/y₁) (g₂:=(x₂-r)/y₂)
    (g₃:=(x₃-r)/y₃) h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  have hc (e y : ℝ) (he : e≠0) : (e/y)/e=1/y := by
    calc (e/y)/e = (e*e⁻¹)*y⁻¹ := by ring
         _ = 1/y := by simp [he]
  simpa only [hc (x₁-r) y₁ h₁,hc (x₂-r) y₂ h₂,hc (x₃-r) y₃ h₃] using h

theorem poleValue_margin (x₁ x₂ x₃ r g₁ g₂ g₃ : ℝ) :
    poleValue x₁ x₂ x₃ r g₁ g₂ g₃-
      ((x₁-r)*(x₂-r))*dd2 x₁ x₂ x₃ g₁ g₂ g₃ =
        g₁-dd1 x₁ x₂ g₁ g₂*(x₁-r) := by
  unfold poleValue
  ring

/-- Exact error product for unscaled values.  Analytic convergence of its
coefficient to K is deliberately not an assumption or conclusion here. -/
theorem baryRoot_true_error {x₁ x₂ x₃ r y₁ y₂ y₃ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/y₃)≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ y₃-r =
      dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃) /
        poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)*
          ((x₁-r)*(x₂-r)*(x₃-r)) := by
  have hp := pole_identity_from_values (y₁:=y₁) (y₂:=y₂) (y₃:=y₃)
    h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  rw [←hp]
  rw [dd2_error_decomposition]
  unfold baryRoot
  field_simp [hH,h₁,h₂,h₃]

/-- Replacing one ordinate by μ times that ordinate changes its reciprocal
by the exact relative increment δ=1/μ-1. -/
theorem rescale_div {a y μ : ℝ} (hy : y≠0) (hμ : μ≠0) :
    a/(μ*y)=a/y+(1/μ-1)*(a/y) := by
  field_simp [hy,hμ]
  ring

theorem dd2_scaled_third {x₁ x₂ x₃ a₁ a₂ a₃ y₁ y₂ y₃ μ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) (hy₃ : y₃≠0) (hμ : μ≠0) :
    dd2 x₁ x₂ x₃ (a₁/y₁) (a₂/y₂) (a₃/(μ*y₃)) =
      dd2 x₁ x₂ x₃ (a₁/y₁) (a₂/y₂) (a₃/y₃)+
        (1/μ-1)*(a₃/y₃)/((x₃-x₁)*(x₃-x₂)) := by
  rw [rescale_div hy₃ hμ,dd2_change_third h₃₂ h₃₁]

theorem poleValue_scaled_third {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0) (hy₃ : y₃≠0) (hμ : μ≠0) :
    poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/(μ*y₃)) =
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)) := by
  unfold poleValue
  rw [dd2_scaled_third h₃₂ h₃₁ hy₃ hμ]
  ring

/-- Exact mixed-ordinate error formula.  This retains the denominator and
therefore remains valid even when two possible leading terms cancel. -/
theorem baryRoot_scaled_error {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0) :
    baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃)-r =
      (dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)/((x₃-x₁)*(x₃-x₂))) /
      (poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)+
        (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)))*
      ((x₁-r)*(x₂-r)*(x₃-r)) := by
  rw [baryRoot_true_error h₂₁ h₃₂ h₃₁ h₁ h₂ h₃ hH,
    dd2_scaled_third h₃₂ h₃₁ hy₃ hμ,
    poleValue_scaled_third h₃₂ h₃₁ hy₃ hμ]

/-- A finite-data crossing lemma: the rational error is a strict fraction of
the error at the opposite endpoint, uniformly over nonnegative scaling δ. -/
theorem scaled_crossing_ratio {E G P g δ π : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) :
    E*(G+δ*g/π)/(P+E*(δ*g/π)) ∈ Set.Ioo (0:ℝ) 1 := by
  have ht : 0≤δ*g/π := div_nonneg (mul_nonneg hδ hg.le) hπ.le
  have hP : 0<P := by nlinarith [mul_pos hE hG]
  have hD : 0<P+E*(δ*g/π) := by nlinarith [mul_nonneg hE.le ht]
  constructor
  · exact div_pos (mul_pos hE (by linarith)) hD
  · apply (div_lt_one hD).2
    nlinarith

theorem scaled_pole_pos {E G P g δ π : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) : 0<P+δ*g*E/π := by
  have ht : 0≤δ*g*E/π :=
    div_nonneg (mul_nonneg (mul_nonneg hδ hg.le) hE.le) hπ.le
  nlinarith [mul_pos hE hG]

/-- Under the local positive-coefficient data conditions, the mixed model's
denominator is nonzero; this condition need not be assumed independently. -/
theorem scaled_denominator_ne_zero {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hE : 0<(x₁-r)*(x₂-r))
    (hG : 0<dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃))
    (hg : 0<(x₃-r)/y₃) (hδ : 0≤1/μ-1)
    (hπ : 0<(x₃-x₁)*(x₃-x₂))
    (hmargin : 0<
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)-
        ((x₁-r)*(x₂-r))*
          dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)) :
    dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0 := by
  have hp := pole_identity_from_values (y₁:=y₁) (y₂:=y₂) (y₃:=μ*y₃)
    h₂₁ h₃₂ h₃₁ h₁ h₂ h₃
  rw [poleValue_scaled_third h₃₂ h₃₁ hy₃ hμ] at hp
  have hpos := scaled_pole_pos hE hG hg hδ hπ hmargin
  intro hz
  rw [hz,zero_mul] at hp
  have he : (1/μ-1)*((x₃-r)/y₃)*((x₁-r)*(x₂-r))/((x₃-x₁)*(x₃-x₂)) =
      (1/μ-1)*((x₃-r)/y₃)*(x₁-r)*(x₂-r)/((x₃-x₁)*(x₃-x₂)) := by ring
  rw [he] at hpos
  linarith

theorem scaled_error_positive_side {E G P g δ π es : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) (hes : 0<es) :
    (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) ∈ Set.Ioo (0:ℝ) es := by
  have hr := scaled_crossing_ratio hE hG hg hδ hπ hmargin
  have he : (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) =
      (E*(G+δ*g/π)/(P+E*(δ*g/π)))*es := by ring
  rw [he]
  exact ⟨mul_pos hr.1 hes,by simpa using mul_lt_mul_of_pos_right hr.2 hes⟩

theorem scaled_error_negative_side {E G P g δ π es : ℝ}
    (hE : 0<E) (hG : 0<G) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hmargin : 0<P-E*G) (hes : es<0) :
    (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) ∈ Set.Ioo es (0:ℝ) := by
  have hr := scaled_crossing_ratio hE hG hg hδ hπ hmargin
  have he : (G+δ*g/π)/(P+E*(δ*g/π))*(E*es) =
      (E*(G+δ*g/π)/(P+E*(δ*g/π)))*es := by ring
  rw [he]
  exact ⟨by simpa using mul_lt_mul_of_neg_right hr.2 hes,mul_neg_of_pos_of_neg hr.1 hes⟩

/-- In a same-side step the scaling cannot grow arbitrarily.  The sign of G
need not be assumed: it follows from the same-side inequality and δ ≥ 0. -/
theorem delta_bound_of_same_side {E G P g δ π : ℝ}
    (hE : 0<E) (hP : 0<P) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π))<0) :
    δ< -π*G/g := by
  have ht : 0≤δ*g/π := div_nonneg (mul_nonneg hδ hg.le) hπ.le
  have hD : 0<P+E*(δ*g/π) := by nlinarith [mul_nonneg hE.le ht]
  have hn : E*(G+δ*g/π)<0 := by
    have h := (div_lt_iff₀ hD).mp hsame
    simpa using h
  have hs : G+δ*g/π<0 := by
    by_contra h
    have hp := mul_nonneg hE.le (le_of_not_gt h)
    linarith
  have hterm : δ*g/π < -G := by linarith
  have hmul := (div_lt_iff₀ hπ).mp hterm
  apply (lt_div_iff₀ hg).2
  nlinarith

theorem negative_coefficient_of_same_side {E G P g δ π : ℝ}
    (hE : 0<E) (hP : 0<P) (hg : 0<g) (hδ : 0≤δ) (hπ : 0<π)
    (hsame : E*(G+δ*g/π)/(P+E*(δ*g/π))<0) : G<0 := by
  have h := delta_bound_of_same_side hE hP hg hδ hπ hsame
  have hp : 0< -π*G/g := lt_of_le_of_lt hδ h
  have hn : 0< -π*G := by
    have hh := (lt_div_iff₀ hg).mp hp
    simpa using hh
  nlinarith

/-- An a priori upper bound on δ gives a positive lower bound for the weight. -/
theorem weight_lower_bound {δ M : ℝ} (hδ : 0≤δ) (hM : δ≤M) :
    0<1/(1+M) ∧ 1/(1+M)≤1/(1+δ) := by
  have hδp : 0<1+δ := by linarith
  have hMp : 0<1+M := by linarith
  exact ⟨one_div_pos.mpr hMp,one_div_le_one_div_of_le hδp (by linarith)⟩

theorem weight_from_delta {μ : ℝ} (hμ : μ≠0) : 1/(1+(1/μ-1))=μ := by
  field_simp [hμ]
  ring

/-- The exact mixed-ordinate root crosses toward x₃ and stops strictly before
x₃.  The hypotheses are finite-data signs, rather than an assumed crossing. -/
theorem baryRoot_scaled_between {x₁ x₂ x₃ r y₁ y₂ y₃ μ : ℝ}
    (h₂₁ : x₂-x₁≠0) (h₃₂ : x₃-x₂≠0) (h₃₁ : x₃-x₁≠0)
    (h₁ : x₁-r≠0) (h₂ : x₂-r≠0) (h₃ : x₃-r≠0)
    (hy₃ : y₃≠0) (hμ : μ≠0)
    (hH : dd2 x₁ x₂ x₃ (1/y₁) (1/y₂) (1/(μ*y₃))≠0)
    (hE : 0<(x₁-r)*(x₂-r))
    (hG : 0<dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃))
    (hg : 0<(x₃-r)/y₃) (hδ : 0≤1/μ-1)
    (hπ : 0<(x₃-x₁)*(x₃-x₂))
    (hmargin : 0<
      poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)-
        ((x₁-r)*(x₂-r))*
          dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)) :
    (r<x₃ → baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃) ∈ Set.Ioo r x₃) ∧
    (x₃<r → baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃) ∈ Set.Ioo x₃ r) := by
  let E := (x₁-r)*(x₂-r)
  let G := dd2 x₁ x₂ x₃ ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)
  let P := poleValue x₁ x₂ x₃ r ((x₁-r)/y₁) ((x₂-r)/y₂) ((x₃-r)/y₃)
  let g := (x₃-r)/y₃
  let δ := 1/μ-1
  let π := (x₃-x₁)*(x₃-x₂)
  have he : baryRoot x₁ x₂ x₃ y₁ y₂ (μ*y₃)-r =
      (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) := by
    rw [baryRoot_scaled_error h₂₁ h₃₂ h₃₁ h₁ h₂ h₃ hy₃ hμ hH]
    dsimp [E,G,P,g,δ,π]
    congr 1
    congr 1
    ring
  constructor
  · intro hs
    have hp := scaled_error_positive_side hE hG hg hδ hπ hmargin (sub_pos.mpr hs)
    change (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) ∈ Set.Ioo 0 (x₃-r) at hp
    rw [←he] at hp
    exact ⟨by linarith [hp.1],by linarith [hp.2]⟩
  · intro hs
    have hp := scaled_error_negative_side hE hG hg hδ hπ hmargin (sub_neg.mpr hs)
    change (G+δ*g/π)/(P+E*(δ*g/π))*(E*(x₃-r)) ∈ Set.Ioo (x₃-r) 0 at hp
    rw [←he] at hp
    exact ⟨by linarith [hp.1],by linarith [hp.2]⟩

/-- The barycentric root agrees with the normalized proposal already used in
the algorithm specification.  Neither quotient is interpreted at a zero
denominator. -/
theorem bary_denominator_eq {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0) :
    dd2 a b p (1/ya) (1/yb) (1/yp) =
      -ModAB.Rational.denominator a b ya yb p yp*yb/
        (ya*yp*(b-a)*(p-b)*(p-a)) := by
  unfold dd2 dd1 ModAB.Rational.denominator ModAB.Rational.normalizedDenominator
  field_simp [hba,hpb,hpa,hya,hyb,hyp]
  ring

theorem bary_denominator_ne_zero {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0) :
    dd2 a b p (1/ya) (1/yb) (1/yp)≠0 := by
  rw [bary_denominator_eq hba hpb hpa hya hyb hyp]
  exact div_ne_zero (mul_ne_zero (neg_ne_zero.mpr hD) hyb)
    (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hya hyp) hba) hpb) hpa)

theorem baryRoot_eq_rawPoint {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0)
    (hH : dd2 a b p (1/ya) (1/yb) (1/yp)≠0) :
    baryRoot a b p ya yb yp = ModAB.Rational.rawPoint a b ya yb p yp := by
  have hD' : ModAB.Rational.normalizedDenominator (b-a) (p-a) (ya/yb) (yp/yb)≠0 := hD
  unfold baryRoot ModAB.Rational.rawPoint ModAB.Rational.normalizedOffset
  field_simp [hH,hD']
  unfold dd2 dd1 ModAB.Rational.normalizedNumerator ModAB.Rational.normalizedDenominator
  field_simp [hba,hpb,hpa,hya,hyb,hyp]
  ring

/-- A source-oriented bridge with just the algorithm's denominator guard. -/
theorem rawPoint_eq_baryRoot {a b p ya yb yp : ℝ}
    (hba : b-a≠0) (hpb : p-b≠0) (hpa : p-a≠0)
    (hya : ya≠0) (hyb : yb≠0) (hyp : yp≠0)
    (hD : ModAB.Rational.denominator a b ya yb p yp≠0) :
    ModAB.Rational.rawPoint a b ya yb p yp = baryRoot a b p ya yb yp := by
  exact (baryRoot_eq_rawPoint hba hpb hpa hya hyb hyp hD
    (bary_denominator_ne_zero hba hpb hpa hya hyb hyp hD)).symm

end
end ModAB.Rational.LocalAlgebra
