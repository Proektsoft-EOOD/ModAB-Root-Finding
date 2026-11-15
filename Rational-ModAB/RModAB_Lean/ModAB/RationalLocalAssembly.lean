import ModAB.RationalLocalSafeguards
import ModAB.RationalOrder
import ModAB.RationalCycleEntry

/-!
Local convergence and the tribonacci Q-order in one theorem. The errors are
not assumed to converge: bounded coefficients and three sufficiently small
seeds imply convergence, and the coefficient limit determines the Q-order.
The displayed recurrence remains an explicit interface to interpolation.
-/
namespace ModAB.RationalLocalAssembly
open Filter Topology
noncomputable section

/-- Three small seeds, a uniformly bounded coefficient and its nonzero limit
suffice for genuine tribonacci Q-order. No error-limit premise is present. -/
theorem bounded_coefficient_recurrence_hasQOrder {e H : ℕ → ℝ} {B ε K : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε) (hq : B * ε ^ 2 < 1)
    (h0 : |e 0| ≤ ε) (h1 : |e 1| ≤ ε) (h2 : |e 2| ≤ ε)
    (hH : ∀ n, |H n| ≤ B) (hK : K ≠ 0)
    (hHlim : Tendsto H atTop (𝓝 K)) (hnz : ∀ n, e n ≠ 0)
    (hrec : ∀ n, e (n + 3) = H n * e (n + 2) * e (n + 1) * e n) :
    ModAB.RationalOrder.HasQOrder e ModAB.RationalOrder.tau
      (|K| ^ ((ModAB.RationalOrder.tau - 1) / 2)) := by
  have hlim := ModAB.RationalLocalSafeguards.bounded_coefficient_recurrence_tendsto
    hB hε hq h0 h1 h2 hH hrec
  exact ModAB.RationalOrder.three_error_recurrence_hasQOrder hnz hK hlim hHlim hrec

/-- The invariant two-step cone and the two coefficient limits imply actual
Q-order three for both subsequences, as well as the unequal consecutive-step
orders. Neither error convergence nor a quotient limit is assumed. -/
theorem cycle_cone_hasQOrders {θ ε C K : ℝ}
    (R : ModAB.Rational.CycleEntry.Run θ ε C)
    (hθ0 : 0 < θ) (hθ1 : θ < 1) (hC : 0 ≤ C) (hε : 0 ≤ ε)
    (hsmall : C * ε ^ 2 ≤ θ) (hK : K ≠ 0)
    (hu : Tendsto R.u atTop (𝓝 K)) (hv : Tendsto R.v atTop (𝓝 K)) :
    ModAB.RationalOrder.HasQOrder R.a 3 |K| ∧
    ModAB.RationalOrder.HasQOrder (fun n => R.s (n + 1)) 3 |K| ∧
    Tendsto (fun n => |R.s (n + 1)| / |R.a n| ^ (9 / 5 : ℝ)) atTop
      (𝓝 (|K| ^ (2 / 5 : ℝ))) ∧
    Tendsto (fun n => |R.a (n + 1)| / |R.s (n + 1)| ^ (5 / 3 : ℝ)) atTop
      (𝓝 (|K| ^ (1 / 3 : ℝ))) := by
  have he := R.errors_tendsto hθ0 hθ1 hC hε hsmall
  have ho := R.orders hθ0 hθ1 hC hε hsmall hK hu hv
  have hκ : 0 < |K| := abs_pos.mpr hK
  have haq := R.latest_cubic_order hθ0 hθ1 hC hε hsmall hK hu hv
  refine ⟨⟨by norm_num, hκ, he.1, ?_⟩,
    ⟨by norm_num, hκ, (tendsto_add_atTop_iff_nat 1).mpr he.2.1, ?_⟩, ho.2⟩
  · simpa only [Real.rpow_ofNat] using haq
  · simpa only [Nat.add_assoc, Nat.reduceAdd, Real.rpow_ofNat] using ho.1

end
end ModAB.RationalLocalAssembly
