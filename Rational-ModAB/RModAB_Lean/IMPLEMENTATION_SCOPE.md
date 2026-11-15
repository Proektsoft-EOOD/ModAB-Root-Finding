# Algorithm reference and the experimental C# guards

The proof release is 0.8.0 (11 November 2026). Its baseline algorithm reference is the immutable commit in `RModAB_Source_Reference.json`. That record is historical evidence, not the future repository address.

The separately supplied C# demonstration adds two guards to the rational proposal:

```csharp
var k = 1.0 / y2;
if (!double.IsFinite(k))
    return Node.SafeSecant(x1, y1, x2, y2);

// The original rational-offset calculation produces t.
if (!double.IsFinite(t))
    return Node.SafeSecant(x1, y1, x2, y2);
```

These guards change floating-point execution when the reciprocal or rational offset is nonfinite, including cases involving missing-history NaN values. They are an experimental numerical change, not a complete remedy for extreme-range arithmetic. The demonstration retains a case that both variants fail to solve within the iteration cap.

The exact-real definition `ModAB.Rational.candidate` already allows a secant fallback when the rational denominator is zero. Its admissibility results and the broader convergence model already contain that mathematical choice. The local asymptotic statements concern nondegenerate exact-real regimes under the paper's explicit hypotheses. The C# guards do not change those Lean statements or proofs.

This package does not prove complete IEEE754 refinement of either executable. Finite-precision `SafeSecant` need not satisfy every exact-real strict-interiority property at all scales. The existing proof-source comment contrasting the old C# exceptional path with the real model refers to the recorded baseline; the proof file is preserved byte-for-byte.

When the implementation is published, record its actual commit and hashes in a new reference. Retain the baseline record, identify whether the guards were adopted, and update public links. The future repository URL remains unset here.
