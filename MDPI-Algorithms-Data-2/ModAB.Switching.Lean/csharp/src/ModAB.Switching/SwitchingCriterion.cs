namespace ModAB.Switching;

public enum SwitchingDecision
{
    InvalidInput,
    Indeterminate,
    CertifiedNoSwitch,
    CertifiedSwitch
}

/// <summary>Computed quantities in the scaled comparison; NaN fields denote invalid input.</summary>
public readonly record struct SwitchingResult(
    SwitchingDecision Decision,
    double Scale,
    double Rho,
    double K,
    double P1,
    double P2,
    double P3,
    double H,
    double Left,
    double Right,
    double Gap)
{
    public bool ShouldSwitch => Decision == SwitchingDecision.CertifiedSwitch;
}

public static class SwitchingCriterion
{
    // 32 * 2^-53 = 2^-48. This is a dimensionless error margin.
    public const double Guard = 3.552713678800500929355621337890625e-15;

    /// <summary>
    /// Certifies the sign of R-L for the supplied doubles treated as exact data.
    /// f1 and f2 must be finite, nonzero and of opposite signs; f3 must be finite.
    /// Indeterminate results do not authorize switching.
    /// </summary>
    /// <remarks>
    /// Assumes round-to-nearest, ties-to-even, and gradual underflow.
    /// Errors in the supplied function values are outside this certificate.
    /// </remarks>
    public static SwitchingResult Evaluate(double f1, double f2, double f3)
    {
        if (!IsAdmissible(f1, f2, f3))
            return new(SwitchingDecision.InvalidInput, double.NaN, double.NaN,
                double.NaN, double.NaN, double.NaN, double.NaN, double.NaN,
                double.NaN, double.NaN, double.NaN);

        double e1 = Math.Abs(f1), e2 = Math.Abs(f2);
        double rho = Binary64.Divide(Math.Min(e1, e2), Math.Max(e1, e2));
        double a = Binary64.Subtract(1.0, rho);
        double b = Binary64.Add(1.0, rho);
        double v = Binary64.Divide(Binary64.Divide(a, 2.0), b);
        double r = Binary64.Subtract(1.0, v);
        double k = Binary64.Multiply(r, r);

        double scale = Math.Max(Math.Max(e1, e2), Math.Abs(f3));
        double p1 = Binary64.Divide(f1, scale);
        double p2 = Binary64.Divide(f2, scale);
        double p3 = Binary64.Divide(f3, scale);
        double h = Binary64.Divide(Binary64.Add(p1, p2), 2.0);
        double left = Math.Abs(Binary64.Subtract(p3, h));
        double term1 = Binary64.Multiply(k, Math.Abs(h));
        double term2 = Binary64.Multiply(k, Math.Abs(p3));
        double right = Binary64.Add(term1, term2);
        double gap = Binary64.Subtract(right, left);
        SwitchingDecision decision = gap > Guard
            ? SwitchingDecision.CertifiedSwitch
            : gap < -Guard
                ? SwitchingDecision.CertifiedNoSwitch
                : SwitchingDecision.Indeterminate;

        return new(decision, scale, rho, k, p1, p2, p3, h, left, right, gap);
    }

    public static bool IsAdmissible(double f1, double f2, double f3) =>
        double.IsFinite(f1) && double.IsFinite(f2) && double.IsFinite(f3) &&
        f1 != 0.0 && f2 != 0.0 && (f1 < 0.0) != (f2 < 0.0);
}
