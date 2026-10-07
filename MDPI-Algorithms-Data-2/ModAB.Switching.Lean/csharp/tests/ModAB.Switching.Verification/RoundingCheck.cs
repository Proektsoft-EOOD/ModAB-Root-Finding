using System.Numerics;

namespace ModAB.Switching.Verification;

/// <summary>Independent integer rounding oracle and per-operation runtime comparison.</summary>
internal static class RoundingCheck
{
    internal static long CheckedOperations { get; private set; }
    internal static long BoundaryOperations { get; private set; }
    private static bool Negative(double x) => (BitConverter.DoubleToUInt64Bits(x) >> 63) != 0;

    // Integer quotient/remainder, with an explicit half-way parity test. No double
    // arithmetic is used to compute the reference result, including subnormal results.
    internal static ulong RoundBits(Rational x, bool negativeZero = false)
    {
        bool negative = x.Sign < 0 || (x.Sign == 0 && negativeZero);
        ulong sign = negative ? 1UL << 63 : 0;
        if (x.Sign == 0) return sign;
        BigInteger n = BigInteger.Abs(x.Numerator), d = x.Denominator;
        int k = checked((int)(n.GetBitLength() - d.GetBitLength()));
        bool below = k >= 0 ? n < (d << k) : (n << -k) < d;
        if (below) --k;
        int exponent = Math.Max(-1074, k - 52);
        BigInteger numerator = exponent < 0 ? n << -exponent : n;
        BigInteger denominator = exponent > 0 ? d << exponent : d;
        BigInteger q = BigInteger.DivRem(numerator, denominator, out BigInteger rem);
        int half = (2 * rem).CompareTo(denominator);
        if (half > 0 || (half == 0 && !q.IsEven)) ++q;
        if (q.IsZero) return sign;
        if (q >= (BigInteger.One << 53)) { q >>= 1; ++exponent; }
        if (exponent > 971) return sign | 0x7ff0000000000000UL;
        if (q < (BigInteger.One << 52)) return sign | (ulong)q;
        ulong storedExponent = checked((ulong)(exponent + 1075));
        return sign | (storedExponent << 52) | ((ulong)q - (1UL << 52));
    }

    private static double Check(double actual, Rational exact, bool negativeZero, string operation)
    {
        ulong expected = RoundBits(exact, negativeZero);
        ulong found = BitConverter.DoubleToUInt64Bits(actual);
        if (found != expected)
            throw new InvalidOperationException($"{operation}: bits {found:X16}, expected {expected:X16}.");
        ++CheckedOperations;
        return actual;
    }

    public static double Add(double x, double y) => Check(Binary64.Add(x, y),
        Rational.FromDouble(x) + Rational.FromDouble(y), Negative(x) && Negative(y), "add");
    public static double Subtract(double x, double y) => Check(Binary64.Subtract(x, y),
        Rational.FromDouble(x) - Rational.FromDouble(y), Negative(x) && !Negative(y), "subtract");
    public static double Multiply(double x, double y) => Check(Binary64.Multiply(x, y),
        Rational.FromDouble(x) * Rational.FromDouble(y), Negative(x) != Negative(y), "multiply");
    public static double Divide(double x, double y) => Check(Binary64.Divide(x, y),
        Rational.FromDouble(x) / Rational.FromDouble(y), Negative(x) != Negative(y), "divide");

    internal static void Verify(double a, double b, double c)
    {
        long start = CheckedOperations;
        var observed = SwitchingCriterion.Evaluate(a, b, c);
        var traced = TracedCriterion.Evaluate(a, b, c);
        if (CheckedOperations - start != 17) throw new InvalidOperationException("Unexpected operation count.");
        if (observed.Decision != traced.Decision) throw new InvalidOperationException("Decision mismatch.");
        double[] x = [observed.Scale, observed.Rho, observed.K, observed.P1, observed.P2,
            observed.P3, observed.H, observed.Left, observed.Right, observed.Gap];
        double[] y = [traced.Scale, traced.Rho, traced.K, traced.P1, traced.P2,
            traced.P3, traced.H, traced.Left, traced.Right, traced.Gap];
        for (int i = 0; i < x.Length; ++i)
            if (BitConverter.DoubleToUInt64Bits(x[i]) != BitConverter.DoubleToUInt64Bits(y[i]))
                throw new InvalidOperationException($"Original/instrumented field {i} differs.");
    }

    internal static void CheckBoundaries()
    {
        long before = CheckedOperations;
        double[] values = [0.0, BitConverter.UInt64BitsToDouble(1UL << 63),
            double.Epsilon, -double.Epsilon, 2 * double.Epsilon, -2 * double.Epsilon,
            Math.BitDecrement(Binary64.MinimumNormal), Binary64.MinimumNormal,
            Math.BitIncrement(Binary64.MinimumNormal), -Binary64.MinimumNormal,
            0.5, Math.BitDecrement(1.0), 1.0, Math.BitIncrement(1.0),
            Math.BitIncrement(Math.BitIncrement(1.0)), 2.0, -1.0, -2.0,
            Binary64.UnitRoundoff, -Binary64.UnitRoundoff, double.MaxValue, -double.MaxValue];
        foreach (double x in values) foreach (double y in values)
        {
            Add(x,y); Subtract(x,y); Multiply(x,y);
            if (y != 0) Divide(x,y);
        }
        // Half-way points and adjacent significands in every finite binade.
        for (int exponent = -1074; exponent <= 1023; ++exponent)
        {
            double p = Math.ScaleB(1.0, exponent);
            Divide(p, 2); Divide(-p, 2);
            if (exponent >= -1021)
            {
                double halfUlp = Math.ScaleB(1.0, exponent - 53);
                Add(p, halfUlp); Add(Math.BitIncrement(p), halfUlp);
                Subtract(-p, halfUlp); Subtract(-Math.BitIncrement(p), halfUlp);
            }
        }
        BoundaryOperations = CheckedOperations - before;
    }
}
