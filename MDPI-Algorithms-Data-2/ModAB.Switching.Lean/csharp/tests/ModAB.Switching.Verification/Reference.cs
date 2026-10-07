namespace ModAB.Switching.Verification;

internal readonly record struct ExactValues(Rational K, Rational H, Rational Left, Rational Right, Rational Gap, Rational Scale)
{
    public static ExactValues Evaluate(double f1, double f2, double f3)
    {
        Rational a = Rational.FromDouble(f1), b = Rational.FromDouble(f2), z = Rational.FromDouble(f3);
        Rational scale = Rational.Max(Rational.Max(Rational.Abs(a), Rational.Abs(b)), Rational.Abs(z));
        Rational y = (a + b) / 2;
        Rational d = Rational.Abs(b - a);
        Rational r = 1 - Rational.Abs(y) / d;
        Rational k = r * r;
        Rational h = y / scale, p3 = z / scale;
        Rational left = Rational.Abs(p3 - h);
        Rational right = k * (Rational.Abs(h) + Rational.Abs(p3));
        return new(k, h, left, right, right - left, scale);
    }
}

internal readonly record struct RawValues(double Left, double Right, double Y, double D, double V, double R, double K)
{
    public bool ShouldSwitch => Left < Right;

    public static RawValues? Evaluate(double f1, double f2, double f3)
    {
        double d = Math.Abs(Binary64.Subtract(f2, f1));
        if (!double.IsFinite(d)) return null;
        double y = Binary64.Divide(Binary64.Add(f1, f2), 2.0);
        double v = Binary64.Divide(Math.Abs(y), d);
        double r = Binary64.Subtract(1.0, v);
        double k = Binary64.Multiply(r, r);
        double left = Math.Abs(Binary64.Subtract(f3, y));
        double right = Binary64.Add(Binary64.Multiply(k, Math.Abs(y)), Binary64.Multiply(k, Math.Abs(f3)));
        return new(left, right, y, d, v, r, k);
    }

    public bool MeetsNormalRangeHypotheses(double a, double b, double z)
    {
        double sum = Binary64.Add(a, b);
        double ky = Binary64.Multiply(K, Math.Abs(Y));
        double kz = Binary64.Multiply(K, Math.Abs(z));
        double[] rounded = [sum, Y, D, V, R, K, Binary64.Subtract(z, Y), ky, kz, Right];
        if ((Y == 0 && sum != 0) || (V == 0 && Y != 0) ||
            (ky == 0 && Y != 0) || (kz == 0 && z != 0) ||
            rounded.Any(t => !double.IsFinite(t) || (t != 0 && Math.Abs(t) < Binary64.MinimumNormal)))
            return false;

        Rational af = Rational.FromDouble(a), bf = Rational.FromDouble(b), zf = Rational.FromDouble(z);
        Rational yf = Rational.FromDouble(Y), df = Rational.FromDouble(D);
        Rational vf = Rational.FromDouble(V), rf = Rational.FromDouble(R), kf = Rational.FromDouble(K);
        Rational[] exact = [af + bf, Rational.FromDouble(sum) / 2, bf - af,
            Rational.Abs(yf) / df, 1 - vf, rf * rf, zf - yf,
            kf * Rational.Abs(yf), kf * Rational.Abs(zf), Rational.FromDouble(ky) + Rational.FromDouble(kz)];
        Rational minimum = Rational.FromDouble(Binary64.MinimumNormal), maximum = Rational.FromDouble(double.MaxValue);
        return exact.All(t => t == 0 || (Rational.Abs(t) >= minimum && Rational.Abs(t) <= maximum));
    }
}
