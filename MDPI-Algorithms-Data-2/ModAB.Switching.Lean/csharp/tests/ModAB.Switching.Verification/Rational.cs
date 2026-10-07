using System.Numerics;

namespace ModAB.Switching.Verification;

/// <summary>Reduced exact rational arithmetic for the independent reference calculation.</summary>
internal readonly struct Rational : IComparable<Rational>, IEquatable<Rational>
{
    public BigInteger Numerator { get; }
    public BigInteger Denominator { get; }
    public static Rational Zero => new(0);
    public static Rational One => new(1);
    public int Sign => Numerator.Sign;

    public Rational(BigInteger numerator) : this(numerator, BigInteger.One) { }

    public Rational(BigInteger numerator, BigInteger denominator)
    {
        if (denominator.IsZero) throw new DivideByZeroException();
        if (denominator.Sign < 0) { numerator = -numerator; denominator = -denominator; }
        BigInteger gcd = BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator), denominator);
        Numerator = numerator / gcd;
        Denominator = denominator / gcd;
    }

    public static Rational FromDouble(double value)
    {
        if (!double.IsFinite(value)) throw new ArgumentOutOfRangeException(nameof(value));
        ulong bits = BitConverter.DoubleToUInt64Bits(value);
        int storedExponent = (int)((bits >> 52) & 0x7ff);
        ulong significand = bits & 0x000fffffffffffffUL;
        int exponent;
        if (storedExponent == 0) exponent = -1074;
        else { significand |= 1UL << 52; exponent = storedExponent - 1023 - 52; }
        BigInteger numerator = significand;
        if ((bits >> 63) != 0) numerator = -numerator;
        return exponent >= 0
            ? new Rational(numerator << exponent)
            : new Rational(numerator, BigInteger.One << -exponent);
    }

    // Used only to display results. All assertions and maxima use exact rationals.
    public double ToApproximateDouble()
    {
        if (Numerator.IsZero) return 0.0;
        BigInteger n = BigInteger.Abs(Numerator);
        int ns = (int)Math.Max(0, n.GetBitLength() - 54);
        int ds = (int)Math.Max(0, Denominator.GetBitLength() - 54);
        return Numerator.Sign * Math.ScaleB((double)(n >> ns) / (double)(Denominator >> ds), ns - ds);
    }

    public static Rational Abs(Rational x) => x.Sign < 0 ? -x : x;
    public static Rational Max(Rational x, Rational y) => x >= y ? x : y;
    public int CompareTo(Rational other) => (Numerator * other.Denominator).CompareTo(other.Numerator * Denominator);
    public bool Equals(Rational other) => Numerator == other.Numerator && Denominator == other.Denominator;
    public override bool Equals(object? other) => other is Rational r && Equals(r);
    public override int GetHashCode() => HashCode.Combine(Numerator, Denominator);
    public override string ToString() => $"{Numerator}/{Denominator}";

    public static implicit operator Rational(int x) => new(x);
    public static Rational operator +(Rational x, Rational y) => new(x.Numerator * y.Denominator + y.Numerator * x.Denominator, x.Denominator * y.Denominator);
    public static Rational operator -(Rational x, Rational y) => new(x.Numerator * y.Denominator - y.Numerator * x.Denominator, x.Denominator * y.Denominator);
    public static Rational operator -(Rational x) => new(-x.Numerator, x.Denominator);
    public static Rational operator *(Rational x, Rational y) => new(x.Numerator * y.Numerator, x.Denominator * y.Denominator);
    public static Rational operator /(Rational x, Rational y) => new(x.Numerator * y.Denominator, x.Denominator * y.Numerator);
    public static bool operator <(Rational x, Rational y) => x.CompareTo(y) < 0;
    public static bool operator >(Rational x, Rational y) => x.CompareTo(y) > 0;
    public static bool operator <=(Rational x, Rational y) => x.CompareTo(y) <= 0;
    public static bool operator >=(Rational x, Rational y) => x.CompareTo(y) >= 0;
    public static bool operator ==(Rational x, Rational y) => x.Equals(y);
    public static bool operator !=(Rational x, Rational y) => !x.Equals(y);
}
