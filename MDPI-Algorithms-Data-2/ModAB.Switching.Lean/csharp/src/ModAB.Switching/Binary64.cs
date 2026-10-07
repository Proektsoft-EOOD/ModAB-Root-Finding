using System.Runtime.CompilerServices;

namespace ModAB.Switching;

/// <summary>Separate binary64 operations in the order used in the error analysis.</summary>
public static class Binary64
{
    public const double UnitRoundoff = 1.1102230246251565404236316680908203125e-16;
    public const double MinimumNormal = 2.2250738585072013830902327173324040642e-308;

    // The method boundaries prevent contraction of a product and an addition.
    [MethodImpl(MethodImplOptions.NoInlining)]
    public static double Add(double x, double y) => x + y;

    [MethodImpl(MethodImplOptions.NoInlining)]
    public static double Subtract(double x, double y) => x - y;

    [MethodImpl(MethodImplOptions.NoInlining)]
    public static double Multiply(double x, double y) => x * y;

    [MethodImpl(MethodImplOptions.NoInlining)]
    public static double Divide(double x, double y) => x / y;

    /// <summary>Checks the rounding and gradual-underflow assumptions on this process.</summary>
    public static bool HasRequiredArithmetic() =>
        Add(1.0, UnitRoundoff) == 1.0 &&
        Add(Math.BitIncrement(1.0), UnitRoundoff) == Math.BitIncrement(Math.BitIncrement(1.0)) &&
        Divide(double.Epsilon, 2.0) == 0.0 &&
        Divide(Multiply(2.0, double.Epsilon), 2.0) == double.Epsilon &&
        Divide(MinimumNormal, 2.0) == Math.ScaleB(1.0, -1023) &&
        Subtract(Multiply(Math.BitIncrement(1.0), Math.BitDecrement(1.0)), 1.0) == 0.0;
}
