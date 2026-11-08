// SgModAB with the scaled switching test of Theorem 3 / Corollary 3 (Section 5.5) instead of the
// unscaled one; everything else is identical to SgModAB. ScaledSwitchGap evaluates the 17 operations
// in the order verified in Lean (ModAB.Switching.Lean/ModAB/Source.lean).
using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{
    public static partial class Solver
    {
        public static double SgModABScaled(Func<double, double> f, 
            double x1, double x2, out ReturnCode returnCode,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(f, x1, x2, aTol, rTol, out Node p1, out Node p2, out var F))
            {
                returnCode = ReturnCode.Invalid;
                return double.NaN;
            }
            returnCode = ReturnCode.Success;
            // p1.Y and p2.Y always hold the true residuals at the endpoints of the bracket.
            if (p1.Y == 0.0)
                return p1.X;
            if (p2.Y == 0.0)
                return p2.X;

            // The Anderson-Björck corrected ordinates, used only on the AB path.
            double f1 = 0.0, f2 = 0.0;
            var isAB = false;
            var sideMoved = 0; // Last endpoint moved in this phase, including the switching bisection.
            const int LEFT = -1, RIGHT = 1;
            var fallbackThreshold = 0.0; // The threshold for switching back to the bisection method.
            var yMin = 0.0;
            const double C = 2.0;
            const int MaxResidualSteps = 3;
            var residualSteps = 0;
            for (var i = 1; i <= MaxIterations; ++i)
            {
                var x3 = isAB
                    ? Node.SafeSecant(p1.X, f1, p2.X, f2)
                    : Node.SafeMidpoint(p1, p2);

                if (p2.X - p1.X <= aTol + rTol * Math.Abs(x3))
                    return x3;

                // If x3 got clamped, reuse the true residual stored at the endpoint.
                var y3 = x3 == p1.X ? p1.Y :
                         x3 == p2.X ? p2.Y :
                         F(x3);

                if (y3 == 0.0)
                    return x3;

                if (double.IsNaN(y3))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                var switchToAB = false;
                if (isAB)
                    yMin = Math.Min(Math.Abs(p1.Y), Math.Abs(p2.Y)); // Best true residual of the bracket.
                else if (double.IsFinite(p1.Y) && double.IsFinite(p2.Y) && double.IsFinite(y3))
                {
                    // Scaled switching test (Theorem 3, Corollary 3). All intermediate values are
                    // bounded by 2, so nothing overflows and subnormal inputs are handled. The switch
                    // is taken only if the computed gap certifies the exact inequality; an uncertified
                    // decision (|gap| <= 2^-48) keeps bisection. Non-finite values keep bisection.
                    switchToAB = ScaledSwitchGap(p1.Y, p2.Y, y3) > SwitchMargin;
                    if (switchToAB)
                    {
                        isAB = true;
                        residualSteps = 0;
                        fallbackThreshold = C * (p2.X - p1.X);
                        f1 = p1.Y;
                        f2 = p2.Y;
                    }
                }
                var p3 = new Node(x3, y3);
                if (isAB) // Anderson-Björck step
                {
                    if (SameSign(p1.Y, y3))
                    {
                        if (sideMoved == LEFT)
                            f2 *= GetABFactor(y3, f1); // Apply Anderson-Björck factor to the right side
                        else
                            sideMoved = LEFT;

                        p1 = p3; f1 = y3;
                    }
                    else
                    {
                        if (sideMoved == RIGHT)
                            f1 *= GetABFactor(y3, f2); // Apply Anderson-Björck factor to the left side
                        else
                            sideMoved = RIGHT;

                        p2 = p3; f2 = y3;
                    }
                    if (p2.X - p1.X > fallbackThreshold)
                    {
                        if (residualSteps >= MaxResidualSteps || 42 * Math.Abs(y3) >= yMin)
                        {
                            isAB = false;
                            sideMoved = 0;
                        }
                        else
                            ++residualSteps;
                    }
                    else
                        residualSteps = 0;

                    fallbackThreshold *= 0.5;
                }
                else // Bisection step
                {
                    if (SameSign(p1.Y, y3))
                        p1 = p3;
                    else
                        p2 = p3;
                }
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }

        // 2^-48 = 32·eps, the certification margin of Corollary 3 (exactly representable).
        private const double SwitchMargin = 3.552713678800501E-15;

        // Scaled difference G^ = R - L of the switching test, evaluated in the order verified in Lean
        // (Source.lean): 17 separately rounded operations; requires finite f1, f2, f3 with f1·f2 < 0.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double ScaledSwitchGap(double f1, double f2, double f3)
        {
            var e1 = Math.Abs(f1);
            var e2 = Math.Abs(f2);
            var rho = Math.Min(e1, e2) / Math.Max(e1, e2);
            var a = 1.0 - rho;
            var b = 1.0 + rho;
            var v = (a / 2.0) / b;
            var r = 1.0 - v;
            var k = r * r;
            var s = Math.Max(Math.Max(e1, e2), Math.Abs(f3));
            var p1 = f1 / s;
            var p2 = f2 / s;
            var p3 = f3 / s;
            var h = (p1 + p2) / 2.0;
            var left = Math.Abs(p3 - h);
            var right = k * Math.Abs(h) + k * Math.Abs(p3);
            return right - left;
        }
    }
}
