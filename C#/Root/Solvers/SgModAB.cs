using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{
    public static partial class Solver
    {
        public static double SgModAB(Func<double, double> f, double x1, double x2, out ReturnCode returnCode,
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
            var sideMoved = 0; // The side that was moved on the previous Anderson-Björck step.
            const int LEFT = -1, RIGHT = 1;
            var fallbackThreshold = 0.0; // The threshold for switching back to the bisection method.
            var yMin = 0.0;
            const double C = 2.0;
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
                else if (double.IsFinite(p2.Y - p1.Y)) // Avoids overflow in the calculations below.
                {
                    var ym = 0.5 * (p1.Y + p2.Y);
                    var r = 1 - Math.Abs(ym / (p2.Y - p1.Y)); // Symmetry factor
                    var k = r * r; // Deviation factor
                    switchToAB = Math.Abs(ym - y3) < k * Math.Abs(ym) + k * Math.Abs(y3);
                    // k·|ym| + k·|y3| cannot overflow; an infinite y3 fails the test.
                }
                var p3 = new Node(x3, y3);
                if (isAB) // Anderson-Björck step
                {
                    if (SameSign(p1.Y, y3))
                    {
                        if (sideMoved == LEFT)
                            f2 *= GetABFactor(y3, p1.Y); // Apply Anderson-Björck factor to the right side
                        else
                            sideMoved = LEFT;

                        p1 = p3; f1 = y3;
                    }
                    else
                    {
                        if (sideMoved == RIGHT)
                            f1 *= GetABFactor(y3, p2.Y); // Apply Anderson-Björck factor to the left side
                        else
                            sideMoved = RIGHT;

                        p2 = p3; f2 = y3;
                    }
                    if (p2.X - p1.X > fallbackThreshold && Math.Abs(y3) > 0.5 * yMin)
                    {
                        isAB = false;
                        sideMoved = 0;
                    }
                    else
                        fallbackThreshold *= 0.5;
                }
                else // Bisection step
                {
                    if (SameSign(p1.Y, y3))
                        p1 = p3;
                    else
                        p2 = p3;

                    if (switchToAB)
                    {
                        isAB = true;
                        fallbackThreshold = C * (p2.X - p1.X);
                        f1 = p1.Y;
                        f2 = p2.Y;
                    }
                }
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double GetABFactor(double y3, double yMoved)
        {
            var m = 1.0 - y3 / yMoved;
            return m > 0.0 ? m : 0.5;
        }
    }
}
