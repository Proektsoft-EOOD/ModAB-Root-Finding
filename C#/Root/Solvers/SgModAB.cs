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
            // The true residuals at the endpoints of the bracket.
            double f1 = p1.Y, f2 = p2.Y;
            if (f1 == 0.0)
                return p1.X;
            if (f2 == 0.0)
                return p2.X;

            var isAB = false;
            var sideMoved = 0; // The side that was moved on the previous Anderson-Björck step.
            const int LEFT = -1, RIGHT = 1;
            var fallbackThreshold = 0.0; // The threshold for switching back to the bisection method.
            var yMin = 0.0;
            const double C = 2.0;
            for (var i = 1; i <= MaxIterations; ++i)
            {
                var x3 = isAB
                    ? Node.SafeSecant(p1, p2)
                    : Node.SafeMidpoint(p1, p2);

                var xTol = aTol + rTol * Math.Abs(x3);
                if (p2.X - p1.X <= xTol)
                    return x3;

                // If x3 got clamped, reuse the true residual stored at the endpoint.
                var y3 = x3 == p1.X ? f1 :
                         x3 == p2.X ? f2 :
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
                    yMin = Math.Min(Math.Abs(f1), Math.Abs(f2)); // Best true residual of the bracket.
                else if (double.IsFinite(f2 - f1)) // Avoids overflow in the calculations below.
                {
                    var ym = 0.5 * (f1 + f2);
                    var r = 1 - Math.Abs(ym / (f2 - f1)); // Symmetry factor
                    var k = r * r; // Deviation factor
                    switchToAB = Math.Abs(ym - y3) < k * Math.Abs(ym) + k * Math.Abs(y3);
                    // k·|ym| + k·|y3| cannot overflow; an infinite y3 fails the test.
                }
                var p3 = new Node(x3, y3);
                if (SameSign(f1, y3))
                {
                    if (sideMoved == LEFT)
                        p2.Y *= GetABFactor(y3, f1);
                    else if (isAB)
                        sideMoved = LEFT;

                    p1 = p3;
                    f1 = y3;
                }
                else
                {
                    if (sideMoved == RIGHT)
                        p1.Y *= GetABFactor(y3, f2);
                    else if (isAB)
                        sideMoved = RIGHT;

                    p2 = p3;
                    f2 = y3;
                }
                if (isAB)
                {
                    if (p2.X - p1.X > fallbackThreshold && Math.Abs(y3) > 0.5 * yMin)
                    {
                        isAB = false;
                        sideMoved = 0;
                        p1.Y = f1;
                        p2.Y = f2;
                    }
                    else
                        fallbackThreshold *= 0.5;
                }
                else if (switchToAB)
                {
                    isAB = true;
                    fallbackThreshold = C * (p2.X - p1.X);
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
