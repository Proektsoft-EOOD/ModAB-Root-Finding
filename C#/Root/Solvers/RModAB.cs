using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{

    public static partial class Solver
    {
        // The Anderson-Björck phase adds two changes:
        //     1. Rational correction: the AB secant is replaced by the secant through the true
        //        residuals, corrected for curvature by a 3-point rational model through the bracket
        //        and the endpoint it most recently dropped. The AB secant remains the fallback when
        //        the correction is unusable or the curvature changes sign within the bracket.
        //     2. Closing step: an AB point is kept at least eps/2 away from the endpoint that was
        //        just evaluated, so that one-sided convergence still closes the bracket.
        public static double RModAB(Func<double, double> f, double x1, double x2, out ReturnCode returnCode,
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
            const double C = 2.0;
            var fallbackThreshold = 0.0; // The threshold for switching back to the bisection method.
            // A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
            // than 1 % lies on a plateau; the first AB step is then not rescaled.
            const double PlateauRatio = 0.99;
            // Each return from AB to bisection halves the margin of the switching test, so that
            // borderline switches are rejected after unsuccessful AB phases.
            var fallbacks = 0;
            var lastMoved = 0; // The endpoint set by the last evaluation: LEFT, RIGHT or 0 for none.
            // The last endpoints dropped on each side.
            double xl = double.NaN, yl = double.NaN, xr = double.NaN, yr = double.NaN;
            for (var i = 1; i <= MaxIterations; ++i)
            {
                var dx = p2.X - p1.X;
                double x3, eps;
                if (isAB)
                {
                    // Rational correction through the endpoint dropped by the last update, or the A&B
                    // secant when the curvature changes sign within the bracket. A missing endpoint is
                    // NaN, which fails both tests, so it needs no check of its own.
                    if (CurvatureChangesSign(xl, yl, p1.X, p1.Y, p2.X, p2.Y, xr, yr))
                        x3 = Node.SafeSecant(p1.X, f1, p2.X, f2);
                    else
                        x3 = lastMoved == LEFT ?
                            RationalStep(p1.X, p1.Y, p2.X, p2.Y, xl, yl, f1, f2) :
                            RationalStep(p1.X, p1.Y, p2.X, p2.Y, xr, yr, f1, f2);

                    eps = aTol + rTol * Math.Abs(x3);
                    // Closing step: when AB converges from one side, its next point lands within eps/2
                    // of the endpoint just evaluated and the bracket does not shrink. A step of eps/2
                    // towards the far side then crosses the root and closes the bracket.
                    var h = 0.5 * eps;
                    if (lastMoved == LEFT)
                    {
                        if (x3 - p1.X < h && p1.X + h < p2.X)
                            x3 = p1.X + h;
                    }
                    else if (p2.X - x3 < h && p2.X - h > p1.X)
                        x3 = p2.X - h;
                }
                else
                {
                    x3 = Node.SafeMidpoint(p1, p2);
                    eps = aTol + rTol * Math.Abs(x3);
                }
                if (dx <= eps)
                    return x3;

                // If x3 got clamped, reuse the true residual stored at the endpoint.
                var y3 = x3 == p1.X ? p1.Y :
                         x3 == p2.X ? p2.Y :
                         F(x3);

                if (!(Math.Abs(y3) > 0.0)) // Exit on zero or NaN
                {
                    if (y3 == 0d)
                        return x3;

                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                var switching = false; // This bisection switches to AB.
                if (!isAB)
                {
                    var dy = p2.Y - p1.Y;
                    if (double.IsFinite(dy)) // Avoids overflow in the calculations below.
                    {
                        var ym = 0.5 * (p1.Y + p2.Y);
                        var r = 1d - Math.Abs(ym / dy); // Symmetry factor
                        var k = r * r; // Deviation factor
                        if (Math.Abs(ym - y3) < Math.ScaleB(k * Math.Abs(ym) + k * Math.Abs(y3), -fallbacks))
                        {   // k·|ym| + k·|y3| cannot overflow; an infinite y3 fails the test.
                            isAB = true;
                            switching = true;
                            fallbackThreshold = C * dx;
                            f1 = p1.Y;
                            f2 = p2.Y;
                        }
                    }
                }
                var p3 = new Node(x3, y3);
                if (isAB) // Anderson-Björck step
                {
                    double y1 = p1.Y, y2 = p2.Y; // True residuals of the bracket before the update.
                    if (SameSign(y1, y3))
                    {
                        if (sideMoved == LEFT)
                            f2 *= GetABFactor(y3, f1); // Apply Anderson-Björck factor to the right side
                        else if (!switching || Math.Abs(y3) <= PlateauRatio * Math.Abs(y1))
                            sideMoved = LEFT; // A switching bisection on a plateau starts the AB phase afresh.

                        if (x3 != p1.X)
                        {
                            xl = p1.X; yl = p1.Y;
                        }
                        p1 = p3; f1 = y3; lastMoved = LEFT;
                    }
                    else
                    {
                        if (sideMoved == RIGHT)
                            f1 *= GetABFactor(y3, f2); // Apply Anderson-Björck factor to the left side
                        else if (!switching || Math.Abs(y3) <= PlateauRatio * Math.Abs(y2))
                            sideMoved = RIGHT; // A switching bisection on a plateau starts the AB phase afresh.

                        if (x3 != p2.X)
                        {
                            xr = p2.X; yr = p2.Y;
                        }
                        p2 = p3; f2 = y3; lastMoved = RIGHT;
                    }
                    dx = p2.X - p1.X;
                    if (dx > fallbackThreshold &&
                        42 * Math.Abs(y3) >= Math.Min(Math.Abs(y1), Math.Abs(y2)))
                    {
                        isAB = false;
                        sideMoved = 0;
                        ++fallbacks;
                    }
                    fallbackThreshold *= 0.5;
                }
                else if (SameSign(p1.Y, y3)) // Bisection step
                {
                    xl = p1.X; yl = p1.Y;
                    p1 = p3; lastMoved = LEFT;
                }
                else
                {
                    xr = p2.X; yr = p2.Y;
                    p2 = p3; lastMoved = RIGHT;
                }
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }

        // Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
        // (x1, f1), (x2, f2) and a third point (x3, f3) outside it, or NaN when it is unusable.
        // The model fits both nearly-quadratic shapes and functions with a pole; for a straight
        // line it gives the secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative
        // to f2, q1 = f1/f2 and q3 = f3/f2 (so the step does not depend on the scale of f), and
        // c = 1 - q1, the root is
        //     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
        // two divisions per step, 1/f2 and this one. The step is used only if it moves the secant
        // point xs = x1 - w*q1/c by at most half the distance to the nearer endpoint (on a plateau
        // the model degenerates towards a step function). Multiplied by c, that test reads
        //     |t*c + w*q1| <= w*min(|q1|, 1)/2,
        // with no division; it also rejects NaN and infinite values.
        // When the step is unusable, including a NaN third point (no endpoint dropped on that
        // side yet), it returns the A&B secant point through (x1, s1), (x2, s2) instead.
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double RationalStep(double x1, double y1, double x2, double y2, double x3, double y3,
            double s1, double s2)
        {
            var k = 1.0 / y2;
            double q1 = y1 * k, q3 = y3 * k;
            var c = 1.0 - q1;
            double w = x2 - x1, u = x3 - x1;
            var t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c);
            var a1 = Math.Abs(q1);
            var x = x1 + t;
            if (Math.Abs(t * c + w * q1) <= 0.5 * w * (a1 < 1.0 ? a1 : 1.0) && x <= x2)
                return x;

            return Node.SafeSecant(x1, s1, x2, s2);
        }

        // True if the second divided differences on the left (xl, x1, x2) and on the right
        // (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
        // one-signed curvature does worse than the secant. Scaled by positive factors, they are
        //     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
        //     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
        // so only their signs are compared, with no division or product of the two.
        // False when xl or xr is NaN (no endpoint dropped on that side yet).
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static bool CurvatureChangesSign(double xl, double yl, double x1, double y1,
            double x2, double y2, double xr, double yr)
        {
            double dy = y2 - y1, w = x2 - x1;
            var pl = dy * (x1 - xl) - (y1 - yl) * w;
            var pr = (yr - y2) * w - dy * (xr - x2);
            return (pl < 0.0 && pr > 0.0) || (pl > 0.0 && pr < 0.0);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double GetABFactor(double y3, double yMoved)
        {
            var m = 1.0 - y3 / yMoved;
            return m > 0.0 ? m : 0.5;
        }
    }
}