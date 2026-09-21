using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{
    public static partial class Solver
    {
        /// <summary>
        /// Finds a zero of <paramref name="f"/> in
        /// [<paramref name="x1"/>, <paramref name="x2"/>] by the
        /// simplified corrected improved modified Anderson-Björck method.
        /// </summary>
        /// <remarks>
        /// Exact-arithmetic assumptions used by the convergence proof:
        /// <list type="bullet">
        /// <item><description>f is continuous on the initial interval;</description></item>
        /// <item><description>the endpoint residuals are finite and have opposite signs;</description></item>
        /// <item><description>aTol and rTol are finite and non-negative.</description></item>
        /// </list>
        ///
        /// The implementation is slightly more permissive than the theorem:
        /// an infinite residual may be retained temporarily. SafeSecant then
        /// returns a midpoint, so the method continues with bisection geometry.
        /// NaN residuals are rejected because they have no usable sign.
        /// </remarks>
        public static double ModABCorr(
            Func<double, double> f,
            double x1, double x2,
            out ReturnCode returnCode,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(f, x1, x2, aTol, rTol, out Node p1, out Node p2, out var F))
            {
                returnCode = ReturnCode.Invalid;
                return double.NaN;
            }
            returnCode = ReturnCode.Success;

            // The true endpoint residuals.
            double y1 = p1.Y;
            if (y1 == 0.0)
                return p1.X;

            double y2 = p2.Y;
            if (y2 == 0.0)
                return p2.X;

            var isAB = false;
            // the endpoint replaced by the preceding AB step:
            //   +1 = left endpoint moved;
            //   -1 = right endpoint moved;
            //    0 = no preceding AB step in the current AB phase.
            var side = 0;
            var fallbackThreshold = 0.0;

            // Consecutive AB steps that failed the width test but were kept
            // because they halved the best residual. Capping them preserves
            // the dyadic worst-case bound: near a root of multiplicity m,
            // halving |f| shrinks the distance only by 2^(-1/m).
            double yMin = 0.0;
            const int MaxResidualSteps = 3;
            var residualSteps = 0;
            const double C = 2.0;

            for (var i = 1; i <= MaxIterations; ++i)
            {
                var x3 = isAB
                    ? Node.SafeSecant(p1, p2)
                    : Node.SafeMidpoint(p1, p2);

                var xTol = aTol + rTol * Math.Abs(x3);
                if (p2.X - p1.X <= xTol)
                    return x3;

                // If roundoff makes the proposal coincide with an endpoint,
                // reuse the true residual already stored there.
                var y3 = x3 == p1.X ? y1 :
                         x3 == p2.X ? y2 :
                         F(x3);

                // Exact zero is checked before NaN because it is more frequent
                if (y3 == 0.0)
                    return x3;

                if (double.IsNaN(y3))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                var switchToAB = false;
                // The switching controller is evaluated only during a genuine
                // bisection step and only from true residuals. This prevents the
                // controller from analysing a surrogate created by previous AB
                // corrections. Non-finite values deliberately disable switching:
                // PassesSwitchingTest rejects a non-finite ordinate itself, and a
                // NaN symmetry factor fails its comparison, so neither needs a
                // guard here.
                if (isAB)
                    yMin = Math.Min(Math.Abs(y1), Math.Abs(y2));
                else
                {
                    var ym = 0.5 * y1 + 0.5 * y2; // avoids overflow in y1 + y2.
                    switchToAB = PassesSwitchingTest(ym, y3, EvaluateSymmetryFactor(y1, y2));
                }
                // Best true residual of the bracket before y3 replaces an endpoint.
                var p3 = new Node(x3, y3);
                if (SameNonzeroSign(y1, y3))
                {
                    if (side == 1)
                        p2.Y = ScalePreservingNonzeroSign(p2.Y, GetABFactor(y3, y1));
                    else if (isAB)
                        side = 1;

                    p1 = p3;
                    y1 = y3;
                }
                else
                {
                    if (side == -1)
                        p1.Y = ScalePreservingNonzeroSign(p1.Y, GetABFactor(y3, y2));
                    else if (isAB)
                        side = -1;

                    p2 = p3;
                    y2 = y3;
                }
                if (isAB)
                {
                    // Fallback if AB fails to reduce the bracket width, unless it still 
                    // halves the best residual (at most MaxResidualSteps times in a row)
                    var fallBack = false;
                    if (p2.X - p1.X > fallbackThreshold)
                    {
                        if (Math.Abs(y3) < 0.5 * yMin && residualSteps < MaxResidualSteps)
                            ++residualSteps;
                        else
                            fallBack = true;
                    }
                    else
                        residualSteps = 0;

                    if (fallBack)
                    {
                        // AB has fallen behind the bisection reference. Return
                        // to bisection and discard phase-dependent interpolation
                        // corrections by restoring the true endpoint residuals.
                        isAB = false;
                        side = 0;
                        p1.Y = y1;
                        p2.Y = y2;
                    }
                    else
                    {
                        // Prepare the next dyadic threshold exactly in binary
                        // arithmetic until the subnormal range is reached.
                        fallbackThreshold *= 0.5;
                    }
                }
                else if (switchToAB)
                {
                    // At the first AB phase p1.Y == y1 and p2.Y == y2. After every
                    // fallback they are restored, and bisection never modifies them. 
                    // Hence a redundant reset is not needed here.
                    isAB = true;
                    residualSteps = 0;
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

        /// <summary>
        /// Returns k=r^2 for the symmetry-sensitive switching criterion.
        /// The calculation is homogeneous in the true endpoint residuals.
        /// </summary>
        private static double EvaluateSymmetryFactor(double y1, double y2)
        { 
            var a = Math.Abs(y1);
            var b = Math.Abs(y2);
            var den = a + b;

            if (double.IsInfinity(den))
            {
                // Infinite true residuals deliberately disable AB switching and
                // keep the controller in bisection mode. NaN is returned rather
                // than an infinity because every exit of PassesSwitchingTest is
                // a "<" comparison, which is false against NaN; an infinity
                // would instead satisfy it and switch. Residuals are never zero
                // here, so only an overflowing sum remains, and halving both
                // restores it without changing the ratio.
                if (double.IsInfinity(a) || double.IsInfinity(b))
                    return double.NaN;

                a *= 0.5;
                b *= 0.5;
                den = a + b;
            }

            // |b-a| <= den, so the quotient lies in [0,1] and halving it after
            // the division avoids forming 2*den, which could overflow.
            var r = 1.0 - Math.Abs(b - a) / den / 2.0;
            return r * r;
        }

        /// <summary>
        /// Tests whether the true midpoint value is sufficiently close to the
        /// midpoint value of the chord through the true endpoint residuals.
        /// </summary>
        private static bool PassesSwitchingTest(double ym, double yf, double symmetryFactor)
        {
            var absYm = Math.Abs(ym);
            var absYf = Math.Abs(yf);
            var sum = absYf + absYm;

            // Fast path. The exact-root case was handled before this method was
            // called, and a non-finite ordinate fails the comparison, which
            // disables switching as intended.
            if (double.IsFinite(sum))
                return Math.Abs(ym - yf) < symmetryFactor * sum;

            // Only reached when the sum overflows. Non-finite values are
            // unsuitable for the linearity comparison.
            if (!double.IsFinite(ym) || !double.IsFinite(yf))
                return false;

            // Normalize both sides of the homogeneous inequality to avoid
            // overflow in subtraction or addition.
            var scale = Math.Max(absYf, absYm);
            var normYm = ym / scale;
            var normYf = yf / scale;
            return Math.Abs(normYm - normYf) < symmetryFactor * (Math.Abs(normYf) + Math.Abs(normYm));
        }

        /// <summary>
        /// Multiplies an auxiliary AB ordinate by a positive factor while
        /// preserving a finite non-zero sign in binary64 arithmetic.
        /// </summary>
        /// <remarks>
        /// This helper is retained in the reference simplified version because
        /// it preserves the finite-precision analogue of the sign invariant and
        /// keeps ordinary regression trajectories unchanged. It acts only on
        /// auxiliary ordinates; an underflowed working value is never accepted
        /// as a root of f.
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double ScalePreservingNonzeroSign(double value, double positiveFactor)
        {
            var scaled = value * positiveFactor;

            if (scaled == 0.0 && value != 0.0)
                return Math.CopySign(double.Epsilon, value);

            if (double.IsInfinity(scaled))
                return Math.CopySign(double.MaxValue, value);

            return scaled;
        }
    }
}
