using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{
    public static partial class Solver
    {
        /// <summary>
        /// Finds a zero of <paramref name="f"/> in
        /// [<paramref name="x1"/>, <paramref name="x2"/>] by the
        /// safeguarded improved modified Anderson-Björck method.
        /// L. Tomov and N. Ganchovski
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
            double x1,
            double x2,
            out ReturnCode returnCode,
            double aTol = 1e-14,
            double rTol = 1e-14)
        {
            if (!Initialize(
                    f, x1, x2, aTol, rTol,
                    out Node p1,
                    out Node p2,
                    out var F))
            {
                returnCode = ReturnCode.Invalid;
                return double.NaN;
            }

            // y1 and y2 are the TRUE endpoint residuals. In AB mode p1.Y and
            // p2.Y may be corrected auxiliary ordinates used only to construct
            // the next secant. Keeping these two roles separate is the central
            // post-publication correction and is required by the proof.
            double y1 = p1.Y;
            double y2 = p2.Y;

            returnCode = ReturnCode.Success;

            // Exact binary64 zero is a valid early exit for the computed
            // machine function. It is checked only on true function values,
            // never on an Anderson-Björck auxiliary ordinate.
            if (y1 == 0.0)
                return p1.X;

            if (y2 == 0.0)
                return p2.X;

            var isAB = false;

            // the endpoint replaced by the preceding AB step:
            //   +1 = left endpoint moved;
            //   -1 = right endpoint moved;
            //    0 = no preceding AB step in the current AB phase.
            var side = 0;

            // The simplified fallback controller is exactly the dyadic rule
            //
            //     T(t) = 32 * Wc / 2^t,
            //
            // where Wc is the width at entry to an AB phase. The first five
            // threshold tests are vacuous because T(t) >= Wc while the bracket
            // width is non-increasing. The first non-trivial test is therefore
            // T(6) = Wc/2. This is C=32, so the corrected convergence theorem
            // uses q=ceil(log2(32))+2=7.
            const int Log2C = 5;
            var uncheckedABSteps = 0;
            var fallbackThreshold = 0.0;

            for (var iteration = 1;
                 iteration <= MaxIterations;
                 ++iteration)
            {
                var x3 = isAB
                    ? Node.SafeSecant(p1, p2)
                    : Node.SafeMidpoint(p1, p2);

                // Form the relative tolerance only from a safeguarded in-bracket
                // point. This avoids false termination caused by a non-finite or
                // wildly out-of-bracket interpolation proposal.
                var xTol = aTol + rTol * Math.Abs(x3);

                // Width-based termination supplies a geometric certificate:
                // for a preserved sign-changing bracket and x3 in the bracket,
                // dist(x3,Z) <= p2.X-p1.X for the zero set Z.
                if (p2.X - p1.X <= xTol)
                    return x3;

                // If roundoff makes the proposal coincide with an endpoint,
                // reuse the TRUE residual already stored there. Reusing p1.Y or
                // p2.Y would be wrong in AB mode because those may be corrected
                // working ordinates rather than values of f.
                var y3 = x3 == p1.X ? y1 :
                         x3 == p2.X ? y2 :
                         F(x3);

                // Exact zero is checked before NaN because it is a frequent and
                // very cheap successful branch; NaN==0.0 is false, so the order
                // does not hide an invalid value.
                if (y3 == 0.0)
                    return x3;

                if (double.IsNaN(y3))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }

                // The switching controller is evaluated only during a genuine
                // bisection step and only from TRUE residuals. Infinite values
                // deliberately disable switching. This prevents the controller
                // from analysing a surrogate created by previous AB corrections.
                var switchToAB = false;
                if (!isAB && double.IsFinite(y3))
                {
                    var symmetryFactor = EvaluateSymmetryFactor(y1, y2);
                    if (double.IsFinite(symmetryFactor))
                    {
                        // This form avoids overflow in y1+y2.
                        var ym = 0.5 * y1 + 0.5 * y2;
                        switchToAB = PassesSwitchingTest(
                            ym,
                            y3,
                            symmetryFactor);
                    }
                }

                // Update the mathematical bracket using y1,y2,y3 only. The
                // corrected p1.Y and p2.Y values are interpolation state and
                // must never decide which endpoint is replaced.
                var p3 = new Node(x3, y3);
                if (SameNonzeroSign(y1, y3))
                {
                    // The new point has the sign of the left endpoint, so the
                    // left endpoint moves. If this happened on the preceding AB
                    // step as well, the right endpoint has remained fixed twice
                    // and receives the Anderson-Björck correction.
                    if (side == 1)
                    {
                        p2.Y = ScalePreservingNonzeroSign(
                            p2.Y,
                            GetABFactor(y3, p1.Y));
                    }
                    else if (isAB)
                    {
                        side = 1;
                    }

                    p1 = p3;
                    y1 = y3;
                }
                else
                {
                    // Since the true endpoint residuals have opposite non-zero
                    // signs and y3 is non-zero, this branch means that y3 has the
                    // sign of the right endpoint.
                    if (side == -1)
                    {
                        p1.Y = ScalePreservingNonzeroSign(
                            p1.Y,
                            GetABFactor(y3, p2.Y));
                    }
                    else if (isAB)
                    {
                        side = -1;
                    }

                    p2 = p3;
                    y2 = y3;
                }

                if (isAB)
                {
                    // Skip the first Log2C=5 vacuous threshold tests. During
                    // them T(t)>=Wc and the current width cannot exceed Wc.
                    if (uncheckedABSteps > 0)
                    {
                        --uncheckedABSteps;
                    }
                    else if (p2.X - p1.X > fallbackThreshold)
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
                    // At the first AB phase p1.Y==y1 and p2.Y==y2. After every
                    // fallback they are restored, and bisection never modifies
                    // them. Hence a redundant reset is not needed here.
                    isAB = true;
                    uncheckedABSteps = Log2C;

                    // First non-trivial threshold for C=32:
                    // T(Log2C+1)=Wc/2, with Wc measured AFTER the switching
                    // bisection has updated the bracket.
                    fallbackThreshold = 0.5 * (p2.X - p1.X);
                }
            }

            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double GetABFactor(double y3, double movingWorkingY)
        {
            var m = 1.0 - y3 / movingWorkingY;

            // The proof needs only a finite positive correction factor. If the
            // standard factor is non-positive or non-finite, 1/2 preserves the
            // working ordinate's sign and is the conventional safe fallback.
            return double.IsFinite(m) && m > 0.0 ? m : 0.5;
        }

        private const double ScaleThreshold = double.MaxValue / 4.0;

        /// <summary>
        /// Returns k=r^2 for the symmetry-sensitive switching criterion.
        /// The calculation is homogeneous in the true endpoint residuals.
        /// </summary>
        private static double EvaluateSymmetryFactor(double y1, double y2)
        {
            var a = Math.Abs(y1);
            var b = Math.Abs(y2);
            var scale = Math.Max(a, b);

            // Infinite true residuals deliberately disable AB switching and
            // keep the controller in bisection mode.
            if (!(scale > 0.0) || !double.IsFinite(scale))
                return double.PositiveInfinity;

            if (scale >= ScaleThreshold)
            {
                a /= scale;
                b /= scale;
            }

            var r = 1.0 - Math.Abs(b - a) / (2.0 * (a + b));
            return r * r;
        }

        /// <summary>
        /// Tests whether the true midpoint value is sufficiently close to the
        /// midpoint value of the chord through the true endpoint residuals.
        /// </summary>
        private static bool PassesSwitchingTest(
            double ym,
            double yf,
            double symmetryFactor)
        {
            var absYm = Math.Abs(ym);
            var absYf = Math.Abs(yf);
            var scale = Math.Max(absYf, absYm);

            // The exact-root case was handled before this method was called.
            // Non-finite values are unsuitable for the linearity comparison.
            if (!(scale > 0.0) || !double.IsFinite(scale))
                return false;

            if (scale < ScaleThreshold)
            {
                return Math.Abs(ym - yf) <
                       symmetryFactor * (absYf + absYm);
            }

            // Normalize both sides of the homogeneous inequality to avoid
            // overflow in subtraction or addition.
            var normYm = ym / scale;
            var normYf = yf / scale;
            return Math.Abs(normYm - normYf) <
                   symmetryFactor *
                   (Math.Abs(normYf) + Math.Abs(normYm));
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
        private static double ScalePreservingNonzeroSign(
            double value,
            double positiveFactor)
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
