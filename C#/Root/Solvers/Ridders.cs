namespace Proektsoft.Root
{
    public static partial class Solver
    {
        /// <summary>
        /// Finds a zero of <paramref name="f"/> in
        /// [<paramref name="x1"/>, <paramref name="x2"/>] using
        /// Ridder's method.
        /// </summary>
        /// <remarks>
        /// C. Ridders, "A new algorithm for computing a single root of a
        /// real continuous function", IEEE Transactions on Circuits and
        /// Systems, 26(11), 979--980, 1979.
        ///
        /// The function must be continuous and the initial endpoint values must
        /// have opposite signs. This implementation terminates either at an
        /// exact binary64 zero or when the preserved bracket satisfies the
        /// requested mixed absolute-relative X tolerance.
        /// </remarks>
        public static double Ridders(
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

            returnCode = ReturnCode.Success;

            if (p1.Y == 0.0)
                return p1.X;

            if (p2.Y == 0.0)
                return p2.X;

            for (var iteration = 1;
                 iteration <= MaxIterations;
                 ++iteration)
            {
                // Ridder's construction first evaluates the midpoint.
                var p3 = new Node(Node.SafeMidpoint(p1, p2), F);

                // Exact binary64 zero is a valid result for the actually
                // evaluated machine function.
                if (p3.Y == 0.0)
                    return p3.X;

                // Reject NaN BEFORE it enters products, square
                // roots, Math.Sign, or bracket updates. Checking only the later
                // p4 value is too late because the proposal may already be NaN.
                if (double.IsNaN(p3.Y))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }

                // Form the Ridders proposal through a normalized helper. The raw expression
                //
                //     p3.Y*p3.Y - p1.Y*p2.Y
                //
                // may overflow even when the normalized Ridders ratio is well
                // defined. SafeRiddersPoint normalizes the residuals, checks the
                // radicand, protects the final arithmetic, and guarantees a
                // finite point in the current bracket. If interpolation is not
                // usable, it returns the already evaluated midpoint p3.X.
                var x4 = SafeRiddersPoint(p1, p2, p3);

                var xTol = aTol + rTol * Math.Abs(x4);

                // Use the whole bracket width rather than the
                // distance between two consecutive Ridders points. Under
                // continuity and preservation of the sign-changing bracket,
                // any returned in-bracket point then satisfies
                //
                //     dist(x4,Z) <= p2.X-p1.X.
                //
                // No new function evaluation is needed once this certificate
                // already meets the requested X tolerance.
                if (p2.X - p1.X <= xTol)
                    return x4;

                // Reuse every value already known. In particular,
                // when the safeguarded proposal is the midpoint, do not evaluate
                // f(p3.X) a second time. The endpoint cases can occur through a
                // final rounding/clamp and likewise reuse their true residuals.
                Node p4;
                if (x4 == p3.X)
                    p4 = p3;
                else if (x4 == p1.X)
                    p4 = p1;
                else if (x4 == p2.X)
                    p4 = p2;
                else
                    p4 = new Node(x4, F);

                if (p4.Y == 0.0)
                    return p4.X;

                // NaN has no usable sign. This check now occurs
                // only after the proposal itself has been made finite, so the
                // user's function is never called with a NaN abscissa.
                if (double.IsNaN(p4.Y))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }

                // Preserve an ORDERED sign-changing bracket. The location of
                // the Ridders point relative to the midpoint is linked to the
                // midpoint sign, so these two compact branches select the
                // tighter valid sub-bracket without reversing p1.X and p2.X.
                if (SameSign(p1.Y, p4.Y))
                {
                    // p4 has the sign of the left endpoint and therefore
                    // replaces it. If the midpoint has the sign of the old
                    // right endpoint, it provides a tighter new right endpoint.
                    p1 = p4;
                    if (SameSign(p2.Y, p3.Y))
                        p2 = p3;
                }
                else
                {
                    // Since p4 is nonzero and the old endpoints have opposite
                    // signs, p4 has the sign of the right endpoint.
                    p2 = p4;
                    if (SameSign(p1.Y, p3.Y))
                        p1 = p3;
                }
            }

            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }

        /// <summary>
        /// Returns a finite Ridders proposal in [p1.X,p2.X]. If the
        /// interpolation formula is unusable in binary64 arithmetic, the
        /// already evaluated midpoint p3.X is returned.
        /// </summary>
        private static double SafeRiddersPoint(
            Node p1,
            Node p2,
            Node p3)
        {
            var scale = Math.Max(
                Math.Abs(p3.Y),
                Math.Max(Math.Abs(p1.Y), Math.Abs(p2.Y)));

            // Infinite residuals still have signs and can be handled by
            // bisection, but they are unsuitable for the Ridders ratio. NaN at
            // p3 is rejected by the caller; this guard also protects against any
            // unexpected non-finite endpoint value.
            if (!(scale > 0.0) || !double.IsFinite(scale))
                return p3.X;

            // Normalisation is homogeneous: the common scale cancels from the
            // Ridders quotient. It keeps every residual magnitude in [0,1], so
            // the square and product below cannot overflow.
            var y1 = p1.Y / scale;
            var y2 = p2.Y / scale;
            var y3 = p3.Y / scale;

            var radicand = y3 * y3 - y1 * y2;
            if (!(radicand > 0.0) || !double.IsFinite(radicand))
                return p3.X;

            // Because p1.Y and p2.Y have opposite signs,
            // radicand >= y3^2 in exact arithmetic and |ratio|<=1. Thus the
            // exact Ridders point lies in the current bracket.
            var ratio = Math.CopySign(1.0, p1.Y) *
                        y3 /
                        Math.Sqrt(radicand);

            if (!double.IsFinite(ratio))
                return p3.X;

            var x = p3.X + (p3.X - p1.X) * ratio;
            if (!double.IsFinite(x))
                return p3.X;

            // In exact arithmetic the proposal is strictly interior. If a
            // last-bit error puts it on or beyond an endpoint, use the already
            // evaluated midpoint instead of returning an endpoint that might
            // leave the bracket unchanged. This guarantees progress whenever
            // the midpoint itself is representably distinct from the endpoints.
            if (x <= p1.X || x >= p2.X)
                return p3.X;

            return x;
        }
    }
}
