namespace Proektsoft.Root
{
    /// <summary>
    /// Stores one abscissa/ordinate pair used by the bracketing solvers.
    /// </summary>
    internal struct Node
    {
        public double X;
        public double Y;

        public Node(double x, double y)
        {
            X = x;
            Y = y;
        }

        public Node(double x, Func<double, double> F)
        {
            X = x;
            Y = F(x);
        }

        private const double ScaleThreshold = double.MaxValue / 4.0;

        /// <summary>
        /// Computes a safeguarded false-position/secant point for an ordered
        /// bracket whose ordinates have opposite non-zero signs.
        /// </summary>
        /// <remarks>
        /// SHARED CORRECTION USED BY modAB3 AND THE OTHER BRACKETING METHODS.
        ///
        /// The textbook formula
        ///
        ///     (x1*y2 - x2*y1) / (y2 - y1)
        ///
        /// may overflow although the mathematical intersection is finite.
        /// With opposite-sign ordinates, writing a=|y1| and b=|y2| gives the
        /// equivalent convex combination
        ///
        ///     x = (b/(a+b))*x1 + (a/(a+b))*x2.
        ///
        /// The previous helper protected the residual arithmetic but did not
        /// verify the final weighted sum. In an extreme binary64 case the two
        /// finite weighted products can still round to a sum of +/-Infinity.
        /// This helper now owns the complete postcondition required by every
        /// caller:
        ///
        ///     the returned point is finite and belongs to [p1.X,p2.X].
        ///
        /// If that postcondition cannot be obtained from the secant geometry,
        /// the safe midpoint is returned. Centralising the safeguard here lets
        /// the numerical methods stay short and prevents different solvers from
        /// implementing slightly different local clamps.
        /// </remarks>
        public static double SafeSecant(Node p1, Node p2)
        {
            var a = Math.Abs(p1.Y);
            var b = Math.Abs(p2.Y);
            var scale = Math.Max(a, b);

            // A zero, NaN, or infinite working magnitude does not provide
            // reliable secant weights. It does NOT indicate a root because the
            // ordinates may be Anderson-Björck auxiliary values. Bisection is
            // the safe and mathematically neutral fallback.
            if (a == 0.0 || b == 0.0 || !double.IsFinite(scale))
                return SafeMidpoint(p1, p2);

            // If the magnitudes are large, normalise before forming a+b.
            // For smaller values, preserving the original arithmetic path keeps
            // ordinary finite cases bitwise identical to the previous version.
            if (scale >= ScaleThreshold)
            {
                a /= scale;
                b /= scale;
            }

            var denominator = a + b;
            if (!(denominator > 0.0) || !double.IsFinite(denominator))
                return SafeMidpoint(p1, p2);

            var w1 = b / denominator;
            var w2 = a / denominator;
            var x = w1 * p1.X + w2 * p2.X;

            // Protect the final weighted sum itself. This also
            // catches NaN, which ordinary comparisons would silently miss.
            if (!double.IsFinite(x))
                return SafeMidpoint(p1, p2);

            // In exact arithmetic the convex combination is strictly inside
            // the bracket. The projection only corrects a possible last-bit
            // excursion caused by rounding. The solvers therefore receive the
            // precise finite-precision analogue of the interiority invariant.
            return x < p1.X ? p1.X :
                   x > p2.X ? p2.X :
                   x;
        }

        /// <summary>
        /// Computes the midpoint without first forming p1.X+p2.X, which can
        /// overflow for finite endpoints of the same sign.
        /// </summary>
        /// <remarks>
        /// For finite ordered endpoints the result lies in the closed bracket,
        /// so an additional clamp in every caller is unnecessary.
        /// </remarks>
        public static double SafeMidpoint(Node p1, Node p2) =>
            0.5 * p1.X + 0.5 * p2.X;
    }
}
