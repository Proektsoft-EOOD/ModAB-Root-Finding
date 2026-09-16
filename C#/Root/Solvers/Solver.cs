namespace Proektsoft.Root
{
    // This class provides methods for solving the equation f(x) = 0 numerically.
    // All methods are of bracketing type. They require a continuous function and
    // an initial interval whose endpoint values have opposite signs.

    public static partial class Solver
    {
        public enum ReturnCode
        {
            Success = 0,
            Invalid = 1,
            MaxIterationsExceeded = 2,
            FalseConvergence = 3
        }

        public const int MaxIterations = 200;

        /// <summary>
        /// Number of function evaluations performed by the latest completed
        /// call. This remains a lightweight diagnostic property rather than a
        /// concurrency-safe result object.
        /// </summary>
        /// <remarks>
        /// CORRECTION: the function delegate itself is no longer stored in a
        /// mutable static field. Consequently, parallel or nested calls cannot
        /// replace the function being solved midway through another algorithm.
        /// EvaluationCount is still a shared diagnostic value; callers that need
        /// per-call counts under concurrency should return the count as part of a
        /// dedicated result structure in a future API revision.
        /// </remarks>
        public static int EvaluationCount { get; private set; }

        private static double EvaluateAndCount(
            Func<double, double> f,
            double x)
        {
            ++EvaluationCount;
            return f(x);
        }

        /// <summary>
        /// Returns true only when both values have the same non-zero sign.
        /// </summary>
        /// <remarks>
        /// The predicate is shared by all bracketing algorithms. Comparisons
        /// with NaN are false, so every solver must reject NaN before using this
        /// predicate to update a bracket.
        /// </remarks>
        private static bool SameNonzeroSign(double x, double y) =>
            (x < 0.0 && y < 0.0) ||
            (x > 0.0 && y > 0.0);

        /// <summary>
        /// Validates common inputs, orders the endpoints, evaluates the two
        /// endpoint residuals, and returns a LOCAL counting wrapper for f.
        /// </summary>
        /// <remarks>
        /// CORRECTION TO THE PREVIOUS SHARED WRAPPER.
        ///
        /// The old implementation stored the active function in
        ///
        ///     private static Func&lt;double,double&gt; F;
        ///
        /// so another or nested solver call could overwrite F while the first
        /// call was still running. The wrapper is now returned through the out
        /// parameter F and is therefore local to the current invocation.
        ///
        /// The null check is deliberately applied to the user delegate f before
        /// the wrapper is created. Checking the wrapper itself would never catch
        /// a null f because a non-null lambda can close over a null reference.
        /// </remarks>
        private static bool Initialize(
            Func<double, double> f,
            double x1,
            double x2,
            double aTol,
            double rTol,
            out Node p1,
            out Node p2,
            out Func<double, double> F)
        {
            EvaluationCount = 0;
            ArgumentNullException.ThrowIfNull(f);

            // The counting wrapper belongs to this one solver call. Numerical
            // results can no longer be corrupted by a different call changing a
            // shared static function reference.
            F = x => EvaluateAndCount(f, x);

            if (!double.IsFinite(aTol) || aTol < 0.0)
                throw new ArgumentOutOfRangeException(
                    nameof(aTol),
                    "The absolute tolerance must be finite and non-negative.");

            if (!double.IsFinite(rTol) || rTol < 0.0)
                throw new ArgumentOutOfRangeException(
                    nameof(rTol),
                    "The relative tolerance must be finite and non-negative.");

            if (double.IsNaN(x1))
                throw new ArgumentOutOfRangeException(
                    nameof(x1),
                    "The left endpoint must not be NaN.");

            if (double.IsNaN(x2))
                throw new ArgumentOutOfRangeException(
                    nameof(x2),
                    "The right endpoint must not be NaN.");

            if (double.IsInfinity(x1) || double.IsInfinity(x2))
            {
                p1 = p2 = default;
                return false;
            }

            if (x1 > x2)
                (x1, x2) = (x2, x1);

            p1 = new Node(x1, F);
            p2 = new Node(x2, F);

            // NaN has no usable sign. Infinite residuals are allowed here: a
            // bisection-capable method can often shrink the bracket until finite
            // values are reached. Each solver decides how aggressively it uses
            // interpolation while such values remain present.
            return !(double.IsNaN(p1.Y) ||
                     double.IsNaN(p2.Y) ||
                     SameNonzeroSign(p1.Y, p2.Y));
        }
    }
}
