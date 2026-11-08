namespace Root.Benchmark
{
    // This class contains the required data
    // to define a test problem for the numerical library
    internal struct Problem
    {
        internal string Name;
        internal Func<double, double> F;
        internal double a;
        internal double b;
        internal double[] Roots;

        // True if at least one finite reference root lies in the initial interval.
        // Cases without such a root are not reference-eligible.
        internal bool HasReference
        {
            get
            {
                if (Roots is null || !double.IsFinite(a) || !double.IsFinite(b))
                    return false;

                double lo = Math.Min(a, b), hi = Math.Max(a, b);
                foreach (var z in Roots)
                    if (double.IsFinite(z) && z >= lo && z <= hi)
                        return true;

                return false;
            }
        }

        // Independent validation of a returned approximation x (MDPI Algorithms
        // revision, Section 7.1, eq. (54)):
        //     min |x - z| <= aTol + rTol*|x|  over the reference roots z in [a, b].
        // x must be finite and lie in [a, b]. A computed zero F(x) == 0 is NOT accepted
        // as a substitute for this test, and the function is never evaluated here, so
        // validation cannot change the evaluation count. The default rTol = 0 keeps
        // legacy two-argument calls absolute-only. Invalid tolerances, a null or empty
        // reference list and non-finite tolerance arithmetic are rejected.
        internal bool IsRoot(double x, double aTol = 1e-14, double rTol = 0.0)
        {
            if (!double.IsFinite(aTol) || aTol < 0.0 || !double.IsFinite(rTol) || rTol < 0.0)
                return false;

            if (Roots is null || !double.IsFinite(a) || !double.IsFinite(b) || !double.IsFinite(x))
                return false;

            double lo = Math.Min(a, b), hi = Math.Max(a, b);
            if (x < lo || x > hi)
                return false;

            var tol = aTol + rTol * Math.Abs(x);
            if (!double.IsFinite(tol))
                return false;

            foreach (var z in Roots)
                if (double.IsFinite(z) && z >= lo && z <= hi && Math.Abs(x - z) <= tol)
                    return true;

            return false;
        }
    }
}
