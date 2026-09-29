namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "f(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the ITP method:
        // Oliveira, I. F. D. and Ricardo H. C. Takahashi.
        // “An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality.”
        // ACM Transactions on Mathematical Software (TOMS) 47 (2021): 1 – 24.
        // https://doi.org/10.1145/3423597
        // f(x) must be continuous and sign(f(x1)) ≠ sign(f(x2))

        public static double ITP(Func<double, double> f, 
            double x1, double x2, 
            out ReturnCode returnCode,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(f, x1, x2, aTol, rTol, out Node p1, out Node p2, out var F))
            {
                returnCode = ReturnCode.Invalid;
                return double.NaN;
            }
            double span = p2.X - p1.X;
            int n0 = 1;
            double k1 = 0.2 / span, k2 = 2d;
            double eps2 = aTol + rTol * span;
            int nb = (int)Math.Ceiling(Math.Log2(span / eps2));
            int nmax = nb + n0;
            returnCode = ReturnCode.Success;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                double xb = Node.SafeMidpoint(p1, p2);
                span = p2.X - p1.X;
                if (span < aTol + rTol * span)
                    return xb;
                
                // Interpolate
                double xf = Node.SafeSecant(p1, p2);
                if (double.IsNaN(xf))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                // Truncate
                double σ = Math.Sign(xb - xf);
                double δ = k1 * Math.Pow(span, k2);
                double xt = δ <= Math.Abs(xb - xf) ? xf + σ * δ : xb;
                // Project
                double r = Math.Max(0, eps2 * Math.Pow(2d, nmax - i) - span / 2d);
                double x = Math.Abs(xt - xb) <= r ? xt : xb - σ * r;
                // Update
                Node p = new(x, F);
                if (SameSign(p.Y, p1.Y))
                    p1 = p;
                else if (SameSign(p.Y, p2.Y))
                    p2 = p;
                else
                    return x;
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }
    }
}
