namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "f(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using Anderson-Bjork's method:
        // Anderson, N., Björck, Å. A new high order method of regula-falsi type
        // for computing a root of an equation.
        // BIT 13, 253–264 (1973). 
        // https://doi.org/10.1007/BF01951936
        // f(x) must be continuous and sign(f(x1)) ≠ sign(f(x2))

        public static double AndersonBjork(Func<double, double> f, 
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
            int side = 0;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = Node.SafeSecant(p1, p2);
                var xTol = aTol + rTol * Math.Abs(x3);
                if (p2.X - p1.X <= xTol)
                    return x3;

                Node p3 = new(x3, F);
                if (p3.Y == 0)
                    return x3;

                if (double.IsNaN(p3.Y))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                if (SameSign(p1.Y, p3.Y))
                {
                    if (side == 1)
                    {
                        double m = 1 - p3.Y / p1.Y;
                        p2.Y *= m <= 0 ? 0.5 : m;
                    }
                    else
                        side = 1;
                    
                    p1 = p3;
                }
                else
                {
                    if (side == -1)
                    {
                        double m = 1 - p3.Y / p2.Y;
                        p1.Y *= m <= 0 ? 0.5 : m;
                    }
                    else
                        side = -1;
                    
                    p2 = p3;
                }
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }
    }
}
