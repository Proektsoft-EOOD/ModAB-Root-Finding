namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "f(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the Illinois method
        // f(x) must be continuous and sign(f(x1)) ≠ sign(f(x2))

        public static double Illinois(Func<double, double> f, 
            double x1, double x2, 
            out ReturnCode returnCode,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(f, x1, x2, aTol, rTol, out Node p1, out Node p2, out var F))
            {
                returnCode = ReturnCode.Invalid;
                return double.NaN;
            }
            int side = 0;
            returnCode = ReturnCode.Success;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = Node.SafeSecant(p1, p2);
                var xTol = aTol + rTol * Math.Abs(x3);
                if (p2.X - p1.X <= xTol)
                    return x3;
                
                Node p3 = new(x3, F);
                if (p3.Y == 0)
                    return p3.X;
                
                if (double.IsNaN(p3.Y))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                if (SameSign(p1.Y, p3.Y))
                {
                    if (side == 1)
                        p2.Y /= 2;

                    side = 1;
                    p1 = p3;
                }
                else
                {
                    if (side == -1)
                        p1.Y /= 2;

                    side = -1;
                    p2 = p3;
                }
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }
    }
}