namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "f(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the bisection method
        // f(x) must be continuous and sign(f(x1)) ≠ sign(f(x2))

        public static double Bisection(Func<double, double> f, 
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
            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = Node.SafeMidpoint(p1, p2);
                var xTol = aTol + rTol * Math.Abs(x3);
                // Check for X-convergence and return the result as a secant step
                if (p2.X - p1.X <= xTol)
                    return Node.SafeSecant(p1, p2);
                
                Node p3 = new(x3, F);
                // Check for Y-convergence and return the result
                if (p3.Y == 0.0)
                    return x3;

                if (double.IsNaN(p3.Y))
                {
                    returnCode = ReturnCode.Invalid;
                    return double.NaN;
                }
                if (SameNonzeroSign(p1.Y, p3.Y))
                    p1 = p3;
                else
                    p2 = p3;
            }
            returnCode = ReturnCode.MaxIterationsExceeded;
            return double.NaN;
        }
    }
}
