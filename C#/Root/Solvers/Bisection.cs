namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the bisection method
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double Bisection(Func<double, double> F, double x1, double x2,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = Node.Mid(p1, p2);
                var eps2 = aTol + rTol * Math.Abs(x3);
                // Check for X-convergence and return the result as a secant step
                if (p2.X - p1.X <= eps2)
                {
                    EvaluationCount = i + 1;
                    return Node.Sec(p1, p2);
                }
                Node p3 = new(x3, F);
                // Check for Y-convergence and return the result
                if (p3.Y == 0.0)
                {
                    EvaluationCount = i + 2;
                    return x3;
                }
                if (Math.Sign(p1.Y) == Math.Sign(p3.Y))
                    p1 = p3;
                else
                    p2 = p3;
            }
            EvaluationCount = MaxIterations + 2;
            return double.NaN;
        }
    }
}
