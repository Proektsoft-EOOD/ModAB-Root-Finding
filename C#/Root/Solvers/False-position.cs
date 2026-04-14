namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using the false-position (regula-falsi) method
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))

        public static double FalsePosition(Func<double, double> F, double x1, double x2,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            var x0 = p1.X;
            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = Node.Sec(p1, p2);
                var eps = (aTol + rTol * Math.Abs(x3)) / 2.0;
                if (Math.Abs(x3 - x0) <= eps)    
                {
                    EvaluationCount = i + 1;
                    return x3;
                }
                x0 = x3;
                Node p3 = new(x3, F);
                if (p3.Y == 0)
                {
                    EvaluationCount = i + 2;
                    return p3.X;
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
