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
        internal bool IsRoot(double x, double tol = 1e-14)
        {
            for (int i = 0, len = Roots.Length; i < len; i++)
                if (Math.Abs(x - Roots[i]) <= tol)
                    return true;

            return F(x) == 0;
        }
    }
}
