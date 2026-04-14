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
    }
}
