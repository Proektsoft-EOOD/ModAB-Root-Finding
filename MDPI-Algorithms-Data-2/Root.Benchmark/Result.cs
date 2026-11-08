using Proektsoft.Root;
namespace Root.Benchmark
{
    internal struct Result
    {
        internal double Root;
        internal double FunctionValue;
        internal int EvaluationCount;
        internal Solver.ReturnCode ReturnCode; // validation/reporting classification
        internal Solver.ReturnCode SolverReturnCode; // unchanged result from the solver
    }
}
