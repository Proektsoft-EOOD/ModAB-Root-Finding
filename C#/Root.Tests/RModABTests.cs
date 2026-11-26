using Proektsoft.Root;

namespace Root.Tests
{
    // Solver.EvaluationCount is a shared static, so the tests in this class must not run in
    // parallel with other tests that call the solvers (xUnit runs the methods of one class in
    // sequence).
    public class RModABTests
    {
        // A skewed monotone function with its root at 0.3: AB steps reduce the residual more than
        // fourfold but less than 42 times without shrinking the bracket enough. They should not
        // fall back to bisection straight away. Ported from the NonlinearSolve.jl ModAB test
        // "Anderson-Bjork keeps a few residual-only steps while the bracket lags".
        private static double Skewed(double x)
        {
            var q = x - 0.3;
            return q < 0.0 ? 0.001 * q : q * (1.0 + 10000.0 * q * q * q * q);
        }

        // The Julia test runs with maxiters = 20: 2 endpoint evaluations and at most 19 loop
        // iterations, one evaluation each. Solver.MaxIterations is a constant, so the C# test
        // checks the same budget on the evaluation count instead.
        private const int EvaluationBudget = 21;

        [Theory]
        [InlineData(1e-6)]  // Float32 tolerance of the Julia test
        [InlineData(1e-12)] // Float64 tolerance of the Julia test
        public void KeepsResidualOnlyStepsWhileBracketLags(double aTol)
        {
            var root = Solver.RModAB(Skewed, 0.0, 1.0, out var returnCode, aTol, 0.0);
            var evaluations = Solver.EvaluationCount;

            Assert.Equal(Solver.ReturnCode.Success, returnCode);
            Assert.True(Math.Abs(root - 0.3) <= aTol, $"|{root:R} - 0.3| > {aTol}");
            Assert.True(evaluations <= EvaluationBudget,
                $"{evaluations} evaluations, budget {EvaluationBudget}");
        }
    }
}
