using Proektsoft.Root;
using BenchmarkDotNet.Attributes;
namespace Root.Benchmark
{
    public class BenchmarkTime
    {
        const double tol = 1e-14;
        private static Problem[] _problems = (
            new[] {
                BenchmarkCount.problems1,
                BenchmarkCount.problems2,
                BenchmarkCount.problems3 
            }).SelectMany(x => x).ToArray();

        [Benchmark]
        public void Bisection()
        {
            foreach (Problem p in _problems)
                Solver.Bisection(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void FalsePosition()
        {
            foreach (Problem p in _problems)
                Solver.FalsePosition(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void Illinois()
        {
            foreach (Problem p in _problems)
                Solver.Illinois(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void AndersonBjork()
        {
            foreach (Problem p in _problems)
                Solver.AndersonBjork(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void ITP()
        {
            foreach (Problem p in _problems)
                Solver.ITP(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void Ridders()
        {
            foreach (Problem p in _problems)
                Solver.Ridders(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void Brent()
        {
            foreach (Problem p in _problems)
                Solver.Brent(p.F, p.a, p.b, tol, tol);
        }

        [Benchmark]
        public void ModAB()
        {
            foreach (Problem p in _problems) 
                Solver.ModAB(p.F, p.a, p.b, tol, tol);
        }
    }
}
