using Proektsoft.Root;
using System.Diagnostics;
namespace Root.Benchmark
{
    internal class BenchmarkCount
    {
        internal static void Run()
        {
            Problem[] problems =
            [.. BenchmarkProblems.Set1, .. BenchmarkProblems.Set2, .. BenchmarkProblems.Set3];
            const double tol = 1e-14;
            const int methodCount = 9;
            var problemCount = problems.Length;
            var results = new Result[problemCount, methodCount];
            var sumAll = new int[methodCount];
            var sumValid = new int[methodCount];
            var average = new double[methodCount];
            var logSum = new double[methodCount];
            var mean = new double[methodCount];
            var variance = new double[methodCount];
            var stddev = new double[methodCount];
            var median = new double[methodCount];
            var max = new int[methodCount];
            var bestAt = new int[methodCount];
            var worstAt = new int[methodCount];
            var successCount = new int[methodCount];
            var invalidCount = new int[methodCount];
            var falseConvergenceCount = new int[methodCount];
            var maxCount = new int[methodCount];
            for (int i = 0; i < problemCount; ++i)
            {
                var p = problems[i];
                Solver.ReturnCode returnCode;
                for (int j = 0; j < methodCount; ++j)
                {
                    int evals = 0;
                    double Counted(double x) { ++evals; return p.F(x); }
                    var root = j switch
                    {
                        0 => Solver.Bisection(Counted, p.a, p.b, out returnCode, tol, tol),
                        1 => Solver.FalsePosition(Counted, p.a, p.b, out returnCode, tol, tol),
                        2 => Solver.Illinois(Counted, p.a, p.b, out returnCode, tol, tol),
                        3 => Solver.AndersonBjork(Counted, p.a, p.b, out returnCode, tol, tol),
                        4 => Solver.ITP(Counted, p.a, p.b, out returnCode, tol, tol),
                        5 => Solver.Ridders(Counted, p.a, p.b, out returnCode, tol, tol),
                        6 => Solver.Brent(Counted, p.a, p.b, out returnCode, tol, tol),
                        7 => Solver.OldModAB(Counted, p.a, p.b, out returnCode, tol, tol),
                        8 => Solver.SgModAB(Counted, p.a, p.b, out returnCode, tol, tol),
                        _ => throw new NotImplementedException()
                    };
                    // Raw number of callback evaluations; validation never overwrites it.
                    var solverReturnCode = returnCode;
                    if (!p.HasReference)
                        returnCode = Solver.ReturnCode.Invalid;
                    else if (returnCode == Solver.ReturnCode.Success && !p.IsRoot(root, tol, tol))
                        returnCode = Solver.ReturnCode.FalseConvergence; // validation status, raw count kept
                    results[i, j] = new Result 
                    { 
                        Root = root, 
                        FunctionValue = double.IsNaN(root) ? double.NaN : p.F(root), 
                        EvaluationCount = evals,
                        ReturnCode = returnCode,
                        SolverReturnCode = solverReturnCode
                    };
                    if (returnCode == Solver.ReturnCode.Success) 
                        successCount[j]++;
                    else if (returnCode == Solver.ReturnCode.Invalid) 
                        invalidCount[j]++;
                    else if (returnCode == Solver.ReturnCode.FalseConvergence) 
                        falseConvergenceCount[j]++;
                    else if (returnCode == Solver.ReturnCode.MaxIterationsExceeded) 
                        maxCount[j]++;

                    sumAll[j] += evals;
                    if (returnCode != Solver.ReturnCode.Invalid)
                    {
                        sumValid[j] += evals;
                        logSum[j] += Math.Log(evals);
                        if (evals > max[j])
                            max[j] = evals;
                    }
                }
            }
            // Penalized score: a run that does not succeed, including one rejected by the
            // validator, is charged a common penalty, one larger than the largest raw count
            // among the reference-eligible cases. It is a scoring convention for ranking only,
            // not an evaluation count; all other rows report raw counts.
            var penalty = 1;
            for (int i = 0; i < problemCount; ++i)
                if (problems[i].HasReference)
                    for (int j = 0; j < methodCount; ++j)
                        penalty = Math.Max(penalty, results[i, j].EvaluationCount + 1);

            int Cost(Result r) =>
                r.ReturnCode == Solver.ReturnCode.Success ? r.EvaluationCount : penalty;

            var score = new int[methodCount];
            for (int i = 0; i < problemCount; ++i)
            {
                var resultRow = Enumerable.Range(0, methodCount).Select(j => results[i, j]).ToArray();
                if (!problems[i].HasReference) // not reference-eligible: excluded from the score
                    continue;

                var rowMinCount = resultRow.Min(Cost);
                var rowMaxCount = resultRow.Max(Cost);
                for (int j = 0; j < methodCount; ++j)
                {
                    var count = Cost(resultRow[j]);
                    score[j] += count;
                    if (count == rowMinCount) bestAt[j]++;
                    if (count == rowMaxCount) worstAt[j]++;
                }
            }
            Console.WriteLine(
@"## Benchmark results (C#)

List of algorithms:  
- bs 	 – Bisection  
- fp	 – False-position  
- ill	 – Illinois  	
- AB	 – Anderson-Björck  
- ITP	 – Interpolate. truncate. project  
- Rid	 – Ridders  
- Brе	 – Brent  
- OldModAB – The old Modified Anderson-Björck from April 2026 MDPI Algorithms + fixes
- SGModab  – Safeguarded Modified Anderson-Björck, revised November 2026"
                );
            for (int k = 0; k < 5; ++k)
            {
                switch (k)
                {
                    case 0: Console.WriteLine("\r\n### Results (roots found)\r\n"); break;
                    case 1: Console.WriteLine("\r\n### Function values\r\n"); break;
                    case 2: Console.WriteLine("\r\n### Validation classifications\r\n"); break;
                    case 3: Console.WriteLine("\r\n### Raw solver return codes\r\n"); break;
                    case 4: Console.WriteLine("\r\n### Evaluation count\r\n"); break;
                }
                Console.WriteLine("|   Func |    bs |    fp |   ill |    ab |   ITP |   rid |    br |OldModAB| SGModab");
                Console.WriteLine("| ------ | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ------");
                for (int i = 0; i < problemCount; ++i)
                {
                    Problem p = problems[i];
                    Console.Write($"|{p.Name.PadLeft(7)} | ");
                    for (int j = 0; j < methodCount; ++j)
                    {
                        var result = results[i, j];
                        if (k == 0)
                            Console.Write(double.IsNaN(result.Root) ? "NaN | " : result.Root);
                        else if (k == 1)
                            Console.Write(double.IsNaN(result.FunctionValue) ? "NaN | " : result.FunctionValue);
                        else if(k == 2)
                            Console.Write($"    {result.ReturnCode.ToString()[0]}");
                        else if (k == 3)
                            Console.Write($"    {result.SolverReturnCode.ToString()[0]}");
                        else
                            Console.Write(result.EvaluationCount.ToString().PadLeft(5));

                        if (j < methodCount - 1) Console.Write(" | ");
                    }
                    Console.WriteLine("");
                }
            }
            for (int j = 0; j < methodCount; ++j)
            {
                var count = problemCount - invalidCount[j];
                average[j] = (double)sumValid[j] / count;
                mean[j] = Math.Exp(logSum[j] / count); // geometric mean
                variance[j] = 0.0;
                for (int i = 0; i < problemCount; ++i)
                {
                    var result = results[i, j];
                    if (result.ReturnCode != Solver.ReturnCode.Invalid)
                    {
                        var diff = result.EvaluationCount - average[j];
                        variance[j] += diff * diff;
                    }
                }
                stddev[j] = Math.Sqrt(variance[j] / count);
                // The median must be taken over column j only - results.Cast<Result>() would
                // flatten the whole problemCount x methodCount matrix
                var sorted = Enumerable.Range(0, problemCount)
                    .Select(i => results[i, j])
                    .Where(r => r.ReturnCode != Solver.ReturnCode.Invalid)
                    .Select(r => r.EvaluationCount)
                    .OrderBy(n => n)
                    .ToArray();
                median[j] = sorted.Length == 0 ? double.NaN
                    : sorted.Length % 2 == 1 ? sorted[sorted.Length / 2]
                    : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);
            }
            static void WriteRow<T>(string label, IEnumerable<T> values, string? format = null)
                where T : struct, IFormattable =>
                Console.WriteLine($"|{label,7} | " + string.Join(" | ",
                    values.Select(v => v.ToString(format, System.Globalization.CultureInfo.InvariantCulture).PadLeft(5))));

            WriteRow("SumAll", sumAll);
            WriteRow("Valid", sumValid);
            WriteRow("Ave", average, "f1");
            WriteRow("Mean", mean, "f1");
            WriteRow("StdDev", stddev, "f1");
            WriteRow("Median", median, "f1");
            WriteRow("Max", max);
            WriteRow("BestAt", bestAt);
            WriteRow("WorstAt", worstAt);
            WriteRow("Succ", successCount);
            WriteRow("Invalid", invalidCount);
            WriteRow("False", falseConvergenceCount);
            WriteRow("MaxIter", maxCount);
            WriteRow("Score", score);
            Console.WriteLine($"\r\nRaw counts include the two endpoint evaluations; validation follows eq. (54) " +
                $"without accepting computed zeros. Score is a penalized score (failed runs charged {penalty}), " +
                "a scoring convention, not an evaluation count.");
            Console.WriteLine("\r\nPress any key to continue...");
            if (!Console.IsInputRedirected) Console.ReadKey();
        }
    }
}
