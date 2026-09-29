using Root.Benchmark;
using BenchmarkDotNet.Running;
//Proektsoft.Root.Solver.ModAB3(x => 1.0 / x - 10.0, 0.0, 1.0);
BenchmarkCount.Run();
BenchmarkRunner.Run<BenchmarkTime>();