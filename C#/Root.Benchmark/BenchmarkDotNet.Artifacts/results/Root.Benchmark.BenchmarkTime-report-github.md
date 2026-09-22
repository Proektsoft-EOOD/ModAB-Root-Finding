```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9457/25H2/2025Update/HudsonValley2)
Intel Core i7-1065G7 CPU 1.30GHz (Max: 1.50GHz), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v4
  DefaultJob : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v4


```
| Method        | Mean      | Error    | StdDev    | Median    |
|-------------- |----------:|---------:|----------:|----------:|
| Bisection     | 100.85 μs | 1.720 μs |  1.609 μs | 100.56 μs |
| FalsePosition | 374.01 μs | 7.375 μs |  9.590 μs | 370.24 μs |
| Illinois      | 150.24 μs | 1.662 μs |  1.554 μs | 149.68 μs |
| AndersonBjork | 214.85 μs | 2.463 μs |  2.057 μs | 214.60 μs |
| ITP           | 269.24 μs | 6.571 μs | 19.063 μs | 263.90 μs |
| Ridders       | 112.38 μs | 2.224 μs |  5.741 μs | 110.93 μs |
| Brent         | 126.07 μs | 1.682 μs |  3.435 μs | 126.29 μs |
| ModAB         |  58.47 μs | 1.124 μs |  1.338 μs |  58.27 μs |
| ModABCorr     |  61.80 μs | 0.894 μs |  0.836 μs |  61.57 μs |
