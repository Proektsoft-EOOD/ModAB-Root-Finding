```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9457/25H2/2025Update/HudsonValley2)
Intel Core i7-1065G7 CPU 1.30GHz (Max: 1.50GHz), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host]     : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v4
  DefaultJob : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v4


```
| Method        | Mean      | Error    | StdDev   |
|-------------- |----------:|---------:|---------:|
| Bisection     |  94.36 μs | 1.872 μs | 1.838 μs |
| FalsePosition | 366.45 μs | 5.343 μs | 4.462 μs |
| Illinois      | 151.03 μs | 1.666 μs | 1.477 μs |
| AndersonBjork | 214.44 μs | 2.019 μs | 1.576 μs |
| ITP           | 239.22 μs | 3.720 μs | 3.480 μs |
| Ridders       | 102.03 μs | 0.947 μs | 0.886 μs |
| Brent         | 119.11 μs | 1.626 μs | 1.521 μs |
| ModAB         |  58.22 μs | 0.791 μs | 0.740 μs |
| SGModab       |  60.62 μs | 1.117 μs | 1.329 μs |
