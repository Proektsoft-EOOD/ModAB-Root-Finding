# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

100 test problems (`f01`–`f100`), `abstol = 1e-14`, `maxiters = 200`.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB
---  | -----: | ----: | -----: | ------: | ---: | ----:
SUM  |   4890 |  6194 |   4026 |   98050 | 2773 |  1944
AVE  |   48,9 |  61,9 |   40,3 |   980,5 | 27,7 |  19,4
REL  |   252% |  319% |   207% |   5044% | 143% |  100%

Failures: `alefeld` returns `NaN` on f43, f44, f45, f78, f91 and f100; `brent` on f100.
`modAB` solves all 100.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB
---  | -----: | ----: | -----: | ------: | ---: | ---: | ----:
SUM  |   4786 |  1944 |   3238 |    1826 | 2515 | 2069 |  1944
AVE  |   47,9 |  19,4 |   32,4 |    18,3 | 25,1 | 20,7 |  19,4
REL  |   246% |  100% |   167% |     94% | 129% | 106% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
`modAB` solves all 100.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop.

Func | bisect   | brent   | ridder  | alefeld   |   ITP   |  modAB
---- | -------: | ------: | ------: | --------: | ------: | ------:
SUM  | 107,5 μs | 97,7 μs | 61,5 μs | 364,2 μs* | 88,2 μs | 34,6 μs
AVE  |  1075 ns |  977 ns |  615 ns |  3642 ns* |  882 ns |  346 ns
REL  |     311% |    282% |    178% |    1053%* |    255% |    100%

\* `alefeld` errors on f43, f44, f45, f78, f91 and f100, so its total covers only the
94 problems it solved and understates the true cost. Every other solver completed all 100.