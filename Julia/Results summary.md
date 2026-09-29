# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `modAB_sg` — the local development version in this machine's `NonlinearSolve.jl` checkout.
- `modAB_rel` — the algorithm exactly as shipped in the registered
  `BracketingNonlinearSolve` **v1.12.7** release (byte-identical to SciML
  `upstream/master` `b15bf2b1`), vendored by `modab_release.jl` under the name
  `ModABRelease` so both can run in one pass.

101 test problems (`f01`–`f101`), `abstol = 1e-14`, `maxiters = 200`. f101 is L. Tomov's
counterexample that needs the `MaxResidualSteps` cap in `modAB_sg`.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | --------: | ----:
SUM  |   4941 |  3118 |   4094 |   72940 | 2834 |      2056 |  2021
AVE  |   48,9 |  30,9 |   40,5 |   722,2 | 28,1 |      20,4 |  20,0
REL  |   244% |  154% |   203% |   3609% | 140% |      102% |  100%

Failures: `alefeld` returns `NaN` on f100. Both modAB versions solve all 101.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4789 |  1950 |   3272 |    1829 | 2544 | 2072 |      2056 |  2021
AVE  |   47,4 |  19,3 |   32,4 |    18,1 | 25,2 | 20,5 |      20,4 |  20,0
REL  |   237% |   96% |   162% |     90% | 126% | 103% |      102% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
Both modAB versions solve all 101.

Across the 101 problems the two versions differ on 33 evaluation counts: the local
version is cheaper on 21, dearer on 12. Its largest gains are on the odd-power
problems and f80 — f94 30→21, f80 26→19, f92 26→22; on f101 it needs 76 against 78.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop. Each entry is the smallest of the medians from the runs of 29.09.2026:
in every run a few entries, scattered over problems and solvers, came out about 10 times
slower because of background load on the machine, while the rest agreed to within a few
percent.

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  modAB_sg
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  |  95,8 μs | 98,7 μs | 52,2 μs | 325,9 μs* | 78,2 μs |   30,9 μs | 30,0 μs
AVE  |   949 ns |  977 ns |  516 ns |  3259 ns* |  775 ns |    306 ns |  297 ns
REL  |     319% |    329% |    174% |    1096%* |    261% |      103% |    100%

\* `alefeld` errors on f100, so its total covers only the 100 problems it solved
and understates the true cost. Every other solver completed all 101.

The local version no longer trades speed for its safeguards. Per function
evaluation the release costs 15,0 ns and the local version 14,9 ns. The local
version needs 1,7% fewer evaluations and is 2,8% faster in wall-clock. The time
difference is close to the run-to-run spread on this machine, so the fair reading
is that the two are equivalent in speed, with the local version ahead on
evaluation count — which is what matters when the objective function is expensive.