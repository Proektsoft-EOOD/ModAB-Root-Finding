# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `modAB_sg` — the local development version in this machine's `NonlinearSolve.jl` checkout:
  branch `modab-residual-cap-f101` on top of SciML `upstream/master` `0055e314`, with the
  `MaxResidualSteps` cap, the refactored A&B update and the unclamped secant with endpoint reuse.
  The first line of each results file records the branch, commit and a hash of `modAB.jl`.
- `modAB_rel` — the algorithm exactly as shipped in the latest registered
  `BracketingNonlinearSolve` **v1.12.8** release (SciML tag `BracketingNonlinearSolve-v1.12.8`,
  `c5f0d24f`), vendored by `modab_release.jl` under the name `ModABRelease` so both can run
  in one pass. v1.12.8 already contains the safeguards of SciML/NonlinearSolve.jl#1326, but not
  the `MaxResidualSteps` cap.

101 test problems (`f01`–`f101`), `abstol = 1e-14`, `maxiters = 200`. f101 is L. Tomov's
counterexample that needs the `MaxResidualSteps` cap in `modAB_sg`.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | --------: | ----:
SUM  |   4941 |  3118 |   4094 |   72940 | 2834 |      2146 |  1995
AVE  |   48,9 |  30,9 |   40,5 |   722,2 | 28,1 |      21,2 |  19,8
REL  |   248% |  156% |   205% |   3656% | 142% |      108% |  100%

Failures: `alefeld` returns `NaN` on f100. `modAB_sg` solves all 101. `modAB_rel` hits
`maxiters = 200` on f101 (201 evaluations); the point it returns is within 1e-41 of the root,
but the solve ends with `ReturnCode.MaxIters`.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4789 |  1950 |   3272 |    1829 | 2544 | 2072 |      2146 |  1995
AVE  |   47,4 |  19,3 |   32,4 |    18,1 | 25,2 | 20,5 |      21,2 |  19,8
REL  |   240% |   98% |   164% |     92% | 128% | 104% |      108% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
`modAB_rel` hits `maxiters` on f101, as above.

Across the 101 problems the two modAB versions differ on 27 evaluation counts, and the
local version is cheaper on all of them. f101 accounts for most of the difference, 201→76.
The other 26 save one evaluation each, mostly because the `FloatingPointLimit` exit reuses
the stored residual instead of evaluating `f(x2)` again.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop. Each entry is the smallest of the medians from two runs on 07.10.2026: in
each run a few entries, scattered over problems and solvers, came out about 10 times
slower because of background load on the machine, while the rest agreed to within a few
percent. The solvers that did not change match the runs of 29.09.2026 to within 1%.

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  modAB_sg
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  |  95,7 μs | 98,5 μs | 52,1 μs | 326,0 μs* | 78,2 μs |   33,5 μs | 28,4 μs
AVE  |   947 ns |  975 ns |  516 ns |  3260 ns* |  774 ns |    332 ns |  281 ns
REL  |     337% |    347% |    184% |    1149%* |    275% |      118% |    100%

\* `alefeld` errors on f100, so its total covers only the 100 problems it solved
and understates the true cost. Every other solver completed all 101.

The local version is 15,4% faster than the v1.12.8 release over all 101 problems, and
6,0% faster without f101. It is faster on 87 problems, within ±3% on 13, and slower only
on f28 (286 ns → 325 ns). Per function evaluation the release costs 15,6 ns and the
local version 14,2 ns: v1.12.8 clamps the secant step and then compares the result with
both endpoints, while the local version makes only the two endpoint comparisons.
