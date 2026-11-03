# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `modAB_sg` — the local development version in this machine's `NonlinearSolve.jl` checkout:
  branch `fix-modab-units-and-stagnation` at `5667c638`, which has the residual test
  `42|y3| >= min(|yl|, |yr|)` in place of the `MaxResidualSteps` cap, the refactored A&B update and
  the unclamped secant with endpoint reuse. On top of that come two uncommitted changes to the switch
  from bisection to A&B: every fallback to bisection halves the margin of the switching test, and a
  switching bisection that reduces the residual of the replaced endpoint by less than 1 % (a plateau)
  does not count as a move of that side, so the first A&B step does not rescale the other ordinate.
  The first line of each results file records the branch, commit and a hash of `modAB.jl`.
- `modAB_rel` — the algorithm exactly as shipped in the latest registered
  `BracketingNonlinearSolve` **v1.12.8** release (SciML tag `BracketingNonlinearSolve-v1.12.8`,
  `c5f0d24f`), vendored by `modab_release.jl` under the name `ModABRelease` so both can run
  in one pass. v1.12.8 already contains the safeguards of SciML/NonlinearSolve.jl#1326, but not
  the fix for f101: its residual test is `|y3| > yMin / 2`, with no cap.

101 test problems (`f01`–`f101`), `abstol = 1e-14`, `maxiters = 200`. f101 is L. Tomov's
counterexample, which the factor of 42 in the residual test of `modAB_sg` sends back to bisection.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | --------: | ----:
SUM  |   4941 |  3118 |   4094 |   72940 | 2834 |      2146 |  1974
AVE  |   48,9 |  30,9 |   40,5 |   722,2 | 28,1 |      21,2 |  19,5
REL  |   250% |  158% |   207% |   3695% | 144% |      109% |  100%

Failures: `alefeld` returns `NaN` on f100. `modAB_sg` solves all 101. `modAB_rel` hits
`maxiters = 200` on f101 (201 evaluations); the point it returns is within 1e-41 of the root,
but the solve ends with `ReturnCode.MaxIters`.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | modAB_sg
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4789 |  1950 |   3272 |    1829 | 2544 | 2072 |      2146 |  1974
AVE  |   47,4 |  19,3 |   32,4 |    18,1 | 25,2 | 20,5 |      21,2 |  19,5
REL  |   243% |   99% |   166% |     93% | 129% | 105% |      109% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
`modAB_rel` hits `maxiters` on f101, as above.

Across the 101 problems the two modAB versions differ on 27 evaluation counts. The local
version is cheaper on 26 of them and needs one more on f76 (40→41). f101 accounts for most of
the difference, 201→56, and f81 saves four, 17→13. The other 24 save one evaluation each,
mostly because the `FloatingPointLimit` exit reuses the stored residual instead of evaluating
`f(x2)` again.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop. Each entry is the smallest of the medians from two runs on 03.11.2026. Each run
had a few problems that came out up to 10 times slower across all solvers because of
background load on the machine: f13, f14, f18, f19 and f26 in the first run and f59, f60,
f67, f68, f74, f81, f90, f91 and f96 in the second. The solvers that did not change match
the runs of 01.11.2026 to within about 1.5%, and `alefeld` to within about 4%.

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  modAB_sg
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  |  96,2 μs | 98,7 μs | 52,2 μs | 329,1 μs* | 78,7 μs |   33,5 μs | 27,8 μs
AVE  |   952 ns |  977 ns |  517 ns |  3291 ns* |  779 ns |    332 ns |  275 ns
REL  |     346% |    355% |    188% |    1184%* |    283% |      121% |    100%

\* `alefeld` errors on f100, so its total covers only the 100 problems it solved
and understates the true cost. Every other solver completed all 101.

The local version is 16,9% faster than the v1.12.8 release over all 101 problems, and
6,6% faster without f101. It is faster on 85 problems, within ±3% on 11, and slower on
f28 (291 ns → 313 ns), f38 (191 ns → 204 ns), f43 (204 ns → 214 ns), f44 (434 ns → 452 ns)
and f72 (222 ns → 232 ns). Per function evaluation the release costs 15,6 ns and the local
version 14,1 ns: v1.12.8 clamps the secant step and then compares the result with
both endpoints, while the local version makes only the two endpoint comparisons.
