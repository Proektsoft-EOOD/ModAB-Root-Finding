# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `RModAB` — the rational modab: `NonlinearSolve.jl/modAB.jl` from this repository, loaded
  over the local checkout (branch `fix-modab-units-and-stagnation` at `78ed8c32`). On top of
  the residual test `42|y3| >= min(|yl|, |yr|)`, the halved switching margin after each
  fallback and the plateau rule of the previous version (`modAB_sg`), it changes the A&B step
  in two ways. A rational correction bends the secant step to follow the curvature of the
  function, using a 3-point model through the bracket and the endpoint it dropped last; the
  plain A&B secant remains the fallback when the curvature changes sign within the bracket or
  the model is unusable. A closing step keeps each A&B point at least `ϵ` away from the endpoint
  just evaluated, so a one-sided convergence still closes the bracket.
  The first line of each results file records the branch, commit and a hash of `modAB.jl`.
- `modAB_rel` — the algorithm exactly as shipped in the latest registered
  `BracketingNonlinearSolve` **v1.12.8** release (SciML tag `BracketingNonlinearSolve-v1.12.8`,
  `c5f0d24f`), vendored by `modab_release.jl` under the name `ModABRelease` so both can run
  in one pass. v1.12.8 already contains the safeguards of SciML/NonlinearSolve.jl#1326, but not
  the fix for f101: its residual test is `|y3| > yMin / 2`, with no cap.

101 test problems (`f01`–`f101`), `abstol = 1e-14`, `maxiters = 200`. f101 is L. Tomov's
counterexample, which the factor of 42 in the residual test of `RModAB` sends back to bisection.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB_rel | RModAB
---  | -----: | ----: | -----: | ------: | ---: | --------: | ----:
SUM  |   4941 |  3118 |   4094 |   72940 | 2834 |      2146 |  1924
AVE  |   48,9 |  30,9 |   40,5 |   722,2 | 28,1 |      21,2 |  19,0
REL  |   257% |  162% |   213% |   3791% | 147% |      112% |  100%

Failures: `alefeld` returns `NaN` on f100. `RModAB` solves all 101. `modAB_rel` hits
`maxiters = 200` on f101 (201 evaluations); the point it returns is within 1e-41 of the root,
but the solve ends with `ReturnCode.MaxIters`.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | RModAB
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4789 |  1950 |   3272 |    1829 | 2544 | 2072 |      2146 |  1924
AVE  |   47,4 |  19,3 |   32,4 |    18,1 | 25,2 | 20,5 |      21,2 |  19,0
REL  |   249% |  101% |   170% |     95% | 132% | 108% |      112% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
`modAB_rel` hits `maxiters` on f101, as above.

Across the 101 problems `RModAB` and `modAB_rel` differ on 65 evaluation counts. `RModAB` is
cheaper on 58 of them and needs one more on 7: f37, f41, f42, f65, f76, f80 and f84. f101
accounts for most of the difference, 201→56.

Against the previous version `modAB_sg` (1974 evaluations on 03.11.2026), `RModAB` saves 50
evaluations, 2.5 %: 44 problems need fewer (32 one fewer, 12 two fewer) and 6 need one more
(f37, f41, f42, f65, f80, f84). The gain is smaller than in the other ports (3.6 % in C, C#
and Python). Those stop on a relative tolerance and this one on the absolute
`abstol = 1e-14`, which is the likely reason.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop. Each entry is the smallest of the medians from two runs on 08.11.2026, on mains
power with `% Processor Performance` mostly between 200 % and 255 %. The solvers that did not
change match the runs of 03.11.2026 to within about 1 % in total (`modAB_rel` 33,6 μs
against 33,5 μs).

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  RModAB
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  |  95,8 μs | 98,7 μs | 52,2 μs | 327,7 μs* | 79,0 μs |   33,6 μs | 28,5 μs
AVE  |   949 ns |  977 ns |  516 ns |  3277 ns* |  782 ns |    333 ns |  282 ns
REL  |     336% |    346% |    183% |    1149%* |    277% |      118% |    100%

\* `alefeld` errors on f100, so its total covers only the 100 problems it solved
and understates the true cost. Every other solver completed all 101.

`RModAB` is 15,1 % faster than the v1.12.8 release over all 101 problems, and 4,8 % faster
without f101. It is faster on 56 problems, within ±3 % on 17 and slower on 28, by 3–42 %; the
largest are f42 (93 ns → 132 ns), f82 (60 ns → 74 ns) and f28 (301 ns → 358 ns). Per
function evaluation the release costs 15,7 ns and `RModAB` 14,8 ns. The previous version
`modAB_sg` took 27,8 μs on 03.11.2026 (14,1 ns per evaluation): the rational correction
saves 2,5 % of the evaluations but makes each step a little dearer, so `RModAB` is 2,5 %
slower than `modAB_sg` on these cheap test functions.
