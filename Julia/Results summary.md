# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `RModAB` — the rational modab: `NonlinearSolve.jl/modAB.jl` from this repository, loaded
  over the local checkout (branch `fix-modab-units-and-stagnation` at `1f0f49fe`). On top of
  the residual test `42|y3| >= min(|yl|, |yr|)`, the halved switching margin after each
  fallback and the plateau rule of the previous version (`modAB_sg`), it changes the A&B step
  in two ways. A rational correction bends the secant step to follow the curvature of the
  function, using a 3-point model through the A&B-corrected bracket and the endpoint it dropped
  last, with the step kept within 0.7 of the distance from the A&B secant point to the nearer
  endpoint; the plain A&B secant is used when the curvature changes sign within the bracket.
  A closing step keeps each A&B point at least `ϵ` away from the endpoint
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
SUM  |   4941 |  3118 |   4094 |   72940 | 2834 |      2146 |  1912
AVE  |   48,9 |  30,9 |   40,5 |   722,2 | 28,1 |      21,2 |  18,9
REL  |   258% |  163% |   214% |   3815% | 148% |      112% |  100%

Failures: `alefeld` returns `NaN` on f100. `RModAB` solves all 101. `modAB_rel` hits
`maxiters = 200` on f101 (201 evaluations); the point it returns is within 1e-41 of the root,
but the solve ends with `ReturnCode.MaxIters`.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | RModAB
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4789 |  1950 |   3272 |    1829 | 2544 | 2072 |      2146 |  1912
AVE  |   47,4 |  19,3 |   32,4 |    18,1 | 25,2 | 20,5 |      21,2 |  18,9
REL  |   250% |  102% |   171% |     96% | 133% | 108% |      112% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
`modAB_rel` hits `maxiters` on f101, as above.

Across the 101 problems `RModAB` and `modAB_rel` differ on 66 evaluation counts. `RModAB` is
cheaper on 61 of them and needs more on 5: one more on f15, f42, f65 and f85, three more on
f100. f101 accounts for most of the difference, 201→56.

Against the previous version `modAB_sg` (1974 evaluations on 03.11.2026), `RModAB` saves 62
evaluations, 3.1 %: 50 problems need fewer (32 one fewer, 14 two fewer, 4 three fewer) and 7
need more (one more on f15, f20, f42, f65 and f81, two more on f85, three more on f100). The
gain is smaller than in C and C# (4.1 %). Those stop on a relative tolerance and this one on
the absolute `abstol = 1e-14`, which is the likely reason. The first `RModAB` of 08.11.2026,
with the rational model through the true residuals and the A&B secant as its fallback, took
1924 evaluations; the A&B-corrected bracket with the 0.7 limit saves 12 more (11 problems need
fewer, 6 need more).

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop. Each entry is the smallest of the medians from two runs on 09.11.2026, on mains
power with `% Processor Performance` between 225 % and 254 %. The solvers that did not change
match the runs of 08.11.2026 to within 0,7 % in total (`modAB_rel` 33,6 μs in both).

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  RModAB
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  |  95,8 μs | 98,7 μs | 52,1 μs | 325,5 μs* | 78,8 μs |   33,6 μs | 28,7 μs
AVE  |   949 ns |  977 ns |  516 ns |  3255 ns* |  780 ns |    332 ns |  285 ns
REL  |     333% |    343% |    181% |    1133%* |    274% |      117% |    100%

\* `alefeld` errors on f100, so its total covers only the 100 problems it solved
and understates the true cost. Every other solver completed all 101.

`RModAB` is 14,4 % faster than the v1.12.8 release over all 101 problems, and 3,7 % faster
without f101. It is faster on 43 problems, within ±3 % on 17 and slower on 41, by 3–56 %; the
largest are f42 (93 ns → 145 ns), f10 (115 ns → 150 ns) and f11 (83 ns → 107 ns). Per
function evaluation the release costs 15,6 ns and `RModAB` 15,0 ns. The previous version
`modAB_sg` took 27,8 μs on 03.11.2026 (14,1 ns per evaluation): the rational correction
saves 3,1 % of the evaluations but makes each step a little dearer, so `RModAB` is 3,4 %
slower than `modAB_sg` on these cheap test functions. The first `RModAB` of 08.11.2026 took
28,5 μs, 0,7 % less, about as much as the unchanged solvers moved between the two days.
