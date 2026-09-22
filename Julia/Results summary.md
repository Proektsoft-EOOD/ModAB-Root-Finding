# Julia benchmark summary

All numbers below come from the benchmark scripts in this directory, run against the
**local `NonlinearSolve.jl` checkout** (`../../NonlinearSolve.jl`) rather than the
registered package. See `Project.toml` / `setup_env.jl` for the wiring.

Two modAB columns are reported side by side in the same run:

- `modAB_SG` — the local development version in this machine's `NonlinearSolve.jl` checkout.
- `modAB_rel` — the algorithm exactly as shipped in the registered
  `BracketingNonlinearSolve` **v1.12.7** release (byte-identical to SciML
  `upstream/master` `b15bf2b1`), vendored by `modab_release.jl` under the name
  `ModABRelease` so both can run in one pass.

100 test problems (`f01`–`f100`), `abstol = 1e-14`, `maxiters = 200`.

## Count of function evaluations — NonlinearSolve.jl solvers

`nonlinearsolve_modab_benchmark.jl` → `BenchmarkResults.txt`

Func | bisect | brent | ridder | alefeld |  ITP | modAB_rel | modAB_SG
---  | -----: | ----: | -----: | ------: | ---: | --------: | ----:
SUM  |   4890 |  6194 |   4026 |   98050 | 2773 |      1978 |  1944
AVE  |   48,9 |  61,9 |   40,3 |   980,5 | 27,7 |      19,8 |  19,4
REL  |   252% |  319% |   207% |   5044% | 143% |      102% |  100%

Failures: `alefeld` returns `NaN` on f43, f44, f45, f78, f91 and f100; `brent` on f100.
Both modAB versions solve all 100.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | modAB_SG
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4786 |  1944 |   3238 |    1826 | 2515 | 2069 |      1978 |  1944
AVE  |   47,9 |  19,4 |   32,4 |    18,3 | 25,1 | 20,7 |      19,8 |  19,4
REL  |   246% |  100% |   167% |     94% | 129% | 106% |      102% |  100%

Failures: `ridder` returns `NaN` on f72 and f73; on f100 `alefeld` and `A42` return `NaN`,
while `brent`, `ridder` and `ITP` return the right endpoint `ℯ`, which is not a root.
Both modAB versions solve all 100.

Across the 100 problems the two versions differ on 32 evaluation counts: the local
version is cheaper on 20, dearer on 12. Its largest gains are on the odd-power
problems and f80 — f94 30→21, f80 26→19, f92 26→22.

## Wall-clock time — NonlinearSolve.jl solvers

`NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl` → `NonlinearSolve.jl/Benchmark Results.txt`

Median of `@belapsed` over `solve()` only, with the problem constructed outside the
timing loop.

Func | bisect   | brent   | ridder  | alefeld   |   ITP   | modAB_rel |  modAB_SG
---- | -------: | ------: | ------: | --------: | ------: | --------: | ------:
SUM  | 111,3 μs | 97,3 μs | 60,7 μs | 374,7 μs* | 91,8 μs |   35,5 μs | 36,8 μs
AVE  |  1113 ns |  973 ns |  607 ns |  3747 ns* |  918 ns |    355 ns |  368 ns
REL  |     302% |    264% |    165% |    1018%* |    249% |       96% |    100%

\* `alefeld` errors on f43, f44, f45, f78, f91 and f100, so its total covers only the
94 problems it solved and understates the true cost. Every other solver completed all 100.

The two modAB versions trade off against each other. Per function evaluation the
release costs 17,9 ns and the local version 18,9 ns, because the local version's
safeguard helpers do more work per iteration. The local version needs 1,7% fewer
evaluations but ends up 3,7% slower in wall-clock. Both figures are close to the
run-to-run spread on this machine, so the fair reading is that the two are
equivalent in speed, with the local version ahead on evaluation count — which is
what matters when the objective function is expensive.