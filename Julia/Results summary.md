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
SUM  |   4890 |  6194 |   4026 |   98050 | 2773 |      1978 |  1945
AVE  |   48,9 |  61,9 |   40,3 |   980,5 | 27,7 |      19,8 |  19,5
REL  |   251% |  318% |   207% |   5041% | 143% |      102% |  100%

Failures: `alefeld` returns `NaN` on f43, f44, f45, f78, f91 and f100; `brent` on f100.
Both modAB versions solve all 100.

## Count of function evaluations — Roots.jl solvers vs. modAB

`roots_modab_benchmark.jl` → `BenchmarkResultsRoots.txt`

Func | bisect | brent | ridder | alefeld |  ITP |  A42 | modAB_rel | modAB_SG
---  | -----: | ----: | -----: | ------: | ---: | ---: | --------: | ----:
SUM  |   4786 |  1944 |   3238 |    1826 | 2515 | 2069 |      1978 |  1945
AVE  |   47,9 |  19,4 |   32,4 |    18,3 | 25,1 | 20,7 |      19,8 |  19,5
REL  |   246% |  100% |   166% |     94% | 129% | 106% |      102% |  100%

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
SUM  | 117,3 μs | 95,9 μs | 58,1 μs | 382,1 μs* | 95,5 μs |   34,2 μs | 32,8 μs
AVE  |  1173 ns |  959 ns |  581 ns |  3821 ns* |  955 ns |    342 ns |  328 ns
REL  |     358% |    292% |    177% |    1165%* |    291% |      104% |    100%

\* `alefeld` errors on f43, f44, f45, f78, f91 and f100, so its total covers only the
94 problems it solved and understates the true cost. Every other solver completed all 100.

The local version no longer trades speed for its safeguards. Per function
evaluation the release costs 17,3 ns and the local version 16,9 ns. The local
version needs 1,7% fewer evaluations and is 4,1% faster in wall-clock. The time
difference is close to the run-to-run spread on this machine, so the fair reading
is that the two are equivalent in speed, with the local version ahead on
evaluation count — which is what matters when the objective function is expensive.