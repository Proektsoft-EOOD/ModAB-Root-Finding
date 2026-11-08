### Timing results

Wall-clock time of the root-fortran solvers over the same problem set as the
[evaluation-count results](../README.md): 157 test functions, 223 cases with all
parameter values and brackets. Each method solves every case with the same settings as
`root_tests`: `atol = 1e-15`, `rtol = 1e-13`, `ftol = 1e-15`, `maxiter = 1000` and
`bisect_on_failure = .true.`. Every timed run reproduced the evaluation counts of the
summary table exactly (modab: 2522).

Times are the sum, over the 223 cases, of the time per call, in µs. For each case all
20 methods are timed back to back (200 calls each, timed with `system_clock`), so slow
drift of the CPU clock affects every method alike. Each build was run 5 times,
interleaved, with 7 passes over the problem set per run: 35 passes per method in total.

* **Median, Mean, StdDev, Min** are taken over the 35 passes.
* **Rel** is the median time relative to modab.
* **Rel range** is the lowest and highest Rel over the 5 runs.

### Solver time

| Method | Evals | Median | Mean | StdDev | Min | Rel | Rel range | ns/eval |
|---|--:|--:|--:|--:|--:|--:|:-:|--:|
| **modab** | **2522** | **69.0** | **69.1** | **2.1** | **65.2** | **1.00** | **1.00–1.00** | **27.4** |
| chandr | 2634 | 70.5 | 70.7 | 2.1 | 67.3 | 1.02 | 1.01–1.03 | 26.8 |
| brenth | 3167 | 80.5 | 81.0 | 2.5 | 77.0 | 1.17 | 1.16–1.19 | 25.4 |
| brentq | 3227 | 81.8 | 81.8 | 2.6 | 76.9 | 1.19 | 1.18–1.21 | 25.3 |
| brent | 3229 | 85.8 | 86.4 | 2.4 | 81.8 | 1.24 | 1.23–1.28 | 26.6 |
| zhang | 3884 | 89.6 | 89.5 | 2.8 | 84.1 | 1.30 | 1.28–1.31 | 23.1 |
| ridders | 3489 | 91.2 | 91.1 | 3.1 | 86.9 | 1.32 | 1.28–1.34 | 26.1 |
| bdqrf | 3503 | 102.6 | 103.0 | 2.8 | 99.5 | 1.49 | 1.47–1.52 | 29.3 |
| abkk | 5080 | 120.1 | 120.5 | 3.9 | 113.4 | 1.74 | 1.73–1.76 | 23.6 |
| toms748 | 3482 | 123.0 | 123.3 | 3.6 | 117.0 | 1.78 | 1.76–1.81 | 35.3 |
| blendtf | 6322 | 128.1 | 129.3 | 4.3 | 123.3 | 1.86 | 1.85–1.90 | 20.3 |
| illinois | 5319 | 146.9 | 147.2 | 4.6 | 139.0 | 2.13 | 2.11–2.16 | 27.6 |
| muller | 4371 | 150.4 | 151.1 | 4.4 | 143.9 | 2.18 | 2.18–2.21 | 34.4 |
| itp | 3925 | 173.7 | 174.7 | 5.4 | 166.1 | 2.52 | 2.53–2.55 | 44.2 |
| barycentr | 8605 | 177.2 | 178.6 | 5.5 | 169.3 | 2.57 | 2.57–2.62 | 20.6 |
| pegasus | 5995 | 199.4 | 199.3 | 5.7 | 191.3 | 2.89 | 2.87–2.91 | 33.3 |
| bisection | 10130 | 220.2 | 220.8 | 5.9 | 211.8 | 3.19 | 3.16–3.26 | 21.7 |
| A&B | 13968 | 434.5 | 435.2 | 11.5 | 414.7 | 6.30 | 6.25–6.52 | 31.1 |
| rbp | 24279 | 801.6 | 806.6 | 19.5 | 773.3 | 11.62 | 11.51–11.91 | 33.0 |
| false pos | 89207 | 2674.6 | 2694.9 | 70.2 | 2589.3 | 38.77 | 38.39–39.77 | 30.0 |

These times are for the solvers themselves, with the call overhead described below
removed. For this measurement, the bisection fallback in `root_solver%solve` was moved
into a separate, non-inlined subroutine; all else is the code in `src/root_module.F90`.
The solvers are allocated and initialised once and called through `solve`.

Most methods take 20–35 ns per evaluation, test function included, so the times
mostly follow the number of evaluations.
Methods with a costlier iteration (`toms748`, `itp`, `muller`) need more time per
evaluation.

### Time through the library API

| Method | Evals | `root_scalar` median | Rel | OO `solve` median | Rel |
|---|--:|--:|--:|--:|--:|
| **modab** | **2522** | **1578** | **1.00** | **765** | **1.00** |
| chandr | 2634 | 1658 | 1.05 | 760 | 0.99 |
| brenth | 3167 | 1579 | 1.00 | 771 | 1.01 |
| brentq | 3227 | 1582 | 1.00 | 768 | 1.00 |
| brent | 3229 | 1578 | 1.00 | 771 | 1.01 |
| zhang | 3884 | 1592 | 1.01 | 781 | 1.02 |
| ridders | 3489 | 1591 | 1.01 | 780 | 1.02 |
| bdqrf | 3503 | 1598 | 1.01 | 790 | 1.03 |
| abkk | 5080 | 1621 | 1.03 | 813 | 1.06 |
| toms748 | 3482 | 1627 | 1.03 | 814 | 1.06 |
| blendtf | 6322 | 1658 | 1.05 | 830 | 1.08 |
| illinois | 5319 | 1657 | 1.05 | 842 | 1.10 |
| muller | 4371 | 1649 | 1.04 | 846 | 1.11 |
| itp | 3925 | 1688 | 1.07 | 874 | 1.14 |
| barycentr | 8605 | 1787 | 1.13 | 884 | 1.16 |
| pegasus | 5995 | 1694 | 1.07 | 894 | 1.17 |
| bisection | 10130 | 1813 | 1.15 | 932 | 1.22 |
| A&B | 13968 | 2043 | 1.29 | 1142 | 1.49 |
| rbp | 24279 | 2328 | 1.48 | 1543 | 2.02 |
| false pos | 89207 | 4340 | 2.75 | 3548 | 4.64 |

These are medians in µs for the code as shipped, called the two ways a user can call it:

* **`root_scalar`**: `call root_scalar(root_method_xxx, f, ax, bx, ...)`.
* **OO `solve`**: a solver object allocated and initialised once, then
  `call solver%solve(ax, bx, ...)` for every case.

Both paths have a fixed cost per call that is far larger than the solvers' own work,
so they hardly tell the methods apart. With gfortran an internal procedure passed as an
argument needs a trampoline on the stack, and on Windows making it executable is a
system call (`__enable_execute_stack`):

* `root_scalar` passes its internal `func_wrapper` to the solver: about 3.8 µs per call.
* `solve` contains a `func_wrapper` that only the bisection fallback uses, but the
  trampoline is set up on every call: about 3.1 µs per call.

Until commit `f754691` modab had a third overhead: `ModAB` and `safe_secant` used
`ieee_arithmetic`, so gfortran saved and restored the floating-point environment around
every call, 1.3 µs per solve. That made modab 5 times slower than chandrupatla in the
solver-time table (361 µs vs 71 µs). With ordered-comparison NaN and infinity tests the
evaluation counts are unchanged and modab is level with chandrupatla.

### Setup

* Intel Core i7-1065G7 @ 1.30 GHz, 16 GB RAM, Windows 11 Home 10.0.26200, on mains power.
  (`% Processor Performance` 193–247 %, mean 226 %, during the runs)
* GNU Fortran (GCC) 16.1.0, fpm 0.13.0, `--profile release`:
  `-O3 -funroll-loops -fPIC -fcoarray=single -fimplicit-none`
* Source: commit `acccfe7` plus the rational modab (RModAB) changes to `ModAB`, 2026-11-08

The timing program is `root_tests` with three changes. `problems` and `func` are moved
into a module, so that passing `func` costs no trampoline. The case list `ns = [...]` is
assigned only when it is requested, instead of allocated on every function evaluation.
The timing loop above replaces the `cpu_time` call around a single solve, which is too
coarse to measure (the `runs/sec` column of `root_tests` shows `Infinity`).
