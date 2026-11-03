### Timing results

Wall-clock time of the root-fortran solvers over the same problem set as the
[evaluation-count results](../README.md): 157 test functions, 223 cases with all
parameter values and brackets. Each method solves every case with the same settings as
`root_tests`: `atol = 1e-15`, `rtol = 1e-13`, `ftol = 1e-15`, `maxiter = 1000` and
`bisect_on_failure = .true.`. Every timed run reproduced the evaluation counts of the
summary table exactly (modab: 2651).

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
| chandr | 2634 | 72.2 | 73.8 | 6.4 | 68.9 | 0.98 | 0.97–0.99 | 27.4 |
| **modab** | **2651** | **73.8** | **75.2** | **6.0** | **70.2** | **1.00** | **1.00–1.00** | **27.9** |
| brenth | 3167 | 83.7 | 84.5 | 6.2 | 80.2 | 1.13 | 1.11–1.14 | 26.4 |
| brentq | 3227 | 84.1 | 85.6 | 7.3 | 80.5 | 1.14 | 1.12–1.16 | 26.1 |
| brent | 3229 | 89.0 | 90.3 | 6.8 | 85.0 | 1.20 | 1.20–1.22 | 27.6 |
| zhang | 3884 | 91.6 | 93.3 | 7.6 | 86.9 | 1.24 | 1.22–1.26 | 23.6 |
| ridders | 3489 | 93.6 | 94.6 | 7.5 | 89.0 | 1.27 | 1.23–1.29 | 26.8 |
| bdqrf | 3503 | 106.0 | 107.8 | 9.1 | 102.2 | 1.44 | 1.42–1.46 | 30.3 |
| abkk | 5080 | 122.9 | 125.2 | 10.6 | 117.9 | 1.66 | 1.64–1.69 | 24.2 |
| toms748 | 3482 | 126.8 | 129.1 | 9.8 | 121.0 | 1.72 | 1.71–1.75 | 36.4 |
| blendtf | 6322 | 131.9 | 134.7 | 11.1 | 126.5 | 1.79 | 1.76–1.81 | 20.9 |
| illinois | 5319 | 150.6 | 153.5 | 10.4 | 145.0 | 2.04 | 2.03–2.07 | 28.3 |
| muller | 4371 | 155.3 | 158.1 | 11.4 | 148.0 | 2.10 | 2.08–2.13 | 35.5 |
| itp | 3925 | 181.3 | 184.4 | 17.3 | 174.9 | 2.46 | 2.43–2.50 | 46.2 |
| barycentr | 8605 | 184.6 | 187.5 | 14.9 | 176.6 | 2.50 | 2.45–2.53 | 21.4 |
| pegasus | 5995 | 205.2 | 208.6 | 12.4 | 197.6 | 2.78 | 2.75–2.82 | 34.2 |
| bisection | 10130 | 226.6 | 231.3 | 20.2 | 219.8 | 3.07 | 3.06–3.09 | 22.4 |
| A&B | 13968 | 451.1 | 462.9 | 63.0 | 434.5 | 6.11 | 5.99–6.21 | 32.3 |
| rbp | 24279 | 830.5 | 849.2 | 65.1 | 814.3 | 11.25 | 11.14–11.44 | 34.2 |
| false pos | 89207 | 2799.2 | 2836.1 | 199.8 | 2708.6 | 37.91 | 37.01–38.71 | 31.4 |

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
| chandr | 2634 | 1611 | 1.00 | 762 | 0.99 |
| **modab** | **2651** | **1609** | **1.00** | **766** | **1.00** |
| brenth | 3167 | 1618 | 1.01 | 774 | 1.01 |
| brentq | 3227 | 1609 | 1.00 | 776 | 1.01 |
| brent | 3229 | 1616 | 1.00 | 776 | 1.01 |
| zhang | 3884 | 1637 | 1.02 | 784 | 1.02 |
| ridders | 3489 | 1631 | 1.01 | 782 | 1.02 |
| bdqrf | 3503 | 1630 | 1.01 | 791 | 1.03 |
| abkk | 5080 | 1653 | 1.03 | 814 | 1.06 |
| toms748 | 3482 | 1664 | 1.03 | 817 | 1.07 |
| blendtf | 6322 | 1663 | 1.03 | 829 | 1.08 |
| illinois | 5319 | 1694 | 1.05 | 838 | 1.09 |
| muller | 4371 | 1692 | 1.05 | 842 | 1.10 |
| itp | 3925 | 1772 | 1.10 | 878 | 1.15 |
| barycentr | 8605 | 1724 | 1.07 | 882 | 1.15 |
| pegasus | 5995 | 1736 | 1.08 | 894 | 1.17 |
| bisection | 10130 | 1765 | 1.10 | 926 | 1.21 |
| A&B | 13968 | 1987 | 1.23 | 1140 | 1.49 |
| rbp | 24279 | 2365 | 1.47 | 1534 | 2.00 |
| false pos | 89207 | 4368 | 2.71 | 3500 | 4.57 |

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

* Intel Core i7-1065G7 @ 1.30 GHz, 16 GB RAM, Windows 11 Home 10.0.26200, on mains power
  (`% Processor Performance` 152–255 %, mean 233 %, during the runs)
* GNU Fortran (GCC) 16.1.0, fpm 0.13.0, `--profile release`:
  `-O3 -funroll-loops -fPIC -fcoarray=single -fimplicit-none`
* Source: commit `f754691`, 2026-11-03

The timing program is `root_tests` with three changes. `problems` and `func` are moved
into a module, so that passing `func` costs no trampoline. The case list `ns = [...]` is
assigned only when it is requested, instead of allocated on every function evaluation.
The timing loop above replaces the `cpu_time` call around a single solve, which is too
coarse to measure (the `runs/sec` column of `root_tests` shows `Infinity`).
