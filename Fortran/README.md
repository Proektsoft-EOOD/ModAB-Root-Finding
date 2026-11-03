### Results summary

Code and benchmarks are based on this [Roots-Fortran](https://github.com/jacobwilliams/roots-fortran) 
by [Jacob Williams](https://github.com/jacobwilliams) fork:   
https://github.com/Proektsoft-EOOD/roots-fortran

Relative number of evaluations to modAB:

<img height="500" alt="image" src="https://github.com/user-attachments/assets/69749408-a623-4b3e-bff2-85c37de666a1" />
 
 Summary of results from the 157 function tests:
Method | bisection | false pos | abkk | barycentr | zhang | ridders | illinois | blendtf | itp | bdqrf | rbp | pegasus | toms748 | brent | brentq | A&B | brenth | muller | modab | chandr
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
Sum | 10130 | 89207 | 5080 | 8605 | 3884 | 3489 | 5319 | 6322 | 3925 | 3503 | 24279 | 5995 | 3482 | 3229 | 3227 | 13968 | 3167 | 4371 | 2651 | 2634
Ave | 45.4 | 400.0 | 22.8 | 38.6 | 17.4 | 15.6 | 23.9 | 28.3 | 17.6 | 15.7 | 108.9 | 26.9 | 15.6 | 14.5 | 14.5 | 62.6 | 14.2 | 19.6 | 11.9 | 11.8
Max | 86 | 1086 | 100 | 334 | 92 | 94 | 175 | 104 | 54 | 92 | 2052 | 204 | 148 | 128 | 128 | 1086 | 128 | 339 | 52 | 54
Mean | 42.9 | 104.9 | 19.9 | 31.2 | 15.8 | 13.8 | 16.4 | 24.3 | 14.8 | 13.7 | 16.0 | 15.8 | 12.2 | 11.7 | 11.7 | 16.5 | 11.4 | 12.5 | 10.9 | 10.8
Best | 10 | 8 | 0 | 3 | 1 | 9 | 11 | 10 | 12 | 9 | 23 | 52 | 58 | 60 | 60 | 66 | 73 | 73 | 42 | 72
Worst | 92 | 108 | 1 | 11 | 0 | 0 | 0 | 0 | 0 | 0 | 10 | 0 | 0 | 0 | 0 | 8 | 0 | 1 | 0 | 0
Median | 47.0 | 57.0 | 20.0 | 30.0 | 16.0 | 14.0 | 12.0 | 26.0 | 13.0 | 14.0 | 12.0 | 11.0 | 11.0 | 11.0 | 11.0 | 10.0 | 10.0 | 10.0 | 11.0 | 10.0
Std Dev | 10.7 | 480.4 | 15.1 | 39.6 | 10.2 | 10.8 | 27.4 | 15.5 | 12.3 | 10.3 | 417.8 | 37.0 | 17.9 | 15.1 | 15.1 | 202.1 | 15.2 | 34.2 | 6.4 | 6.8
Std Error | 0.72 | 32.17 | 1.01 | 2.65 | 0.68 | 0.73 | 1.83 | 1.04 | 0.82 | 0.69 | 27.98 | 2.48 | 1.20 | 1.01 | 1.01 | 13.53 | 1.02 | 2.29 | 0.43 | 0.46
Rel | 3.82 | 33.65 | 1.92 | 3.25 | 1.47 | 1.32 | 2.01 | 2.38 | 1.48 | 1.32 | 9.16 | 2.26 | 1.31 | 1.22 | 1.22 | 5.27 | 1.19 | 1.65 | 1.00 | 0.99

Detailed results are available here:
[Test_results.pdf](/Fortran/results/Test_results.pdf)

### Timing results

Wall-clock time of each solver over the 223 cases, median of 35 passes over the whole
set, with the per-call overhead of the library API removed (see below):

| Method | Evals | Time, µs | Rel |
|---|--:|--:|--:|
| chandr | 2634 | 72.2 | 0.98 |
| **modab** | **2651** | **73.8** | **1.00** |
| brenth | 3167 | 83.7 | 1.13 |
| brentq | 3227 | 84.1 | 1.14 |
| brent | 3229 | 89.0 | 1.20 |
| zhang | 3884 | 91.6 | 1.24 |
| ridders | 3489 | 93.6 | 1.27 |
| bdqrf | 3503 | 106.0 | 1.44 |
| abkk | 5080 | 122.9 | 1.66 |
| toms748 | 3482 | 126.8 | 1.72 |
| blendtf | 6322 | 131.9 | 1.79 |
| illinois | 5319 | 150.6 | 2.04 |
| muller | 4371 | 155.3 | 2.10 |
| itp | 3925 | 181.3 | 2.46 |
| barycentr | 8605 | 184.6 | 2.50 |
| pegasus | 5995 | 205.2 | 2.78 |
| bisection | 10130 | 226.6 | 3.07 |
| A&B | 13968 | 451.1 | 6.11 |
| rbp | 24279 | 830.5 | 11.25 |
| false pos | 89207 | 2799.2 | 37.91 |

Rel is relative to modab. ModAB and Chandrupatla are the fastest, within 2 %, because
they need the fewest evaluations, at a similar cost of about 28 ns per evaluation.

Through the library API every call also pays a fixed 3–7 µs (stack trampolines for
internal procedures; a system call on Windows), which is 10–20 times modab's own run
time and makes the methods hard to tell apart: 1609 µs per pass through `root_scalar`
and 766 µs through `solver%solve` for modab, against 1611 µs and 762 µs for Chandrupatla.

Intel Core i7-1065G7, Windows 11, gfortran 16.1, fpm `--profile release`.
Details and method: [Timing_results.md](/Fortran/results/Timing_results.md)
