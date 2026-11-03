A fast and robust root-finding library written in C for Python, using the Modified Anderson-Bjork method:
Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
Improvements to the Modified Anderson–Björck(modAB) Root-Finding Algorithm.
Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
It finds the root of a single nonlinear equation `f(x) = 0` within the specified interval `[x1, x2]`.

### Installation

```bash
pip install pymodab
```

### Usage

```python
import math
from pymodab import find_root, get_evaluation_count

# Find the root of cos(x) - x = 0 in [0, 1]
root = find_root(lambda x: math.cos(x) - x, 0, 1, 1e-3, 1e-3, 10)
print(f"Root: {root}")  # 0.7390851332086904

# Get the number of function evaluations
print(f"Evaluations: {get_evaluation_count()}")
print(f"Error:       {math.cos(root) - root}")

# Using default tolerances
root = find_root(lambda x: x**2 - 2, 1, 2)
print(f"sqrt(2) = {root}")  # 1.414213562373095

# Get the number of function evaluations
print(f"Evaluations: {get_evaluation_count()}")
print(f"Error:       {root**2 - 2}")
```

### API

#### `find_root(f, x1, x2, atol=1e-14, rtol=1e-14, max_iter=200)`

Find the root of `f(x) = 0` within the interval `[x1, x2]`.

**Parameters:**
- `f`: A continuous function of one variable, or the address (`int`) of a compiled C function `double f(double)`
- `x1`, `x2`: Bracket interval endpoints (must satisfy `f(x1) * f(x2) < 0`)
- `atol`: Absolute tolerance (default: 1e-14)
- `rtol`: Relative tolerance (default: 1e-14)
- `max_iter`: Maximum iterations (default: 200)

**Returns:** The root, or `NaN` if not found. Exceptions raised by `f` are propagated to the caller. This holds for the native extension only: with the ctypes fallback (`pymodab.NATIVE` is `False`), an exception in `f` is only printed to stderr, the search goes on with an undefined value of `f(x)`, and the result cannot be trusted; it may be `NaN` or a wrong root.

Since version 1.0.6, `find_root` is a native CPython extension that calls `f` directly, which is about 2x faster than the previous ctypes wrapper. Version 1.0.7 ships it as a wheel for Windows, Linux (x86_64 and aarch64, glibc and musl) and macOS (Intel and Apple Silicon); one wheel serves every CPython from 3.8 up. On any other platform pymodab falls back to ctypes automatically; `pymodab.NATIVE` tells which one is used.

Version 1.0.8 improves the fallback from secant to bisection steps: a secant step that does not shrink the bracket fast enough is still kept if it at least halves the best residual of the bracket. This avoids needless bisection steps when the secant converges from one side, and falls back sooner when it stalls. `get_evaluation_count()` now returns the exact number of calls to `f`. Version 1.0.9 ships that work as wheels for every supported platform.

Version 1.0.10 adds overflow and NaN safeguards to the arithmetic. The midpoint and the secant point are computed without forming `x1 + x2` or `x1*y2 - x2*y1`, so they stay finite and inside the bracket even for endpoints or residuals near the floating-point limits. A `NaN` from `f` stops the search with `NaN` instead of corrupting the bracket. The Anderson-Björck scaling never lets an auxiliary ordinate underflow to zero. The new secant formula rounds differently, so single problems may need a few evaluations more or less; the 100 benchmark problems below need 1895 in total instead of 1897.

Version 1.1.0 simplifies these safeguards without weakening them. The test that switches from bisection to secant steps is now computed inline as `|ym - y3| < k·|ym| + k·|y3|`, which cannot overflow, and it is skipped while `f(x2) - f(x1)` overflows, so the search keeps bisecting until the residuals shrink. The bracket is updated from the sign of the true residuals instead of the Anderson-Björck auxiliary ordinates, so an auxiliary ordinate that underflows to zero can no longer steer the update, and the extra underflow guard of 1.0.10 is no longer needed. The 100 benchmark problems give the same roots and evaluation counts as with 1.0.10.

Version 1.1.1 separates the bisection and Anderson-Björck steps in the solver loop. A bisection step now updates only the bracket and the true residuals. The Anderson-Björck corrections and the fallback test run only on secant steps, and the auxiliary ordinates start from the true residuals whenever the search switches to secant steps. It also caps the residual exception of 1.0.8: a secant step that does not shrink the bracket fast enough is kept on the strength of halving the residual for at most three consecutive steps (`MAX_RESIDUAL_STEPS`), after which the search falls back to bisection anyway. Without the cap, `f(x) = -x*x` for `x < 0` and `7*x*x` for `x >= 0` on `[-π, 1]`, a counterexample by L. Tomov, keeps halving the residual while the bracket never gets narrower than 1, and the search runs out of iterations; with the cap it converges in 76 evaluations. This function is now problem f101 of the benchmark. The other 100 problems give the same roots and evaluation counts as with 1.1.0; all 101 need 1971 evaluations in total.

Version 1.1.2 streamlines the solver loop. The bracket width and the difference of the true residuals are computed once per step and reused. The exits on a zero and on a `NaN` residual share a single test on the hot path. The best true residual of the bracket, which the fallback test needs, is now computed only when a secant step fails to shrink the bracket fast enough, instead of on every secant step, and the test is written as `2·|y3| >= ymin` instead of `|y3| >= 0.5·ymin`. The 101 benchmark problems give the same roots and evaluation counts as with 1.1.1.

Version 1.1.3 replaces the capped residual exception of 1.1.1 with a stricter test and removes the `MAX_RESIDUAL_STEPS` counter. A secant step that does not shrink the bracket fast enough is now kept only if it reduces the residual more than 42 times, that is `42·|y3| < ymin`, where `ymin` is the smaller true residual of the bracket. The factor comes from 2^(1074/200) ≈ 41.4: at this rate a residual that starts near 1 underflows to zero within the default limit of 200 iterations, because the smallest positive double is 2^-1074. A search that stays on secant steps this way therefore ends by itself, without a counter. The factor is a fixed constant and does not follow `max_iter`, so raising the limit does not weaken the test. On the benchmark, problem f101 now needs 65 evaluations instead of 76, and the other 100 problems give the same roots and evaluation counts as with 1.1.2; all 101 need 1960 evaluations in total instead of 1971.

Version 1.1.4 makes the switch from bisection to secant steps more careful in two ways. First, every fallback from secant steps to bisection halves the margin of the switching test, so after n fallbacks it reads `|ym - y3| < (k·|ym| + k·|y3|)·2^-n`. A function on which the secant steps already failed must look more linear before it gets another try, and borderline switches are rejected. Second, if the bisection step that switches to secant steps reduces the residual of the endpoint it replaces by less than 1 %, the function is on a plateau there. That step then does not count as a move of that side, so the next secant step does not apply the Anderson-Björck correction to the opposite ordinate. On the benchmark, f81 needs 13 evaluations instead of 16 and f101 needs 56 instead of 65, while f76, f92 and f95 need one more each; the other 96 problems give the same evaluation counts as with 1.1.3. All 101 need 1951 evaluations in total instead of 1960.

For maximum speed, pass a compiled function, e.g. from [numba](https://numba.pydata.org/). It is then called from C with no Python overhead:

```python
import math
from numba import cfunc
from pymodab import find_root

@cfunc("float64(float64)")
def f(x):
    return math.cos(x) - x

root = find_root(f.address, 0, 1)
```

#### `get_evaluation_count()`

Returns the number of function evaluations from the last root-finding call.

### Algorithm

Modified Anderson-Björck's method is a new robust and efficient bracketing root-finding algorithm. It combines bisection with Anderson-Björk's method to achieve both fast performance and worst-case optimality.

#### References:

Ganchovski N.; Traykov A. Modified Anderson-Björck's method for solving non-linear equations in structural mechanics. IOP Conference Series: Materials Science and Engineering 2023, 1276 (1) 012010, IOP Publishing.  
https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010/pdf

Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. 
https://doi.org/10.3390/a19050332


### License

MIT License

### Benchmark results

The modAB algorithm is benchmarked against the available algorithms in Python/SciPy in respect to number of evaluations and execution times:
* `bisect`- Bisection method
* `brentq` - Brent’s method (van Wijngaarden–Dekker–Brent, 1973)
* `brenth` - Brent–Dekker variant (hyperbolic extrapolation variant, 1975)
* `ridder` - Ridder’s method (1979)
* `toms748` - Alefeld–Potra–Shi method (1995 - TOMS Algorithm 748)
* `chandr` - Chandrupatla's method (1997) - `scipy.optimize.elementwise.find_root`
* `cybrentq` - Cython implementation of brentq by Gledis Caushaj
* `modAB_SG` - Safeguarded Modified Anderson Bjork's method (Ganchovski & Traykov, 2023; improved 2026 by L.Tomov and N.Ganchovski)
* `modAB_ct` - the old version of pymodab 1.0.5 implemented with ctypes

#### Function evaluations

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**modAB_SG**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4883|     2810|     2769|     3425|     2112|     2810|     1993| <mark>**1951**</mark>| 
   AVG|       48|       28|       27|       34|       21|       28|       20| <mark>**19**</mark>| 
MEDIAN|       49|       13|       12|       16|       12|       13|       13| <mark>**13**</mark>| 
   MIN|        3|        4|        4|        3|        3|        4|        3| <mark>**3**</mark>| 
   MAX|       53|      102|      102|      202|       70|      102|       78| <mark>**56**</mark>| 
FACTOR|   2.503x|   1.440x|   1.419x|   1.756x|   1.083x|   1.440x|   1.022x| <mark>**1.000x**</mark>| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**modAB_SG**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1414.31|   841.66|   817.12|   982.57| 58199.46|   118.47|   173.30| <mark>**87.38**</mark>| 
   AVG|  14.0031|   8.3333|   8.0903|   9.7284| 576.2322|   1.1730|   1.7158| <mark>**0.8651**</mark>| 
MEDIAN|  13.8669|   4.4447|   4.1354|   5.4038| 341.9733|   0.6342|   1.3461| <mark>**0.6657**</mark>| 
   MIN|   1.7627|   1.5331|   2.0402|   1.4579|  74.3486|   0.2054|   0.5725| <mark>**0.2388**</mark>| 
   MAX|  39.6453|  30.6501|  33.2314|  54.1355|1943.2977|   6.1242|   5.0479| <mark>**2.7570**</mark>| 
FACTOR|  16.186x|   9.633x|   9.352x|  11.245x| 666.080x|   1.356x|   1.983x| <mark>**1.000x**</mark>| 

#### Notes:  

Last Run on: 03.11.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 from PyPI - the previous implementation with ctypes  
modAB_SG: pymodab version 1.1.4 - the latest implementation of the safeguarded modab as native C extension  