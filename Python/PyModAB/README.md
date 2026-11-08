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

`find_root` is a native CPython extension that calls `f` directly. Wheels are provided for Windows, Linux (x86_64 and aarch64, glibc and musl) and macOS (Intel and Apple Silicon); one wheel serves every CPython from 3.8 up. On any other platform pymodab falls back to a ctypes wrapper automatically; `pymodab.NATIVE` tells which one is used.

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

### Version history

Version 1.2.0 is the rational modab (RModAB). It changes the secant steps in two ways. First, a rational correction: the secant step through the true residuals is corrected for the curvature of the function by the root of the 3-point model `f(x) = (x - r)/(p·x + q)` through the bracket and the endpoint it dropped last. The correction is capped at half the distance from the secant point to the nearer endpoint, and it is skipped, in favour of the Anderson-Björck secant, when the curvature changes sign within the bracket or the model is degenerate. The step takes two divisions, the same as the secant step, and the curvature test none, so a step costs little more than before. Second, a closing step: when the secant steps converge from one side, the next point is kept at least half the tolerance away from the endpoint just evaluated, so it crosses the root and closes the bracket. On the benchmark, 55 problems need fewer evaluations (39 one fewer, 15 two fewer, f35 three fewer), f37 and f80 need one more each, and the median drops from 13 to 11. All 101 need 1881 evaluations in total instead of 1951, and the benchmark runs about 6 % faster relative to cybrentq.

Earlier versions:
- **1.1.4** – Each fallback to bisection halves the margin of the switching test, and a switching bisection on a plateau does not count as a side move; 1951 evaluations on the benchmark.
- **1.1.3** – A slow secant step is kept only if it cuts the residual more than 42 times (2^(1074/200) ≈ 41.4), replacing the step cap of 1.1.1; f101 needs 65 evaluations instead of 76.
- **1.1.2** – Streamlined solver loop with the same roots and evaluation counts.
- **1.1.1** – Separate bisection and Anderson-Björck steps, and a cap on slow secant steps, which fixes L. Tomov's counterexample, now problem f101.
- **1.1.0** – Simpler overflow-safe switching test; the bracket is updated from the true residuals.
- **1.0.10** – Overflow and NaN safeguards: the midpoint and secant point stay finite and inside the bracket, and a `NaN` from `f` stops the search.
- **1.0.8 / 1.0.9** – A slow secant step is kept if it halves the residual; `get_evaluation_count()` counts every call to `f`.
- **1.0.6 / 1.0.7** – Native CPython extension, about 2x faster than the ctypes wrapper, with wheels for all major platforms.

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
* `RModAB` - Rational Modified Anderson Bjork's method (Ganchovski & Traykov, 2023; improved 2026 by L.Tomov and N.Ganchovski; rational correction and closing step added in November 2026)
* `modAB_ct` - the old version of pymodab 1.0.5 implemented with ctypes

#### Function evaluations

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**RModAB**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4883|     2810|     2769|     3425|     2112|     2810|     1993| <mark>**1881**</mark>| 
   AVG|       48|       28|       27|       34|       21|       28|       20| <mark>**19**</mark>| 
MEDIAN|       49|       13|       12|       16|       12|       13|       13| <mark>**11**</mark>| 
   MIN|        3|        4|        4|        3|        3|        4|        3| <mark>**3**</mark>| 
   MAX|       53|      102|      102|      202|       70|      102|       78| <mark>**56**</mark>| 
FACTOR|   2.596x|   1.494x|   1.472x|   1.821x|   1.123x|   1.494x|   1.060x| <mark>**1.000x**</mark>| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**RModAB**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1377.01|   835.83|   819.48|   981.71| 58154.94|   116.80|   175.41| <mark>**85.42**</mark>| 
   AVG|  13.6338|   8.2755|   8.1136|   9.7199| 575.7915|   1.1564|   1.7368| <mark>**0.8457**</mark>| 
MEDIAN|  13.7756|   4.2841|   4.2259|   5.4295| 331.5251|   0.6410|   1.3860| <mark>**0.6582**</mark>| 
   MIN|   1.4452|   1.7892|   2.0950|   1.3739|  66.0735|   0.2271|   0.6304| <mark>**0.1844**</mark>| 
   MAX|  18.0682|  32.1124|  32.3811|  55.2827| 2097.4978|   6.1896|   5.5730| <mark>**3.0084**</mark>| 
FACTOR|  16.121x|   9.785x|   9.594x|  11.493x| 680.824x|   1.367x|   2.054x| <mark>**1.000x**</mark>| 

#### Notes:  

Last Run on: 08.11.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 from PyPI - the previous implementation with ctypes  
RModAB: pymodab version 1.2.0 - the latest implementation, rational modab (RModAB), as native C extension  