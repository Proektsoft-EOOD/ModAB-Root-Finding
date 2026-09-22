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

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| modAB_SG|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4832|     2733|     2693|     3357|     2042|     2733|     1915|     1895| 
   AVG|       48|       27|       27|       34|       20|       27|       19|       19| 
MEDIAN|       49|       13|       12|       16|       12|       13|       12|       12| 
   MIN|        3|        4|        4|        3|        3|        4|        3|        3| 
   MAX|       53|      102|      102|      202|       58|      102|       55|       55| 
FACTOR|   2.550x|   1.442x|   1.421x|   1.772x|   1.078x|   1.442x|   1.011x|   1.000x| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| modAB_SG|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1461.96|   865.40|   846.89|  1047.16| 60698.42|   122.77|   174.44|    92.97| 
   AVG|  14.6196|   8.6540|   8.4689|  10.4716| 606.9842|   1.2277|   1.7444|   0.9297| 
MEDIAN|  14.4848|   4.6459|   4.5338|   5.9401| 387.4567|   0.6761|   1.3868|   0.7057| 
   MIN|   1.4859|   1.6424|   1.8570|   1.7509|  75.4128|   0.2147|   0.5570|   0.1720| 
   MAX|  26.4949|  31.4656|  32.3145|  67.3989|1722.7349|   6.7877|   4.7697|   3.0769| 
FACTOR|  15.725x|   9.308x|   9.109x|  11.263x| 652.883x|   1.321x|   1.876x|   1.000x| 

#### Notes:  

Last Run on: 22.09.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 from PyPI - the previous implementation with ctypes  
modAB_SG: pymodab version 1.0.10 - the latest implementation of the safeguarded modab as native C extension  