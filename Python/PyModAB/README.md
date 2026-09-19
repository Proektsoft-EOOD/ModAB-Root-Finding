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

**Returns:** The root, or `NaN` if not found. Exceptions raised by `f` are propagated to the caller.

Since version 1.0.6, `find_root` is a native CPython extension that calls `f` directly, which is about 2x faster than the previous ctypes wrapper. Version 1.0.7 ships it as a wheel for Windows, Linux (x86_64 and aarch64, glibc and musl) and macOS (Intel and Apple Silicon); one wheel serves every CPython from 3.8 up. On any other platform pymodab falls back to ctypes automatically; `pymodab.NATIVE` tells which one is used.

Version 1.0.8 improves the fallback from secant to bisection steps: a secant step that does not shrink the bracket fast enough is still kept if it at least halves the best residual, up to three times in a row. This avoids needless bisection steps when the secant converges from one side, and falls back sooner when it stalls. `get_evaluation_count()` now returns the exact number of calls to `f`.

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
* `modAB` - Modified Anderson Bjork's method (Ganchovski & Traykov, 2023; improved 2026)
* `modAB_ct` - the old version of pymodab 1.0.5 implemented with ctypes

#### Function evaluations

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4832|     2733|     2693|     3357|     2042|     2733|     1915|     1897| 
   AVG|       48|       27|       27|       34|       20|       27|       19|       19| 
MEDIAN|       49|       13|       12|       16|       12|       13|       12|       12| 
   MIN|        3|        4|        4|        3|        3|        4|        3|        3| 
   MAX|       53|      102|      102|      202|       58|      102|       55|       55| 
FACTOR|   2.547x|   1.441x|   1.420x|   1.770x|   1.076x|   1.441x|   1.009x|   1.000x| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1313.69|   780.79|   758.26|   908.26| 52391.22|   107.38|   156.39|    83.91| 
  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
   AVG|  13.1369|   7.8079|   7.5826|   9.0826| 523.9122|   1.0738|   1.5639|   0.8391| 
MEDIAN|  13.2611|   3.9732|   3.8331|   5.0113| 320.0253|   0.5897|   1.2359|   0.6279| 
   MIN|   1.5521|   1.5666|   1.8513|   1.3097|  92.8237|   0.2102|   0.7051|   0.2053| 
   MAX|  28.8559|  30.3630|  29.9153|  52.3420| 1524.1777|   5.5486|   4.4191|   2.9652| 
FACTOR|  15.656x|   9.305x|   9.037x|  10.825x| 624.393x|   1.280x|   1.864x|   1.000x| 

#### Notes:  

Last Run on: 20.09.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 from PyPI - the previous implementation with ctypes  
modAB: pymodab version 1.0.8 - the latest implementation as native C extension  