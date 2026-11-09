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

Version 1.2.1 refines the rational step of 1.2.0. The 3-point model `f(x) = (x - r)/(p·x + q)` now goes through the Anderson-Björck-corrected bracket instead of the true residuals, with the true residual of the dropped endpoint as its third point. A step that would move the Anderson-Björck secant point by more than 0.7 of its distance to the nearer endpoint is no longer replaced by the secant point but clamped to that limit, on the side the model points to. The step still takes two divisions. On the benchmark, 11 problems need fewer evaluations (7 one fewer, f48 and f76 two fewer, f41 three fewer, f80 four fewer) and 4 need more (f13 and f15 one more, f85 and f100 three more), so all 101 need 1871 evaluations instead of 1881. On 200,000 random functions from ten families the new step saves 1.1 % of the evaluations, with no failures. In the timing benchmark below, pymodab is 1.34 times as fast as cybrentq and 9.9 times as fast as brentq.

Earlier versions:
- **1.2.0** – The rational modab (RModAB): a rational correction bends the secant step to follow the curvature of the function, and a closing step keeps each new point at least half the tolerance away from the endpoint just evaluated; 1881 evaluations on the benchmark.
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
   SUM|     4883|     2810|     2769|     3425|     2112|     2810|     1993| <mark>**1871**</mark>| 
   AVG|       48|       28|       27|       34|       21|       28|       20| <mark>**19**</mark>| 
MEDIAN|       49|       13|       12|       16|       12|       13|       13| <mark>**11**</mark>| 
   MIN|        3|        4|        4|        3|        3|        4|        3| <mark>**3**</mark>| 
   MAX|       53|      102|      102|      202|       70|      102|       78| <mark>**56**</mark>| 
FACTOR|   2.610x|   1.502x|   1.480x|   1.831x|   1.129x|   1.502x|   1.065x| <mark>**1.000x**</mark>| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**RModAB**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1480.45|   935.32|   887.77|  1023.61| 57058.75|   113.52|   172.83| <mark>**84.39**</mark>| 
   AVG|  14.6579|   9.2606|   8.7898|  10.1348| 564.9381|   1.1240|   1.7112|<mark>**0.8355**</mark>| 
MEDIAN|  14.3018|   5.1766|   4.7206|   5.7367| 334.1344|   0.6116|   1.2809|<mark>**0.6183**</mark>| 
   MIN|   1.8050|   1.5092|   1.7821|   1.5801|  73.4375|   0.2257|   0.6471|<mark>**0.1898**</mark>| 
   MAX|  48.4630|  42.4292|  36.0671|  53.3243|1928.2951|   6.6491|   5.1035|<mark>**2.8232**</mark>| 
FACTOR|  17.543x|  11.083x|  10.520x|  12.130x| 676.144x|   1.345x|   2.048x|<mark>**1.000x**</mark>|
#### Notes:  

Last Run on: 09.11.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 from PyPI - the previous implementation with ctypes  
RModAB: pymodab version 1.2.1 - the latest implementation, rational modab (RModAB), as native C extension  