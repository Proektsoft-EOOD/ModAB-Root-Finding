## PyModAB root-finding library

A fast and robust root-finding library for Python, using the Modified Anderson-Bjork method (Ganchovski & Traykov, 2023), written in C.
It finds the root of a single nonlinear equation `f(x) = 0` within the specified interval `[x1, x2]`.  
Works in Windows, Linux and Mac OS.

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
- `f`: A continuous function of one variable
- `x1`, `x2`: Bracket interval endpoints (must satisfy `f(x1) * f(x2) < 0`)
- `atol`: Absolute tolerance (default: 1e-14)
- `rtol`: Relative tolerance (default: 1e-14)
- `max_iter`: Maximum iterations (default: 200)

**Returns:** The root, or `NaN` if not found.

#### `get_evaluation_count()`

Returns the number of function evaluations from the last root-finding call.

### Algorithm

Modified Anderson-Björck's method is a new robust and efficient bracketing root-finding algorithm. It combines bisection with Anderson-Björk's method to achieve both fast performance and worst-case optimality.

#### References:

Ganchovski N.; Traykov A. Modified Anderson-Björck's method for solving non-linear equations in structural mechanics. IOP Conference Series: Materials Science and Engineering 2023, 1276 (1) 012010, IOP Publishing.  
https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010/pdf

Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements of the Modified Anderson-Björck (modAB) Root-Finding Algorithm. Preprints 2026, 2026032190.  
https://www.preprints.org/manuscript/202603.2190

### License

MIT License

### Benchmark results

The modAB algorithm is benchmarked against the available algorithms in Python/SciPy in respect to number of evaluations and execution times:
* `bisect` - Bisection method
* `brentq` - Brent’s method (van Wijngaarden–Dekker–Brent, 1973)
* `brenth` - Brent–Dekker variant (hyperbolic extrapolation variant, 1975)
* `ridder` - Ridder’s method (1979)
* `toms748` - Alefeld–Potra–Shi method (1995 - TOMS Algorithm 748)
* `chandr` - Chandrupatla's method (1997) - `scipy.optimize.elementwise.find_root`
* `modAB` - Modified Anderson Bjork's method (Ganchovski & Traykov, 2023)

#### Function evaluations

  Func|   bisect|   brentq|   brenth|   ridder|   chandr|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4411|     2548|     2512|     3158|     1873|     1723|
   AVG|       48|       28|       27|       34|       20|       19|
MEDIAN|       49|       12|       12|       16|       12|       12|
   MIN|        3|        4|        4|        4|        3|        3|
   MAX|       53|      102|      102|      202|       58|       55|
FACTOR|   2.560x|   1.479x|   1.458x|   1.833x|   1.087x|   1.000x|

#### Execution times  (ms per problem, 100 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  2842.38|  2077.80|  1464.31|  1854.27| 100429.98|   287.89|
   AVG|  30.8954|  22.5848|  15.9164|  20.1551| 1091.6302|   3.1293|
MEDIAN|  29.4968|   9.9226|   9.2430|  11.4640| 702.8896|   2.7059|
   MIN|   4.4171|   4.0385|   3.7629|   2.3574| 173.5578|   1.1338|
   MAX|  60.1204| 292.9884|  81.5352| 141.5226| 3600.9297|   9.2801|
FACTOR|   9.873x|   7.217x|   5.086x|   6.441x| 348.847x|   1.000x|

#### Notes:

Last Run on: 02.04.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
numpy Version: 2.4.4  
scipy Version: 1.17.1  
pymodab Version: 1.0.3  

The complete source code to reproduce the above benchmarks is available in [RootBenchmark.py](RootBenchmark.py).  
Detailed benchmark results are listed in [BenchmarkResults.md](BenchmarkResults.md).  
