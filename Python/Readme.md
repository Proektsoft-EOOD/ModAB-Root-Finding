## PyModAB root-finding library

A fast and robust root-finding library for Python, using the Modified Anderson-Bjork method (Ganchovski & Traykov, 2023; improved 2026), written in C.
It finds the root of a single nonlinear equation `f(x) = 0` within the specified interval `[x1, x2]`.  
Works in Windows, Linux and Mac OS.

### 💾 Installation

```bash
pip install pymodab
```
Download stats: https://pepy.tech/projects/pymodab

### 🛠️ Usage

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

### 🔌API

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

### 📕 Algorithm

Modified Anderson-Björck's method is a new robust and efficient bracketing root-finding algorithm. It combines bisection with Anderson-Björk's method to achieve both fast performance and worst-case optimality.

#### References:

Ganchovski N.; Traykov A. Modified Anderson-Björck's method for solving non-linear equations in structural mechanics. IOP Conference Series: Materials Science and Engineering 2023, 1276 (1) 012010, IOP Publishing.  
https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010/pdf

Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. 
https://doi.org/10.3390/a19050332

### 📄 License

MIT License

### ⏱ Benchmark results

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
   SUM|     4464|     2571|     2534|     3177|     1891|     2571|     1746|     1746| 
   AVG|       48|       28|       27|       34|       20|       28|       19|       19| 
MEDIAN|       49|       12|       12|       16|       12|       12|       12|       12| 
   MIN|        3|        4|        4|        3|        3|        4|        3|        3| 
   MAX|       53|      102|      102|      202|       58|      102|       55|       55| 
FACTOR|   2.557x|   1.473x|   1.451x|   1.820x|   1.083x|   1.473x|   1.000x|   1.000x| 

#### Execution times  (ms per problem, 100 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1359.38|   946.80|   795.71|   928.54| 50164.55|   103.09|   145.35|    75.77| 
  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
   AVG|  14.6170|  10.1806|   8.5561|   9.9843| 539.4037|   1.1085|   1.5629|   0.8147| 
MEDIAN|  14.1612|   4.7324|   4.5651|   5.9825| 307.2908|   0.5875|   1.2580|   0.6201| 
   MIN|   2.1026|   2.5171|   1.7910|   1.4065|  84.3124|   0.2047|   0.5414|   0.1632| 
   MAX|  55.4387|  78.7136|  32.4967|  52.4362|1574.8965|   5.4855|   4.2678|   2.8956| 
FACTOR|  17.941x|  12.496x|  10.502x|  12.255x| 662.085x|   1.361x|   1.918x|   1.000x| 

#### Notes:

Last Run on: 18.09.2026  
Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz (1.50 GHz) with 16.0 GB RAM  
Windows 11 Home  
Python Version: 3.14.7  
numpy Version: 2.4.6  
scipy Version: 1.18.0  
cybrentq Version: 0.1.5 - by Gledis Caushaj (https://github.com/gledi-ai/cybrentq)  
modAB_ct: pymodab version 1.0.5 - the previous implementation with ctypes  
modAB: pymodab version 1.0.7 - the latest implementation as native C extension  

The complete source code to reproduce the above benchmarks is available in [RootBenchmark.py](RootBenchmark.py).  
Detailed benchmark results are listed in [BenchmarkResults.md](BenchmarkResults.md).  
