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
   SUM|     4832|     2733|     2693|     3357|     2042|     2733|     1915|     1897| 
   AVG|       48|       27|       27|       34|       20|       27|       19|       19| 
MEDIAN|       49|       13|       12|       16|       12|       13|       12|       12| 
   MIN|        3|        4|        4|        3|        3|        4|        3|        3| 
   MAX|       53|      102|      102|      202|       58|      102|       55|       55| 
FACTOR|   2.547x|   1.441x|   1.420x|   1.770x|   1.076x|   1.441x|   1.009x|   1.000x| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1351.21|   807.07|   787.88|   933.13| 53784.09|   113.74|   163.98|    85.48| 
  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct|    modAB|
   AVG|  13.5121|   8.0707|   7.8788|   9.3313| 537.8409|   1.1374|   1.6398|   0.8548| 
MEDIAN|  13.6084|   4.0090|   4.1190|   5.2633| 322.7620|   0.6304|   1.2578|   0.6489| 
   MIN|   1.7064|   1.5499|   1.7578|   1.3562|  66.9817|   0.2363|   0.5673|   0.1710| 
   MAX|  18.8409|  30.7261|  29.6728|  53.3488| 1530.3598|   5.8249|   4.4991|   2.9119| 
FACTOR|  15.807x|   9.442x|   9.217x|  10.916x| 629.199x|   1.331x|   1.918x|   1.000x| 

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

The complete source code to reproduce the above benchmarks is available in [RootBenchmarkSciPy.py](RootBenchmarkSciPy.py).  
Detailed benchmark results are listed in [BenchmarkResultsSciPy.md](BenchmarkResultsSciPy.md).  
A similar benchmark against the [PyRoot](https://github.com/SimpleArt/pyroot) library is available in [RootBenchmarkPyRoot.py](RootBenchmarkPyRoot.py), with results in [BenchmarkResultsPyRoot.md](BenchmarkResultsPyRoot.md).  
