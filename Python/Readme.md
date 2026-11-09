## PyModAB root-finding library

[![PyPI Downloads](https://static.pepy.tech/personalized-badge/pymodab?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/pymodab)

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
   SUM|  1391.32|   835.48|   858.20|  1032.86| 57706.69|   113.17|   175.57| <mark>**84.23**</mark>| 
   AVG|  13.7755|   8.2721|   8.4970|  10.2264| 571.3533|   1.1204|   1.7383| <mark>**0.8339**</mark>| 
MEDIAN|  13.6447|   4.5126|   4.1708|   5.3115| 342.9902|   0.6397|   1.3843| <mark>**0.6200**</mark>| 
   MIN|   1.4040|   1.6441|   1.8657|   1.3201|  66.3280|   0.1948|   0.5811| <mark>**0.1695**</mark>| 
   MAX|  26.0941|  35.8400|  85.4255|  89.4150| 1952.4837|   5.7567|   4.9504| <mark>**2.8539**</mark>| 
FACTOR|  16.519x|   9.919x|  10.189x|  12.263x| 685.129x|   1.344x|   2.084x| <mark>**1.000x**</mark>| 

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

The complete source code to reproduce the above benchmarks is available in [RootBenchmarkSciPy.py](RootBenchmarkSciPy.py).  
Detailed benchmark results are listed in [BenchmarkResultsSciPy.md](BenchmarkResultsSciPy.md).  
A similar benchmark against the [PyRoot](https://github.com/SimpleArt/pyroot) library is available in [RootBenchmarkPyRoot.py](RootBenchmarkPyRoot.py), with results in [BenchmarkResultsPyRoot.md](BenchmarkResultsPyRoot.md).  
