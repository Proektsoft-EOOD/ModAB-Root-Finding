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
   SUM|     4883|     2810|     2769|     3425|     2112|     2810|     1993| <mark>**1881**</mark>| 
   AVG|       48|       28|       27|       34|       21|       28|       20| <mark>**19**</mark>| 
MEDIAN|       49|       13|       12|       16|       12|       13|       13| <mark>**11**</mark>| 
   MIN|        3|        4|        4|        3|        3|        4|        3| <mark>**3**</mark>| 
   MAX|       53|      102|      102|      202|       70|      102|       78| <mark>**56**</mark>| 
FACTOR|   2.596x|   1.494x|   1.472x|   1.821x|   1.123x|   1.494x|   1.060x| <mark>**1.000x**</mark>| 


#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**RModAB**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1458.11|   901.30|   891.77|  1006.99| 57970.14|   115.41|   174.74| <mark>**85.14**</mark>| 
   AVG|  14.4367|   8.9238|   8.8294|   9.9702| 573.9617|   1.1426|   1.7301| <mark>**0.8430**</mark>| 
MEDIAN|  14.5688|   4.8678|   4.6496|   5.7422| 344.0286|   0.6270|   1.3687| <mark>**0.6272**</mark>| 
   MIN|   1.6549|   1.8202|   2.0263|   1.7382|  76.1830|   0.2046|   0.6398| <mark>**0.1824**</mark>| 
   MAX|  19.6980|  33.9563|  60.6919|  60.1620| 1910.0981|   5.6508|   5.5485| <mark>**2.8541**</mark>| 
FACTOR|  17.126x|  10.586x|  10.474x|  11.828x| 680.890x|   1.356x|   2.052x| <mark>**1.000x**</mark>| 

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

The complete source code to reproduce the above benchmarks is available in [RootBenchmarkSciPy.py](RootBenchmarkSciPy.py).  
Detailed benchmark results are listed in [BenchmarkResultsSciPy.md](BenchmarkResultsSciPy.md).  
A similar benchmark against the [PyRoot](https://github.com/SimpleArt/pyroot) library is available in [RootBenchmarkPyRoot.py](RootBenchmarkPyRoot.py), with results in [BenchmarkResultsPyRoot.md](BenchmarkResultsPyRoot.md).  
