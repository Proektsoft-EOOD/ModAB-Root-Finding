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
* `modAB_SG` - Safeguarded Modified Anderson Bjork's method (Ganchovski & Traykov, 2023; improved 2026 by L.Tomov and N.Ganchovski)
* `modAB_ct` - the old version of pymodab 1.0.5 implemented with ctypes

#### Function evaluations

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**modAB_SG**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|     4832|     2733|     2693|     3357|     2042|     2733|     1915| <mark>**1895**</mark>| 
   AVG|       48|       27|       27|       34|       20|       27|       19| <mark>**19**</mark>| 
MEDIAN|       49|       13|       12|       16|       12|       13|       12| <mark>**12**</mark>| 
   MIN|        3|        4|        4|        3|        3|        4|        3| <mark>**3**</mark>| 
   MAX|       53|      102|      102|      202|       58|      102|       55| <mark>**55**</mark>| 
FACTOR|   2.550x|   1.442x|   1.421x|   1.772x|   1.078x|   1.442x|   1.011x| <mark>**1.000x**</mark>| 

#### Execution times  (ms per problem, 200 iterations)

  Func|   bisect|   brentq|   brenth|   ridder|   chandr| cybrentq| modAB_ct| <mark>**modAB_SG**</mark>|
----- | ------: | ------: | ------: | ------: | ------: | ------: | ------: | ------: |
   SUM|  1461.96|   865.40|   846.89|  1047.16| 60698.42|   122.77|   174.44| <mark>**92.97**</mark>| 
   AVG|  14.6196|   8.6540|   8.4689|  10.4716| 606.9842|   1.2277|   1.7444| <mark>**0.9297**</mark>| 
MEDIAN|  14.4848|   4.6459|   4.5338|   5.9401| 387.4567|   0.6761|   1.3868| <mark>**0.7057**</mark>| 
   MIN|   1.4859|   1.6424|   1.8570|   1.7509|  75.4128|   0.2147|   0.5570| <mark>**0.1720**</mark>| 
   MAX|  26.4949|  31.4656|  32.3145|  67.3989|1722.7349|   6.7877|   4.7697| <mark>**3.0769**</mark>| 
FACTOR|  15.725x|   9.308x|   9.109x|  11.263x| 652.883x|   1.321x|   1.876x| <mark>**1.000x**</mark>| 

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

The complete source code to reproduce the above benchmarks is available in [RootBenchmarkSciPy.py](RootBenchmarkSciPy.py).  
Detailed benchmark results are listed in [BenchmarkResultsSciPy.md](BenchmarkResultsSciPy.md).  
A similar benchmark against the [PyRoot](https://github.com/SimpleArt/pyroot) library is available in [RootBenchmarkPyRoot.py](RootBenchmarkPyRoot.py), with results in [BenchmarkResultsPyRoot.md](BenchmarkResultsPyRoot.md).  
