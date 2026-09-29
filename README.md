
# ModAB Root-Finding Algorithm

This repository provides efficient implementations of the **Modified Anderson-Björck (ModAB)** bracketing root-finding algorithm [1] for solving nonlinear equations of the form `f(x) = 0`.

Available in **C**, **C#**, **C++ (ROOT.CERN framework)**, **Python**, **Julia**, **Zig**, **Rust**, **Java**, **TypeScript**, and **Excel/VBA**.   

[![PyPI Downloads](https://static.pepy.tech/personalized-badge/pymodab?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/pymodab)

## ✨ Features

- **Guaranteed convergence** - Bracketing approach ensures the root is always found.
- **Great performance** - Fewer function evaluations than Brent, Ridders, and other popular methods. Being simple and lightweight, it has very little computational overhead per iteration.
- **Worst-case optimality** - Retains the worst-case optimality of Bisection for the hard cases.
- **Cross-platform** - Python package that works on Windows, Linux, and macOS.
- **Multiple languages** - Use in your preferred environment.
- **Extensive testing and benchmarking** - The algorithm is benchmark against a broad set of test functions of many different types.

## 📕 Algorithm Overview

The ModAB algorithm combines the reliability of bisection with the speed of the secant method:

1. **Bisection phase** - Starts with bisection to ensure stability
2. **Linearity detection** - Monitors when the function behavior is close enough to linear
3. **False-position acceleration** - Switches to false-position method when conditions are favorable
4. **Anderson-Björck correction** - Applies A&B corrections to prevent stalling
5. **Adaptive fallback** - Returns to bisection if progress slows

This hybrid approach achieves superlinear convergence while maintaining the worst-case optimality of bisection.

## ⏱️ Benchmarks

&emsp;&emsp;[For 101 functions test set](Test%20functions.pdf)

<img height="300" alt="image" src="https://github.com/user-attachments/assets/435c4fb9-7f79-4a58-a33b-23a84d82311a" />

### Times
| Method        | Mean      | Error    | StdDev   |
|-------------- |----------:|---------:|---------:|
| Bisection     | 105.77 us | 1.800 us | 3.246 us |
| FalsePosition | 379.02 us | 4.851 us | 4.538 us |
| Illinois      | 157.61 us | 1.998 us | 1.668 us |
| AndersonBjork | 225.40 us | 1.217 us | 1.016 us |
| ITP           | 246.71 us | 3.727 us | 3.304 us |
| Ridders       | 110.09 us | 2.127 us | 2.364 us |
| Brent         | 125.04 us | 2.383 us | 2.112 us |
| ModAB         |  59.34 us | 0.728 us | 0.645 us |
| SGModab       |  61.32 us | 1.221 us | 1.588 us |
| <mark>**SGModab**</mark>|  <mark>61.32 us</mark> | <mark>1.221 us</mark> |  <mark>1.588 us</mark> |

*BenchmarkDotNet v0.15.8, .NET 10.0.12 x64 RyuJIT x86-64-v4  
Windows 11 (10.0.26200.9457/25H2/2025Update/HudsonValley2)  
Intel Core i7-1065G7 CPU 1.30GHz 8 logical and 4 physical cores + 16 GB RAM

### Number of Evaluations
|Func | bs | fp | ill | AB | ITP | Rid | Br | ModAB | <mark>**SGModAB**</mark>
| --: | --: | --: | --: | --: | --: | --: | --: | --: | --: 
| SumAll |  4882 | 11192 |  4087 |  5267 |  3337 |  3354 |  3141 |  1992 |  <mark>2000</mark>
|  Valid |  4475 | 10715 |  3486 |  4666 |  2925 |  2897 |  2699 |  1705 |  <mark>1715</mark>
|    Ave |  48.1 | 115.2 |  37.5 |  50.2 |  31.5 |  31.2 |  29.0 |  18.3 |  <mark>18.4</mark>
|   Mean |  46.5 |  68.4 |  22.1 |  22.6 |  24.8 |  22.2 |  17.2 |  14.7 |  <mark>14.9</mark>
| StdDev |   7.6 |  88.2 |  49.4 |  69.6 |  19.1 |  33.7 |  39.4 |  14.7 |  <mark>15.0</mark>
| Median |  49.0 | 202.0 |  15.0 |  13.0 |  23.0 |  17.0 |  12.0 |  12.0 |  <mark>12.0</mark>
|    Max |    53 |   202 |   202 |   202 |    55 |   202 |   142 |    78 |  <mark>&emsp;87</mark>
| BestAt |    14 |     6 |     2 |    26 |     7 |     5 |    35 |    44 |  <mark>&emsp;42</mark>
|WorstAt |    38 |    50 |     5 |    13 |     2 |     3 |     0 |     1 |  <mark>&emsp;&ensp;1</mark>
|   Succ |    93 |    46 |    88 |    80 |    93 |    91 |    93 |    93 |  <mark>&emsp;93</mark>
|Invalid |     8 |     8 |     8 |     8 |     8 |     8 |     8 |     8 |  <mark>&emsp;&ensp;8</mark>
|  False |     0 |    19 |     0 |     0 |     0 |     2 |     0 |     0 |  <mark>&emsp;&ensp;0</mark>
|MaxIter |     0 |    28 |     5 |    13 |     0 |     0 |     0 |     0 |  <mark>&emsp;&ensp;0</mark>


[Detailed results](Benchmark%20results%20Safeguarded.md)

## 📄 License

MIT License - see the [LICENSE](C/LICENSE) file for details.

## 🧮 Implementations in other software/libraries

Calcpad - 					https://calcpad.eu   							- C#  
Root-Fortran - 				https://github.com/jacobwilliams/roots-fortran	- Fortran  
ROOT.CERN - 				https://github.com/root-project/root			- C++  
SCiML/NonlinearSolve.jl - 	https://github.com/SciML/NonlinearSolve.jl		- Julia  
JuliaMath/Roots.jl -		https://github.com/JuliaMath/Roots.jl			- Julia  
MultiFloats.jl - 			https://github.com/dzhang314/MultiFloats.jl		- Julia  
mpmath - 					https://pypi.org/project/mpmath/				- Python (https://github.com/mpmath/mpmath)  
PyModAB - 					https://pypi.org/project/pymodab/				- Python/C

## 📖 References

1. Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332

2. Ganchovski N. Structural Analysis by Functional Modeling in the Cloud.  PhD Thesis **2025**, UACEG, Sofia

3. Ganchovski, N.; Traykov, A. (2023). "Modified Anderson-Björck's method for solving non-linear equations in structural mechanics." *IOP Conference Series: Materials Science and Engineering*, 1276(1), 012010. [DOI: 10.1088/1757-899X/1276/1/012010](https://doi.org/10.1088/1757-899X/1276/1/012010)

## 📚 Citations

1. Galvão, Henrique & Silva, Valdelírio. (2024). Classes de Métodos Numéricos não Convencionais para Determinação de Raízes de Funções. [DOI: 10.13140/RG.2.2.24088.57606](https://doi.org/10.13140/RG.2.2.24088.57606)

2. Błonka, Adrian. (2025). Optimization of Network Tied-Arch Bridges with Metaheuristic and Gradient-based Algorithms. https://www.researchgate.net/publication/405267708

3. Yarndley, Jack & Evans, Adam & Zhou, Xingyu & Wijayatunga, Minduli & Armellin, Roberto. (2026). Exoplanetary Tour Design with Solar Sails: TheAntipodes Results in the GTOC13 Problem. [DOI: 10.48550/arXiv.2607.10150](https://doi.org/10.48550/arXiv.2607.10150)



## 🐝Contributing

Contributions are welcome. You can open an issue or submit a pull request.
