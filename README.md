
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
| Method        | Mean      | Error    | StdDev   | Median    |
|-------------- |----------:|---------:|---------:|----------:|
| Bisection     |  95.50 us | 0.806 us | 0.673 us |  95.38 us |
| FalsePosition | 370.83 us | 7.344 us | 7.542 us | 374.59 us |
| Illinois      | 160.40 us | 3.161 us | 3.640 us | 160.53 us |
| AndersonBjork | 225.74 us | 2.502 us | 2.218 us | 225.50 us |
| ITP           | 244.70 us | 2.759 us | 2.445 us | 245.08 us |
| Ridders       | 108.69 us | 1.210 us | 1.132 us | 108.14 us |
| Brent         | 126.79 us | 2.502 us | 2.457 us | 126.16 us |
| ModAB         |  59.62 us | 0.614 us | 0.544 us |  59.74 us |
| **SGModab**   | **60.73 us** | **1.209 us** | **2.576 us** | **59.64 us** |

*BenchmarkDotNet v0.15.8, .NET 10.0.12 x64 RyuJIT x86-64-v4  
Windows 11 (10.0.26200.9457/25H2/2025Update/HudsonValley2)  
Intel Core i7-1065G7 CPU 1.30GHz 8 logical and 4 physical cores + 16 GB RAM

### Number of Evaluations
|Func | bs | fp | ill | AB | ITP | Rid | Br | ModAB | **SGModAB**
| --: | --: | --: | --: | --: | --: | --: | --: | --: | --: 
| SumAll |  4835 | 11192 |  4087 |  5267 |  3337 |  3354 |  3141 |  1992 |  **1973**
|  Valid |  4428 | 10715 |  3486 |  4666 |  2925 |  2897 |  2699 |  1705 |  **1686**
|    Ave |  47.6 | 115.2 |  37.5 |  50.2 |  31.5 |  31.2 |  29.0 |  18.3 |  **18.1**
|   Mean |  45.2 |  68.4 |  22.1 |  22.6 |  24.8 |  22.2 |  17.2 |  14.7 |  **14.6**
| StdDev |   8.9 |  88.2 |  49.4 |  69.6 |  19.1 |  33.7 |  39.4 |  14.7 |  **14.6**
| Median |  49.0 | 202.0 |  15.0 |  13.0 |  23.0 |  17.0 |  12.0 |  12.0 |  **12.0**
|    Max |    53 |   202 |   202 |   202 |    55 |   202 |   142 |    78 |    **77**
| BestAt |    14 |     6 |     2 |    27 |     7 |     5 |    36 |    46 |    **52**
|WorstAt |    38 |    50 |     5 |    13 |     2 |     3 |     0 |     1 |     **1**
|   Succ |    93 |    46 |    88 |    80 |    93 |    91 |    93 |    93 |    **93**
|Invalid |     8 |     8 |     8 |     8 |     8 |     8 |     8 |     8 |     **8**
|  False |     0 |    19 |     0 |     0 |     0 |     2 |     0 |     0 |     **0**
|MaxIter |     0 |    28 |     5 |    13 |     0 |     0 |     0 |     0 |     **0**


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
