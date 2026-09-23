
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

&emsp;&emsp;[For 100 functions test set](Test%20functions.pdf)

<img height="300" alt="image" src="https://github.com/user-attachments/assets/435c4fb9-7f79-4a58-a33b-23a84d82311a" />

### Times
| Method        | Mean      | Error    | StdDev    | Median    |
|-------------- |----------:|---------:|----------:|----------:|
| Bisection     |  94.99 us | 1.845 us |  1.636 us |  95.06 us |
| FalsePosition | 399.32 us | 7.982 us | 18.657 us | 399.70 us |
| Illinois      | 159.51 us | 2.934 us |  7.932 us | 157.14 us |
| AndersonBjork | 219.17 us | 1.645 us |  1.284 us | 219.52 us |
| ITP           | 246.22 us | 2.900 us |  2.712 us | 246.88 us |
| Ridders       | 107.08 us | 2.082 us |  2.314 us | 106.71 us |
| Brent         | 123.07 us | 1.724 us |  1.528 us | 123.31 us |
| ModAB         |  59.39 us | 1.169 us |  1.251 us |  58.96 us |
| <mark>**SGModab**</mark>|  <mark>60.26 us</mark> | <mark>1.193 us</mark> |  <mark>1.116 us</mark> |  <mark>59.99 us</mark> |

*BenchmarkDotNet v0.15.8, .NET 10.0.12 x64 RyuJIT x86-64-v4  
Windows 11 (10.0.26200.9457/25H2/2025Update/HudsonValley2)  
Intel Core i7-1065G7 CPU 1.30GHz 8 logical and 4 physical cores + 16 GB RAM

### Number of Evaluations
|Func | bs | fp | ill | AB | ITP | Rid | Br | ModAB | <mark>**SGModAB**</mark>
| --: | --: | --: | --: | --: | --: | --: | --: | --: | --: 
|    Sum |  4424 | 10513 |  3284 |  4464 |  2872 |  2828 |  2622 |  1627 | <mark>1628</mark>
|    Ave |  48.1 | 114.3 |  35.7 |  48.5 |  31.2 |  30.7 |  28.5 |  17.7 | <mark>17.7</mark>
|   Mean |  46.4 |  67.6 |  21.6 |  22.0 |  24.6 |  22.0 |  16.9 |  14.5 | <mark>14.6</mark>
| StdDev |   7.7 |  88.2 |  46.6 |  68.1 |  19.1 |  33.7 |  39.3 |  13.4 | <mark>13.3</mark>
| Median |  49.0 | 133.0 |  15.0 |  13.0 |  23.0 |  17.0 |  12.0 |  12.0 | <mark>12.0</mark>
|    Max |    53 |   202 |   202 |   202 |    55 |   202 |   142 |    55 | <mark>&emsp;55</mark>
| BestAt |    13 |     6 |     2 |    26 |     7 |     5 |    35 |    44 | <mark>&emsp;42</mark>
|WorstAt |    38 |    49 |     4 |    12 |     2 |     3 |     0 |     1 | <mark>&emsp;&ensp;1</mark>
|   Succ |    92 |    46 |    88 |    80 |    92 |    90 |    92 |    92 | <mark>&emsp;92</mark>
|Invalid |     8 |     8 |     8 |     8 |     8 |     8 |     8 |     8 | <mark>&emsp;&ensp;8</mark>
|  False |     0 |    19 |     0 |     0 |     0 |     2 |     0 |     0 | <mark>&emsp;&ensp;0</mark>
|MaxIter |     0 |    27 |     4 |    12 |     0 |     0 |     0 |     0 | <mark>&emsp;&ensp;0</mark>

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
