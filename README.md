
# ModAB Root-Finding Algorithm

This repository provides efficient implementations of the **Modified Anderson-Björck (ModAB)** bracketing root-finding algorithm [1] for solving nonlinear equations of the form `f(x) = 0`.

Available in **C**, **C#**, **C++ (ROOT.CERN framework)**, **Python**, **Julia**, **Zig**, **Rust**, **Java**, **TypeScript**, and **Excel/VBA**.

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

&emsp;&emsp;[For 92 functions test set](/MDPI-Algorithms-Data/Test%20Functions.pdf) + 1: f93(x) = x^7 - 0.001

<img height="300" alt="image" src="https://github.com/user-attachments/assets/cd8b63e9-15d6-4cb0-863a-ccdd445de80e" />

### Times
| Method        | Mean      | Error    | StdDev   |
|-------------- |----------:|---------:|---------:|
| Bisection     |  87.17 us | 1.229 us | 1.150 us |
| FalsePosition | 328.11 us | 4.240 us | 3.966 us |
| Illinois      | 104.43 us | 1.507 us | 1.410 us |
| AndersonBjork | 128.49 us | 1.380 us | 1.223 us |
| ITP           | 196.66 us | 3.165 us | 2.961 us |
| Ridders       |  65.79 us | 1.077 us | 0.899 us |
| Brent         | 118.80 us | 2.198 us | 2.056 us |
| ModAB         |  49.07 us | 0.620 us | 0.549 us |

*BenchmarkDotNet v0.15.8, .NET 10.0.5 x64 RyuJIT x86-64-v4  
Windows 11 (10.0.26200.8037/25H2/2025Update/HudsonValley2)  
Intel Core i7-1065G7 CPU 1.30GHz 8 logical and 4 physical cores + 16 GB RAM

### Number of Evaluations
Func | bs | fp | ill | AB | ITP | Rid | Br | ModAB
-- | -- | -- | -- | -- | -- | -- | -- | --
Sum | 4463 | 8334 | 2950 | 3297 | 2853 | 2278 | 2903 | 1764
Ave | 47.9 | 88.4 | 31.6 | 33.6 | 30.4 | 24.5 | 31.3 | 18.9
Rel | 2.53 | 4.67 | 1.67 | 1.78 | 1.61 | 1.30 | 1.65 | 1.00
Max | 53 | 202 | 202 | 202 | 55 | 84 | 142 | 55
Mean | 46.4 | 52.7 | 19.4 | 18.0 | 24.1 | 19.1 | 18.3 | 15.1
Best | 13 | 8 | 10 | 52 | 6 | 5 | 22 | 38
Worst | 38 | 43 | 6 | 7 | 3 | 1 | 1 | 1
Median | 49.0 | 49.0 | 13.0 | 11.0 | 22.0 | 16.0 | 12.0 | 12.0
Std Dev | 7.6 | 79.4 | 40.8 | 51.0 | 19.0 | 20.4 | 40.0 | 15.1
Std Error | 0.79 | 8.19 | 4.25 | 5.00 | 1.97 | 2.12 | 4.17 | 1.57

&emsp;&emsp;[Detailed results](Benchmark%20results.md)

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

## 📚 References

1. Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332

2. Ganchovski N. Structural Analysis by Functional Modeling in the Cloud.  PhD Thesis **2025**, UACEG, Sofia

3. Ganchovski, N.; Traykov, A. (2023). "Modified Anderson-Björck's method for solving non-linear equations in structural mechanics." *IOP Conference Series: Materials Science and Engineering*, 1276(1), 012010. [DOI: 10.1088/1757-899X/1276/1/012010](https://doi.org/10.1088/1757-899X/1276/1/012010)

## 🐝Contributing

Contributions are welcome. You can open an issue or submit a pull request.
