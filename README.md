
# ModAB Root-Finding Algorithm

This repository provides efficient implementations of the **Modified Anderson-Björck (ModAB)** bracketing root-finding algorithm for solving nonlinear equations of the form `f(x) = 0`.

Available in **C**, **C#**, **Python**, **Julia**, **Zig**, **Rust** and **C++ (ROOT.CERN framework)**.

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

&emsp;&emsp;[For 92 functions test set](/MDPI-Algorithms-Data/Test%20Functions.pdf)

<img height="300" alt="image" src="https://github.com/user-attachments/assets/cd8b63e9-15d6-4cb0-863a-ccdd445de80e" />

### Times
Algorithm | Mean | Error | StdDev | Relative
-- | -- | -- | -- | --
Bisection | 86.57 us | 0.832 us | 0.738 us | 1,75
FalsePosition | 321.77 us | 3.561 us | 2.974 us | 6,51
Illinois | 104.75 us | 1.353 us | 1.265 us | 2,12
Anderson Björck | 122.19 us | 1.205 us | 1.127 us | 2,47
ITP | 193.37 us | 1.604 us | 1.339 us | 3,92
Ridders | 66.61 us | 1.275 us | 1.309 us | 1,35
Brent | 119.56 us | 1.213 us | 1.134 us | 2,42
modAB | 51.58 us | 1.013 us | 1.352 us | 1,04

*BenchmarkDotNet v0.15.8, .NET 10.0.5 x64 RyuJIT x86-64-v4  
Windows 11 (10.0.26200.8037/25H2/2025Update/HudsonValley2)  
Intel Core i7-1065G7 CPU 1.30GHz 8 logical and 4 physical cores + 16 GB RAM

### Number of Evaluations
Func | bs | fp | ill | AB | ITP | Rid | Br | ModAB
-- | -- | -- | -- | -- | -- | -- | -- | --
Sum | 4410 | 8132 | 2907 | 3095 | 2798 | 2256 | 2880 | 1741
Ave | 47.9 | 88.4 | 31.6 | 33.6 | 30.4 | 24.5 | 31.3 | 18.9
Rel | 2.53 | 4.67 | 1.67 | 1.78 | 1.61 | 1.30 | 1.65 | 1.00
Max | 53 | 202 | 202 | 202 | 55 | 84 | 142 | 55
Mean | 46.3 | 52.0 | 19.2 | 17.5 | 23.9 | 19.1 | 18.3 | 15.0
Best | 13 | 8 | 10 | 52 | 6 | 5 | 22 | 38
Worst | 38 | 43 | 6 | 7 | 3 | 1 | 1 | 1
Median | 49.0 | 48.5 | 13.0 | 11.0 | 22.0 | 16.0 | 12.0 | 12.0
Std Dev | 7.6 | 79.0 | 41.0 | 48.2 | 19.0 | 20.5 | 40.2 | 15.1
Std Error | 0.80 | 8.23 | 4.27 | 5.03 | 1.98 | 2.13 | 4.19 | 1.58

&emsp;&emsp;[Detailed results](Benchmark%20results.md)

## 📄 License

MIT License - see the [LICENSE](C/LICENSE) file for details.

## 🧮 Implementations in other software/libraries

Calcpad - 					https://calcpad.eu   - C#  
Root-Fortran - 				https://github.com/jacobwilliams/roots-fortran	- Fortran  
ROOT.CERN - 				https://github.com/root-project/root			- C++  
SCiML/NonlinearSolve.jl - 	https://github.com/SciML/NonlinearSolve.jl		- Julia  
JuliaMath/Roots.jl -		https://github.com/JuliaMath/Roots.jl			- Julia  
MultiFloats.jl - 			https://github.com/dzhang314/MultiFloats.jl		- Julia  
PyModAB - 			https://pypi.org/project/pymodab/		- Python/C

## 📚 References

1. Ganchovski N. Structural Analysis by Functional Modeling in the Cloud.  PhD Thesis **2025**, UACEG, Sofia

2. Ganchovski, N.; Traykov, A. (2023). "Modified Anderson-Björck's method for solving non-linear equations in structural mechanics." *IOP Conference Series: Materials Science and Engineering*, 1276(1), 012010. [DOI: 10.1088/1757-899X/1276/1/012010](https://doi.org/10.1088/1757-899X/1276/1/012010)

3. Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332

## 🐝Contributing

Contributions are welcome. You can open an issue or submit a pull request.
