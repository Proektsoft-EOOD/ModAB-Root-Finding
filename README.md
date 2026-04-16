
# ModAB Root-Finding Algorithm

This repository provides efficient implementations of the **Modified Anderson-Björck (ModAB)** bracketing root-finding algorithm for solving nonlinear equations of the form `f(x) = 0`.

Available in **C**, **C#**,  **Python**, **Julia**, **Zig** and **C++ (ROOT.CERN framework)**.

## Features

- **Guaranteed convergence** - Bracketing approach ensures the root is always found.
- **Great performance** - Fewer function evaluations than Brent, Ridders, and other popular methods. Being simple and lightweight, it has very little computational overhead per iteration.
- **Worst-case optimality** - Retains the worst-case optimality of Bisection for the hard cases.
- **Cross-platform** - Python package that works on Windows, Linux, and macOS.
- **Multiple languages** - Use in your preferred environment.
- **Extensive testing and benchmarking** - The algorithm is benchmark against a broad set of test functions of many different types.

## Algorithm Overview

The ModAB algorithm combines the reliability of bisection with the speed of the secant method:

1. **Bisection phase** - Starts with bisection to ensure stability
2. **Linearity detection** - Monitors when the function behavior is close enough to linear
3. **False-position acceleration** - Switches to false-position method when conditions are favorable
4. **Anderson-Björck correction** - Applies A&B corrections to prevent stalling
5. **Adaptive fallback** - Returns to bisection if progress slows

This hybrid approach achieves superlinear convergence while maintaining the worst-case optimality of bisection.

## License

MIT License - see the [LICENSE](C/LICENSE) file for details.

## Implementations in other software/libraries

Calcpad - 					https://calcpad.eu   - C#  
Root-Fortran - 				https://github.com/jacobwilliams/roots-fortran	- Fortran  
ROOT.CERN - 				https://github.com/root-project/root			- C++  
SCiML/NonlinearSolve.jl - 	https://github.com/SciML/NonlinearSolve.jl		- Julia  
JuliaMath/Roots.jl -		https://github.com/JuliaMath/Roots.jl			- Julia  
MultiFloats.jl - 			https://github.com/dzhang314/MultiFloats.jl		- Julia  

## References

1. Ganchovski N. Structural Analysis by Functional Modelling in the Cloud.  PhD Thesis **2025**, UACEG, Sofia

2. Ganchovski, N.; Traykov, A. (2023). "Modified Anderson-Björck's method for solving non-linear equations in structural mechanics." *IOP Conference Series: Materials Science and Engineering*, 1276(1), 012010. [DOI: 10.1088/1757-899X/1276/1/012010](https://doi.org/10.1088/1757-899X/1276/1/012010)

3. Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. (2026). "Improvements of the Modified Anderson-Björck (modAB) Root-Finding Algorithm." *Preprints*, 2026032190. [DOI: 10.20944/preprints202603.2190.v1](https://doi.org/10.20944/preprints202603.2190.v1)

## Contributing

Contributions are welcome. You can open an issue or submit a pull request.
