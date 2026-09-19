# Cython ModAB Implementation

## Building

To build the Cython extension:

```bash
python setup.py build_ext --inplace
```

This will generate a compiled extension module (`.pyd` on Windows, `.so` on Linux/Mac).

## Usage

```python
import cymodab

# Example: Find root of x^2 - 2 = 0 in interval [0, 2]
def f(x):
    return x * x - 2

root = cymodab.modAB_root(f, 0.0, 2.0)
print(f"Root: {root}")  # Should be approximately sqrt(2) ≈ 1.414
```

## Optimizations Applied

1. **C-level types**: All variables use `double` and `int` C types
2. **IEEE-strict optimization**: Compiled with `/O2 /fp:precise` (MSVC) or `-O3 -ffp-contract=off` (GCC/Clang). Fast-math is not used: the solver relies on NaN/inf semantics, and reassociation changes the iterates
3. **No Python overhead**:
   - Uses `libc.math.fabs` instead of Python's `abs()`
   - Uses `NAN` constant from C math library
   - Direct C comparisons and arithmetic
4. **Compiler directives**:
   - `boundscheck=False`: Disables array bounds checking
   - `wraparound=False`: Disables negative indexing
   - `cdivision=True`: Uses C division (faster, no zero-check)
   - `initializedcheck=False`: Assumes variables are initialized
5. **Inline functions**: Helper functions (`c_max`, `c_min`, `c_clamp`) are inlined
6. **Integer bisection flag**: Uses `int` instead of `bool` for better performance

## Performance

Expected speedup: **10-50x** compared to pure Python, depending on:
- Function evaluation cost (if `f()` is Python, it will dominate)
- Number of iterations
- Compiler optimizations available

For maximum performance, consider implementing the target function `f()` in Cython as well.

## Annotation

The build process creates `cymodab.html` showing which lines are pure C (white) vs Python interaction (yellow). Review this to identify any remaining bottlenecks.
