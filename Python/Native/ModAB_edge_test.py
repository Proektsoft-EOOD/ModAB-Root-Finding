"""
Edge-case tests for the shared modAB safeguards.

The 100-problem benchmark suite never produces a non-finite or overflowing
residual, so the overflow/NaN branches of same_nonzero_sign, safe_midpoint,
safe_secant, symmetry_factor and passes_switching_test are covered here.

Expected values come from the C# reference in C#/Root/Node.cs and
C#/Root/Solvers/ModABCorr.cs; every language port is held to the same table.
"""
import math
import sys

INF = math.inf
NAN = math.nan

from ModAB import (modAB_root, same_nonzero_sign, safe_midpoint, safe_secant,
                   symmetry_factor, passes_switching_test)

passed = 0
total = 0


def ck(name, i, got, want):
    global passed, total
    total += 1
    if got == want or (isinstance(got, float) and isinstance(want, float)
                       and math.isnan(got) and math.isnan(want)):
        passed += 1
        return
    print(f"FAIL {name}[{i}]: got {got!r} want {want!r}")


# --- same_nonzero_sign ---
ck("same_nonzero_sign", 0, bool(same_nonzero_sign(1.0, 2.0)), True)
ck("same_nonzero_sign", 1, bool(same_nonzero_sign(-1.0, -2.0)), True)
ck("same_nonzero_sign", 2, bool(same_nonzero_sign(1.0, -2.0)), False)
ck("same_nonzero_sign", 3, bool(same_nonzero_sign(-1.0, 2.0)), False)
ck("same_nonzero_sign", 4, bool(same_nonzero_sign(0.0, 1.0)), False)
ck("same_nonzero_sign", 5, bool(same_nonzero_sign(1.0, 0.0)), False)
ck("same_nonzero_sign", 6, bool(same_nonzero_sign(0.0, 0.0)), False)
ck("same_nonzero_sign", 7, bool(same_nonzero_sign(-0.0, -1.0)), False)
ck("same_nonzero_sign", 8, bool(same_nonzero_sign(NAN, 1.0)), False)
ck("same_nonzero_sign", 9, bool(same_nonzero_sign(1.0, NAN)), False)
ck("same_nonzero_sign", 10, bool(same_nonzero_sign(NAN, NAN)), False)
ck("same_nonzero_sign", 11, bool(same_nonzero_sign(INF, 1.0)), True)
ck("same_nonzero_sign", 12, bool(same_nonzero_sign(-INF, -1.0)), True)
ck("same_nonzero_sign", 13, bool(same_nonzero_sign(INF, -INF)), False)

# --- safe_midpoint ---
ck("safe_midpoint", 0, safe_midpoint(2.0, 4.0), 3.0)
ck("safe_midpoint", 1, safe_midpoint(1e+308, 1e+308), 1e+308)
ck("safe_midpoint", 2, safe_midpoint(-1e+308, 1e+308), 0.0)
ck("safe_midpoint", 3, safe_midpoint(1e+308, 1.7e+308), 1.35e+308)
ck("safe_midpoint", 4, safe_midpoint(-1.7e+308, -1e+308), -1.35e+308)
ck("safe_midpoint", 5, safe_midpoint(0.0, 1.0), 0.5)
ck("safe_midpoint", 6, safe_midpoint(-1.0, 1.0), 0.0)

# --- safe_secant ---
ck("safe_secant", 0, safe_secant(0.0, -1.0, 1.0, 1.0), 0.5)
ck("safe_secant", 1, safe_secant(0.0, -1.0, 1.0, 3.0), 0.25)
ck("safe_secant", 2, safe_secant(0.0, -1e+308, 1.0, 1e+308), 0.5)
ck("safe_secant", 3, safe_secant(0.0, -1.7e+308, 1.0, 1.7e+308), 0.5)
ck("safe_secant", 4, safe_secant(0.0, INF, 1.0, -1.0), 0.5)
ck("safe_secant", 5, safe_secant(0.0, -1.0, 1.0, INF), 0.5)
ck("safe_secant", 6, safe_secant(0.0, 0.0, 1.0, 0.0), 0.5)
ck("safe_secant", 7, safe_secant(0.0, NAN, 1.0, 1.0), 0.5)
ck("safe_secant", 8, safe_secant(0.0, -1e-300, 1.0, 1e+300), 0.0)
ck("safe_secant", 9, safe_secant(0.0, -1e+300, 1.0, 1e-300), 1.0)
ck("safe_secant", 10, safe_secant(1e+308, -1.0, 1.7e+308, 1.0), 1.35e+308)

# --- symmetry_factor ---
ck("symmetry_factor", 0, symmetry_factor(-1.0, 1.0), 1.0)
ck("symmetry_factor", 1, symmetry_factor(-1.0, 3.0), 0.5625)
ck("symmetry_factor", 2, symmetry_factor(-3.0, 1.0), 0.5625)
ck("symmetry_factor", 3, symmetry_factor(-1e+308, 1e+308), 1.0)
ck("symmetry_factor", 4, symmetry_factor(-1.7e+308, 1.0), 0.25)
ck("symmetry_factor", 5, symmetry_factor(INF, -1.0), NAN)
ck("symmetry_factor", 6, symmetry_factor(-1.0, INF), NAN)
ck("symmetry_factor", 7, symmetry_factor(NAN, 1.0), NAN)
ck("symmetry_factor", 8, symmetry_factor(-1e-300, 1e+300), 0.25)

# --- passes_switching_test ---
ck("passes_switching_test", 0, bool(passes_switching_test(1.0, 1.0, 0.5)), True)
ck("passes_switching_test", 1, bool(passes_switching_test(1.0, -1.0, 0.5)), False)
ck("passes_switching_test", 2, bool(passes_switching_test(1.0, 0.5, 0.5)), True)
ck("passes_switching_test", 3, bool(passes_switching_test(1.0, 0.5, 0.1)), False)
ck("passes_switching_test", 4, bool(passes_switching_test(1e+308, 1e+308, 0.5)), True)
ck("passes_switching_test", 5, bool(passes_switching_test(1.7e+308, -1.7e+308, 0.5)), False)
ck("passes_switching_test", 6, bool(passes_switching_test(1e+308, 1.6e+308, 0.5)), True)
ck("passes_switching_test", 7, bool(passes_switching_test(INF, 1.0, 0.5)), False)
ck("passes_switching_test", 8, bool(passes_switching_test(1.0, INF, 0.5)), False)
ck("passes_switching_test", 9, bool(passes_switching_test(NAN, 1.0, 0.5)), False)
ck("passes_switching_test", 10, bool(passes_switching_test(1.0, 1.0, NAN)), False)
ck("passes_switching_test", 11, bool(passes_switching_test(0.0, 0.0, 0.5)), False)

# --- solver-level regression: x1 + x2 overflows, this returned inf before ---
total += 1
_root = modAB_root(lambda x: x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200)
if abs(_root - 1.2e308) <= 1e-13 * 1.2e308:
    passed += 1
else:
    print(f"FAIL overflowing bracket: got {_root!r} want 1.2e308")

print(f"Python edge-cases: {passed}/{total} "
      f"{'PASS' if passed == total else 'FAIL'}")
sys.exit(0 if passed == total else 1)
