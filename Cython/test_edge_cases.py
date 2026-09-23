"""
Edge-case tests for the shared modAB safeguards.

The 100-problem benchmark suite never produces a non-finite or overflowing
residual, so the overflow/NaN branches of same_sign, safe_midpoint and
safe_secant, and the overflow cases of the switching test, are covered here.

Expected values come from the C# reference in C#/Root/Node.cs and
C#/Root/Solvers/SgModAB.cs; every language port is held to the same table.
"""
import math
import sys

INF = math.inf
NAN = math.nan

import cymodab

# The safeguards are `cdef inline` for speed; these are the test hooks.
modAB_root = cymodab.modAB_root
same_sign = cymodab.t_same_sign
safe_midpoint = cymodab.t_safe_midpoint
safe_secant = cymodab.t_safe_secant

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


# --- same_sign ---
ck("same_sign", 0, bool(same_sign(1.0, 2.0)), True)
ck("same_sign", 1, bool(same_sign(-1.0, -2.0)), True)
ck("same_sign", 2, bool(same_sign(1.0, -2.0)), False)
ck("same_sign", 3, bool(same_sign(-1.0, 2.0)), False)
ck("same_sign", 4, bool(same_sign(0.0, 1.0)), False)
ck("same_sign", 5, bool(same_sign(1.0, 0.0)), False)
ck("same_sign", 6, bool(same_sign(0.0, 0.0)), False)
ck("same_sign", 7, bool(same_sign(-0.0, -1.0)), False)
ck("same_sign", 8, bool(same_sign(NAN, 1.0)), False)
ck("same_sign", 9, bool(same_sign(1.0, NAN)), False)
ck("same_sign", 10, bool(same_sign(NAN, NAN)), False)
ck("same_sign", 11, bool(same_sign(INF, 1.0)), True)
ck("same_sign", 12, bool(same_sign(-INF, -1.0)), True)
ck("same_sign", 13, bool(same_sign(INF, -INF)), False)

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

# --- solver-level regression: x1 + x2 overflows, this returned inf before ---
total += 1
_root = modAB_root(lambda x: x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200)
if abs(_root - 1.2e308) <= 1e-13 * 1.2e308:
    passed += 1
else:
    print(f"FAIL overflowing bracket: got {_root!r} want 1.2e308")

# --- solver-level: overflow and infinite residuals in the switching test ---
_M = 1.7976931348623157e308


def ck_root(name, f, a, b, want, rtol=1e-13):
    global passed, total
    total += 1
    got = modAB_root(f, a, b, 0.0, 1e-14, 0.0, 200)
    if abs(got - want) <= rtol * max(abs(want), 1.0):
        passed += 1
    else:
        print(f"FAIL {name}: got {got!r} want {want!r}")


# |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite
ck_root("hump", lambda x: _M * (-0.2 + 1.2 * (x + 1)) if x <= 0 else _M * (1 - 0.2 * x),
        -1.0, 1.0, -5.0 / 6.0)
# f2 - f1 overflows: switching is disabled until the residuals shrink
ck_root("f2-f1 overflow", lambda x: 1.7e308 * math.tanh(10 * (x - 0.3)), -1.0, 1.0, 0.3)
# infinite residuals at one or both ends of the bracket
ck_root("inf left", lambda x: -INF if x < -0.5 else x - 0.1, -1.0, 1.0, 0.1)
ck_root("inf both", lambda x: -INF if x < -0.5 else (INF if x > 0.9 else x - 0.1),
        -1.0, 1.0, 0.1)
# subnormal residuals: AB corrections underflow towards zero
ck_root("subnormal", lambda x: 1e-300 * (x ** 3 - 0.2), -1.0, 1.0, 0.2 ** (1.0 / 3.0))

print(f"Cython edge-cases: {passed}/{total} "
      f"{'PASS' if passed == total else 'FAIL'}")
sys.exit(0 if passed == total else 1)
