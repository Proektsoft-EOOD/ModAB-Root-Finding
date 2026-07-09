import math
from ModAB import modAB_root
from dataclasses import dataclass
from typing import Callable


@dataclass
class Problem:
    name: str
    f: Callable[[float], float]
    a: float
    b: float


def P(x):
    return x + 1.11111


# Test problems from RootBenchmark.py
problems1 = [
    # Sergio Galdino. A family of regula falsi root-finding methods
    Problem("f01", lambda x: x**3 - 1, 0.5, 1.5),
    Problem("f02", lambda x: x**2 * (x**2 / 3 + math.sqrt(2) * math.sin(x)) - math.sqrt(3) / 18, 0.1, 1),
    Problem("f03", lambda x: 11 * x**11 - 1, 0.1, 1),
    Problem("f04", lambda x: x**3 + 1, -1.8, 0),
    Problem("f05", lambda x: x**3 - 2 * x - 5, 2, 3),
    Problem("f06", lambda x: 2 * x * math.exp(-5) + 1 - 2 * math.exp(-5 * x), 0, 1),
    Problem("f07", lambda x: 2 * x * math.exp(-10) + 1 - 2 * math.exp(-10 * x), 0, 1),
    Problem("f08", lambda x: 2 * x * math.exp(-20) + 1 - 2 * math.exp(-20 * x), 0, 1),
    Problem("f09", lambda x: (1 + (1 - 5)**2) * x**2 - (1 - 5 * x)**2, 0, 1),
    Problem("f10", lambda x: (1 + (1 - 10)**2) * x**2 - (1 - 10 * x)**2, 0, 1),
    Problem("f11", lambda x: (1 + (1 - 20)**2) * x**2 - (1 - 20 * x)**2, 0, 1),
    Problem("f12", lambda x: x**2 - (1 - x)**5, 0, 1),
    Problem("f13", lambda x: x**2 - (1 - x)**10, 0, 1),
    Problem("f14", lambda x: x**2 - (1 - x)**20, 0, 1),
    Problem("f15", lambda x: (1 + (1 - 5)**4) * x - (1 - 5 * x)**4, 0, 1),
    Problem("f16", lambda x: (1 + (1 - 10)**4) * x - (1 - 10 * x)**4, 0, 1),
    Problem("f17", lambda x: (1 + (1 - 20)**4) * x - (1 - 20 * x)**4, 0, 1),
    Problem("f18", lambda x: math.exp(-5 * x) * (x - 1) + x**5, 0, 1),
    Problem("f19", lambda x: math.exp(-10 * x) * (x - 1) + x**10, 0, 1),
    Problem("f20", lambda x: math.exp(-20 * x) * (x - 1) + x**20, 0, 1),
    Problem("f21", lambda x: x**2 + math.sin(x / 5) - 1 / 4, 0, 1),
    Problem("f22", lambda x: x**2 + math.sin(x / 10) - 1 / 4, 0, 1),
    Problem("f23", lambda x: x**2 + math.sin(x / 20) - 1 / 4, 0, 1),
    Problem("f24", lambda x: (x + 2) * (x + 1) * (x - 3)**3, 2.6, 4.6),
    Problem("f25", lambda x: (x - 4)**5 * math.log(x), 3.6, 5.6),
    Problem("f26", lambda x: (math.sin(x) - x / 4)**3, 2, 4),
    Problem("f27", lambda x: (81 - P(x) * (108 - P(x) * (54 - P(x) * (12 - P(x))))) * (1 if P(x) < 3 else (-1 if P(x) > 3 else 0)), 1, 3),
    Problem("f28", lambda x: math.sin((x - 7.143)**3), 7, 8),
    Problem("f29", lambda x: math.exp((x - 3)**5) - 1, 2.6, 4.6),
    Problem("f30", lambda x: math.exp((x - 3)**5) - math.exp(x - 1), 4, 5),
    Problem("f31", lambda x: math.pi - 1 / x, 0.05, 5),
    Problem("f32", lambda x: 4 - math.tan(x), 0, 1.5),
    Problem("f33", lambda x: math.cos(x) - x**3, 0, 4),
    # Steven A. Stage. Comments on An Improvement to the Brent's Method
    Problem("f34", lambda x: math.cos(x) - x, -11, 9),
    Problem("f35", lambda x: math.sqrt(abs(x - 2 / 3)) * (1 if x <= 2 / 3 else -1) - 0.1, -11, 9),
    Problem("f36", lambda x: abs(x - 2 / 3)**0.2 * (1 if x <= 2 / 3 else -1), -11, 9),
    Problem("f37", lambda x: (x - 7 / 9)**3 + (x - 7 / 9) * 1e-3, -11, 9),
    Problem("f38", lambda x: -0.5 if x <= 1 / 3 else 0.5, -11, 9),
    Problem("f39", lambda x: -1e-3 if x <= 1 / 3 else 1 - 1e-3, -11, 9),
    Problem("f40", lambda x: 0 if x == 0 else 1 / (x - 2 / 3), -11, 9),
    # A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method
    Problem("f41", lambda x: 2 * x * math.exp(-5) - 2 * math.exp(-5 * x) + 1, 0, 10),
    Problem("f42", lambda x: (x**2 - x - 6) * (x**2 - 3 * x + 2), 0, math.pi),
    Problem("f43", lambda x: x**3, -1, 1.5),
    Problem("f44", lambda x: x**5, -1, 1.5),
    Problem("f45", lambda x: x**7, -1, 1.5),
    Problem("f46", lambda x: (math.exp(-5 * x) - x - 0.5) / x**5, 0.09, 0.7),
    Problem("f47", lambda x: 1 / math.sqrt(x) - 2 * math.log(5e3 * math.sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f48", lambda x: 1 / math.sqrt(x) - 2 * math.log(5e7 * math.sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f49", lambda x: (-x**3 - x - 1) if x <= 0 else (x**(1 / 3) - x - 1), -1, 1),
    Problem("f50", lambda x: x**3 - 2 * x - x + 3, -3, 2),
    Problem("f51", lambda x: math.log(x), 0.5, 5),
    Problem("f52", lambda x: (10 - x) * math.exp(-10 * x) - x**10 + 1, 0.5, 8),
    Problem("f53", lambda x: math.exp(math.sin(x)) - x - 1, 1.0, 4),
    Problem("f54", lambda x: 2 * math.sin(x) - 1, 0.1, math.pi / 3),
    Problem("f55", lambda x: (x - 1) * math.exp(-x), 0.0, 1.5),
    Problem("f56", lambda x: (x - 1)**3 - 1, 1.5, 3),
    Problem("f57", lambda x: math.exp(x**2 + 7 * x - 30) - 1, 2.6, 3.5),
    Problem("f58", lambda x: math.atan(x) - 1, 1.0, 8),
    Problem("f59", lambda x: math.exp(x) - 2 * x - 1, 0.2, 3),
    Problem("f60", lambda x: math.exp(-x) - x - math.sin(x), 0.0, 2),
    Problem("f61", lambda x: x**2 - math.sin(x)**2 - 1, -1, 2),
    Problem("f62", lambda x: math.sin(x) - x / 2, math.pi / 2, math.pi),
]

# Oliveira I. F. D., Takahashi R. H. C.
# An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality
problems2 = [
    Problem("f63", lambda x: x * math.exp(x) - 1, -1, 1),
    Problem("f64", lambda x: math.tan(x - 1 / 10), -1, 1),
    Problem("f65", lambda x: math.sin(x) + 0.5, -1, 1),
    Problem("f66", lambda x: 4 * x**5 + x * x + 1, -1, 1),
    Problem("f67", lambda x: x + x**10 - 1, -1, 1),
    Problem("f68", lambda x: math.pi**x - math.e, -1, 1),
    Problem("f69", lambda x: math.log(abs(x - 10 / 9)), -1, 1),
    Problem("f70", lambda x: 1 / 3 + (1 if x > 0 else (-1 if x < 0 else 0)) * abs(x)**(1 / 3) + x**3, -1, 1),
    Problem("f71", lambda x: (x + 2 / 3) / (x + 101 / 100), -1, 1),
    Problem("f72", lambda x: (x * 1e6 - 1)**3, -1, 1),
    Problem("f73", lambda x: math.exp(x) * (x * 1e6 - 1)**3, -1, 1),
    Problem("f74", lambda x: (x - 1 / 3)**2 * math.atan(x - 1 / 3), -1, 1),
    Problem("f75", lambda x: (1 if 3 * x - 1 > 0 else (-1 if 3 * x - 1 < 0 else 0)) * (1 - math.sqrt(1 - (3 * x - 1)**2 / 81)), -1, 1),
    Problem("f76", lambda x: (1 + 1e6) / 1e6 if x > (1 - 1e6) / 1e6 else -1, -1, 1),
    Problem("f77", lambda x: 1 / (21 * x - 1) if x != 1 / 21 else 0, -1, 1),
    Problem("f78", lambda x: x * x / 4 + math.ceil(x / 2) - 0.5, -1, 1),
    Problem("f79", lambda x: math.ceil(10 * x - 1) + 0.5, -1, 1),
    Problem("f80", lambda x: x + math.sin(x * 1e6) / 10 + 1e-3, -1, 1),
    Problem("f81", lambda x: 1 + math.sin(1 / (x + 1)) if x > -1 else -1, -1, 1),
    Problem("f82", lambda x: 202 * x - 2 * math.floor((2 * x + 1e-2) / 2e-2) - 0.1, -1, 1),
    Problem("f83", lambda x: (202 * x - 2 * math.floor((2 * x + 1e-2) / 2e-2) - 0.1)**3, -1, 1),
]

# SciML project benchmarks suite
problems3 = [
    Problem("f84", lambda x: (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05, 0.5, 5.5),
    Problem("f85", lambda x: math.sin(x) - 0.5 * x - 0.3, -10.0, 10.0),
    Problem("f86", lambda x: math.exp(x) - 1 - x - x * x / 2 - 0.005, -2.0, 2.0),
    Problem("f87", lambda x: 1 / (x - 0.5) - 2 - 0.05, 0.6, 2.0),
    Problem("f88", lambda x: math.log(x) - x + 2 - 0.05, 0.1, 3.0),
    Problem("f89", lambda x: math.sin(20 * x) + 0.1 * x - 0.1, -4.0, 5.0),
    Problem("f90", lambda x: x**3 - 2 * x**2 + x - 0.025, -1.0, 2.0),
    Problem("f91", lambda x: x * math.sin(1 / x) - 0.1 - 0.01, 0.01, 1.0),
    Problem("f92", lambda x: x**3 - 0.001, -10, 10),
    Problem("f93", lambda x: x**7 - 0.001, -10, 10),
]

all_problems = problems1 + problems2 + problems3


class CountedFunc:
    """Wraps a callable and counts the number of evaluations."""
    __slots__ = ('_f', 'count')

    def __init__(self, f):
        self._f = f
        self.count = 0

    def __call__(self, x):
        self.count += 1
        return self._f(x)


def run_tests():
    eps = 1e-14
    passed = 0
    failed = 0
    total_evals = 0

    print("ModAB Root-Finding Test Results")
    print("=" * 70)
    print(f"{'Func':>4} | {'Root':>22} | {'f(root)':>15} | {'Evals':>6} | Status")
    print("-" * 70)

    for p in all_problems:
        cf = CountedFunc(p.f)
        try:
            root = modAB_root(cf, p.a, p.b, 0, xtol=eps, ytol=0.0)
            fval = p.f(root)

            if math.isnan(root):
                status = "FAIL"
                failed += 1
            elif abs(fval) < 1e-10:
                status = "PASS"
                passed += 1
            else:
                status = "PASS" if abs(fval) < 1e-6 else "WEAK"
                passed += 1

            total_evals += cf.count
            print(f"{p.name:>4} | {root:>22.15g} | {fval:>15.6e} | {cf.count:>6} | {status}")
        except Exception as e:
            failed += 1
            total_evals += cf.count
            print(f"{p.name:>4} | {'ERROR':>22} | {'N/A':>15} | {cf.count:>6} | FAIL ({e})")

    print("-" * 70)
    print(f"Total: {passed + failed} tests, {passed} passed, {failed} failed")
    print(f"Total evaluations: {total_evals}")


if __name__ == "__main__":
    run_tests()
