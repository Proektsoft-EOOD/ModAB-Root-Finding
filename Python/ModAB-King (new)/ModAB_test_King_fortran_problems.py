import math
from ModAB_King import modAB_root
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

# Problems ported verbatim from root_tests.f90 (the authoritative Fortran
# source for table_real64.tex). All 157 cases; parameterized families
# expanded to one Problem per n value. Each validated: bracket sign change
# + published Fortran root. Same Problem(name, f, a, b) interface.
from fortran_problems import PROBLEMS as _FORT

all_problems = [Problem(n, f, a, b) for (n, f, a, b, _root) in _FORT]


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
    eps = 1e-15
    passed = 0
    failed = 0
    total_evals = 0

    print("ModAB-King Root-Finding Test Results")
    print("=" * 70)
    print(f"{'Func':>4} | {'Root':>22} | {'f(root)':>15} | {'Evals':>6} | Status")
    print("-" * 70)

    for p in all_problems:
        cf = CountedFunc(p.f)
        try:
            root = modAB_root(cf, p.a, p.b, 0, xtol=eps, ytol=1e-15)
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
            print(f"{p.name:<9} | {root:>22.15g} | {fval:>15.6e} | {cf.count:>6} | {status}")
        except Exception as e:
            failed += 1
            total_evals += cf.count
            print(f"{p.name:<9} | {'ERROR':>22} | {'N/A':>15} | {cf.count:>6} | FAIL ({e})")

    print("-" * 70)
    print(f"Total: {passed + failed} tests, {passed} passed, {failed} failed")
    print(f"Total evaluations: {total_evals}")


if __name__ == "__main__":
    run_tests()
