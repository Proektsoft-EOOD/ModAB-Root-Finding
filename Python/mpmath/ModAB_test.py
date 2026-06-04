import math
import time
from dataclasses import dataclass
from typing import Callable
from mpmath import mp, mpf, findroot


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
    Problem("f02", lambda x: x**2 * (x**2 / 3 + mp.sqrt(2) * mp.sin(x)) - mp.sqrt(3) / 18, 0.1, 1),
    Problem("f03", lambda x: 11 * x**11 - 1, 0.1, 1),
    Problem("f04", lambda x: x**3 + 1, -1.8, 0),
    Problem("f05", lambda x: x**3 - 2 * x - 5, 2, 3),
    Problem("f06", lambda x: 2 * x * mp.exp(-5) + 1 - 2 * mp.exp(-5 * x), 0, 1),
    Problem("f07", lambda x: 2 * x * mp.exp(-10) + 1 - 2 * mp.exp(-10 * x), 0, 1),
    Problem("f08", lambda x: 2 * x * mp.exp(-20) + 1 - 2 * mp.exp(-20 * x), 0, 1),
    Problem("f09", lambda x: (1 + (1 - 5)**2) * x**2 - (1 - 5 * x)**2, 0, 1),
    Problem("f10", lambda x: (1 + (1 - 10)**2) * x**2 - (1 - 10 * x)**2, 0, 1),
    Problem("f11", lambda x: (1 + (1 - 20)**2) * x**2 - (1 - 20 * x)**2, 0, 1),
    Problem("f12", lambda x: x**2 - (1 - x)**5, 0, 1),
    Problem("f13", lambda x: x**2 - (1 - x)**10, 0, 1),
    Problem("f14", lambda x: x**2 - (1 - x)**20, 0, 1),
    Problem("f15", lambda x: (1 + (1 - 5)**4) * x - (1 - 5 * x)**4, 0, 1),
    Problem("f16", lambda x: (1 + (1 - 10)**4) * x - (1 - 10 * x)**4, 0, 1),
    Problem("f17", lambda x: (1 + (1 - 20)**4) * x - (1 - 20 * x)**4, 0, 1),
    Problem("f18", lambda x: mp.exp(-5 * x) * (x - 1) + x**5, 0, 1),
    Problem("f19", lambda x: mp.exp(-10 * x) * (x - 1) + x**10, 0, 1),
    Problem("f20", lambda x: mp.exp(-20 * x) * (x - 1) + x**20, 0, 1),
    Problem("f21", lambda x: x**2 + mp.sin(x / 5) - mpf('0.25'), 0, 1),
    Problem("f22", lambda x: x**2 + mp.sin(x / 10) - mpf('0.25'), 0, 1),
    Problem("f23", lambda x: x**2 + mp.sin(x / 20) - mpf('0.25'), 0, 1),
    Problem("f24", lambda x: (x + 2) * (x + 1) * (x - 3)**3, 2.6, 4.6),
    Problem("f25", lambda x: (x - 4)**5 * mp.log(x), 3.6, 5.6),
    Problem("f26", lambda x: (mp.sin(x) - x / 4)**3, 2, 4),
    Problem("f27", lambda x: (81 - P(x) * (108 - P(x) * (54 - P(x) * (12 - P(x))))) * (1 if P(x) < 3 else (-1 if P(x) > 3 else 0)), 1, 3),
    Problem("f28", lambda x: mp.sin((x - 7.143)**3), 7, 8),
    Problem("f29", lambda x: mp.exp((x - 3)**5) - 1, 2.6, 4.6),
    Problem("f30", lambda x: mp.exp((x - 3)**5) - mp.exp(x - 1), 4, 5),
    Problem("f31", lambda x: mp.pi - 1 / x, 0.05, 5),
    Problem("f32", lambda x: 4 - mp.tan(x), 0, 1.5),
    Problem("f33", lambda x: mp.cos(x) - x**3, 0, 4),
    # Steven A. Stage. Comments on An Improvement to the Brent's Method
    Problem("f34", lambda x: mp.cos(x) - x, -11, 9),
    Problem("f35", lambda x: mp.sqrt(abs(x - mpf('2') / 3)) * (1 if x <= mpf('2') / 3 else -1) - mpf('0.1'), -11, 9),
    Problem("f36", lambda x: abs(x - mpf('2') / 3)**mpf('0.2') * (1 if x <= mpf('2') / 3 else -1), -11, 9),
    Problem("f37", lambda x: (x - mpf('7') / 9)**3 + (x - mpf('7') / 9) * mpf('1e-3'), -11, 9),
    Problem("f38", lambda x: mpf('-0.5') if x <= mpf('1') / 3 else mpf('0.5'), -11, 9),
    Problem("f39", lambda x: mpf('-1e-3') if x <= mpf('1') / 3 else 1 - mpf('1e-3'), -11, 9),
    Problem("f40", lambda x: mpf('0') if x == 0 else 1 / (x - mpf('2') / 3), -11, 9),
    # A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method
    Problem("f41", lambda x: 2 * x * mp.exp(-5) - 2 * mp.exp(-5 * x) + 1, 0, 10),
    Problem("f42", lambda x: (x**2 - x - 6) * (x**2 - 3 * x + 2), 0, float(mp.pi)),
    Problem("f43", lambda x: x**3, -1, 1.5),
    Problem("f44", lambda x: x**5, -1, 1.5),
    Problem("f45", lambda x: x**7, -1, 1.5),
    Problem("f46", lambda x: (mp.exp(-5 * x) - x - mpf('0.5')) / x**5, 0.09, 0.7),
    Problem("f47", lambda x: 1 / mp.sqrt(x) - 2 * mp.log(5e3 * mp.sqrt(x)) + mpf('0.8'), 0.0005, 0.5),
    Problem("f48", lambda x: 1 / mp.sqrt(x) - 2 * mp.log(5e7 * mp.sqrt(x)) + mpf('0.8'), 0.0005, 0.5),
    Problem("f49", lambda x: (-x**3 - x - 1) if x <= 0 else (x**mpf('0.333333333333333') - x - 1), -1, 1),
    Problem("f50", lambda x: x**3 - 2 * x - x + 3, -3, 2),
    Problem("f51", lambda x: mp.log(x), 0.5, 5),
    Problem("f52", lambda x: (10 - x) * mp.exp(-10 * x) - x**10 + 1, 0.5, 8),
    Problem("f53", lambda x: mp.exp(mp.sin(x)) - x - 1, 1.0, 4),
    Problem("f54", lambda x: 2 * mp.sin(x) - 1, 0.1, float(mp.pi) / 3),
    Problem("f55", lambda x: (x - 1) * mp.exp(-x), 0.0, 1.5),
    Problem("f56", lambda x: (x - 1)**3 - 1, 1.5, 3),
    Problem("f57", lambda x: mp.exp(x**2 + 7 * x - 30) - 1, 2.6, 3.5),
    Problem("f58", lambda x: mp.atan(x) - 1, 1.0, 8),
    Problem("f59", lambda x: mp.exp(x) - 2 * x - 1, 0.2, 3),
    Problem("f60", lambda x: mp.exp(-x) - x - mp.sin(x), 0.0, 2),
    Problem("f61", lambda x: x**2 - mp.sin(x)**2 - 1, -1, 2),
    Problem("f62", lambda x: mp.sin(x) - x / 2, float(mp.pi) / 2, float(mp.pi)),
]

# Oliveira I. F. D., Takahashi R. H. C.
# An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality
problems2 = [
    Problem("f63", lambda x: x * mp.exp(x) - 1, -1, 1),
    Problem("f64", lambda x: mp.tan(x - mpf('0.1')), -1, 1),
    Problem("f65", lambda x: mp.sin(x) + mpf('0.5'), -1, 1),
    Problem("f66", lambda x: 4 * x**5 + x * x + 1, -1, 1),
    Problem("f67", lambda x: x + x**10 - 1, -1, 1),
    Problem("f68", lambda x: mp.pi**x - mp.e, -1, 1),
    Problem("f69", lambda x: mp.log(abs(x - mpf('10') / 9)), -1, 1),
    Problem("f70", lambda x: mpf('1') / 3 + (1 if x > 0 else (-1 if x < 0 else 0)) * abs(x)**mpf('0.333333333333333') + x**3, -1, 1),
    Problem("f71", lambda x: (x + mpf('2') / 3) / (x + mpf('101') / 100), -1, 1),
    Problem("f72", lambda x: (x * mpf('1e6') - 1)**3, -1, 1),
    Problem("f73", lambda x: mp.exp(x) * (x * mpf('1e6') - 1)**3, -1, 1),
    Problem("f74", lambda x: (x - mpf('1') / 3)**2 * mp.atan(x - mpf('1') / 3), -1, 1),
    Problem("f75", lambda x: (1 if 3 * x - 1 > 0 else (-1 if 3 * x - 1 < 0 else 0)) * (1 - mp.sqrt(1 - (3 * x - 1)**2 / 81)), -1, 1),
    Problem("f76", lambda x: (1 + mpf('1e6')) / mpf('1e6') if x > (1 - mpf('1e6')) / mpf('1e6') else mpf('-1'), -1, 1),
    Problem("f77", lambda x: 1 / (21 * x - 1) if x != mpf('1') / 21 else mpf('0'), -1, 1),
    Problem("f78", lambda x: x * x / 4 + math.ceil(float(x) / 2) - mpf('0.5'), -1, 1),
    Problem("f79", lambda x: math.ceil(10 * float(x) - 1) + mpf('0.5'), -1, 1),
    Problem("f80", lambda x: x + mp.sin(x * mpf('1e6')) / 10 + mpf('1e-3'), -1, 1),
    Problem("f81", lambda x: 1 + mp.sin(1 / (x + 1)) if x > -1 else mpf('-1'), -1, 1),
    Problem("f82", lambda x: 202 * x - 2 * math.floor((2 * float(x) + 1e-2) / 2e-2) - mpf('0.1'), -1, 1),
    Problem("f83", lambda x: (202 * x - 2 * math.floor((2 * float(x) + 1e-2) / 2e-2) - mpf('0.1'))**3, -1, 1),
]

# SciML project benchmarks suite
problems3 = [
    Problem("f84", lambda x: (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - mpf('0.05'), 0.5, 5.5),
    Problem("f85", lambda x: mp.sin(x) - mpf('0.5') * x - mpf('0.3'), -10.0, 10.0),
    Problem("f86", lambda x: mp.exp(x) - 1 - x - x * x / 2 - mpf('0.005'), -2.0, 2.0),
    Problem("f87", lambda x: 1 / (x - mpf('0.5')) - 2 - mpf('0.05'), 0.6, 2.0),
    Problem("f88", lambda x: mp.log(x) - x + 2 - mpf('0.05'), 0.1, 3.0),
    Problem("f89", lambda x: mp.sin(20 * x) + mpf('0.1') * x - mpf('0.1'), -4.0, 5.0),
    Problem("f90", lambda x: x**3 - 2 * x**2 + x - mpf('0.025'), -1.0, 2.0),
    Problem("f91", lambda x: x * mp.sin(1 / x) - mpf('0.1') - mpf('0.01'), 0.01, 1.0),
    Problem("f92", lambda x: x**3 - mpf('0.001'), -10, 10),
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

    def reset(self):
        self.count = 0


# Bracketing solvers that can use intervals [a, b]
BRACKETING_SOLVERS = ['bisect', 'illinois', 'pegasus', 'anderson', 'ridder', 'modAB']


def run_single_solver(solver_name, problems, tol=1e-20):
    """Run a single solver on all problems and return results."""
    results = []

    for p in problems:
        cf = CountedFunc(p.f)
        try:
            root = findroot(cf, (mpf(p.a), mpf(p.b)), solver=solver_name, tol=tol, verify=False, maxsteps=400)
            fval = float(p.f(root))

            if math.isnan(float(root)):
                status = "FAIL"
            elif abs(fval) < 1e-18:
                status = "PASS"
            elif abs(fval) < 1e-12:
                status = "WEAK"
            else:
                status = "FAIL"

            results.append({
                'name': p.name,
                'root': float(root),
                'fval': fval,
                'evals': cf.count,
                'status': status
            })
        except Exception as e:
            results.append({
                'name': p.name,
                'root': None,
                'fval': None,
                'evals': cf.count,
                'status': "FAIL",
                'error': str(e)
            })

    return results


def run_benchmark():
    """Run benchmark comparing all bracketing solvers."""
    mp.dps = 30  # Set precision

    print("mpmath Root-Finding Benchmark")
    print("=" * 100)
    print(f"Precision: {mp.dps} decimal places")
    print(f"Number of test problems: {len(all_problems)}")
    print()

    # Collect results for all solvers
    all_results = {}
    all_times = {}
    for solver in BRACKETING_SOLVERS:
        print(f"Running {solver}...", end=" ", flush=True)
        start_time = time.perf_counter()
        all_results[solver] = run_single_solver(solver, all_problems)
        elapsed = time.perf_counter() - start_time
        all_times[solver] = elapsed
        passed = sum(1 for r in all_results[solver] if r['status'] in ['PASS', 'WEAK'])
        total_evals = sum(r['evals'] for r in all_results[solver])
        print(f"done ({passed}/{len(all_problems)} passed, {total_evals} evals, {elapsed:.3f}s)")

    print()

    # Print summary table
    print("=" * 100)
    print("SUMMARY: Function evaluations by solver")
    print("=" * 100)

    # Header
    header = f"{'Func':>4} |"
    for solver in BRACKETING_SOLVERS:
        header += f" {solver:>10} |"
    print(header)
    print("-" * 100)

    # Data rows
    for i, p in enumerate(all_problems):
        row = f"{p.name:>4} |"
        for solver in BRACKETING_SOLVERS:
            result = all_results[solver][i]
            status_letter = ' ' if result['status'] == 'PASS' else result['status'][0]  # space, W, or F
            row += f" {result['evals']:>9}{status_letter} |"
        print(row)

    print("-" * 100)

    # Totals row
    totals_row = "Total|"
    for solver in BRACKETING_SOLVERS:
        total_evals = sum(r['evals'] for r in all_results[solver])
        totals_row += f" {total_evals:>10} |"
    print(totals_row)

    # Pass count row
    pass_row = "Pass |"
    for solver in BRACKETING_SOLVERS:
        passed = sum(1 for r in all_results[solver] if r['status'] in ['PASS', 'WEAK'])
        pass_row += f" {passed:>10} |"
    print(pass_row)

    # Fail count row
    fail_row = "Fail |"
    for solver in BRACKETING_SOLVERS:
        failed = sum(1 for r in all_results[solver] if r['status'] == 'FAIL')
        fail_row += f" {failed:>10} |"
    print(fail_row)

    # Time row
    time_row = "Time |"
    for solver in BRACKETING_SOLVERS:
        time_row += f" {all_times[solver]:>9.3f}s |"
    print(time_row)

    print("=" * 100)

    # Print average evaluations for passed problems
    print()
    print("Average evaluations per solve:")
    for solver in BRACKETING_SOLVERS:
            avg = sum(r['evals'] for r in all_results[solver]) / len(all_results[solver])
            print(f"  {solver:>10}: {avg:.2f}")



def run_detailed_test(solver_name='modAB'):
    """Run detailed test for a single solver."""
    mp.dps = 30

    print(f"Detailed Test Results for: {solver_name}")
    print("=" * 70)
    print(f"{'Func':>4} | {'Root':>22} | {'f(root)':>15} | Evals")
    print("-" * 70)

    passed = 0
    failed = 0
    total_evals = 0

    for p in all_problems:
        cf = CountedFunc(p.f)
        try:
            root = findroot(cf, (mpf(p.a), mpf(p.b)), solver=solver_name, verify=False, maxsteps=400)
            fval = float(p.f(root))

            if math.isnan(float(root)):
                status = "F"
                failed += 1
            elif abs(fval) < 1e-18:
                status = " "
                passed += 1
            elif abs(fval) < 1e-12:
                status = "W"
                passed += 1
            else:
                status = "F"
                failed += 1

            total_evals += cf.count
            print(f"{p.name:>4} | {float(root):>22.15g} | {fval:>15.6e} | {cf.count:>5}{status}")
        except Exception as e:
            failed += 1
            total_evals += cf.count
            print(f"{p.name:>4} | {'ERROR':>22} | {'N/A':>15} | {cf.count:>5}F ({e})")

    print("-" * 70)
    print(f"Total: {passed + failed} tests, {passed} passed, {failed} failed")
    print(f"Total evaluations: {total_evals}")
    if passed > 0:
        print(f"Average evaluations per pass: {total_evals / passed:.2f}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        if sys.argv[1] == "benchmark":
            run_benchmark()
        elif sys.argv[1] in BRACKETING_SOLVERS:
            run_detailed_test(sys.argv[1])
        else:
            print(f"Unknown argument: {sys.argv[1]}")
            print(f"Usage: python {sys.argv[0]} [benchmark | {' | '.join(BRACKETING_SOLVERS)}]")
    else:
        # Default: run benchmark
        run_benchmark()
