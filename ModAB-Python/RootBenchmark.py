import math
import time
from pymodab import find_root as pymodab_find_root
from dataclasses import dataclass
from typing import Callable
import numpy as np
from scipy.optimize import bisect as sp_bisect, brentq, brenth, ridder as sp_ridder, elementwise

# Function-call counting wrapper
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

def make_scipy_solver(scipy_func, name):
    """Create a wrapper that matches the (f, a, b, precision) signature."""
    def solver(f, a, b, precision=1e-14):
        try:
            root, info = scipy_func(f, a, b, xtol=precision, rtol=precision,
                                    maxiter=100, full_output=True, disp=False)
            return root
        except (ValueError, RuntimeError):
            return float('nan')
    solver.__name__ = name
    return solver

def wrap_find_root():
    """Create a wrapper for find_root that matches the required signature."""
    def solver(f, a, b, precision=1e-14):
        g = np.vectorize(f, otypes=[np.float64])
        try:
            tolerances = dict(xatol=precision, xrtol=precision, fatol=0, frtol=0)
            res = elementwise.find_root(g, (a, b), tolerances=tolerances)
            return res.x
        except (ValueError, RuntimeError):
            return float('nan')
    solver.__name__ = 'sp_chandrupatla'
    return solver

scipy_bisect = make_scipy_solver(sp_bisect, "sp_bisect")
scipy_brentq = make_scipy_solver(brentq,    "sp_brentq")
scipy_brenth = make_scipy_solver(brenth,    "sp_brenth")
scipy_ridder = make_scipy_solver(sp_ridder, "sp_ridder")
scipy_chandrupatla = wrap_find_root()

def mod_ab(f, a, b, precision=1e-14):
    """Wrapper for pymodab.find_root that matches the benchmark signature."""
    try:
        return pymodab_find_root(f, a, b, atol=precision, rtol=precision)
    except (ValueError, RuntimeError):
        return float('nan')

# Problem definition
@dataclass
class Problem:
    name: str
    f: Callable[[float], float]
    a: float
    b: float
    value: float = 0.0


def P(x):
    return x + 1.11111

# Test problems
problems1 = [
    # Sérgio Galdino. A family of regula falsi root-finding methods
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
    # A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method for the Numerical Solution of a Single Nonlinear Equation
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
]

all_problems = problems1 + problems2 + problems3

# Solver table — modAB must be last for the relative-factor column
solvers = [
    ("bisect", scipy_bisect),
    ("brentq", scipy_brentq),
    ("brenth", scipy_brenth),
    ("ridder", scipy_ridder),
    ("chandr", scipy_chandrupatla),
    (" modAB", mod_ab),
]

MODAB_IDX = len(solvers) - 1  # index of modAB in the solvers list


def _fmt_stat(val, width=8):
    """Format a float stat value for printing."""
    if math.isnan(val):
        return f"{'N/A':>{width}}"
    return f"{val:>{width}.2f}"


def _print_stats_block(data, col_w, solver_names, *, integer_vals=False):
    """
    Print avg / median / min / max / factor-vs-modAB summary rows for one metric.

    data  : list of lists  data[problem_idx][solver_idx]  (NaN for errors)
    """
    n_solvers = len(solver_names)
    stats = {}
    for stat in ("avg", "median", "min", "max"):
        row = []
        for j in range(n_solvers):
            vals = [data[i][j] for i in range(len(data)) if not math.isnan(data[i][j])]
            if not vals:
                row.append(float('nan'))
            elif stat == "avg":
                row.append(sum(vals) / len(vals))
            elif stat == "median":
                s = sorted(vals)
                mid = len(s) // 2
                row.append((s[mid - 1] + s[mid]) / 2 if len(s) % 2 == 0 else s[mid])
            elif stat == "min":
                row.append(min(vals))
            elif stat == "max":
                row.append(max(vals))
        stats[stat] = row

    fmt = (lambda v: f"{int(round(v)):>{col_w}}") if integer_vals else (lambda v: f"{v:>{col_w}.4f}")

    for stat_name, row in stats.items():
        line = f"{stat_name.upper():>6}| "
        for v in row:
            line += (f"{'N/A':>{col_w}}| " if math.isnan(v) else fmt(v) + "| ")
        print(line)

    # Relative factor vs modAB (modAB = 1.00x baseline)
    modab_avg = stats["avg"][MODAB_IDX]
    line = f"{'FACTOR':>6}| "
    for j in range(n_solvers):
        v = stats["avg"][j]
        if math.isnan(v) or math.isnan(modab_avg) or modab_avg == 0:
            line += f"{'N/A':>{col_w}}| "
        else:
            factor = v / modab_avg
            line += f"{factor:>{col_w - 1}.3f}x| "
    print(line)
    print()

def run():
    eps = 1e-14
    col_w = 22  # column width for results / function-value tables
    cnt_w = 8   # column width for count / time tables

    print("Results (root values)")
    print("===")
    header = f"{'Func':>4}| " + "| ".join(f"{name:>{col_w}}" for name, _ in solvers)
    print(header)
    print("------ | ------: | ------: | ------: | ------: | ------: | ------: |")
    for p in all_problems:
        line = f"{p.name:>4}| "
        for name, solver in solvers:
            cf = CountedFunc(p.f)
            try:
                result = solver(cf, p.a, p.b, eps)
                line += f"{result:>{col_w}.15g}| "
            except Exception:
                line += f"{'ERR':>{col_w}}| "
        print(line)
    print()

    print("Function values  f(root)")
    print("===")
    header = f"{'Func':>4}| " + "| ".join(f"{name:>{col_w}}" for name, _ in solvers)
    print(header)
    print("------ | ------: | ------: | ------: | ------: | ------: | ------: |")
    fval_data = []
    for p in all_problems:
        line = f"{p.name:>4}| "
        row = []
        for name, solver in solvers:
            cf = CountedFunc(p.f)
            try:
                result = solver(cf, p.a, p.b, eps)
                fv = abs(cf(result))
                row.append(fv)
                line += f"{cf(result):>{col_w}.15g}| "
            except Exception:
                row.append(float('nan'))
                line += f"{'ERR':>{col_w}}| "
        fval_data.append(row)
        print(line)
    print()

    print("Function evaluations")
    print("===")
    header = f"{'Func':>4}| " + "| ".join(f"{name:>{cnt_w}}" for name, _ in solvers)
    print(header)
    print("------ | ------: | ------: | ------: | ------: | ------: | ------: |")
    count_data = []
    totals = [0] * len(solvers)
    for p in all_problems:
        line = f"{p.name:>4}| "
        row = []
        for j, (name, solver) in enumerate(solvers):
            cf = CountedFunc(p.f)
            try:
                solver(cf, p.a, p.b, eps)
                totals[j] += cf.count
                row.append(float(cf.count))
                line += f"{cf.count:>{cnt_w}}| "
            except Exception:
                row.append(float('nan'))
                line += f"{'ERR':>{cnt_w}}| "
        count_data.append(row)
        print(line)

    # Totals row
    line = f"{'SUM':>4}| "
    for t in totals:
        line += f"{t:>{cnt_w}}| "
    print(line)

    # Stats block for eval counts
    _print_stats_block(count_data, cnt_w, [n for n, _ in solvers], integer_vals=True)

    print("Execution times  (ms per problem, 100 iterations)")
    print("===")
    iterations = 200
    header = f"{'Func':>4}| " + "| ".join(f"{name:>{cnt_w}}" for name, _ in solvers)
    print(header)
    print("-- | -- | -- | -- | -- | -- ")
    time_data = []
    total_time = [0.0] * len(solvers)
    for p in all_problems:
        line = f"{p.name:>4}| "
        row = []
        for j, (name, solver) in enumerate(solvers):
            start = time.perf_counter()
            for _ in range(iterations):
                try:
                    solver(p.f, p.a, p.b, eps)
                except Exception:
                    pass
            elapsed = (time.perf_counter() - start) * 1000  # ms
            total_time[j] += elapsed
            row.append(elapsed)
            line += f"{elapsed:>{cnt_w}.2f}| "
        time_data.append(row)
        print(line)

    # Totals row
    line = f"{'SUM':>4}| "
    for t in total_time:
        line += f"{t:>{cnt_w}.2f}| "
    print(line)

    # Stats block for timing
    _print_stats_block(time_data, cnt_w, [n for n, _ in solvers], integer_vals=False)


if __name__ == "__main__":
    run()