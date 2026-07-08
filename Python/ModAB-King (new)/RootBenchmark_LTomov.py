import math
import time
from ModAB import modAB_root as _modab_root
from ModAB_King import modAB_root as _king_root
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
    """Wrapper for ModAB.modAB_root that matches the benchmark signature."""
    try:
        return _modab_root(f, a, b, 0, xtol=precision, ytol=0.0)
    except (ValueError, RuntimeError):
        return float('nan')


def mod_ab_king(f, a, b, precision=1e-14):
    """Wrapper for ModAB_King.modAB_root that matches the benchmark signature."""
    try:
        return _king_root(f, a, b, 0, xtol=precision, ytol=0.0)
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


# ---- Test problems: 18 functions from chandrupatla_vs_ab_test_table.docx ----
# d09 uses [-1,0] (docx [-1,1] has no sign change: two roots of |x|-0.2);
# d18 uses [0.01,0.2] (docx [0,0.2] has f(0)=0). Both per the docx note to
# shift intervals lacking a sign change. All others verbatim from the table.

def _f06(x):
    s = 1.0 if x > 0 else (-1.0 if x < 0 else 0.0)
    return s * abs(x) ** (1.0 / 3.0)

def _f11(x):
    return (x + 0.3) if x < 0 else (2.0 * x - 0.3)

all_problems = [
    Problem("d01", lambda x: math.cos(x) - x,                    0.0, 1.5),
    Problem("d02", lambda x: math.exp(x) - 3.0 * x,              0.0, 1.5),
    Problem("d03", lambda x: x**3 - x - 1.0,                     0.5, 2.0),
    Problem("d04", lambda x: math.log(x) - 0.5,                  1.0, 3.0),
    Problem("d05", lambda x: math.sqrt(x) - 0.1,                 0.0, 1.5),
    Problem("d06", _f06,                                        -1.0, 1.5),
    Problem("d07", lambda x: (x - 1.0)**3,                       0.0, 3.0),
    Problem("d08", lambda x: (x - 1.0)**5,                       0.0, 3.0),
    Problem("d09", lambda x: abs(x) - 0.2,                      -1.1, 0.1),
    Problem("d10", lambda x: max(x, 0.0) - 0.2,                 -0.6, 1.6),
    Problem("d11", _f11,                                        -0.6, 1.6),
    Problem("d12", lambda x: math.tanh(10.0 * x),               -0.6, 1.6),
    Problem("d13", lambda x: 1.0 - math.exp(-50.0 * x) - 1.0e-6, 0.0, 1.0),
    Problem("d14", lambda x: x**50 - 1.0e-6,                     0.0, 1.0),
    Problem("d15", lambda x: math.exp(x) - 1.0,                -20.0, 1.0),
    Problem("d16", lambda x: math.log(1.0 + x) - 0.5,            0.0, 2.0),
    Problem("d17", lambda x: math.sin(20.0 * x),                 0.1, 0.2),
    Problem("d18", lambda x: math.sin(20.0 * x) + 0.1 * x,      0.01, 0.2),
]

# Solver table — modAB is the relative-factor reference; King added last
solvers = [
    ("bisect", scipy_bisect),
    ("brentq", scipy_brentq),
    ("brenth", scipy_brenth),
    ("ridder", scipy_ridder),
    ("chandr", scipy_chandrupatla),
    (" modAB", mod_ab),
    ("  King", mod_ab_king),
]

MODAB_IDX = len(solvers) - 2  # index of modAB (King is last; modAB stays the ref)


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
    print(header)
    # Stats block for timing
    _print_stats_block(time_data, cnt_w, [n for n, _ in solvers], integer_vals=False)


if __name__ == "__main__":
    run()