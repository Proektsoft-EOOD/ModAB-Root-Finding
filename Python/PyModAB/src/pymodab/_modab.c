/*
 * Native CPython extension for PyModAB.
 *
 * Implements the modAB algorithm directly as a Python extension module, so that
 * Python callables are invoked through the C API instead of a ctypes/libffi
 * trampoline, and exceptions raised in f(x) propagate to the caller.
 * A raw C function pointer "double f(double)" (e.g. numba cfunc .address) can
 * also be passed as an int, in which case no Python code runs during the solve.
 *
 * Built against the Limited API (abi3), so one binary per platform serves all
 * Python versions >= 3.8.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>

static int evaluation_count = 0;

/* Evaluator: returns f(x). On error sets *err = 1 (Python exception is set). */
typedef double (*eval_fn)(void *ctx, double x, int *err);

static double eval_python(void *ctx, double x, int *err) {
    PyObject *xo = PyFloat_FromDouble(x);
    if (xo == NULL) {
        *err = 1;
        return 0.0;
    }
#if !defined(Py_LIMITED_API) || Py_LIMITED_API + 0 >= 0x030C0000
    PyObject *yo = PyObject_Vectorcall((PyObject *)ctx, &xo, 1, NULL);
#else
    PyObject *yo = PyObject_CallFunctionObjArgs((PyObject *)ctx, xo, NULL);
#endif
    Py_DECREF(xo);
    if (yo == NULL) {
        *err = 1;
        return 0.0;
    }
    double y = PyFloat_AsDouble(yo);
    Py_DECREF(yo);
    if (y == -1.0 && PyErr_Occurred())
        *err = 1;
    return y;
}

static double eval_cfunc(void *ctx, double x, int *err) {
    (void)err;
    return ((double (*)(double))ctx)(x);
}

/* -------------------------------------------------------------------------
 * Shared safeguards
 *
 * These helpers carry the overflow/NaN postconditions that the bracketing
 * solver relies on. They mirror C/src/ModAB.c and the reference implementation
 * in C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs.
 * ------------------------------------------------------------------------- */

/* True only when both values have the same non-zero sign. Comparisons with NaN
   are false, so a caller must reject NaN before using this to update a bracket. */
static inline int same_sign(double x, double y) {
    return (x < 0.0 && y < 0.0) || (x > 0.0 && y > 0.0);
}

/* Midpoint without first forming x1 + x2, which can overflow for finite
   endpoints of the same sign. For finite ordered endpoints the result lies in
   the closed bracket, so callers need no extra clamp. */
static inline double safe_midpoint(double x1, double x2) {
    return 0.5 * x1 + 0.5 * x2;
}

/* Safeguarded false-position/secant point for an ordered bracket [x1, x2].

   The textbook formula (x1*y2 - x2*y1) / (y2 - y1) may overflow although the
   mathematical intersection is finite. With opposite-sign ordinates, writing
   a = |y1| and b = |y2| gives the equivalent convex combination

       x = (b/(a+b))*x1 + (a/(a+b))*x2,

   whose weights lie in [0, 1] and sum to 1. This helper owns the complete
   postcondition every caller needs: the returned point is finite and lies in
   [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint is
   returned instead. */
static inline double safe_secant(double x1, double y1, double x2, double y2) {
    double a = fabs(y1);
    double b = fabs(y2);
    double den = a + b;
    double x;

    /* One test on the denominator covers every unusable case: a NaN ordinate
       propagates into it, two zero ordinates make it zero, and an infinite
       ordinate or an overflowing sum makes it infinite. A zero magnitude does
       NOT indicate a root here, because the ordinates may be Anderson-Bjorck
       auxiliary values, so bisection is the safe and neutral fallback. */
    if (!(den > 0.0))
        return safe_midpoint(x1, x2);

    if (isinf(den)) {
        /* An infinite ordinate carries no usable slope. Otherwise a + b merely
           overflowed, and halving both restores it without changing the ratio
           that defines the weights. */
        if (isinf(a) || isinf(b))
            return safe_midpoint(x1, x2);
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    x = (b / den) * x1 + (a / den) * x2;

    /* In exact arithmetic the convex combination is strictly inside the
       bracket; the projection only corrects a possible last-ulp excursion. */
    return x < x1 ? x1 : (x > x2 ? x2 : x);
}

/* The Anderson-Bjorck contraction factor for the ordinate that did not move. */
static inline double ab_factor(double y3, double y_moved) {
    double m = 1.0 - y3 / y_moved;
    return m > 0.0 ? m : 0.5;
}

/* Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
   (x1, f1), (x2, f2) and a third point (x3, f3) outside it, or NaN when it is unusable.
   The model fits both nearly-quadratic shapes and functions with a pole; for a straight
   line it gives the secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative
   to f2, q1 = f1/f2 and q3 = f3/f2 (so the step does not depend on the scale of f), and
   c = 1 - q1, the root is
       r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
   two divisions per step, 1/f2 and this one. The step is used only if it moves the secant
   point xs = x1 - w*q1/c by at most half the distance to the nearer endpoint (on a plateau
   the model degenerates towards a step function). Multiplied by c, that test reads
       |t*c + w*q1| <= w*min(|q1|, 1)/2,
   with no division; it also rejects NaN and infinite values.
   When the step is unusable, including a NaN third point (no endpoint dropped on that
   side yet), it returns the A&B secant point through (x1, y1), (x2, y2) instead. */
static inline double rational_step(double x1, double f1, double x2, double f2,
                                   double x3, double f3, double y1, double y2) {
    double k = 1.0 / f2;
    double q1 = f1 * k, q3 = f3 * k;
    double c = 1.0 - q1;
    double w = x2 - x1, u = x3 - x1;
    double t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c);
    double a1 = fabs(q1);
    double x = x1 + t;
    if (fabs(t * c + w * q1) <= 0.5 * w * (a1 < 1.0 ? a1 : 1.0) && x <= x2)
        return x;
    return safe_secant(x1, y1, x2, y2);
}

/* True if the second divided differences on the left (xl, x1, x2) and on the right
   (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
   one-signed curvature does worse than the secant. Scaled by positive factors, they are
       pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
       pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
   so only their signs are compared, with no division or product of the two.
   False when xl or xr is NaN (no endpoint dropped on that side yet). */
static inline int curvature_changes_sign(double xl, double fl, double x1, double f1,
                                         double x2, double f2, double xr, double fr) {
    double df = f2 - f1, w = x2 - x1;
    double pl = df * (x1 - xl) - (f1 - fl) * w;
    double pr = (fr - f2) * w - df * (xr - x2);
    return (pl < 0.0 && pr > 0.0) || (pl > 0.0 && pr < 0.0);
}

/* Same algorithm as C/src/ModAB.c, with an abort path for evaluator errors. */
/* The A&B phase adds two changes:
   1. Rational correction: the A&B secant is replaced by the secant through the true
      residuals, corrected for curvature by a 3-point rational model through the bracket
      and the endpoint it most recently dropped. The A&B secant remains the fallback when
      the correction is unusable or the curvature changes sign within the bracket.
   2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
      just evaluated, so that one-sided convergence still closes the bracket. */
/* Every call to f is counted here, so the count is exact on every exit path. */
#define EVAL(x) (++evaluation_count, f(ctx, (x), err)); if (*err) return NAN

static double modab_core(eval_fn f, void *ctx, double x1, double x2,
                         double aTol, double rTol, int maxIter, int *err) {
    evaluation_count = 0;

    if (x1 > x2) {
        double temp = x1; x1 = x2; x2 = temp;
    }

    double y1 = EVAL(x1);
    if (y1 == 0.0)
        return x1;

    double y2 = EVAL(x2);
    if (y2 == 0.0)
        return x2;

    /* NaN has no usable sign, and same_sign is false for it, so it
       must be rejected before the predicate is used to update a bracket. */
    if (isnan(y1) || isnan(y2) || same_sign(y1, y2))
        return NAN;

    int bisection = 1;
    int side = 0;
    double threshold = x2 - x1;
    double f1 = y1, f2 = y2; /* True residuals, kept unmodified by A&B corrections */
    const double C = 2.0;
    /* A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
       than 1 % lies on a plateau; the first A&B step is then not rescaled. */
    const double PLATEAU_RATIO = 0.99;
    /* Each return from A&B to bisection halves the margin of the switching test, so that
       borderline switches are rejected after unsuccessful A&B phases. */
    double margin = 1.0;
    int last = 0; /* Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none */
    double xl = NAN, fl = NAN, xr = NAN, fr = NAN; /* Last endpoints dropped on each side */
    for (int i = 1; i <= maxIter; ++i) {
        double dx = x2 - x1; /* Bracket width */
        double x3, epsx;
        if (bisection) {
            x3 = safe_midpoint(x1, x2);
            epsx = aTol + rTol * fabs(x3);
        } else {
            /* Rational correction through the endpoint dropped by the last update, or the A&B
               secant when the curvature changes sign within the bracket. A missing endpoint is
               NaN, which fails both tests, so it needs no check of its own. */
            if (curvature_changes_sign(xl, fl, x1, f1, x2, f2, xr, fr))
                x3 = safe_secant(x1, y1, x2, y2);
            else
                x3 = last == 1 ? rational_step(x1, f1, x2, f2, xl, fl, y1, y2)
                               : rational_step(x1, f1, x2, f2, xr, fr, y1, y2);
            epsx = aTol + rTol * fabs(x3);
            /* Closing step: when A&B converges from one side, its next point lands within epsx/2
               of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
               towards the far side then crosses the root and closes the bracket. */
            double h = 0.5 * epsx;
            if (last == 1) {
                if (x3 - x1 < h && x1 + h < x2) x3 = x1 + h;
            } else if (x2 - x3 < h && x2 - h > x1) {
                x3 = x2 - h;
            }
        }
        /* Check for x-convergence */
        if (dx <= epsx) {
            return x3;
        }

        double y3;
        int switching = 0; /* This bisection switches to A&B */
        if (bisection) {
            y3 = EVAL(x3);
            double dy = f2 - f1;
            if (isfinite(dy)) { /* Avoids overflow in the calculations below */
                double ym = (f1 + f2) * 0.5; /* Chord ordinate at midpoint; f1, f2 have opposite signs */
                double r = 1.0 - fabs(ym / dy); /* Symmetry factor */
                double k = r * r; /* Deviation factor */
                /* k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test. */
                if (fabs(ym - y3) < (k * fabs(ym) + k * fabs(y3)) * margin) {
                    bisection = 0;
                    switching = 1;
                    threshold = C * dx;
                    y1 = f1; y2 = f2; /* A&B starts from the true residuals */
                }
            }
        } else {
            /* If x3 got clamped, reuse the true residual stored at the endpoint. */
            if (x3 == x1) {
                y3 = f1;
            } else if (x3 == x2) {
                y3 = f2;
            } else {
                y3 = EVAL(x3);
            }
            threshold *= 0.5;
        }

        if (!(fabs(y3) > 0.0)) { /* Exit on zero or NaN with a single branch on the hot path */
            if (y3 == 0.0)
                return x3;
            /* A NaN residual has no usable sign, so the bracket cannot be updated. */
            return NAN;
        }

        if (bisection) { /* Bisection step: only the true residuals are tracked */
            if (same_sign(f1, y3)) {
                xl = x1; fl = f1;
                x1 = x3; f1 = y3; last = 1;
            } else {
                xr = x2; fr = f2;
                x2 = x3; f2 = y3; last = -1;
            }
        } else { /* Anderson-Bjorck step */
            double yl = f1, yr = f2; /* True residuals of the bracket before the update */
            if (same_sign(f1, y3)) {
                if (side == 1)
                    y2 *= ab_factor(y3, y1); /* Anderson-Bjorck correction */
                else if (!switching || fabs(y3) <= PLATEAU_RATIO * fabs(yl))
                    side = 1; /* A switching bisection on a plateau starts the A&B phase afresh */
                if (x3 != x1) { xl = x1; fl = f1; }
                x1 = x3; f1 = y1 = y3; last = 1;
            } else {
                if (side == -1)
                    y1 *= ab_factor(y3, y2); /* Anderson-Bjorck correction */
                else if (!switching || fabs(y3) <= PLATEAU_RATIO * fabs(yr))
                    side = -1; /* A switching bisection on a plateau starts the A&B phase afresh */
                if (x3 != x2) { xr = x2; fr = f2; }
                x2 = x3; f2 = y2 = y3; last = -1;
            }
            /* Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
               (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1) */
            if (x2 - x1 > threshold && 42.0 * fabs(y3) >= fmin(fabs(yl), fabs(yr))) {
                bisection = 1;
                side = 0;
                margin *= 0.5;
            }
        }
    }
    return NAN;
}

#undef EVAL

static PyObject *py_find_root(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"f", "x1", "x2", "atol", "rtol", "max_iter", NULL};
    PyObject *f;
    double x1, x2, atol = 1e-14, rtol = 1e-14;
    int max_iter = 200;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odd|ddi:find_root", kwlist,
                                     &f, &x1, &x2, &atol, &rtol, &max_iter))
        return NULL;

    int err = 0;
    double root;
    if (PyLong_Check(f)) {
        /* Raw address of a C function "double f(double)" */
        void *p = PyLong_AsVoidPtr(f);
        if (p == NULL) {
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_ValueError, "function address must not be NULL");
            return NULL;
        }
        Py_BEGIN_ALLOW_THREADS
        root = modab_core(eval_cfunc, p, x1, x2, atol, rtol, max_iter, &err);
        Py_END_ALLOW_THREADS
    } else if (PyCallable_Check(f)) {
        root = modab_core(eval_python, f, x1, x2, atol, rtol, max_iter, &err);
        if (err)
            return NULL;
    } else {
        PyErr_SetString(PyExc_TypeError,
                        "f must be a callable or the int address of a C function double(double)");
        return NULL;
    }
    return PyFloat_FromDouble(root);
}

static PyObject *py_get_evaluation_count(PyObject *self, PyObject *unused) {
    (void)self; (void)unused;
    return PyLong_FromLong(evaluation_count);
}

static PyMethodDef modab_methods[] = {
    {"find_root", (PyCFunction)(void (*)(void))py_find_root, METH_VARARGS | METH_KEYWORDS,
     "find_root(f, x1, x2, atol=1e-14, rtol=1e-14, max_iter=200)\n--\n\n"
     "Find the root of f(x) = 0 within the interval [x1, x2], using an improved\n"
     "version of the Modified Anderson-Bjork method.\n\n"
     "Parameters\n"
     "----------\n"
     "f : callable or int\n"
     "    A continuous function of one variable, or the address (int) of a\n"
     "    compiled C function double f(double), e.g. numba.cfunc(...).address.\n"
     "    A compiled function is called directly, without any Python overhead.\n"
     "x1, x2 : float\n"
     "    Bracket interval endpoints; f(x1) and f(x2) must have opposite signs.\n"
     "atol, rtol : float, optional\n"
     "    Absolute and relative tolerances (default: 1e-14).\n"
     "max_iter : int, optional\n"
     "    Maximum number of iterations (default: 200).\n\n"
     "Returns\n"
     "-------\n"
     "float\n"
     "    The root, or NaN if no root is found or f(x1) and f(x2) have the same sign.\n"
     "    Exceptions raised by f are propagated to the caller."},
    {"get_evaluation_count", py_get_evaluation_count, METH_NOARGS,
     "get_evaluation_count()\n--\n\n"
     "Number of function evaluations from the last find_root call."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modab_module = {
    PyModuleDef_HEAD_INIT, "_modab", "Native modAB root finder.", -1, modab_methods
};

PyMODINIT_FUNC PyInit__modab(void) {
    return PyModule_Create(&modab_module);
}
