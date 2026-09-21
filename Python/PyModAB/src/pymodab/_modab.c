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

static inline int same_sign(double a, double b) {
    return (a > 0) == (b > 0);
}

static inline double clamp(double d, double min, double max) {
    return d < min ? min : (d > max ? max : d);
}

/* Same algorithm as C/src/ModAB.c, with an abort path for evaluator errors. */
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

    if (same_sign(y1, y2))
        return NAN;

    int bisection = 1;
    int side = 0;
    double threshold = x2 - x1;
    double f1 = y1, f2 = y2, ymin = 0.0;
    const double C = 2.0;
    for (int i = 1; i <= maxIter; ++i) {
        double x3 = bisection ? 0.5 * (x1 + x2) : (x1 * y2 - y1 * x2) / (y2 - y1);
        double eps = aTol + rTol * fabs(x3);
        if (x2 - x1 <= eps) {
            return bisection ? x3 : clamp(x3, x1, x2);
        }

        double y3;
        if (bisection) {
            y3 = EVAL(x3);
            double ym = 0.5 * (f1 + f2);
            double r = 1.0 - fabs(ym / (f2 - f1));
            double k = r * r;
            if (fabs(ym - y3) < k * (fabs(y3) + fabs(ym))) {
                bisection = 0;
                threshold = (x2 - x1) * C;
            }
        } else {
            if (x3 <= x1) {
                x3 = x1; y3 = f1;
            } else if (x3 >= x2) {
                x3 = x2; y3 = f2;
            } else {
                y3 = EVAL(x3);
            }
            threshold *= 0.5;
            // Best true residual of the bracket BEFORE y3 replaces an endpoint.
            ymin = fmin(fabs(f1), fabs(f2));
        }

        if (y3 == 0.0)
            return x3;

        if (same_sign(y1, y3)) {
            if (side == 1) {
                double m = 1.0 - y3 / y1;
                y2 = (m > 0.0) ? y2 * m : y2 * 0.5;
            } else if (!bisection) {
                side = 1;
            }
            x1 = x3; f1 = y1 = y3;
        } else {
            if (side == -1) {
                double m = 1.0 - y3 / y2;
                y1 = (m > 0.0) ? y1 * m : y1 * 0.5;
            } else if (!bisection) {
                side = -1;
            }
            x2 = x3; f2 = y2 = y3;
        }

        // Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if (!bisection && x2 - x1 > threshold && fabs(y3) > 0.5 * ymin) {
            bisection = 1;
            side = 0;
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
