# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

from libc.math cimport fabs, NAN, isnan, isinf, isfinite

ctypedef double (*func_type)(double) nogil

cdef inline double c_max(double a, double b) noexcept nogil:
    return a if a > b else b

cdef inline double c_min(double a, double b) noexcept nogil:
    return a if a < b else b

cdef inline double c_clamp(double x, double xmin, double xmax) noexcept nogil:
    return c_max(xmin, c_min(x, xmax))


# ---------------------------------------------------------------------------
# Shared safeguards
#
# These helpers carry the overflow/NaN postconditions that the bracketing
# solver relies on. They are the Cython form of the reference implementation in
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs, so a change to the
# numerics belongs in one place rather than inline in the solver.
# ---------------------------------------------------------------------------

cdef inline bint same_sign(double x, double y) noexcept nogil:
    return (x < 0.0 and y < 0.0) or (x > 0.0 and y > 0.0)


cdef inline double safe_midpoint(double x1, double x2) noexcept nogil:
    return 0.5 * x1 + 0.5 * x2


cdef inline double safe_secant(double x1, double y1, double x2, double y2) noexcept nogil:
    cdef double a = fabs(y1)
    cdef double b = fabs(y2)
    cdef double den = a + b
    cdef double x
    if not den > 0.0:
        return safe_midpoint(x1, x2)

    if isinf(den):
        if isinf(a) or isinf(b):
            return safe_midpoint(x1, x2)
        a *= 0.5
        b *= 0.5
        den = a + b

    x = (b / den) * x1 + (a / den) * x2
    return c_clamp(x, x1, x2)


cdef inline double ab_factor(double y3, double y_moved) noexcept nogil:
    cdef double m = 1.0 - y3 / y_moved
    return m if m > 0.0 else 0.5


cpdef double modAB_root(object f, double x1, double x2, double y=0.0,
                         double xtol=1e-14, double ytol=0.0, int maxiter=200):
    """
    Finds the root of "F(x) = 0" within the interval [x1, x2]
    with the specified precisions - absolute: aTol and relative: rTol,
    using an improved version of the modified Anderson Bjork's method:
        Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
        Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
        Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
    Additional fixes proposed by L. Tomov are applied in this version:
        1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
        2. The original function values y1 and y2 (without A&B corrections)
           are stored for later use in bisection fallback
    The overflow- and NaN-safe form of the interpolation lives in the shared
    safeguards above; the switching test is written so that it cannot overflow.
    F(x) must be continuous and sign(F(x1)) != sign(F(x2))
    """
    cdef double epsy, y1, y2, f1, f2, x3, epsx, y3, ay3, ym, r, k, threshold, ymin, dx, dy, yl, yr
    cdef double C = 2.0  # Threshold safety factor
    cdef int MAX_RESIDUAL_STEPS = 3  # Max consecutive A&B steps kept by the residual test alone
    cdef int side, bisection, residual_steps, _
    if x2 < x1:
        x1, x2 = x2, x1

    epsy = ytol * c_max(fabs(y), 1.0)
    y1 = f(x1) - y
    if fabs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if fabs(y2) <= epsy:
        return x2

    # NaN has no usable sign, and same_sign is false for it, so it must
    # be rejected before the predicate is used to update a bracket.
    if isnan(y1) or isnan(y2) or same_sign(y1, y2):
        return NAN  # No sign change - no root guaranteed

    f1 = y1
    f2 = y2  # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = 1
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    residual_steps = 0
    for _ in range(maxiter):
        if bisection:
            x3 = safe_midpoint(x1, x2)
        else:
            x3 = safe_secant(x1, y1, x2, y2)

        dx = x2 - x1  # Bracket width
        epsx = xtol * c_max(fabs(x3), 1.0)
        if dx <= epsx:  # x-convergence check
            return x3

        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            dy = f2 - f1
            if isfinite(dy):  # Avoids overflow in the calculations below
                ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint; f1, f2 have opposite signs
                r = 1.0 - fabs(ym / dy)  # Symmetry factor
                k = r * r  # Deviation factor
                # k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if fabs(ym - y3) < k * fabs(ym) + k * fabs(y3):
                    bisection = 0
                    threshold = C * dx  # Safety factor: skips two AB steps before the first fallback
                    residual_steps = 0
                    y1 = f1  # A&B starts from the true residuals
                    y2 = f2
        else:
            # If x3 got clamped, reuse the true residual stored at the endpoint.
            if x3 == x1:
                y3 = f1
            elif x3 == x2:
                y3 = f2
            else:
                y3 = f(x3) - y

            threshold *= 0.5

        ay3 = fabs(y3)
        if not ay3 > epsy:  # Exit on y-convergence or NaN with a single branch on the hot path
            if ay3 <= epsy:
                return x3
            # A NaN residual has no usable sign, so the bracket cannot be updated.
            return NAN

        if bisection:
            if same_sign(f1, y3):  # Same sign check
                x1, f1 = x3, y3
            else:
                x2, f2 = x3, y3
        else:
            yl, yr = f1, f2  # True residuals of the bracket before the update
            if same_sign(f1, y3):  # Same sign check
                if side == 1:
                    y2 *= ab_factor(y3, y1)
                else:
                    side = 1
                x1, y1, f1 = x3, y3, y3  # Also store the unmodified y1 value to be used for bisection fallback
            else:
                if side == -1:
                    y1 *= ab_factor(y3, y2)
                else:
                    side = -1
                x2, y2, f2 = x3, y3, y3  # Also store the unmodified y2 value to be used for bisection fallback

            # Fallback if AB fails to reduce the bracket width, unless it still halves the residual,
            # but for no more than MAX_RESIDUAL_STEPS consecutive steps
            if x2 - x1 > threshold:
                ymin = c_min(fabs(yl), fabs(yr))  # Best true residual of the bracket
                if residual_steps >= MAX_RESIDUAL_STEPS or 2.0 * ay3 >= ymin:
                    bisection = 1
                    side = 0
                else:
                    residual_steps += 1
            else:
                residual_steps = 0
    return NAN


# ---------------------------------------------------------------------------
# Test hooks
#
# The safeguards above are `cdef inline ... nogil` so they stay inlined in the
# hot loop and are therefore invisible to Python. These thin `cpdef` wrappers
# exist solely so the edge-case suite can exercise the overflow/NaN paths
# directly; the solver never calls them.
# ---------------------------------------------------------------------------

cpdef bint t_same_sign(double x, double y):
    return same_sign(x, y)

cpdef double t_safe_midpoint(double x1, double x2):
    return safe_midpoint(x1, x2)

cpdef double t_safe_secant(double x1, double y1, double x2, double y2):
    return safe_secant(x1, y1, x2, y2)
