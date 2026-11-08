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
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs, so a change to the
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


# Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
# (x1, f1), (x2, f2) and a third point (x3, f3) outside it, or NaN when it is unusable.
# The model fits both nearly-quadratic shapes and functions with a pole; for a straight
# line it gives the secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative
# to f2, q1 = f1/f2 and q3 = f3/f2 (so the step does not depend on the scale of f), and
# c = 1 - q1, the root is
#     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
# two divisions per step, 1/f2 and this one. The step is used only if it moves the secant
# point xs = x1 - w*q1/c by at most half the distance to the nearer endpoint (on a plateau
# the model degenerates towards a step function). Multiplied by c, that test reads
#     |t*c + w*q1| <= w*min(|q1|, 1)/2,
# with no division; it also rejects NaN and infinite values.
# When the step is unusable, including a NaN third point (no endpoint dropped on that
# side yet), it returns the A&B secant point through (x1, y1), (x2, y2) instead.
cdef inline double rational_step(double x1, double f1, double x2, double f2,
                                 double x3, double f3, double y1, double y2) noexcept nogil:
    cdef double k = 1.0 / f2
    cdef double q1 = f1 * k
    cdef double q3 = f3 * k
    cdef double c = 1.0 - q1
    cdef double w = x2 - x1
    cdef double u = x3 - x1
    cdef double t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c)
    cdef double a1 = fabs(q1)
    cdef double x = x1 + t
    if fabs(t * c + w * q1) <= 0.5 * w * (a1 if a1 < 1.0 else 1.0) and x <= x2:
        return x
    return safe_secant(x1, y1, x2, y2)


# True if the second divided differences on the left (xl, x1, x2) and on the right
# (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
# one-signed curvature does worse than the secant. Scaled by positive factors, they are
#     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
#     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
# so only their signs are compared, with no division or product of the two.
# False when xl or xr is NaN (no endpoint dropped on that side yet).
cdef inline bint curvature_changes_sign(double xl, double fl, double x1, double f1,
                                        double x2, double f2, double xr, double fr) noexcept nogil:
    cdef double df = f2 - f1
    cdef double w = x2 - x1
    cdef double pl = df * (x1 - xl) - (f1 - fl) * w
    cdef double pr = (fr - f2) * w - df * (xr - x2)
    return (pl < 0.0 and pr > 0.0) or (pl > 0.0 and pr < 0.0)


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
    The A&B phase adds two changes:
        1. Rational correction: the A&B secant is replaced by the secant through the true
           residuals, corrected for curvature by a 3-point rational model through the bracket
           and the endpoint it most recently dropped. The A&B secant remains the fallback when
           the correction is unusable or the curvature changes sign within the bracket.
        2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
           just evaluated, so that one-sided convergence still closes the bracket.
    F(x) must be continuous and sign(F(x1)) != sign(F(x2))
    """
    cdef double epsy, y1, y2, f1, f2, x3, epsx, y3, ay3, ym, r, k, threshold, dx, dy, yl, yr, h, xp, fp
    cdef double C = 2.0  # Threshold safety factor
    cdef int side, bisection, switching, _
    cdef int last = 0  # Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none
    # Last endpoints dropped on the left and on the right
    cdef double xl = NAN, fl = NAN, xr = NAN, fr = NAN
    # A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    # than 1 % lies on a plateau; the first A&B step is then not rescaled.
    cdef double PLATEAU_RATIO = 0.99
    # Each return from A&B to bisection halves the margin of the switching test, so that
    # borderline switches are rejected after unsuccessful A&B phases.
    cdef double margin = 1.0
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
    for _ in range(maxiter):
        dx = x2 - x1  # Bracket width
        if bisection:
            x3 = safe_midpoint(x1, x2)
            epsx = xtol * c_max(fabs(x3), 1.0)
        else:
            # Rational correction through the endpoint dropped by the last update, or the A&B
            # secant when the curvature changes sign within the bracket. A missing endpoint is
            # NaN, which fails both tests, so it needs no check of its own.
            if curvature_changes_sign(xl, fl, x1, f1, x2, f2, xr, fr):
                x3 = safe_secant(x1, y1, x2, y2)
            elif last == 1:
                x3 = rational_step(x1, f1, x2, f2, xl, fl, y1, y2)
            else:
                x3 = rational_step(x1, f1, x2, f2, xr, fr, y1, y2)
            epsx = xtol * c_max(fabs(x3), 1.0)
            # Closing step: when A&B converges from one side, its next point lands within epsx/2
            # of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
            # towards the far side then crosses the root and closes the bracket.
            h = 0.5 * epsx
            if last == 1:
                if x3 - x1 < h and x1 + h < x2:
                    x3 = x1 + h
            elif x2 - x3 < h and x2 - h > x1:
                x3 = x2 - h

        if dx <= epsx:  # x-convergence check
            return x3

        switching = 0  # This bisection switches to A&B
        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            dy = f2 - f1
            if isfinite(dy):  # Avoids overflow in the calculations below
                ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint; f1, f2 have opposite signs
                r = 1.0 - fabs(ym / dy)  # Symmetry factor
                k = r * r  # Deviation factor
                # k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if fabs(ym - y3) < (k * fabs(ym) + k * fabs(y3)) * margin:
                    bisection = 0
                    switching = 1
                    threshold = C * dx  # Safety factor: skips two AB steps before the first fallback
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
                xl, fl = x1, f1
                x1, f1, last = x3, y3, 1
            else:
                xr, fr = x2, f2
                x2, f2, last = x3, y3, -1
        else:
            yl, yr = f1, f2  # True residuals of the bracket before the update
            if same_sign(f1, y3):  # Same sign check
                if side == 1:
                    y2 *= ab_factor(y3, y1)
                elif not switching or ay3 <= PLATEAU_RATIO * fabs(yl):
                    side = 1  # A switching bisection on a plateau starts the A&B phase afresh
                if x3 != x1:
                    xl, fl = x1, f1
                x1, y1, f1, last = x3, y3, y3, 1  # Also store the unmodified y1 value to be used for bisection fallback
            else:
                if side == -1:
                    y1 *= ab_factor(y3, y2)
                elif not switching or ay3 <= PLATEAU_RATIO * fabs(yr):
                    side = -1  # A switching bisection on a plateau starts the A&B phase afresh
                if x3 != x2:
                    xr, fr = x2, f2
                x2, y2, f2, last = x3, y3, y3, -1  # Also store the unmodified y2 value to be used for bisection fallback

            # Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            # (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if x2 - x1 > threshold and 42.0 * ay3 >= c_min(fabs(yl), fabs(yr)):
                bisection = 1
                side = 0
                margin *= 0.5
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
