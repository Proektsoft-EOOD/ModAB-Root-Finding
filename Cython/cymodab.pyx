# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

from libc.math cimport fabs, NAN, isnan

ctypedef double (*func_type)(double) nogil

cdef inline double c_max(double a, double b) noexcept nogil:
    return a if a > b else b

cdef inline double c_min(double a, double b) noexcept nogil:
    return a if a < b else b

cdef inline double c_clamp(double x, double xmin, double xmax) noexcept nogil:
    return c_max(xmin, c_min(x, xmax))

cpdef double modAB_root(object f, double x1, double x2, double y=0.0,
                         double xtol=1e-14, double ytol=0.0, int maxiter=200):
    """
    Finds the root of "F(x) = 0" within the interval [x1, x2]
    with the specified precisions - absolute: aTol and relative: rTol,
    using an improved version of the modified Anderson Bjork's method:
        Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
        Improvements to the Modified Anderson–Björck(modAB) Root-Finding Algorithm.
        Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
    Additional fixes proposed by L. Tomov are applied in this version:
        1. The secant point is clamped to the interval [p1.X, p2.X] before the X-convergence exit
        2. The original function values y1 and y2 (without A&B corrections)
           are stored for later use in bisection fallback
    F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))
    """
    cdef double epsy, y1, y2, f1, f2, x3, epsx, y3, ym, dy, r, k, m, threshold, ymin
    cdef int side, bisection, _
    cdef bint same_sign

    if x2 < x1:
        x1, x2 = x2, x1

    epsy = ytol * c_max(fabs(y), 1.0)
    y1 = f(x1) - y
    if fabs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if fabs(y2) <= epsy:
        return x2

    if (y1 > 0.0) == (y2 > 0.0):
        return NAN  # No sign change - no root guaranteed

    f1 = y1
    f2 = y2  # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = 1
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    ymin = 0.0
    for _ in range(maxiter):
        if bisection:
            x3 = (x1 + x2) * 0.5
        else:
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)

        epsx = xtol * c_max(fabs(x3), 1.0)
        if x2 - x1 <= epsx:  # x-convergence check
            if bisection:
                return x3
            else:
                return c_clamp(x3, x1, x2)  # Clamp the secant value

        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint
            dy = f2 - f1
            r = 1.0 - fabs(ym / dy)  # Symmetry factor
            k = r * r  # Deviation factor
            if fabs(ym - y3) < k * (fabs(y3) + fabs(ym)):
                bisection = 0
                threshold = (x2 - x1) * 2.0  # Safety factor: skips two AB steps before the first fallback
        else:
            if x3 <= x1:
                x3 = x1
                y3 = f1
            elif x3 >= x2:
                x3 = x2
                y3 = f2
            else:
                y3 = f(x3) - y

            threshold *= 0.5
            # Best true residual of the bracket BEFORE y3 replaces an endpoint.
            ymin = c_min(fabs(f1), fabs(f2))

        if fabs(y3) <= epsy:  # y-convergence check
            return x3

        if (y1 > 0.0) == (y3 > 0.0):  # Same sign check
            if side == 1:
                m = 1.0 - y3 / y1
                if m > 0.0:
                    y2 = y2 * m
                else:
                    y2 = y2 * 0.5
            elif not bisection:
                side = 1
            x1 = x3
            y1 = y3
            f1 = y3  # Also store the unmodified y1 value to be used for bisection fallback
        else:
            if side == -1:
                m = 1.0 - y3 / y2
                if m > 0.0:
                    y1 = y1 * m
                else:
                    y1 = y1 * 0.5
            elif not bisection:
                side = -1
            x2 = x3
            y2 = y3
            f2 = y3  # Also store the unmodified y2 value to be used for bisection fallback

        # Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if not bisection and x2 - x1 > threshold and abs(y3) > 0.5 * ymin:
            bisection = True
            side = 0
    return NAN
