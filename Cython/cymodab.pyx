# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

from libc.math cimport fabs, copysign, NAN, isnan, isinf, isfinite
from libc.float cimport DBL_MAX

ctypedef double (*func_type)(double) nogil

# Smallest positive denormal; libc.float exposes DBL_TRUE_MIN only on C11
# headers, so it is spelled out here to keep the extension portable.
cdef double DBL_TRUE_MIN = 5e-324

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
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,ModABCorr}.cs, so a change to the
# numerics belongs in one place rather than inline in the solver.
# ---------------------------------------------------------------------------

cdef inline bint same_nonzero_sign(double x, double y) noexcept nogil:
    """True only when both values have the same non-zero sign. Comparisons with
    NaN are false, so a caller must reject NaN before using this predicate to
    update a bracket."""
    return (x < 0.0 and y < 0.0) or (x > 0.0 and y > 0.0)


cdef inline double safe_midpoint(double x1, double x2) noexcept nogil:
    """Midpoint without first forming x1 + x2, which can overflow for finite
    endpoints of the same sign. For finite ordered endpoints the result lies in
    the closed bracket, so callers need no extra clamp."""
    return 0.5 * x1 + 0.5 * x2


cdef inline double safe_secant(double x1, double y1, double x2, double y2) noexcept nogil:
    """Safeguarded false-position/secant point for an ordered bracket [x1, x2].

    The textbook formula (x1*y2 - x2*y1) / (y2 - y1) may overflow although the
    mathematical intersection is finite. With opposite-sign ordinates, writing
    a = |y1| and b = |y2| gives the equivalent convex combination

        x = (b/(a+b))*x1 + (a/(a+b))*x2,

    whose weights lie in [0, 1] and sum to 1. This helper owns the complete
    postcondition every caller needs: the returned point is finite and lies in
    [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint
    is returned instead."""
    cdef double a = fabs(y1)
    cdef double b = fabs(y2)
    cdef double den = a + b
    cdef double x

    # One test on the denominator covers every unusable case: a NaN ordinate
    # propagates into it, two zero ordinates make it zero, and an infinite
    # ordinate or an overflowing sum makes it infinite. A zero magnitude does
    # NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    # auxiliary values, so bisection is the safe and neutral fallback.
    if not den > 0.0:
        return safe_midpoint(x1, x2)

    if isinf(den):
        # An infinite ordinate carries no usable slope. Otherwise a + b merely
        # overflowed, and halving both restores it without changing the ratio
        # that defines the weights.
        if isinf(a) or isinf(b):
            return safe_midpoint(x1, x2)
        a *= 0.5
        b *= 0.5
        den = a + b

    x = (b / den) * x1 + (a / den) * x2

    # In exact arithmetic the convex combination is strictly inside the bracket;
    # the projection only corrects a possible last-ulp excursion.
    return c_clamp(x, x1, x2)


cdef inline double symmetry_factor(double y1, double y2) noexcept nogil:
    """Returns k = r^2 for the symmetry-sensitive switching criterion. The
    calculation is homogeneous in the true endpoint residuals."""
    cdef double a = fabs(y1)
    cdef double b = fabs(y2)
    cdef double den = a + b
    cdef double r

    if isinf(den):
        # Infinite true residuals deliberately disable switching and keep the
        # controller in bisection mode. NaN is returned rather than an infinity
        # because every exit of passes_switching_test is a "<" comparison,
        # which is false against NaN; an infinity would instead satisfy it and
        # switch. Residuals are never zero here, so only an overflowing sum
        # remains, and halving both restores it without changing the ratio.
        if isinf(a) or isinf(b):
            return NAN
        a *= 0.5
        b *= 0.5
        den = a + b

    # |b-a| <= den, so the quotient lies in [0,1]; halving after the division
    # avoids forming 2*den, which could overflow.
    r = 1.0 - fabs(b - a) / den / 2.0
    return r * r


cdef inline bint passes_switching_test(double ym, double yf, double symmetry) noexcept nogil:
    """Tests whether the true midpoint value yf is close enough to the midpoint
    value ym of the chord through the true endpoint residuals."""
    cdef double abs_ym = fabs(ym)
    cdef double abs_yf = fabs(yf)
    cdef double total = abs_yf + abs_ym
    cdef double scale, norm_ym, norm_yf

    # Fast path. The exact-root case is handled before this is called, and a
    # non-finite ordinate or a NaN symmetry factor fails the comparison, which
    # disables switching as intended.
    if isfinite(total):
        return fabs(ym - yf) < symmetry * total

    # Only reached when the sum overflows. Non-finite values are unsuitable for
    # the linearity comparison.
    if not (isfinite(ym) and isfinite(yf)):
        return False

    # Normalize both sides of the homogeneous inequality to avoid overflow.
    scale = c_max(abs_yf, abs_ym)
    norm_ym = ym / scale
    norm_yf = yf / scale
    return fabs(norm_ym - norm_yf) < symmetry * (fabs(norm_yf) + fabs(norm_ym))


cdef inline double ab_factor(double y3, double y_moved) noexcept nogil:
    """The Anderson-Bjorck contraction factor for the ordinate that did not move."""
    cdef double m = 1.0 - y3 / y_moved
    return m if m > 0.0 else 0.5


cdef inline double scale_preserving_nonzero_sign(double value, double positive_factor) noexcept nogil:
    """Multiplies an auxiliary Anderson-Bjorck ordinate by a positive factor
    while preserving a finite non-zero sign in binary64 arithmetic. This keeps
    same_nonzero_sign sound: an auxiliary ordinate that underflowed to zero
    would otherwise silently change which branch of the bracket update is
    taken. It acts only on auxiliary ordinates; an underflowed working value is
    never accepted as a root of f."""
    cdef double scaled = value * positive_factor

    if scaled == 0.0 and value != 0.0:
        return copysign(DBL_TRUE_MIN, value)

    if isinf(scaled):
        return copysign(DBL_MAX, value)

    return scaled


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
    The overflow- and NaN-safe forms of the interpolation and switching
    arithmetic live in the shared safeguards above.
    F(x) must be continuous and sign(F(x1)) != sign(F(x2))
    """
    cdef double epsy, y1, y2, f1, f2, x3, epsx, y3, ym, threshold, ymin
    cdef int side, bisection, _

    if x2 < x1:
        x1, x2 = x2, x1

    epsy = ytol * c_max(fabs(y), 1.0)
    y1 = f(x1) - y
    if fabs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if fabs(y2) <= epsy:
        return x2

    # NaN has no usable sign, and same_nonzero_sign is false for it, so it must
    # be rejected before the predicate is used to update a bracket.
    if isnan(y1) or isnan(y2) or same_nonzero_sign(y1, y2):
        return NAN  # No sign change - no root guaranteed

    f1 = y1
    f2 = y2  # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = 1
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    ymin = 0.0
    for _ in range(maxiter):
        # safe_secant already returns a point inside [x1, x2], so the separate
        # clamp on the convergence exit is no longer needed.
        if bisection:
            x3 = safe_midpoint(x1, x2)
        else:
            x3 = safe_secant(x1, y1, x2, y2)

        epsx = xtol * c_max(fabs(x3), 1.0)
        if x2 - x1 <= epsx:  # x-convergence check
            return x3

        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            ym = safe_midpoint(f1, f2)  # Ordinate of chord at midpoint
            if passes_switching_test(ym, y3, symmetry_factor(f1, f2)):
                bisection = 0
                threshold = (x2 - x1) * 2.0  # Safety factor: skips two AB steps before the first fallback
        else:
            # If rounding makes the proposal coincide with an endpoint, reuse
            # the true residual already stored there.
            if x3 == x1:
                y3 = f1
            elif x3 == x2:
                y3 = f2
            else:
                y3 = f(x3) - y

            threshold *= 0.5
            # Best true residual of the bracket BEFORE y3 replaces an endpoint.
            ymin = c_min(fabs(f1), fabs(f2))

        if fabs(y3) <= epsy:  # y-convergence check
            return x3

        # A NaN residual has no usable sign, so the bracket cannot be updated.
        if isnan(y3):
            return NAN

        if same_nonzero_sign(y1, y3):  # Same sign check
            if side == 1:
                y2 = scale_preserving_nonzero_sign(y2, ab_factor(y3, y1))
            elif not bisection:
                side = 1
            x1 = x3
            y1 = y3
            f1 = y3  # Also store the unmodified y1 value to be used for bisection fallback
        else:
            if side == -1:
                y1 = scale_preserving_nonzero_sign(y1, ab_factor(y3, y2))
            elif not bisection:
                side = -1
            x2 = x3
            y2 = y3
            f2 = y3  # Also store the unmodified y2 value to be used for bisection fallback

        # Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if not bisection and x2 - x1 > threshold and fabs(y3) > 0.5 * ymin:
            bisection = 1
            side = 0
    return NAN


# ---------------------------------------------------------------------------
# Test hooks
#
# The safeguards above are `cdef inline ... nogil` so they stay inlined in the
# hot loop and are therefore invisible to Python. These thin `cpdef` wrappers
# exist solely so the edge-case suite can exercise the overflow/NaN paths
# directly; the solver never calls them.
# ---------------------------------------------------------------------------

cpdef bint t_same_nonzero_sign(double x, double y):
    return same_nonzero_sign(x, y)

cpdef double t_safe_midpoint(double x1, double x2):
    return safe_midpoint(x1, x2)

cpdef double t_safe_secant(double x1, double y1, double x2, double y2):
    return safe_secant(x1, y1, x2, y2)

cpdef double t_symmetry_factor(double y1, double y2):
    return symmetry_factor(y1, y2)

cpdef bint t_passes_switching_test(double ym, double yf, double symmetry):
    return passes_switching_test(ym, yf, symmetry)
