import math

_NAN = float('nan')
_MAX_FLOAT = 1.7976931348623157e308   # sys.float_info.max
_MIN_SUBNORMAL = 5e-324               # smallest positive denormal


# ---------------------------------------------------------------------------
# Shared safeguards
#
# These helpers carry the overflow/NaN postconditions that the bracketing
# solver relies on. They are the Python form of the reference implementation in
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,ModABCorr}.cs, so a change to the
# numerics belongs in one place rather than inline in the solver.
# ---------------------------------------------------------------------------

def same_nonzero_sign(x, y):
    """True only when both values have the same non-zero sign.

    Comparisons with NaN are false, so a caller must reject NaN before using
    this predicate to update a bracket.
    """
    return (x < 0.0 and y < 0.0) or (x > 0.0 and y > 0.0)


def safe_midpoint(x1, x2):
    """Midpoint without first forming x1 + x2, which can overflow for finite
    endpoints of the same sign.

    For finite ordered endpoints the result lies in the closed bracket, so
    callers need no extra clamp.
    """
    return 0.5 * x1 + 0.5 * x2


def safe_secant(x1, y1, x2, y2):
    """Safeguarded false-position/secant point for an ordered bracket [x1, x2].

    The textbook formula ``(x1*y2 - x2*y1) / (y2 - y1)`` may overflow although
    the mathematical intersection is finite. With opposite-sign ordinates,
    writing ``a = |y1|`` and ``b = |y2|`` gives the equivalent convex
    combination::

        x = (b/(a+b))*x1 + (a/(a+b))*x2

    whose weights lie in [0, 1] and sum to 1. This helper owns the complete
    postcondition every caller needs: the returned point is finite and lies in
    [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint
    is returned instead.
    """
    a = abs(y1)
    b = abs(y2)
    den = a + b

    # One test on the denominator covers every unusable case: a NaN ordinate
    # propagates into it, two zero ordinates make it zero, and an infinite
    # ordinate or an overflowing sum makes it infinite. A zero magnitude does
    # NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    # auxiliary values, so bisection is the safe and neutral fallback.
    if not den > 0.0:
        return safe_midpoint(x1, x2)

    if math.isinf(den):
        # An infinite ordinate carries no usable slope. Otherwise a + b merely
        # overflowed, and halving both restores it without changing the ratio
        # that defines the weights.
        if math.isinf(a) or math.isinf(b):
            return safe_midpoint(x1, x2)
        a *= 0.5
        b *= 0.5
        den = a + b

    x = (b / den) * x1 + (a / den) * x2

    # In exact arithmetic the convex combination is strictly inside the
    # bracket; the projection only corrects a possible last-ulp excursion.
    return x1 if x < x1 else (x2 if x > x2 else x)


def symmetry_factor(y1, y2):
    """Return k = r^2 for the symmetry-sensitive switching criterion.

    The calculation is homogeneous in the true endpoint residuals.
    """
    a = abs(y1)
    b = abs(y2)
    den = a + b

    if math.isinf(den):
        # Infinite true residuals deliberately disable switching and keep the
        # controller in bisection mode. NaN is returned rather than an infinity
        # because every exit of passes_switching_test is a "<" comparison,
        # which is false against NaN; an infinity would instead satisfy it and
        # switch. Residuals are never zero here, so only an overflowing sum
        # remains, and halving both restores it without changing the ratio.
        if math.isinf(a) or math.isinf(b):
            return _NAN
        a *= 0.5
        b *= 0.5
        den = a + b

    # |b-a| <= den, so the quotient lies in [0,1]; halving after the division
    # avoids forming 2*den, which could overflow.
    r = 1.0 - abs(b - a) / den / 2.0
    return r * r


def passes_switching_test(ym, yf, symmetry):
    """True when the true midpoint value ``yf`` is close enough to the midpoint
    value ``ym`` of the chord through the true endpoint residuals.
    """
    abs_ym = abs(ym)
    abs_yf = abs(yf)
    total = abs_yf + abs_ym

    # Fast path. The exact-root case is handled before this is called, and a
    # non-finite ordinate or a NaN symmetry factor fails the comparison, which
    # disables switching as intended.
    if math.isfinite(total):
        return abs(ym - yf) < symmetry * total

    # Only reached when the sum overflows. Non-finite values are unsuitable for
    # the linearity comparison.
    if not (math.isfinite(ym) and math.isfinite(yf)):
        return False

    # Normalize both sides of the homogeneous inequality to avoid overflow.
    scale = max(abs_yf, abs_ym)
    norm_ym = ym / scale
    norm_yf = yf / scale
    return abs(norm_ym - norm_yf) < symmetry * (abs(norm_yf) + abs(norm_ym))


def ab_factor(y3, y_moved):
    """The Anderson-Bjorck contraction factor for the ordinate that did not move."""
    m = 1.0 - y3 / y_moved
    return m if m > 0.0 else 0.5


def scale_preserving_nonzero_sign(value, positive_factor):
    """Multiply an auxiliary Anderson-Bjorck ordinate by a positive factor while
    preserving a finite non-zero sign in binary64 arithmetic.

    This keeps :func:`same_nonzero_sign` sound: an auxiliary ordinate that
    underflowed to zero would otherwise silently change which branch of the
    bracket update is taken. It acts only on auxiliary ordinates; an
    underflowed working value is never accepted as a root of f.
    """
    scaled = value * positive_factor

    if scaled == 0.0 and value != 0.0:
        return math.copysign(_MIN_SUBNORMAL, value)

    if math.isinf(scaled):
        return math.copysign(_MAX_FLOAT, value)

    return scaled


def modAB_root(f, x1, x2, y, xtol=1e-14, ytol=0.0, maxiter=200):
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
    if x2 < x1:
        x1, x2 = x2, x1

    epsy = ytol * max(abs(y), 1)
    y1 = f(x1) - y
    if abs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if abs(y2) <= epsy:
        return x2

    # NaN has no usable sign, and same_nonzero_sign is false for it, so it must
    # be rejected before the predicate is used to update a bracket.
    if math.isnan(y1) or math.isnan(y2) or same_nonzero_sign(y1, y2):
        return _NAN  # No sign change - no root guaranteed

    f1, f2 = y1, y2  # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = True
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    C = 2  # Threshold safety factor
    # Best residual of the bracket, refreshed before each AB step. An AB step that
    # fails the width test is still kept if it at least halved this residual.
    ymin = 0
    for _ in range(maxiter):
        # safe_secant already returns a point inside [x1, x2], so the separate
        # clamp on the convergence exit is no longer needed.
        x3 = safe_midpoint(x1, x2) if bisection else safe_secant(x1, y1, x2, y2)
        epsx = xtol * max(abs(x3), 1)
        if x2 - x1 <= epsx:  # x-convergence check
            return x3

        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            ym = safe_midpoint(f1, f2)  # Ordinate of chord at midpoint
            if passes_switching_test(ym, y3, symmetry_factor(f1, f2)):
                bisection = False
                threshold = (x2 - x1) * C  # Safety factor: skips two AB steps before the first fallback
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
            ymin = min(abs(f1), abs(f2))

        if abs(y3) <= epsy:  # y-convergence check
            return x3

        # A NaN residual has no usable sign, so the bracket cannot be updated.
        if math.isnan(y3):
            return _NAN

        if same_nonzero_sign(y1, y3):  # Same sign check
            if side == 1:
                y2 = scale_preserving_nonzero_sign(y2, ab_factor(y3, y1))
            elif not bisection:
                side = 1
            x1, y1, f1 = x3, y3, y3  # Also store the unmodified y1 value to be used for bisection fallback
        else:
            if side == -1:
                y1 = scale_preserving_nonzero_sign(y1, ab_factor(y3, y2))
            elif not bisection:
                side = -1
            x2, y2, f2 = x3, y3, y3  # Also store the unmodified y2 value to be used for bisection fallback

        # Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if not bisection and x2 - x1 > threshold and abs(y3) > 0.5 * ymin:
            bisection = True
            side = 0

    return _NAN
