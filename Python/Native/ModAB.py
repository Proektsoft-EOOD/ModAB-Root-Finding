import math

_NAN = float('nan')


# ---------------------------------------------------------------------------
# Shared safeguards
#
# These helpers carry the overflow/NaN postconditions that the bracketing
# solver relies on. They are the Python form of the reference implementation in
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs, so a change to the
# numerics belongs in one place rather than inline in the solver.
# ---------------------------------------------------------------------------

def same_sign(x, y):
    return (x < 0.0 and y < 0.0) or (x > 0.0 and y > 0.0)


def safe_midpoint(x1, x2):
    return 0.5 * x1 + 0.5 * x2


def safe_secant(x1, y1, x2, y2):
    a = abs(y1)
    b = abs(y2)
    den = a + b
    if not den > 0.0:
        return safe_midpoint(x1, x2)

    if math.isinf(den):
        if math.isinf(a) or math.isinf(b):
            return safe_midpoint(x1, x2)
        a *= 0.5
        b *= 0.5
        den = a + b

    x = (b / den) * x1 + (a / den) * x2
    return x1 if x < x1 else (x2 if x > x2 else x)


def ab_factor(y3, y_moved):
    m = 1.0 - y3 / y_moved
    return m if m > 0.0 else 0.5


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
    The overflow- and NaN-safe form of the interpolation lives in the shared
    safeguards above; the switching test is written so that it cannot overflow.
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

    # NaN has no usable sign, and same_sign is false for it, so it must
    # be rejected before the predicate is used to update a bracket.
    if math.isnan(y1) or math.isnan(y2) or same_sign(y1, y2):
        return _NAN  # No sign change - no root guaranteed

    f1, f2 = y1, y2  # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = True
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    C = 2  # Threshold safety factor
    # Best residual of the bracket.
    ymin = 0
    for _ in range(maxiter):
        x3 = safe_midpoint(x1, x2) if bisection else safe_secant(x1, y1, x2, y2)
        epsx = xtol * max(abs(x3), 1)
        if x2 - x1 <= epsx:  # x-convergence check
            return x3

        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            if math.isfinite(f2 - f1):  # Avoids overflow in the calculations below
                ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint; f1, f2 have opposite signs
                r = 1.0 - abs(ym / (f2 - f1))  # Symmetry factor
                k = r * r  # Deviation factor
                # k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if abs(ym - y3) < k * abs(ym) + k * abs(y3):
                    bisection = False
                    threshold = C * (x2 - x1)  # Safety factor: skips two AB steps before the first fallback
                    y1, y2 = f1, f2  # A&B starts from the true residuals
        else:
            # If x3 got clamped, reuse the true residual stored at the endpoint.
            if x3 == x1:
                y3 = f1
            elif x3 == x2:
                y3 = f2
            else:
                y3 = f(x3) - y

            threshold *= 0.5
            ymin = min(abs(f1), abs(f2)) # Best true residual of the bracket.

        if abs(y3) <= epsy:  # y-convergence check
            return x3

        # A NaN residual has no usable sign, so the bracket cannot be updated.
        if math.isnan(y3):
            return _NAN
        if bisection:
            if same_sign(f1, y3):  # Same sign check
                x1, f1 = x3, y3,
            else:
                x2, f2 = x3, y3
        else:
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

            # Fallback if AB fails to reduce the bracket width, unless it still halves the residual
            if x2 - x1 > threshold and abs(y3) > 0.5 * ymin:
                bisection = True
                side = 0

    return _NAN
