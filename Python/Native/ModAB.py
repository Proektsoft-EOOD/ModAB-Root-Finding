import math

_NAN = float('nan')
_INF = float('inf')


# ---------------------------------------------------------------------------
# Shared safeguards
#
# These helpers carry the overflow/NaN postconditions that the bracketing
# solver relies on. They are the Python form of the reference implementation in
# C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs, so a change to the
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
    The A&B phase adds two changes:
        1. Rational correction: the A&B secant is corrected for curvature by a 3-point rational
           model through the A&B-corrected bracket and the endpoint it most recently dropped,
           with the step kept within 0.7 of the distance from the A&B secant point to the
           nearer endpoint. The plain A&B secant is used when the curvature changes sign
           within the bracket.
        2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
           just evaluated, so that one-sided convergence still closes the bracket.

    Rational correction. The linear-fractional model  f(x) = (x - r) / (p*x + q)  through
    the A&B-corrected bracket (x1, y1), (x2, y2) and the dropped endpoint (xp, fp) outside it,
    with its true residual, fits both nearly-quadratic shapes and functions with a pole; for a
    straight line it gives the secant point. With w = x2 - x1, u = xp - x1 and the residuals
    relative to y2, q1 = y1/y2 and qp = fp/y2 (so the step does not depend on the scale of f),
    and c = 1 - q1 > 1, its root is
        r = x1 + t,   t = w*u*q1*(qp - 1) / (w*(qp - q1) - u*qp*c),
    two divisions per step, 1/y2 and this one. The A&B secant point is xs = x1 - w*q1/c and its
    distance to the nearer endpoint is d = w*min(|q1|, 1)/c. The step is kept within 0.7*d of
    xs (on a plateau the model degenerates towards a step function), so it stays at least 0.3*d
    inside the bracket. Multiplied by c, the test reads
        |t*c + w*q1| <= 0.7*w*min(|q1|, 1),
    with no division. A step outside it is clamped to xs +- 0.7*d on the side the model points
    to; this did better than the A&B secant in randomized-bracket tests and on 200,000 random
    functions, where factors 0.6 to 0.9 differed little. A NaN step, or one that rounds onto
    x2, gives x2, whose residual is already stored. The curvature changes sign
    within the bracket when the second divided differences on the left (xl, x1, x2) and on
    the right (x1, x2, xr) have opposite signs. Scaled by positive factors, they are
        pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*w,   pr = (fr - f2)*w - (f2 - f1)*(xr - x2),
    and only their signs are compared.
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
    # A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    # than 1 % lies on a plateau; the first A&B step is then not rescaled.
    PLATEAU_RATIO = 0.99
    # Each return from A&B to bisection halves the margin of the switching test, so that
    # borderline switches are rejected after unsuccessful A&B phases.
    margin = 1.0
    last = 0  # Endpoint set by the last evaluation: 1 = x1, -1 = x2
    xl = fl = xr = fr = _NAN  # Last endpoints dropped on the left and on the right
    for _ in range(maxiter):
        dx = x2 - x1  # Bracket width
        if bisection:
            x3 = safe_midpoint(x1, x2)
            epsx = xtol * max(abs(x3), 1)
        else:
            # Rational correction through the endpoint dropped by the last update, or the A&B
            # secant when the curvature changes sign within the bracket. A missing endpoint is
            # NaN, which fails the curvature test and makes the rational step NaN.
            df = f2 - f1
            pl = df * (x1 - xl) - (f1 - fl) * dx
            pr = (fr - f2) * dx - df * (xr - x2)
            if pl < 0.0 < pr or pr < 0.0 < pl:
                x3 = safe_secant(x1, y1, x2, y2)
            else:
                if last == 1:
                    xp, fp = xl, fl
                else:
                    xp, fp = xr, fr
                # Rational step through the A&B-corrected bracket (x1, y1), (x2, y2) and the
                # dropped endpoint (xp, fp), clamped to 0.7 of the distance from the A&B secant
                # point to the nearer endpoint; a NaN step, or one that rounds onto x2, gives x2.
                # Python raises where C gives inf or NaN, so those cases are spelled out: a zero
                # y2 makes the C step NaN, and a zero denominator gives an infinite or NaN t.
                x3 = x2
                if y2 != 0.0:
                    k = 1.0 / y2  # Residuals relative to y2, so the step does not depend on the scale of f
                    q1 = y1 * k
                    qp = fp * k
                    c = 1.0 - q1
                    u = xp - x1
                    num = dx * u * q1 * (qp - 1.0)
                    den = dx * (qp - q1) - u * qp * c
                    if den != 0.0:
                        t = num / den
                    elif num != 0.0:
                        t = math.copysign(_INF, num) * math.copysign(1.0, den)
                    else:
                        t = _NAN
                    a = abs(q1)
                    lim = 0.7 * (dx * a if a < 1.0 else dx)  # c * 0.7 * d
                    dev = t * c + dx * q1  # c * (x - xs)
                    if abs(dev) > lim:  # False for NaN: then t is NaN and the test below gives x2
                        t = (math.copysign(lim, dev) - dx * q1) / c
                    t += x1
                    if t < x2:
                        x3 = t
            epsx = xtol * max(abs(x3), 1)
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

        switching = False  # This bisection switches to A&B
        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            dy = f2 - f1
            if math.isfinite(dy):  # Avoids overflow in the calculations below
                ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint; f1, f2 have opposite signs
                r = 1.0 - abs(ym / dy)  # Symmetry factor
                k = r * r  # Deviation factor
                # k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if abs(ym - y3) < (k * abs(ym) + k * abs(y3)) * margin:
                    bisection = False
                    switching = True
                    threshold = C * dx  # Safety factor: skips two AB steps before the first fallback
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

        ay3 = abs(y3)
        if not ay3 > epsy:  # Exit on y-convergence or NaN with a single branch on the hot path
            if ay3 <= epsy:
                return x3
            # A NaN residual has no usable sign, so the bracket cannot be updated.
            return _NAN
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
                elif not switching or ay3 <= PLATEAU_RATIO * abs(yl):
                    side = 1  # A switching bisection on a plateau starts the A&B phase afresh
                if x3 != x1:
                    xl, fl = x1, f1
                x1, y1, f1, last = x3, y3, y3, 1  # Also store the unmodified y1 value to be used for bisection fallback
            else:
                if side == -1:
                    y1 *= ab_factor(y3, y2)
                elif not switching or ay3 <= PLATEAU_RATIO * abs(yr):
                    side = -1  # A switching bisection on a plateau starts the A&B phase afresh
                if x3 != x2:
                    xr, fr = x2, f2
                x2, y2, f2, last = x3, y3, y3, -1  # Also store the unmodified y2 value to be used for bisection fallback

            # Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            # (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if x2 - x1 > threshold and 42 * ay3 >= min(abs(yl), abs(yr)):
                bisection = True
                side = 0
                margin *= 0.5

    return _NAN
