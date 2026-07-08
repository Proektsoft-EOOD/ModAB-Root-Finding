import math
def modAB_root(f, x1, x2, y, xtol=1e-14, ytol=0.0, maxiter=200):
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
    if x2 < x1:
        x1, x2 = x2, x1 
        
    epsy = ytol * max(abs(y), 1)
    y1 = f(x1) - y
    if abs(y1) <= epsy:
        return x1

    y2 = f(x2) - y
    if abs(y2) <= epsy:
        return x2   

    if ((y1 > 0) == (y2 > 0)):
        return float('nan') # No sign change - no root guaranteed

    f1, f2 = y1, y2 # Values for symmetry check kept unmodified by A&B corrections
    side = 0
    bisection = True
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    for _ in range(maxiter):
        x3 = (x1 + x2) * 0.5 if bisection else (x1 * y2 - y1 * x2) / (y2 - y1)
        epsx = xtol * max(abs(x3), 1)
        if x2 - x1 <= epsx: # x-convergence check
            return x3 if bisection else max(x1, min(x3, x2)) # Clamp the secant value
        
        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            ym = (f1 + f2) * 0.5  # Ordinate of chord at midpoint
            dy = f2 - f1
            r = 1 - abs(ym / dy)  # Symmetry factor
            k = r * r             # Deviation factor
            if abs(ym - y3) < k * (abs(y3) + abs(ym)):
                bisection = False
                threshold = (x2 - x1) * 16  # Safety factor: 4 bisection iterations = 2^4
        else:
            if x3 <= x1:
                x3, y3 = x1, f1
            elif x3 >= x2:
                x3, y3 = x2, f2
            else:
                y3 = f(x3) - y
                
            threshold *= 0.5

        if abs(y3) <= epsy: # y-convergence check
            return x3

        if (y1 > 0) == (y3 > 0):  # Same sign check
            if side == 1:
                m = 1 - y3 / y1
                y2 = y2 * m if m > 0 else y2 * 0.5
            elif not bisection:
                side = 1
            x1, y1, f1 = x3, y3, y3 # Also store the unmodified y1 value to be used for bisection fallback
        else:
            if side == -1:
                m = 1 - y3 / y2
                y1 = y1 * m if m > 0 else y1 * 0.5
            elif not bisection:
                side = -1
            x2, y2, f2 = x3, y3, y3 # Akso store the unmodified y2 value to be used for bisection fallback

        if x2 - x1 > threshold:  # AB failed to shrink the interval enough
            bisection = True
            side = 0

    return float('nan')
