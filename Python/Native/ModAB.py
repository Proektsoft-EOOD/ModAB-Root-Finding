import math
def modAB_root(f, x1, x2, y, xtol=1e-14, ytol=0.0, maxiter=200):
    """
    Finds the root of f(x) = y within [x1, x2] using
    modified Anderson-Björk method (Ganchovski, Traykov, 2023).
    f(x) must be continuous and sign(f(x1) - y) ≠ sign(f(x2) - y).
    It laso includes the recent changes from the following paper:
    Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A. 
    Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332. 
    https://doi.org/10.3390/a19050332
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

    side = 0
    bisection = True
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    for _ in range(maxiter):
        x3 = (x1 + x2) * 0.5 if bisection else (x1 * y2 - y1 * x2) / (y2 - y1)
        epsx = xtol * max(abs(x3), 1)
        if x2 - x1 <= epsx: # x-convergence check
            return x3
        
        if bisection:
            y3 = f(x3) - y  # Function value at midpoint
            ym = (y1 + y2) * 0.5  # Ordinate of chord at midpoint
            dy = y2 - y1
            r = 1 - abs(ym / dy)  # Symmetry factor
            k = r * r             # Deviation factor
            if abs(ym - y3) < k * (abs(y3) + abs(ym)):
                bisection = False
                threshold = (x2 - x1) * 8  # Safety factor: 4 bisection iterations = 2^4
        else:
            if x3 <= x1:
                x3, y3 = x1, y1
            elif x3 >= x2:
                x3, y3 = x2, y2
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
            x1, y1 = x3, y3
        else:
            if side == -1:
                m = 1 - y3 / y2
                y1 = y1 * m if m > 0 else y1 * 0.5
            elif not bisection:
                side = -1
            x2, y2 = x3, y3

        if x2 - x1 > threshold:  # AB failed to shrink the interval enough
            bisection = True
            side = 0

    return float('nan')
