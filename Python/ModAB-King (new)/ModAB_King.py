import math

def modAB_root(f, x1, x2, y, xtol=1e-14, ytol=0.0, maxiter=200):
    """
    King-order variant of the modified Anderson-Björk method.
    It is identical to the original ModAB except for the lines marked "#King (new)".
    The only change: in the interpolation phase, attempt an explicit
    3-point Anderson-Björck step (Newton on the parabola through the last
    three evaluated points). It is used only if it stays inside the bracket
    and makes at least Brent-grade progress; otherwise the algorithm falls
    back to the original modAB step exactly. This lifts the asymptotic
    order from 5**(1/3) ~ 1.71 to the tribonacci constant ~ 1.839 while
    preserving the bracketing safety net and worst-case factor-1/2 floor.
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
        return float('nan')  # No sign change - no root guaranteed

    side = 0
    bisection = True
    threshold = x2 - x1                      # Threshold to fall back to bisection if AB fails
    xb, fb, xc, fc, nh = x1, y1, x2, y2, 2   #King (new): last evaluated points
    last_span = abs(x2 - x1)                 #King (new): previous bracketing move
    for _ in range(maxiter):
        used_ab = False                      #King (new)
        if bisection:
            x3 = (x1 + x2) * 0.5
        else:
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if nh >= 3:                                              #King (new)
                if xc != xb and xb != xa and xc != xa:               #King (new)
                    fbc = (fc - fb) / (xc - xb)                      #King (new)
                    fab = (fb - fa) / (xb - xa)                      #King (new)
                    fabc = (fbc - fab) / (xc - xa)                   #King (new)
                    gp = fbc + fabc * (xc - xb)                      #King (new)
                    if gp != 0.0:                                    #King (new)
                        x_ab = xc - fc / gp                          #King (new)
                        if (x1 < x_ab < x2 and                       #King (new)
                                abs(x_ab - xc) < 0.5 * last_span):   #King (new)
                            x3 = x_ab                                #King (new)
                            used_ab = True                           #King (new)

        epsx = xtol * max(abs(x3), 1)
        if x2 - x1 <= epsx:  # x-convergence check
            return x3

        if bisection:
            y3 = f(x3) - y
            ym = (y1 + y2) * 0.5
            dy = y2 - y1
            r = 1 - abs(ym / dy)
            k = r * r
            if abs(ym - y3) < k * (abs(y3) + abs(ym)):
                bisection = False
                threshold = (x2 - x1) * 16
        else:
            if x3 <= x1:
                x3, y3 = x1, y1
            elif x3 >= x2:
                x3, y3 = x2, y2
            else:
                y3 = f(x3) - y
            threshold *= 0.5

        if abs(y3) <= epsy:  # y-convergence check
            return x3

        if not (x3 == x1 or x3 == x2):   #King (new): record evaluated point
            xa, fa, xb, fb, xc, fc = xb, fb, xc, fc, x3, y3   #King (new)
            if nh < 3: nh += 1           #King (new)

        if (y1 > 0) == (y3 > 0):  # Same sign check
            span = abs(x3 - x1)                                  #King (new)
            if side == 1 and not used_ab:                      #King (new - added: and not used_ab)
                m = 1 - y3 / y1
                y2 = y2 * m if m > 0 else y2 * 0.5
            elif not bisection:
                side = 1
            x1, y1 = x3, y3
        else:
            span = abs(x3 - x2)                                  #King (new)
            if side == -1 and not used_ab:                     #King (new - added: and not used_ab)
                m = 1 - y3 / y2
                y1 = y1 * m if m > 0 else y1 * 0.5
            elif not bisection:
                side = -1
            x2, y2 = x3, y3

        if span > 0:                       #King (new - track move size for the guard)
            last_span = span               #King (new)

        if x2 - x1 > threshold:  # AB failed to shrink the interval enough
            bisection = True
            side = 0

    return float('nan')