import java.util.function.DoubleUnaryOperator;
/*
* Finds the root of "F(x) = 0" within the interval [x1, x2]
* with the specified precisions - absolute: aTol and relative: rTol,
* using an improved version of the modified Anderson Bjork's method:
*     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
*     Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
*     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
* Additional fixes proposed by L. Tomov are applied in this version:
*     1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
*     2. The original function values y1 and y2 (without A&B corrections)
*        are stored for later use in bisection fallback
 */
public class ModAB {

    // -----------------------------------------------------------------------
    // Shared safeguards
    //
    // These helpers carry the overflow/NaN postconditions that the bracketing
    // solver relies on. They are the Java form of the reference implementation
    // in C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs, so a change
    // to the numerics belongs in one place rather than inline in the solver.
    // -----------------------------------------------------------------------

    /**
     * Returns true only when both values have the same non-zero sign.
     * <p>
     * Comparisons with NaN are false, so a caller must reject NaN before using
     * this predicate to update a bracket.
     */
    public static boolean sameSign(double x, double y) {
        return (x < 0.0 && y < 0.0) || (x > 0.0 && y > 0.0);
    }

    /**
     * Midpoint without first forming x1 + x2, which can overflow for finite
     * endpoints of the same sign. For finite ordered endpoints the result lies
     * in the closed bracket, so callers need no extra clamp.
     */
    public static double safeMidpoint(double x1, double x2) {
        return 0.5 * x1 + 0.5 * x2;
    }

    /**
     * Safeguarded false-position/secant point for an ordered bracket [x1, x2].
     * <p>
     * The textbook formula {@code (x1*y2 - x2*y1) / (y2 - y1)} may overflow
     * although the mathematical intersection is finite. With opposite-sign
     * ordinates, writing a = |y1| and b = |y2| gives the equivalent convex
     * combination
     * <pre>
     *     x = (b/(a+b))*x1 + (a/(a+b))*x2,
     * </pre>
     * whose weights lie in [0, 1] and sum to 1. This helper owns the complete
     * postcondition every caller needs: the returned point is finite and lies
     * in [x1, x2]. Where the secant geometry cannot deliver that, the safe
     * midpoint is returned instead.
     */
    public static double safeSecant(double x1, double y1, double x2, double y2) {
        double a = Math.abs(y1);
        double b = Math.abs(y2);
        double den = a + b;

        // One test on the denominator covers every unusable case: a NaN
        // ordinate propagates into it, two zero ordinates make it zero, and an
        // infinite ordinate or an overflowing sum makes it infinite. A zero
        // magnitude does NOT indicate a root here, because the ordinates may be
        // Anderson-Bjorck auxiliary values, so bisection is the safe and
        // mathematically neutral fallback.
        if (!(den > 0.0)) {
            return safeMidpoint(x1, x2);
        }

        if (Double.isInfinite(den)) {
            // An infinite ordinate carries no usable slope. Otherwise a + b
            // merely overflowed, and halving both restores it without changing
            // the ratio that defines the weights.
            if (Double.isInfinite(a) || Double.isInfinite(b)) {
                return safeMidpoint(x1, x2);
            }
            a *= 0.5;
            b *= 0.5;
            den = a + b;
        }

        double x = (b / den) * x1 + (a / den) * x2;

        // In exact arithmetic the convex combination is strictly inside the
        // bracket; the projection only corrects a possible last-ulp excursion.
        return x < x1 ? x1 : (x > x2 ? x2 : x);
    }

    /** The Anderson-Bjorck contraction factor for the ordinate that did not move. */
    private static double abFactor(double y3, double yMoved) {
        double m = 1.0 - y3 / yMoved;
        return m > 0.0 ? m : 0.5;
    }

    /**
     * Finds the root of f(x) = y within [x1, x2] using modified Anderson-Björk method.
     * f(x) must be continuous and sign(f(x1) - y) != sign(f(x2) - y).
     * <p>
     * The overflow- and NaN-safe form of the interpolation lives in the shared
     * safeguards above; the switching test is written so that it cannot overflow.
     *
     * @param f       The function to find the root of
     * @param x1      Left boundary of the interval
     * @param x2      Right boundary of the interval
     * @param y       Target value (finds x where f(x) = y)
     * @param xtol    Tolerance for x convergence (default 1e-14)
     * @param ytol    Tolerance for y convergence (default 0.0)
     * @param maxiter Maximum number of iterations (default 200)
     * @return The root, or NaN if not found within maxiter
     */
    public static double modABRoot(DoubleUnaryOperator f, double x1, double x2,
                                    double y, double xtol, double ytol, int maxiter) {
        if (x2 < x1) {
            double temp = x1;
            x1 = x2;
            x2 = temp;
        }
        double epsy = ytol * Math.max(Math.abs(y), 1);
        double y1 = f.applyAsDouble(x1) - y;
        if (Math.abs(y1) <= epsy) {
            return x1;
        }
        double y2 = f.applyAsDouble(x2) - y;
        if (Math.abs(y2) <= epsy) {
            return x2;
        }
        // NaN has no usable sign, and sameSign is false for it, so it
        // must be rejected before the predicate is used to update a bracket.
        if (Double.isNaN(y1) || Double.isNaN(y2) || sameSign(y1, y2)) {
            return Double.NaN; // No sign change - no root guaranteed
        }
        double f1 = y1, f2 = y2; // True residuals, kept unmodified by A&B corrections
        double ymin = 0.0; // Best true residual of the bracket
        int side = 0;
        boolean bisection = true;
        double threshold = x2 - x1;
        for (int i = 0; i < maxiter; i++) {
            double x3 = bisection ? safeMidpoint(x1, x2) : safeSecant(x1, y1, x2, y2);
            double epsx = xtol * Math.max(Math.abs(x3), 1);
            if (x2 - x1 <= epsx) { // x-convergence check
                return x3;
            }
            double y3;
            if (bisection) {
                y3 = f.applyAsDouble(x3) - y;
                if (Double.isFinite(f2 - f1)) { // Avoids overflow in the calculations below
                    double ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                    double r = 1.0 - Math.abs(ym / (f2 - f1)); // Symmetry factor
                    double k = r * r; // Deviation factor
                    // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                    if (Math.abs(ym - y3) < k * Math.abs(ym) + k * Math.abs(y3)) {
                        bisection = false;
                        threshold = 2.0 * (x2 - x1);
                    }
                }
            } else {
                // If x3 got clamped, reuse the true residual stored at the endpoint.
                if (x3 == x1) {
                    y3 = f1;
                } else if (x3 == x2) {
                    y3 = f2;
                } else {
                    y3 = f.applyAsDouble(x3) - y;
                }
                threshold *= 0.5;
                ymin = Math.min(Math.abs(f1), Math.abs(f2));
            }

            if (Math.abs(y3) <= epsy) {
                return x3;
            }

            // A NaN residual has no usable sign, so the bracket cannot be updated.
            if (Double.isNaN(y3)) {
                return Double.NaN;
            }

            if (sameSign(f1, y3)) {
                if (side == 1) {
                    y2 *= abFactor(y3, y1);
                } else if (!bisection) {
                    side = 1;
                }
                x1 = x3;
                y1 = y3;
                f1 = y3;
            } else {
                if (side == -1) {
                    y1 *= abFactor(y3, y2);
                } else if (!bisection) {
                    side = -1;
                }
                x2 = x3;
                y2 = y3;
                f2 = y3;
            }
            if (!bisection && x2 - x1 > threshold && Math.abs(y3) > 0.5 * ymin) {
                bisection = true;
                side = 0;
            }
        }
        return Double.NaN;
    }

    /**
     * Overloaded method with default parameters.
     */
    public static double modABRoot(DoubleUnaryOperator f, double x1, double x2, double y) {
        return modABRoot(f, x1, x2, y, 1e-14, 0.0, 200);
    }

    /**
     * Overloaded method for finding zero (y=0) with default parameters.
     */
    public static double modABRoot(DoubleUnaryOperator f, double x1, double x2) {
        return modABRoot(f, x1, x2, 0, 1e-14, 0.0, 200);
    }
}
