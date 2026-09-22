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
    // in C#/Root/Node.cs and C#/Root/Solvers/{Solver,ModABCorr}.cs, so a change
    // to the numerics belongs in one place rather than inline in the solver.
    // -----------------------------------------------------------------------

    /** Smallest positive subnormal double (2^-1074), the analogue of C#'s double.Epsilon. */
    private static final double MIN_SUBNORMAL = Double.MIN_VALUE;

    /**
     * Returns true only when both values have the same non-zero sign.
     * <p>
     * Comparisons with NaN are false, so a caller must reject NaN before using
     * this predicate to update a bracket.
     */
    public static boolean sameNonzeroSign(double x, double y) {
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

    /**
     * Returns k = r^2 for the symmetry-sensitive switching criterion. The
     * calculation is homogeneous in the true endpoint residuals.
     */
    public static double symmetryFactor(double y1, double y2) {
        double a = Math.abs(y1);
        double b = Math.abs(y2);
        double den = a + b;

        if (Double.isInfinite(den)) {
            // Infinite true residuals deliberately disable switching and keep
            // the controller in bisection mode. NaN is returned rather than an
            // infinity because every exit of passesSwitchingTest is a "<"
            // comparison, which is false against NaN; an infinity would instead
            // satisfy it and switch. Residuals are never zero here, so only an
            // overflowing sum remains, and halving both restores it without
            // changing the ratio.
            if (Double.isInfinite(a) || Double.isInfinite(b)) {
                return Double.NaN;
            }
            a *= 0.5;
            b *= 0.5;
            den = a + b;
        }

        // |b-a| <= den, so the quotient lies in [0,1]; halving after the
        // division avoids forming 2*den, which could overflow.
        double r = 1.0 - Math.abs(b - a) / den / 2.0;
        return r * r;
    }

    /**
     * Tests whether the true midpoint value yf is close enough to the midpoint
     * value ym of the chord through the true endpoint residuals.
     */
    public static boolean passesSwitchingTest(double ym, double yf, double symmetry) {
        double absYm = Math.abs(ym);
        double absYf = Math.abs(yf);
        double sum = absYf + absYm;

        // Fast path. The exact-root case is handled before this is called, and
        // a non-finite ordinate or a NaN symmetry factor fails the comparison,
        // which disables switching as intended.
        if (Double.isFinite(sum)) {
            return Math.abs(ym - yf) < symmetry * sum;
        }

        // Only reached when the sum overflows. Non-finite values are unsuitable
        // for the linearity comparison.
        if (!Double.isFinite(ym) || !Double.isFinite(yf)) {
            return false;
        }

        // Normalize both sides of the homogeneous inequality to avoid overflow.
        double scale = Math.max(absYf, absYm);
        double normYm = ym / scale;
        double normYf = yf / scale;
        return Math.abs(normYm - normYf) < symmetry * (Math.abs(normYf) + Math.abs(normYm));
    }

    /** The Anderson-Bjorck contraction factor for the ordinate that did not move. */
    private static double abFactor(double y3, double yMoved) {
        double m = 1.0 - y3 / yMoved;
        return m > 0.0 ? m : 0.5;
    }

    /**
     * Multiplies an auxiliary Anderson-Bjorck ordinate by a positive factor
     * while preserving a finite non-zero sign in binary64 arithmetic. This
     * keeps {@link #sameNonzeroSign} sound: an auxiliary ordinate that
     * underflowed to zero would otherwise silently change which branch of the
     * bracket update is taken. It acts only on auxiliary ordinates; an
     * underflowed working value is never accepted as a root of f.
     */
    private static double scalePreservingNonzeroSign(double value, double positiveFactor) {
        double scaled = value * positiveFactor;

        if (scaled == 0.0 && value != 0.0) {
            return Math.copySign(MIN_SUBNORMAL, value);
        }

        if (Double.isInfinite(scaled)) {
            return Math.copySign(Double.MAX_VALUE, value);
        }

        return scaled;
    }

    /**
     * Finds the root of f(x) = y within [x1, x2] using modified Anderson-Björk method.
     * f(x) must be continuous and sign(f(x1) - y) != sign(f(x2) - y).
     * <p>
     * The overflow- and NaN-safe forms of the interpolation and switching
     * arithmetic live in the shared safeguards above.
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
        // NaN has no usable sign, and sameNonzeroSign is false for it, so it
        // must be rejected before the predicate is used to update a bracket.
        if (Double.isNaN(y1) || Double.isNaN(y2) || sameNonzeroSign(y1, y2)) {
            return Double.NaN; // No sign change - no root guaranteed
        }
        double f1 = y1, f2 = y2, ymin = 0.0;
        int side = 0;
        boolean bisection = true;
        double threshold = x2 - x1;
        for (int i = 0; i < maxiter; i++) {
            // safeSecant already returns a point inside [x1, x2], so the
            // separate clamp on the convergence exit is no longer needed.
            double x3 = bisection ? safeMidpoint(x1, x2) : safeSecant(x1, y1, x2, y2);
            double epsx = xtol * Math.max(Math.abs(x3), 1);
            if (x2 - x1 <= epsx) { // x-convergence check
                return x3;
            }
            double y3;
            if (bisection) {
                y3 = f.applyAsDouble(x3) - y;
                double ym = safeMidpoint(f1, f2);
                if (passesSwitchingTest(ym, y3, symmetryFactor(f1, f2))) {
                    bisection = false;
                    threshold = (x2 - x1) * 2.0;
                }
            } else {
                // If rounding makes the proposal coincide with an endpoint,
                // reuse the true residual already stored there.
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

            if (sameNonzeroSign(y1, y3)) {
                if (side == 1) {
                    y2 = scalePreservingNonzeroSign(y2, abFactor(y3, y1));
                } else if (!bisection) {
                    side = 1;
                }
                x1 = x3;
                y1 = y3;
                f1 = y3;
            } else {
                if (side == -1) {
                    y1 = scalePreservingNonzeroSign(y1, abFactor(y3, y2));
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
