import java.util.function.DoubleUnaryOperator;

/*
 * Modified Anderson-Björk root-finding algorithm.
 * Based on: Ganchovski, Traykov (2023) and improvements from:
 * Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
 * Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm.
 * Algorithms 2026, 19, 332.
 * https://doi.org/10.3390/a19050332
 */
public class ModAB {

    /**
     * Finds the root of f(x) = y within [x1, x2] using modified Anderson-Björk method.
     * f(x) must be continuous and sign(f(x1) - y) ≠ sign(f(x2) - y).
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
        if ((y1 > 0) == (y2 > 0)) {
            return Double.NaN; // No sign change - no root guaranteed
        }
        int side = 0;
        boolean bisection = true;
        double threshold = x2 - x1;
        for (int i = 0; i < maxiter; i++) {
            double x3 = bisection ? (x1 + x2) * 0.5 : (x1 * y2 - y1 * x2) / (y2 - y1);
            double epsx = xtol * Math.max(Math.abs(x3), 1);
            if (x2 - x1 <= epsx) {
                return x3;
            }
            double y3;
            if (bisection) {
                y3 = f.applyAsDouble(x3) - y;
                double ym = (y1 + y2) * 0.5;
                double dy = y2 - y1;
                double r = 1 - Math.abs(ym / dy);
                double k = r * r;
                if (Math.abs(ym - y3) < k * (Math.abs(y3) + Math.abs(ym))) {
                    bisection = false;
                    threshold = (x2 - x1) * 8;
                }
            } else {
                if (x3 <= x1) {
                    x3 = x1;
                    y3 = y1;
                } else if (x3 >= x2) {
                    x3 = x2;
                    y3 = y2;
                } else {
                    y3 = f.applyAsDouble(x3) - y;
                }
                threshold *= 0.5;
            }
            if (Math.abs(y3) <= epsy) {
                return x3;
            }
            if ((y1 > 0) == (y3 > 0)) {
                if (side == 1) {
                    double m = 1 - y3 / y1;
                    y2 = m > 0 ? y2 * m : y2 * 0.5;
                } else if (!bisection) {
                    side = 1;
                }
                x1 = x3;
                y1 = y3;
            } else {
                if (side == -1) {
                    double m = 1 - y3 / y2;
                    y1 = m > 0 ? y1 * m : y1 * 0.5;
                } else if (!bisection) {
                    side = -1;
                }
                x2 = x3;
                y2 = y3;
            }
            if (x2 - x1 > threshold) {
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
