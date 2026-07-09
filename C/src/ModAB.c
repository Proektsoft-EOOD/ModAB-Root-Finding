#include <math.h>
#include <stdbool.h>

// Cross-platform export macro
#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT __attribute__((visibility("default")))
#endif

// Static variable to store evaluation count
static int evaluation_count = 0;

// Function to get the evaluation count
EXPORT int get_evaluation_count(void) {
    return evaluation_count;
}

// Inline sign check - returns true if same sign (both positive or both negative)
static inline bool same_sign(double a, double b) {
    return (a > 0) == (b > 0);
}

static inline double clamp(double d, double min, double max) {
  return d < min ? min : (d > max ? max : d);
}

// Finds the root of "F(x) = 0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method:
//     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//     Improvements to the Modified Anderson–Björck(modAB) Root-Finding Algorithm.
//     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
// Additional fixes proposed by L. Tomov are included in this version:
//     1. The secant point is clamped to the interval [p1.X, p2.X] before the X-convergence exit
//     2. The original function values y1 and y2 (without A&B corrections) 
//        are stored for later use in bisection fallback
// F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))
EXPORT double modAB_find_root(double (*f)(double), double x1, double x2, double aTol, double rTol, int maxIter) {
    evaluation_count = 0;

    // Ensure x1 < x2
    if (x1 > x2) {
        double temp = x1; x1 = x2; x2 = temp;
    }

    double y1 = f(x1);
    if (y1 == 0.0)
        return x1;

    double y2 = f(x2);
    if (y2 == 0.0)
        return x2;

    if (same_sign(y1, y2))
        return NAN; // No sign change - no root guaranteed

    bool bisection = true;
    int side = 0; // -1 for left, 1 for right, 0 for none
    double threshold = x2 - x1;
    double f1 = y1, f2 = y2;
    const double C = 16.0; // Safety factor

    for (int i = 1; i <= maxIter; ++i) {
        double x3 = bisection ? 0.5 * (x1 + x2) : (x1 * y2 - y1 * x2) / (y2 - y1);
        // Check for x-convergence
        double eps = aTol + rTol * fabs(x3);
        if (x2 - x1 <= eps) {
            evaluation_count = i + 1; // Saves one function evaluation if satisfied
            return bisection ? x3 : clamp(x3, x1, x2);
        }
        
        double y3;
        if (bisection) {
            y3 = f(x3);
            double ym = 0.5 * (f1 + f2);
            double r = 1.0 - fabs(ym / (f2 - f1)); // Symmetry factor
            double k = r * r; // Deviation factor
            // Check if function is close enough to straight line
            if (fabs(ym - y3) < k * (fabs(y3) + fabs(ym))) {
                bisection = false;
                threshold = (x2 - x1) * C;
            }
        } else {
            // Clamp secant point to interval to handle floating-point errors
            if (x3 <= x1) {
                x3 = x1; y3 = f1;
            } else if (x3 >= x2) {
                x3 = x2; y3 = f2;
            } else {
                y3 = f(x3);
            }
            threshold *= 0.5;
        }

        // Check for y-convergence
        if (y3 == 0.0) {
            evaluation_count = i + 2;
            return x3;
        }

        if (same_sign(y1, y3)) {
            if (side == 1) { // Anderson-Bjork correction
                double m = 1.0 - y3 / y1;
                y2 = (m > 0.0) ? y2 * m : y2 * 0.5;
            } else if (!bisection) {
                side = 1;
            }
            x1 = x3; f1 = y1 = y3;
        } else {
            if (side == -1) { // Anderson-Bjork correction
                double m = 1.0 - y3 / y2;
                y1 = (m > 0.0) ? y1 * m : y1 * 0.5;
            } else if (!bisection) {
                side = -1;
            }
            x2 = x3; f2 = y2 = y3;
        }

        if (x2 - x1 > threshold) {
            bisection = true;
            side = 0;
        }
    }
    evaluation_count = maxIter + 2;
    return NAN;
}