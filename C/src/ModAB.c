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

// ---------------------------------------------------------------------------
// Shared safeguards
//
// These helpers carry the overflow/NaN postconditions that the bracketing
// solvers rely on. They are the C form of the reference implementation in
// C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs, so a change to the
// numerics belongs in one place per language rather than inline in the solver.
// ---------------------------------------------------------------------------

// Returns true only when both values have the same non-zero sign.
// Comparisons with NaN are false, so a caller must reject NaN before using this
// predicate to update a bracket.
static inline bool same_sign(double x, double y) {
    return (x < 0.0 && y < 0.0) || (x > 0.0 && y > 0.0);
}

// Midpoint without first forming x1 + x2, which can overflow for finite
// endpoints of the same sign. For finite ordered endpoints the result lies in
// the closed bracket, so callers need no extra clamp.
static inline double safe_midpoint(double x1, double x2) {
    return 0.5 * x1 + 0.5 * x2;
}

// Safeguarded false-position/secant point for an ordered bracket [x1, x2].
//
// The textbook formula (x1*y2 - x2*y1) / (y2 - y1) may overflow although the
// mathematical intersection is finite. With opposite-sign ordinates, writing
// a = |y1| and b = |y2| gives the equivalent convex combination
//
//     x = (b/(a+b))*x1 + (a/(a+b))*x2,
//
// whose weights lie in [0, 1] and sum to 1. This helper owns the complete
// postcondition every caller needs: the returned point is finite and lies in
// [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint is
// returned instead.
static inline double safe_secant(double x1, double y1, double x2, double y2) {
    double a = fabs(y1);
    double b = fabs(y2);
    double den = a + b;

    // One test on the denominator covers every unusable case: a NaN ordinate
    // propagates into it, two zero ordinates make it zero, and an infinite
    // ordinate or an overflowing sum makes it infinite. A zero magnitude does
    // NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    // auxiliary values, so bisection is the safe and neutral fallback.
    if (!(den > 0.0))
        return safe_midpoint(x1, x2);

    if (isinf(den)) {
        // An infinite ordinate carries no usable slope. Otherwise a + b merely
        // overflowed, and halving both restores it without changing the ratio
        // that defines the weights.
        if (isinf(a) || isinf(b))
            return safe_midpoint(x1, x2);
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    double x = (b / den) * x1 + (a / den) * x2;

    // In exact arithmetic the convex combination is strictly inside the
    // bracket; the projection only corrects a possible last-ulp excursion.
    return x < x1 ? x1 : (x > x2 ? x2 : x);
}

// The Anderson-Bjorck contraction factor for the ordinate that did not move.
static inline double ab_factor(double y3, double y_moved) {
    double m = 1.0 - y3 / y_moved;
    return m > 0.0 ? m : 0.5;
}

// Evaluates f(x) and counts the call, so the count is exact on every exit path
static inline double eval(double (*f)(double), double x) {
    ++evaluation_count;
    return f(x);
}

static inline void swap(double *a, double *b) {
    double c = *a;
    *a = *b;
    *b = c;
}

// Finds the root of "F(x) = 0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method:
//     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//     Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
//     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
// Additional fixes proposed by L. Tomov are included in this version:
//     1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
//     2. The original function values y1 and y2 (without A&B corrections)
//        are stored for later use in bisection fallback
// The overflow- and NaN-safe form of the interpolation lives in the shared
// safeguards above; the switching test is written so that it cannot overflow.
// F(x) must be continuous and sign(F(x1)) != sign(F(x2))
EXPORT double modAB_find_root(double (*f)(double), double x1, double x2, double aTol, double rTol, int maxIter) {
    evaluation_count = 0;
    if (x1 > x2) // Ensure x1 < x2
        swap(&x1, &x2);

    double y1 = eval(f, x1);
    if (y1 == 0.0) return x1;
    double y2 = eval(f, x2);
    if (y2 == 0.0) return x2;

    // NaN has no usable sign, and same_sign is false for it, so it must
    // be rejected before the predicate is used to update a bracket.
    if (isnan(y1) || isnan(y2) || same_sign(y1, y2))
        return NAN; // No sign change - no root guaranteed

    bool bisection = true;
    int side = 0; // -1 for left, 1 for right, 0 for none
    double threshold = x2 - x1; // Bisection fallback threshold
    double f1 = y1, f2 = y2; // True residuals, kept unmodified by A&B corrections
    double ymin = 0.0; // Best true residual of the bracket
    const double C = 2.0; // Safety factor
    for (int i = 1; i <= maxIter; ++i) {
        double x3 = bisection ? safe_midpoint(x1, x2) : safe_secant(x1, y1, x2, y2);
        // Check for x-convergence
        if (x2 - x1 <= aTol + rTol * fabs(x3)) return x3;
        double y3;
        if (bisection) {
            y3 = eval(f, x3);
            if (isfinite(f2 - f1)) { // Avoids overflow in the calculations below
                double ym = (f1 + f2) * 0.5; // Ordinate of chord at midpoint; f1, f2 have opposite signs
                double r = 1.0 - fabs(ym / (f2 - f1)); // Symmetry factor
                double k = r * r; // Deviation factor
                // Check if function is close enough to straight line.
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (fabs(ym - y3) < k * fabs(ym) + k * fabs(y3)) {
                    bisection = false;
                    threshold = C * (x2 - x1);
                    y1 = f1; y2 = f2; // A&B starts from the true residuals
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            y3 = (x3 == x1) ? f1 : 
                 (x3 == x2) ? f2 : 
                 eval(f, x3);
            threshold *= 0.5;
            ymin = fmin(fabs(f1), fabs(f2));
        }
        if (y3 == 0.0) return x3; // Check for y-convergence
        if (isnan(y3)) return NAN; // A NaN residual has no usable sign, so the bracket cannot be updated.
        if (bisection)
        {
            if (same_sign(f1, y3)) {
                x1 = x3; f1 = y3;
            } else {
                x2 = x3; f2 = y3;
            }
        }
        else
        {
            if (same_sign(f1, y3)) {
                if (side == 1) {
                    y2 *= ab_factor(y3, y1); // Anderson-Bjork correction
                } else {
                    side = 1;
                }
                x1 = x3; f1 = y1 = y3;
            } else {
                if (side == -1) {
                    y1 *= ab_factor(y3, y2); // Anderson-Bjork correction
                } else {
                    side = -1;
                }
                x2 = x3; f2 = y2 = y3;
            }
            // Fallback if AB fails to reduce the bracket width, unless it still halves the residual
            if (x2 - x1 > threshold && fabs(y3) > 0.5 * ymin) {
                bisection = true;
                side = 0;
            }
        }
    }
    return NAN;
}
