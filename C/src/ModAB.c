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
// C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs, so a change to the
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

// Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
// (x1, f1), (x2, f2) and a third point (x3, f3) outside it, or NaN when it is unusable.
// The model fits both nearly-quadratic shapes and functions with a pole; for a straight
// line it gives the secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative
// to f2, q1 = f1/f2 and q3 = f3/f2 (so the step does not depend on the scale of f), and
// c = 1 - q1, the root is
//     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
// two divisions per step, 1/f2 and this one. The step is used only if it moves the secant
// point xs = x1 - w*q1/c by at most half the distance to the nearer endpoint (on a plateau
// the model degenerates towards a step function). Multiplied by c, that test reads
//     |t*c + w*q1| <= w*min(|q1|, 1)/2,
// with no division; it also rejects NaN and infinite values.
// When the step is unusable, including a NaN third point (no endpoint dropped on that
// side yet), it returns the A&B secant point through (x1, y1), (x2, y2) instead.
static inline double rational_step(double x1, double f1, double x2, double f2,
                                   double x3, double f3, double y1, double y2) {
    double k = 1.0 / f2;
    double q1 = f1 * k, q3 = f3 * k;
    double c = 1.0 - q1;
    double w = x2 - x1, u = x3 - x1;
    double t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c);
    double a1 = fabs(q1);
    double x = x1 + t;
    if (fabs(t * c + w * q1) <= 0.5 * w * (a1 < 1.0 ? a1 : 1.0) && x <= x2)
        return x;
    return safe_secant(x1, y1, x2, y2);
}

// True if the second divided differences on the left (xl, x1, x2) and on the right
// (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
// one-signed curvature does worse than the secant. Scaled by positive factors, they are
//     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
//     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
// so only their signs are compared, with no division or product of the two.
// False when xl or xr is NaN (no endpoint dropped on that side yet).
static inline int curvature_changes_sign(double xl, double fl, double x1, double f1,
                                         double x2, double f2, double xr, double fr) {
    double df = f2 - f1, w = x2 - x1;
    double pl = df * (x1 - xl) - (f1 - fl) * w;
    double pr = (fr - f2) * w - df * (xr - x2);
    return (pl < 0.0 && pr > 0.0) || (pl > 0.0 && pr < 0.0);
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
// The A&B phase adds two changes:
//     1. Rational correction: the A&B secant is replaced by the secant through the true
//        residuals, corrected for curvature by a 3-point rational model through the bracket
//        and the endpoint it most recently dropped. The A&B secant remains the fallback when
//        the correction is unusable or the curvature changes sign within the bracket.
//     2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
//        just evaluated, so that one-sided convergence still closes the bracket.
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
    const double C = 2.0; // Safety factor
    // A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    // than 1 % lies on a plateau; the first A&B step is then not rescaled.
    const double PLATEAU_RATIO = 0.99;
    // Each return from A&B to bisection halves the margin of the switching test, so that
    // borderline switches are rejected after unsuccessful A&B phases.
    double margin = 1.0;
    int last = 0; // Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none
    double xl = NAN, fl = NAN, xr = NAN, fr = NAN; // Last endpoints dropped on each side
    for (int i = 1; i <= maxIter; ++i) {
        double dx = x2 - x1; // Bracket width
        double x3, epsx;
        if (bisection) {
            x3 = safe_midpoint(x1, x2);
            epsx = aTol + rTol * fabs(x3);
        } else {
            // Rational correction through the endpoint dropped by the last update, or the A&B
            // secant when the curvature changes sign within the bracket. A missing endpoint is
            // NaN, which fails both tests, so it needs no check of its own.
            if (curvature_changes_sign(xl, fl, x1, f1, x2, f2, xr, fr))
                x3 = safe_secant(x1, y1, x2, y2);
            else
                x3 = last == 1 ? rational_step(x1, f1, x2, f2, xl, fl, y1, y2)
                               : rational_step(x1, f1, x2, f2, xr, fr, y1, y2);
            epsx = aTol + rTol * fabs(x3);
            // Closing step: when A&B converges from one side, its next point lands within epsx/2
            // of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
            // towards the far side then crosses the root and closes the bracket.
            double h = 0.5 * epsx;
            if (last == 1) {
                if (x3 - x1 < h && x1 + h < x2) x3 = x1 + h;
            } else if (x2 - x3 < h && x2 - h > x1) {
                x3 = x2 - h;
            }
        }
        // Check for x-convergence
        if (dx <= epsx) return x3;
        double y3;
        bool switching = false; // This bisection switches to A&B
        if (bisection) {
            y3 = eval(f, x3);
            double dy = f2 - f1;
            if (isfinite(dy)) { // Avoids overflow in the calculations below
                double ym = (f1 + f2) * 0.5; // Ordinate of chord at midpoint; f1, f2 have opposite signs
                double r = 1.0 - fabs(ym / dy); // Symmetry factor
                double k = r * r; // Deviation factor
                // Check if function is close enough to straight line.
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (fabs(ym - y3) < (k * fabs(ym) + k * fabs(y3)) * margin) {
                    bisection = false;
                    switching = true;
                    threshold = C * dx;
                    y1 = f1; y2 = f2; // A&B starts from the true residuals
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            y3 = (x3 == x1) ? f1 : 
                 (x3 == x2) ? f2 : 
                 eval(f, x3);
            threshold *= 0.5;
        }
        if (!(fabs(y3) > 0.0)) { // Exit on zero or NaN with a single branch on the hot path
            if (y3 == 0.0) return x3; // Check for y-convergence
            return NAN; // A NaN residual has no usable sign, so the bracket cannot be updated.
        }
        if (bisection)
        {
            if (same_sign(f1, y3)) {
                xl = x1; fl = f1;
                x1 = x3; f1 = y3; last = 1;
            } else {
                xr = x2; fr = f2;
                x2 = x3; f2 = y3; last = -1;
            }
        }
        else
        {
            double yl = f1, yr = f2; // True residuals of the bracket before the update
            if (same_sign(f1, y3)) {
                if (side == 1) {
                    y2 *= ab_factor(y3, y1); // Anderson-Bjork correction
                } else if (!switching || fabs(y3) <= PLATEAU_RATIO * fabs(yl)) {
                    side = 1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 != x1) { xl = x1; fl = f1; }
                x1 = x3; f1 = y1 = y3; last = 1;
            } else {
                if (side == -1) {
                    y1 *= ab_factor(y3, y2); // Anderson-Bjork correction
                } else if (!switching || fabs(y3) <= PLATEAU_RATIO * fabs(yr)) {
                    side = -1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 != x2) { xr = x2; fr = f2; }
                x2 = x3; f2 = y2 = y3; last = -1;
            }
            // Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            // (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if (x2 - x1 > threshold && 42.0 * fabs(y3) >= fmin(fabs(yl), fabs(yr))) {
                bisection = true;
                side = 0;
                margin *= 0.5;
            }
        }
    }
    return NAN;
}
