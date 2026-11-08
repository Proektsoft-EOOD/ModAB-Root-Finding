const std = @import("std");

// ---------------------------------------------------------------------------
// Shared safeguards
//
// These helpers carry the overflow/NaN postconditions that the bracketing
// solver relies on. They are the Zig form of the reference implementation in
// C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs, so a change to the
// numerics belongs in one place rather than inline in the solver.
// ---------------------------------------------------------------------------

/// Returns true only when both values have the same non-zero sign.
///
/// Comparisons with NaN are false, so a caller must reject NaN before using
/// this predicate to update a bracket.
pub inline fn sameSign(x: f64, y: f64) bool {
    return (x < 0.0 and y < 0.0) or (x > 0.0 and y > 0.0);
}

/// Midpoint without first forming x1 + x2, which can overflow for finite
/// endpoints of the same sign. For finite ordered endpoints the result lies in
/// the closed bracket, so callers need no extra clamp.
pub inline fn safeMidpoint(x1: f64, x2: f64) f64 {
    return 0.5 * x1 + 0.5 * x2;
}

/// Safeguarded false-position/secant point for an ordered bracket [x1, x2].
///
/// The textbook formula (x1*y2 - x2*y1) / (y2 - y1) may overflow although the
/// mathematical intersection is finite. With opposite-sign ordinates, writing
/// a = |y1| and b = |y2| gives the equivalent convex combination
///
///     x = (b/(a+b))*x1 + (a/(a+b))*x2,
///
/// whose weights lie in [0, 1] and sum to 1. This helper owns the complete
/// postcondition every caller needs: the returned point is finite and lies in
/// [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint
/// is returned instead.
pub inline fn safeSecant(x1: f64, y1: f64, x2: f64, y2: f64) f64 {
    var a = @abs(y1);
    var b = @abs(y2);
    var den = a + b;

    // One test on the denominator covers every unusable case: a NaN ordinate
    // propagates into it, two zero ordinates make it zero, and an infinite
    // ordinate or an overflowing sum makes it infinite. A zero magnitude does
    // NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    // auxiliary values, so bisection is the safe and neutral fallback.
    if (!(den > 0.0)) return safeMidpoint(x1, x2);

    if (std.math.isInf(den)) {
        // An infinite ordinate carries no usable slope. Otherwise a + b merely
        // overflowed, and halving both restores it without changing the ratio
        // that defines the weights.
        if (std.math.isInf(a) or std.math.isInf(b)) return safeMidpoint(x1, x2);
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    const x = (b / den) * x1 + (a / den) * x2;

    // In exact arithmetic the convex combination is strictly inside the
    // bracket; the projection only corrects a possible last-ulp excursion.
    return std.math.clamp(x, x1, x2);
}

/// The Anderson-Bjorck contraction factor for the ordinate that did not move.
inline fn abFactor(y3: f64, y_moved: f64) f64 {
    const m = 1.0 - y3 / y_moved;
    return if (m > 0.0) m else 0.5;
}

/// Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
/// (x1, f1), (x2, f2) and a third point (x3, f3) outside it, or NaN when it is unusable.
/// The model fits both nearly-quadratic shapes and functions with a pole; for a straight
/// line it gives the secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative
/// to f2, q1 = f1/f2 and q3 = f3/f2 (so the step does not depend on the scale of f), and
/// c = 1 - q1, the root is
///     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
/// two divisions per step, 1/f2 and this one. The step is used only if it moves the secant
/// point xs = x1 - w*q1/c by at most half the distance to the nearer endpoint (on a plateau
/// the model degenerates towards a step function). Multiplied by c, that test reads
///     |t*c + w*q1| <= w*min(|q1|, 1)/2,
/// with no division; it also rejects NaN and infinite values.
/// When the step is unusable, including a NaN third point (no endpoint dropped on that
/// side yet), it returns the A&B secant point through (x1, y1), (x2, y2) instead.
pub inline fn rationalStep(x1: f64, f1: f64, x2: f64, f2: f64, x3: f64, f3: f64, y1: f64, y2: f64) f64 {
    const k = 1.0 / f2;
    const q1 = f1 * k;
    const q3 = f3 * k;
    const c = 1.0 - q1;
    const w = x2 - x1;
    const u = x3 - x1;
    const t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c);
    const a1 = @abs(q1);
    const x = x1 + t;
    if (@abs(t * c + w * q1) <= 0.5 * w * (if (a1 < 1.0) a1 else 1.0) and x <= x2) {
        return x;
    }
    return safeSecant(x1, y1, x2, y2);
}

/// True if the second divided differences on the left (xl, x1, x2) and on the right
/// (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
/// one-signed curvature does worse than the secant. Scaled by positive factors, they are
///     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
///     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
/// so only their signs are compared, with no division or product of the two.
/// False when xl or xr is NaN (no endpoint dropped on that side yet).
pub inline fn curvatureChangesSign(xl: f64, fl: f64, x1: f64, f1: f64, x2: f64, f2: f64, xr: f64, fr: f64) bool {
    const df = f2 - f1;
    const w = x2 - x1;
    const pl = df * (x1 - xl) - (f1 - fl) * w;
    const pr = (fr - f2) * w - df * (xr - x2);
    return (pl < 0.0 and pr > 0.0) or (pl > 0.0 and pr < 0.0);
}

// Finds the root of "F(x) = y0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method:
//     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//     Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
//     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
// Zig translation developed by @Ramsyana (https://github.com/ramsyana/Zig-Math-Algorithms)
// Additional fixes proposed by L. Tomov are applied in this version:
//     1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
//     2. The original function values y1 and y2 (without A&B corrections)
//        are stored for later use in bisection fallback
// The overflow- and NaN-safe form of the interpolation lives in the shared
// safeguards above; the switching test is written so that it cannot overflow.
// The A&B phase adds two changes:
//     1. Rational correction: the A&B secant is replaced by the secant through the true
//        residuals, corrected for curvature by a 3-point rational model through the bracket
//        and the endpoint it most recently dropped. The A&B secant remains the fallback when
//        the correction is unusable or the curvature changes sign within the bracket.
//     2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
//        just evaluated, so that one-sided convergence still closes the bracket.
// F(x) must be continuous and sign(F(x1)) != sign(F(x2))
pub fn modAB(F: *const fn (f64) f64, x1_: f64, x2_: f64, y0: f64, xtol: f64, ytol: f64, maxiter: usize) f64 {
    var x1 = @min(x1_, x2_);
    var x2 = @max(x1_, x2_);

    const epsy = ytol * @max(@abs(y0), 1.0);
    var y1 = F(x1) - y0;
    if (@abs(y1) <= epsy) return x1;
    var y2 = F(x2) - y0;
    if (@abs(y2) <= epsy) return x2;

    // NaN has no usable sign, and sameSign is false for it, so it must
    // be rejected before the predicate is used to update a bracket.
    if (std.math.isNan(y1) or std.math.isNan(y2) or sameSign(y1, y2)) {
        return std.math.nan(f64); // No sign change - no root guaranteed
    }
    var f1 = y1; // Values for the symmetry check, kept unmodified by A&B corrections
    var f2 = y2;
    var bisecting = true;
    var side: i32 = 0;
    var threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    const C: f64 = 2.0; // Threshold safety factor
    // A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    // than 1 % lies on a plateau; the first A&B step is then not rescaled.
    const PLATEAU_RATIO: f64 = 0.99;
    // Each return from A&B to bisection halves the margin of the switching test, so that
    // borderline switches are rejected after unsuccessful A&B phases.
    var margin: f64 = 1.0;
    var last: i32 = 0; // Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none
    const nan = std.math.nan(f64);
    var xl = nan; // Last endpoints dropped on each side
    var fl = nan;
    var xr = nan;
    var fr = nan;
    var i: usize = 0;

    while (i < maxiter) : (i += 1) {
        const dx = x2 - x1; // Bracket width
        var x3: f64 = undefined;
        var epsx: f64 = undefined;
        if (bisecting) {
            x3 = safeMidpoint(x1, x2);
            epsx = xtol * @max(@abs(x3), 1.0);
        } else {
            // Rational correction through the endpoint dropped by the last update, or the A&B
            // secant when the curvature changes sign within the bracket. A missing endpoint is
            // NaN, which fails both tests, so it needs no check of its own.
            if (curvatureChangesSign(xl, fl, x1, f1, x2, f2, xr, fr)) {
                x3 = safeSecant(x1, y1, x2, y2);
            } else if (last == 1) {
                x3 = rationalStep(x1, f1, x2, f2, xl, fl, y1, y2);
            } else {
                x3 = rationalStep(x1, f1, x2, f2, xr, fr, y1, y2);
            }
            epsx = xtol * @max(@abs(x3), 1.0);
            // Closing step: when A&B converges from one side, its next point lands within epsx/2
            // of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
            // towards the far side then crosses the root and closes the bracket.
            const h = 0.5 * epsx;
            if (last == 1) {
                if (x3 - x1 < h and x1 + h < x2) x3 = x1 + h;
            } else if (x2 - x3 < h and x2 - h > x1) {
                x3 = x2 - h;
            }
        }
        if (dx <= epsx) { // x-convergence check
            return x3;
        }

        var y3: f64 = undefined;
        var switching = false; // This bisection switches to A&B
        if (bisecting) {
            y3 = F(x3) - y0; // Function value at midpoint
            const dy = f2 - f1;
            if (std.math.isFinite(dy)) { // Avoids overflow in the calculations below
                const ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                const r = 1.0 - @abs(ym / dy); // Symmetry factor
                const k = r * r; // Deviation factor
                // Check if function is close enough to straight line.
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (@abs(ym - y3) < (k * @abs(ym) + k * @abs(y3)) * margin) {
                    bisecting = false;
                    switching = true;
                    threshold = C * dx;
                    y1 = f1; // A&B starts from the true residuals
                    y2 = f2;
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            if (x3 == x1) {
                y3 = f1;
            } else if (x3 == x2) {
                y3 = f2;
            } else {
                y3 = F(x3) - y0;
            }
            threshold *= 0.5;
        }

        const ay3 = @abs(y3);
        if (!(ay3 > epsy)) { // Exit on y-convergence or NaN with a single branch on the hot path
            if (ay3 <= epsy) {
                return x3;
            }
            // A NaN residual has no usable sign, so the bracket cannot be updated.
            return std.math.nan(f64);
        }
        if (bisecting) {
            if (sameSign(f1, y3)) { // Same sign check
                xl = x1;
                fl = f1;
                x1 = x3;
                f1 = y3; // Also store the unmodified y1 value to be used for bisection fallback
                last = 1;
            } else {
                xr = x2;
                fr = f2;
                x2 = x3;
                f2 = y3; // Also store the unmodified y2 value to be used for bisection fallback
                last = -1;
            }
        } else {
            const yl = f1; // True residuals of the bracket before the update
            const yr = f2;
            if (sameSign(f1, y3)) { // Same sign check
                if (side == 1) { // Anderson-Bjork correction
                    y2 *= abFactor(y3, y1);
                } else if (!switching or ay3 <= PLATEAU_RATIO * @abs(yl)) {
                    side = 1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 != x1) {
                    xl = x1;
                    fl = f1;
                }
                x1 = x3;
                y1 = y3;
                f1 = y3; // Also store the unmodified y1 value to be used for bisection fallback
                last = 1;
            } else {
                if (side == -1) { // Anderson-Bjork correction
                    y1 *= abFactor(y3, y2);
                } else if (!switching or ay3 <= PLATEAU_RATIO * @abs(yr)) {
                    side = -1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 != x2) {
                    xr = x2;
                    fr = f2;
                }
                x2 = x3;
                y2 = y3;
                f2 = y3; // Also store the unmodified y2 value to be used for bisection fallback
                last = -1;
            }
            // Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            // (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if (x2 - x1 > threshold and 42.0 * ay3 >= @min(@abs(yl), @abs(yr))) {
                bisecting = true;
                side = 0;
                margin *= 0.5;
            }
        }
    }
    return std.math.nan(f64);
}
