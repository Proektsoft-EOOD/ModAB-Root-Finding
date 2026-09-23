const std = @import("std");

// ---------------------------------------------------------------------------
// Shared safeguards
//
// These helpers carry the overflow/NaN postconditions that the bracketing
// solver relies on. They are the Zig form of the reference implementation in
// C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs, so a change to the
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
    var ymin: f64 = 0.0; // Best true residual of the bracket
    var i: usize = 0;

    while (i < maxiter) : (i += 1) {
        const x3 = if (bisecting) safeMidpoint(x1, x2) else safeSecant(x1, y1, x2, y2);
        const epsx = xtol * @max(@abs(x3), 1.0);
        if (x2 - x1 <= epsx) { // x-convergence check
            return x3;
        }

        var y3: f64 = undefined;
        if (bisecting) {
            y3 = F(x3) - y0; // Function value at midpoint
            if (std.math.isFinite(f2 - f1)) { // Avoids overflow in the calculations below
                const ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                const r = 1.0 - @abs(ym / (f2 - f1)); // Symmetry factor
                const k = r * r; // Deviation factor
                // Check if function is close enough to straight line.
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (@abs(ym - y3) < k * @abs(ym) + k * @abs(y3)) {
                    bisecting = false;
                    threshold = C * (x2 - x1);
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
            ymin = @min(@abs(f1), @abs(f2));
        }

        if (@abs(y3) <= epsy) { // y-convergence check
            return x3;
        }

        // A NaN residual has no usable sign, so the bracket cannot be updated.
        if (std.math.isNan(y3)) {
            return std.math.nan(f64);
        }

        if (sameSign(f1, y3)) { // Same sign check
            if (side == 1) { // Anderson-Bjork correction
                y2 *= abFactor(y3, y1);
            } else if (!bisecting) {
                side = 1;
            }
            x1 = x3;
            y1 = y3;
            f1 = y3; // Also store the unmodified y1 value to be used for bisection fallback
        } else {
            if (side == -1) { // Anderson-Bjork correction
                y1 *= abFactor(y3, y2);
            } else if (!bisecting) {
                side = -1;
            }
            x2 = x3;
            y2 = y3;
            f2 = y3; // Also store the unmodified y2 value to be used for bisection fallback
        }

        // Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if (!bisecting and x2 - x1 > threshold and @abs(y3) > 0.5 * ymin) {
            bisecting = true;
            side = 0;
        }
    }
    return std.math.nan(f64);
}
