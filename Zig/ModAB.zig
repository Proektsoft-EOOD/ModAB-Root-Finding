const std = @import("std");

// ---------------------------------------------------------------------------
// Shared safeguards
//
// These helpers carry the overflow/NaN postconditions that the bracketing
// solver relies on. They are the Zig form of the reference implementation in
// C#/Root/Node.cs and C#/Root/Solvers/{Solver,ModABCorr}.cs, so a change to the
// numerics belongs in one place rather than inline in the solver.
// ---------------------------------------------------------------------------

/// Smallest positive subnormal f64 (2^-1074), the analogue of C#'s double.Epsilon.
const MIN_SUBNORMAL: f64 = std.math.floatTrueMin(f64);

/// Largest finite f64, the analogue of C#'s double.MaxValue.
const MAX_DOUBLE: f64 = std.math.floatMax(f64);

/// Returns true only when both values have the same non-zero sign.
///
/// Comparisons with NaN are false, so a caller must reject NaN before using
/// this predicate to update a bracket.
pub inline fn sameNonzeroSign(x: f64, y: f64) bool {
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

/// Returns k = r^2 for the symmetry-sensitive switching criterion. The
/// calculation is homogeneous in the true endpoint residuals.
pub inline fn symmetryFactor(y1: f64, y2: f64) f64 {
    var a = @abs(y1);
    var b = @abs(y2);
    var den = a + b;

    if (std.math.isInf(den)) {
        // Infinite true residuals deliberately disable switching and keep the
        // controller in bisection mode. NaN is returned rather than an infinity
        // because every exit of passesSwitchingTest is a "<" comparison, which
        // is false against NaN; an infinity would instead satisfy it and
        // switch. Residuals are never zero here, so only an overflowing sum
        // remains, and halving both restores it without changing the ratio.
        if (std.math.isInf(a) or std.math.isInf(b)) return std.math.nan(f64);
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    // |b-a| <= den, so the quotient lies in [0,1]; halving after the division
    // avoids forming 2*den, which could overflow.
    const r = 1.0 - @abs(b - a) / den / 2.0;
    return r * r;
}

/// Tests whether the true midpoint value yf is close enough to the midpoint
/// value ym of the chord through the true endpoint residuals.
pub inline fn passesSwitchingTest(ym: f64, yf: f64, symmetry: f64) bool {
    const abs_ym = @abs(ym);
    const abs_yf = @abs(yf);
    const sum = abs_yf + abs_ym;

    // Fast path. The exact-root case is handled before this is called, and a
    // non-finite ordinate or a NaN symmetry factor fails the comparison, which
    // disables switching as intended.
    if (std.math.isFinite(sum)) return @abs(ym - yf) < symmetry * sum;

    // Only reached when the sum overflows. Non-finite values are unsuitable for
    // the linearity comparison.
    if (!std.math.isFinite(ym) or !std.math.isFinite(yf)) return false;

    // Normalize both sides of the homogeneous inequality to avoid overflow.
    const scale = @max(abs_yf, abs_ym);
    const norm_ym = ym / scale;
    const norm_yf = yf / scale;
    return @abs(norm_ym - norm_yf) < symmetry * (@abs(norm_yf) + @abs(norm_ym));
}

/// The Anderson-Bjorck contraction factor for the ordinate that did not move.
inline fn abFactor(y3: f64, y_moved: f64) f64 {
    const m = 1.0 - y3 / y_moved;
    return if (m > 0.0) m else 0.5;
}

/// Multiplies an auxiliary Anderson-Bjorck ordinate by a positive factor while
/// preserving a finite non-zero sign in binary64 arithmetic. This keeps
/// sameNonzeroSign sound: an auxiliary ordinate that underflowed to zero would
/// otherwise silently change which branch of the bracket update is taken. It
/// acts only on auxiliary ordinates; an underflowed working value is never
/// accepted as a root of f.
inline fn scalePreservingNonzeroSign(value: f64, positive_factor: f64) f64 {
    const scaled = value * positive_factor;

    if (scaled == 0.0 and value != 0.0) return std.math.copysign(MIN_SUBNORMAL, value);

    if (std.math.isInf(scaled)) return std.math.copysign(MAX_DOUBLE, value);

    return scaled;
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
// The overflow- and NaN-safe forms of the interpolation and switching
// arithmetic live in the shared safeguards above.
// F(x) must be continuous and sign(F(x1)) != sign(F(x2))
pub fn modAB(F: *const fn (f64) f64, x1_: f64, x2_: f64, y0: f64, xtol: f64, ytol: f64, maxiter: usize) f64 {
    var x1 = @min(x1_, x2_);
    var x2 = @max(x1_, x2_);

    const epsy = ytol * @max(@abs(y0), 1.0);
    var y1 = F(x1) - y0;
    if (@abs(y1) <= epsy) return x1;
    var y2 = F(x2) - y0;
    if (@abs(y2) <= epsy) return x2;

    // NaN has no usable sign, and sameNonzeroSign is false for it, so it must
    // be rejected before the predicate is used to update a bracket.
    if (std.math.isNan(y1) or std.math.isNan(y2) or sameNonzeroSign(y1, y2)) {
        return std.math.nan(f64); // No sign change - no root guaranteed
    }
    var f1 = y1; // Values for the symmetry check, kept unmodified by A&B corrections
    var f2 = y2;
    var bisecting = true;
    var side: i32 = 0;
    var threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    const C: f64 = 2.0; // Threshold safety factor
    // Best residual of the bracket, refreshed before each AB step. An AB step that
    // fails the width test is still kept if it at least halved this residual.
    var ymin: f64 = 0.0;
    var i: usize = 0;

    while (i < maxiter) : (i += 1) {
        // safeSecant already returns a point inside [x1, x2], so the separate
        // clamp on the convergence exit is no longer needed.
        const x3 = if (bisecting) safeMidpoint(x1, x2) else safeSecant(x1, y1, x2, y2);
        const epsx = xtol * @max(@abs(x3), 1.0);
        if (x2 - x1 <= epsx) { // x-convergence check
            return x3;
        }

        var y3: f64 = undefined;
        if (bisecting) {
            y3 = F(x3) - y0; // Function value at midpoint
            const ym = safeMidpoint(f1, f2); // Ordinate of chord at midpoint
            // Check if function is close enough to straight line
            if (passesSwitchingTest(ym, y3, symmetryFactor(f1, f2))) {
                bisecting = false;
                threshold = (x2 - x1) * C;
            }
        } else {
            // If rounding makes the proposal coincide with an endpoint, reuse
            // the true residual already stored there.
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

        if (sameNonzeroSign(y1, y3)) { // Same sign check
            if (side == 1) { // Anderson-Bjork correction
                y2 = scalePreservingNonzeroSign(y2, abFactor(y3, y1));
            } else if (!bisecting) {
                side = 1;
            }
            x1 = x3;
            y1 = y3;
            f1 = y3; // Also store the unmodified y1 value to be used for bisection fallback
        } else {
            if (side == -1) { // Anderson-Bjork correction
                y1 = scalePreservingNonzeroSign(y1, abFactor(y3, y2));
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
