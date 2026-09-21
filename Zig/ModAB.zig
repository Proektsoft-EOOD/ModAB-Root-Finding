const std = @import("std");
// Finds the root of "F(x) = y0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method:
//     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//     Improvements to the Modified Anderson–Björck(modAB) Root-Finding Algorithm.
//     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
// Zig translation developed by @Ramsyana (https://github.com/ramsyana/Zig-Math-Algorithms)
// Additional fixes proposed by L. Tomov are applied in this version:
//     1. The secant point is clamped to the interval [p1.X, p2.X] before the X-convergence exit
//     2. The original function values y1 and y2 (without A&B corrections)
//        are stored for later use in bisection fallback
// F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))
pub fn modAB(F: *const fn (f64) f64, x1_: f64, x2_: f64, y0: f64, xtol: f64, ytol: f64, maxiter: usize) f64 {
    var x1 = @min(x1_, x2_);
    var x2 = @max(x1_, x2_);

    const epsy = ytol * @max(@abs(y0), 1.0);
    var y1 = F(x1) - y0;
    if (@abs(y1) <= epsy) return x1;
    var y2 = F(x2) - y0;
    if (@abs(y2) <= epsy) return x2;

    if ((y1 > 0.0) == (y2 > 0.0)) {
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
        var x3 = if (bisecting) (x1 + x2) * 0.5 else (x1 * y2 - y1 * x2) / (y2 - y1);
        const epsx = xtol * @max(@abs(x3), 1.0);
        if (x2 - x1 <= epsx) { // x-convergence check
            return if (bisecting) x3 else std.math.clamp(x3, x1, x2); // Clamp the secant value
        }

        var y3: f64 = undefined;
        if (bisecting) {
            y3 = F(x3) - y0; // Function value at midpoint
            const ym = (f1 + f2) * 0.5; // Ordinate of chord at midpoint
            const dy = f2 - f1;
            const r = 1.0 - @abs(ym / dy); // Symmetry factor
            const k = r * r; // Deviation factor
            // Check if function is close enough to straight line
            if (@abs(ym - y3) < k * (@abs(y3) + @abs(ym))) {
                bisecting = false;
                threshold = (x2 - x1) * C;
            }
        } else {
            // Clamp secant point to interval to handle floating-point errors
            if (x3 <= x1) {
                x3 = x1;
                y3 = f1;
            } else if (x3 >= x2) {
                x3 = x2;
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

        if ((y1 > 0.0) == (y3 > 0.0)) { // Same sign check
            if (side == 1) { // Anderson-Bjork correction
                const m = 1.0 - y3 / y1;
                y2 = if (m > 0.0) y2 * m else y2 * 0.5;
            } else if (!bisecting) {
                side = 1;
            }
            x1 = x3;
            y1 = y3;
            f1 = y3; // Also store the unmodified y1 value to be used for bisection fallback
        } else {
            if (side == -1) { // Anderson-Bjork correction
                const m = 1.0 - y3 / y2;
                y1 = if (m > 0.0) y1 * m else y1 * 0.5;
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
