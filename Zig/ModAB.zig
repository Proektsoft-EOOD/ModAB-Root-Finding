const std = @import("std");
// Finds the root of "F(x) = 0" within the interval [x1, x2]
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
pub fn modAB(F: *const fn (f64) f64, x1_: f64, x2_: f64, y0: f64, precision: f64) f64 {
    var x1 = @min(x1_, x2_);
    var x2 = @max(x1_, x2_);
    var y1 = F(x1) - y0;
    if (y1 == 0.0) return x1;
    var y2 = F(x2) - y0;
    if (y2 == 0.0) return x2;

    if (std.math.sign(y1) == std.math.sign(y2)) {
        return std.math.nan(f64);
    }
    var f1 = y1;
    var f2 = y2;
    var bisecting = true;
    var side: i32 = 0;
    const epsilon = if (precision > 0.0) precision else 1e-14;
    var threshold = x2 - x1;
    const C: f64 = 16.0;
    var i: usize = 0;
    const max_iterations: usize = 200;

    while (i < max_iterations) : (i += 1) {
        var x3: f64 = undefined;
        var y3: f64 = undefined;
        if (bisecting) {
            x3 = (x1 + x2) / 2.0;
            y3 = F(x3) - y0;
            const ym = (f1 + f2) / 2.0;
            const dy = f2 - f1;
            const r = 1.0 - @abs(ym / dy);
            const k = r * r;
            if (@abs(ym - y3) < k * (@abs(ym) + @abs(y3))) {
                bisecting = false;
                threshold = (x2 - x1) * C;
            }
        } else {
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1);
            if (x3 <= x1) {
                x3 = x1;
                y3 = f1;
            } else if (x3 >= x2) {
                x3 = x2;
                y3 = f2;
            } else {
                y3 = F(x3) - y0;
            }
            threshold /= 2.0;
        }

        if (y3 == 0.0) {
            return x3;
        }

        if ((x2 - x1) < 2.0 * epsilon) {
            return x3;
        }

        if (std.math.sign(y1) == std.math.sign(y3)) {
            if (side == 1) {
                const m = 1.0 - y3 / y1;
                if (m <= 0.0) {
                    y2 /= 2.0;
                } else {
                    y2 *= m;
                }
            } else if (!bisecting) {
                side = 1;
            }
            x1 = x3;
            y1 = y3;
            f1 = y3;
        } else {
            if (side == -1) {
                const m = 1.0 - y3 / y2;
                if (m <= 0.0) {
                    y1 /= 2.0;
                } else {
                    y1 *= m;
                }
            } else if (!bisecting) {
                side = -1;
            }
            x2 = x3;
            y2 = y3;
            f2 = y3;
        }

        if ((x2 - x1) > threshold) {
            bisecting = true;
            side = 0;
        }
    }
    return std.math.nan(f64);
}
