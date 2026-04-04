
const std = @import("std");

/// Finds the root of `F(x) = y0` within the interval [x1, x2] with the specified precision,
/// using the modified Anderson-Bjork method (Ganchovski, Traykov).
/// `F(x)` must be continuous and sign(F(x1) - y0) ≠ sign(F(x2) - y0).
/// Returns the root if found, or NaN on failure (e.g., max iterations exceeded or invalid bracket).
/// Zig translation developed by @Ramsyana (https://github.com/ramsyana/Zig-Math-Algorithms)
pub fn modAB(F: *const fn (f64) f64, x1_: f64, x2_: f64, y0: f64, precision: f64) f64 {
    var x1 = @min(x1_, x2_);
    var x2 = @max(x1_, x2_);
    var y1 = F(x1) - y0;
    var y2 = F(x2) - y0;

    if (y1 == 0.0) return x1;
    if (y2 == 0.0) return x2;

    if (std.math.sign(y1) == std.math.sign(y2)) {
        return std.math.nan(f64);
    }

    var bisecting = true;
    var side: i32 = 0;
    const epsilon = if (precision > 0.0) precision else 1e-14;
    var threshold = x2 - x1;
    const C: f64 = 16.0;
    var i: usize = 0;
    const max_iterations: usize = 1000;

    while (i < max_iterations) : (i += 1) {
        var x3: f64 = undefined;
        var y3: f64 = undefined;

        if (bisecting) {
            x3 = (x1 + x2) / 2.0;
            y3 = F(x3) - y0;
            const ym = (y1 + y2) / 2.0;
            const r = 1.0 - @abs(ym / (y2 - y1));
            const k = r * r;
            if (@abs(ym - y3) < k * (@abs(ym) + @abs(y3))) {
                bisecting = false;
                threshold = (x2 - x1) * C;
            }
        } else {
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1);
            if (x3 <= x1) {
                x3 = x1;
                y3 = y1;
            } else if (x3 >= x2) {
                x3 = x2;
                y3 = y2;
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
        }

        if ((x2 - x1) > threshold) {
            bisecting = true;
            side = 0;
        }
    }

    return std.math.nan(f64);
}

pub fn main() void {
    const cos = struct {
        fn call(x: f64) f64 { return std.math.cos(x); }
    }.call;
    const result = modAB(cos, 0.0, 2.0, 0.0, 1e-14);
    std.debug.print("cos(x)=0 root: {d:.15}\n", .{result});
}