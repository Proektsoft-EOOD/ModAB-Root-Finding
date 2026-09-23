//! Edge-case tests for the shared modAB safeguards.
//!
//! The 100-problem benchmark suite never produces a non-finite or overflowing
//! residual, so the overflow/NaN branches of same_sign, safe_midpoint and
//! safe_secant, and the overflow cases of the switching test, are covered here.
//!
//! Expected values come from the C# reference in C#/Root/Node.cs and
//! C#/Root/Solvers/SgModAB.cs; every language port is held to the same table.
//!
//! Run (from the Zig directory):
//!     zig run ModAB_edge_test.zig
const std = @import("std");
const M = @import("ModAB.zig");

const inf = std.math.inf(f64);
const nan = std.math.nan(f64);

var passed: usize = 0;
var total: usize = 0;

fn ck(name: []const u8, i: usize, got: f64, want: f64) void {
    total += 1;
    if ((std.math.isNan(got) and std.math.isNan(want)) or got == want) {
        passed += 1;
        return;
    }
    std.debug.print("FAIL {s}[{d}]: got {d} want {d}\n", .{ name, i, got, want });
}

fn ckb(name: []const u8, i: usize, got: bool, want: bool) void {
    total += 1;
    if (got == want) {
        passed += 1;
        return;
    }
    std.debug.print("FAIL {s}[{d}]: got {} want {}\n", .{ name, i, got, want });
}

/// f(x) = x*1e-308 - 1.2; the root is 1.2e308 and x1 + x2 overflows.
fn fOverflow(x: f64) f64 {
    return x * 1e-308 - 1.2;
}

// Solver-level cases for overflow and infinite residuals in the switching test.
const max_f64 = std.math.floatMax(f64);

/// |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite.
fn fHump(x: f64) f64 {
    return if (x <= 0) max_f64 * (-0.2 + 1.2 * (x + 1)) else max_f64 * (1 - 0.2 * x);
}
/// f2 - f1 overflows: switching is disabled until the residuals shrink.
fn fDiffOverflow(x: f64) f64 {
    return 1.7e308 * std.math.tanh(10 * (x - 0.3));
}
/// Infinite residuals at one or both ends of the bracket.
fn fInfLeft(x: f64) f64 {
    return if (x < -0.5) -inf else x - 0.1;
}
fn fInfBoth(x: f64) f64 {
    return if (x < -0.5) -inf else if (x > 0.9) inf else x - 0.1;
}
/// Subnormal residuals: AB corrections underflow towards zero.
fn fSubnormal(x: f64) f64 {
    return 1e-300 * (x * x * x - 0.2);
}

fn ckRoot(name: []const u8, f: *const fn (f64) f64, a: f64, b: f64, want: f64) void {
    total += 1;
    const got = M.modAB(f, a, b, 0.0, 1e-14, 0.0, 200);
    if (@abs(got - want) <= 1e-13 * @max(@abs(want), 1.0)) {
        passed += 1;
        return;
    }
    std.debug.print("FAIL {s}: got {d} want {d}\n", .{ name, got, want });
}

pub fn main() void {
    // --- same_sign ---
    ckb("same_sign", 0, M.sameSign(1.0, 2.0), true);
    ckb("same_sign", 1, M.sameSign(-1.0, -2.0), true);
    ckb("same_sign", 2, M.sameSign(1.0, -2.0), false);
    ckb("same_sign", 3, M.sameSign(-1.0, 2.0), false);
    ckb("same_sign", 4, M.sameSign(0.0, 1.0), false);
    ckb("same_sign", 5, M.sameSign(1.0, 0.0), false);
    ckb("same_sign", 6, M.sameSign(0.0, 0.0), false);
    ckb("same_sign", 7, M.sameSign(-0.0, -1.0), false);
    ckb("same_sign", 8, M.sameSign(nan, 1.0), false);
    ckb("same_sign", 9, M.sameSign(1.0, nan), false);
    ckb("same_sign", 10, M.sameSign(nan, nan), false);
    ckb("same_sign", 11, M.sameSign(inf, 1.0), true);
    ckb("same_sign", 12, M.sameSign(-inf, -1.0), true);
    ckb("same_sign", 13, M.sameSign(inf, -inf), false);

    // --- safe_midpoint ---
    ck("safe_midpoint", 0, M.safeMidpoint(2.0, 4.0), 3.0);
    ck("safe_midpoint", 1, M.safeMidpoint(1e+308, 1e+308), 1e+308);
    ck("safe_midpoint", 2, M.safeMidpoint(-1e+308, 1e+308), 0.0);
    ck("safe_midpoint", 3, M.safeMidpoint(1e+308, 1.7e+308), 1.35e+308);
    ck("safe_midpoint", 4, M.safeMidpoint(-1.7e+308, -1e+308), -1.35e+308);
    ck("safe_midpoint", 5, M.safeMidpoint(0.0, 1.0), 0.5);
    ck("safe_midpoint", 6, M.safeMidpoint(-1.0, 1.0), 0.0);

    // --- safe_secant ---
    ck("safe_secant", 0, M.safeSecant(0.0, -1.0, 1.0, 1.0), 0.5);
    ck("safe_secant", 1, M.safeSecant(0.0, -1.0, 1.0, 3.0), 0.25);
    ck("safe_secant", 2, M.safeSecant(0.0, -1e+308, 1.0, 1e+308), 0.5);
    ck("safe_secant", 3, M.safeSecant(0.0, -1.7e+308, 1.0, 1.7e+308), 0.5);
    ck("safe_secant", 4, M.safeSecant(0.0, inf, 1.0, -1.0), 0.5);
    ck("safe_secant", 5, M.safeSecant(0.0, -1.0, 1.0, inf), 0.5);
    ck("safe_secant", 6, M.safeSecant(0.0, 0.0, 1.0, 0.0), 0.5);
    ck("safe_secant", 7, M.safeSecant(0.0, nan, 1.0, 1.0), 0.5);
    ck("safe_secant", 8, M.safeSecant(0.0, -1e-300, 1.0, 1e+300), 0.0);
    ck("safe_secant", 9, M.safeSecant(0.0, -1e+300, 1.0, 1e-300), 1.0);
    ck("safe_secant", 10, M.safeSecant(1e+308, -1.0, 1.7e+308, 1.0), 1.35e+308);

    // Solver-level regression: this returned inf before safeMidpoint.
    total += 1;
    const root = M.modAB(fOverflow, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
    if (@abs(root - 1.2e308) <= 1e-13 * 1.2e308) {
        passed += 1;
    } else {
        std.debug.print("FAIL overflowing bracket: got {d} want 1.2e308\n", .{root});
    }

    ckRoot("hump", fHump, -1.0, 1.0, -5.0 / 6.0);
    ckRoot("f2-f1 overflow", fDiffOverflow, -1.0, 1.0, 0.3);
    ckRoot("inf left", fInfLeft, -1.0, 1.0, 0.1);
    ckRoot("inf both", fInfBoth, -1.0, 1.0, 0.1);
    ckRoot("subnormal", fSubnormal, -1.0, 1.0, std.math.cbrt(@as(f64, 0.2)));

    std.debug.print("Zig edge-cases: {d}/{d} {s}\n", .{
        passed, total, if (passed == total) "PASS" else "FAIL",
    });
    if (passed != total) std.process.exit(1);
}
