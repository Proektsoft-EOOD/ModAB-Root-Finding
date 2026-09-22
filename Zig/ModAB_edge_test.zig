//! Edge-case tests for the shared modAB safeguards.
//!
//! The 100-problem benchmark suite never produces a non-finite or overflowing
//! residual, so the overflow/NaN branches of same_nonzero_sign, safe_midpoint,
//! safe_secant, symmetry_factor and passes_switching_test are covered here.
//!
//! Expected values come from the C# reference in C#/Root/Node.cs and
//! C#/Root/Solvers/ModABCorr.cs; every language port is held to the same table.
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

pub fn main() void {
    // --- same_nonzero_sign ---
    ckb("same_nonzero_sign", 0, M.sameNonzeroSign(1.0, 2.0), true);
    ckb("same_nonzero_sign", 1, M.sameNonzeroSign(-1.0, -2.0), true);
    ckb("same_nonzero_sign", 2, M.sameNonzeroSign(1.0, -2.0), false);
    ckb("same_nonzero_sign", 3, M.sameNonzeroSign(-1.0, 2.0), false);
    ckb("same_nonzero_sign", 4, M.sameNonzeroSign(0.0, 1.0), false);
    ckb("same_nonzero_sign", 5, M.sameNonzeroSign(1.0, 0.0), false);
    ckb("same_nonzero_sign", 6, M.sameNonzeroSign(0.0, 0.0), false);
    ckb("same_nonzero_sign", 7, M.sameNonzeroSign(-0.0, -1.0), false);
    ckb("same_nonzero_sign", 8, M.sameNonzeroSign(nan, 1.0), false);
    ckb("same_nonzero_sign", 9, M.sameNonzeroSign(1.0, nan), false);
    ckb("same_nonzero_sign", 10, M.sameNonzeroSign(nan, nan), false);
    ckb("same_nonzero_sign", 11, M.sameNonzeroSign(inf, 1.0), true);
    ckb("same_nonzero_sign", 12, M.sameNonzeroSign(-inf, -1.0), true);
    ckb("same_nonzero_sign", 13, M.sameNonzeroSign(inf, -inf), false);

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

    // --- symmetry_factor ---
    ck("symmetry_factor", 0, M.symmetryFactor(-1.0, 1.0), 1.0);
    ck("symmetry_factor", 1, M.symmetryFactor(-1.0, 3.0), 0.5625);
    ck("symmetry_factor", 2, M.symmetryFactor(-3.0, 1.0), 0.5625);
    ck("symmetry_factor", 3, M.symmetryFactor(-1e+308, 1e+308), 1.0);
    ck("symmetry_factor", 4, M.symmetryFactor(-1.7e+308, 1.0), 0.25);
    ck("symmetry_factor", 5, M.symmetryFactor(inf, -1.0), nan);
    ck("symmetry_factor", 6, M.symmetryFactor(-1.0, inf), nan);
    ck("symmetry_factor", 7, M.symmetryFactor(nan, 1.0), nan);
    ck("symmetry_factor", 8, M.symmetryFactor(-1e-300, 1e+300), 0.25);

    // --- passes_switching_test ---
    ckb("passes_switching_test", 0, M.passesSwitchingTest(1.0, 1.0, 0.5), true);
    ckb("passes_switching_test", 1, M.passesSwitchingTest(1.0, -1.0, 0.5), false);
    ckb("passes_switching_test", 2, M.passesSwitchingTest(1.0, 0.5, 0.5), true);
    ckb("passes_switching_test", 3, M.passesSwitchingTest(1.0, 0.5, 0.1), false);
    ckb("passes_switching_test", 4, M.passesSwitchingTest(1e+308, 1e+308, 0.5), true);
    ckb("passes_switching_test", 5, M.passesSwitchingTest(1.7e+308, -1.7e+308, 0.5), false);
    ckb("passes_switching_test", 6, M.passesSwitchingTest(1e+308, 1.6e+308, 0.5), true);
    ckb("passes_switching_test", 7, M.passesSwitchingTest(inf, 1.0, 0.5), false);
    ckb("passes_switching_test", 8, M.passesSwitchingTest(1.0, inf, 0.5), false);
    ckb("passes_switching_test", 9, M.passesSwitchingTest(nan, 1.0, 0.5), false);
    ckb("passes_switching_test", 10, M.passesSwitchingTest(1.0, 1.0, nan), false);
    ckb("passes_switching_test", 11, M.passesSwitchingTest(0.0, 0.0, 0.5), false);

    // Solver-level regression: this returned inf before safeMidpoint.
    total += 1;
    const root = M.modAB(fOverflow, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
    if (@abs(root - 1.2e308) <= 1e-13 * 1.2e308) {
        passed += 1;
    } else {
        std.debug.print("FAIL overflowing bracket: got {d} want 1.2e308\n", .{root});
    }

    std.debug.print("Zig edge-cases: {d}/{d} {s}\n", .{
        passed, total, if (passed == total) "PASS" else "FAIL",
    });
    if (passed != total) std.process.exit(1);
}
