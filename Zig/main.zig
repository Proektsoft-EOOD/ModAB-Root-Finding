//! ModAB Root-Finding Test Suite
//!
//! Tests the modified Anderson-Bjork algorithm against 92 benchmark functions
//! from various academic sources:
//!   - Sergio Galdino: Regula falsi methods (f01-f33)
//!   - Steven A. Stage: Brent's method improvements (f34-f40)
//!   - Swift & Lindfield: Continuation methods (f41-f62)
//!   - Oliveira & Takahashi: Bisection enhancement (f63-f83)
//!   - SciML project: Benchmark suite (f84-f92)

const std = @import("std");
const modAB = @import("ModAB.zig").modAB;

const Problem = struct {
    name: []const u8,
    f: *const fn (f64) f64,
    a: f64,
    b: f64,
};

// Helper function P(x) = x + 1.11111
fn pp(x: f64) f64 {
    return x + 1.11111;
}

// Sign function
fn sgn(x: f64) f64 {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

// ============================================================================
// Test functions f01-f33: Sergio Galdino (Regula falsi methods)
// ============================================================================

fn f01(x: f64) f64 { return std.math.pow(f64, x, 3) - 1.0; }
fn f02(x: f64) f64 { return x * x * (x * x / 3.0 + @sqrt(2.0) * @sin(x)) - @sqrt(3.0) / 18.0; }
fn f03(x: f64) f64 { return 11.0 * std.math.pow(f64, x, 11) - 1.0; }
fn f04(x: f64) f64 { return std.math.pow(f64, x, 3) + 1.0; }
fn f05(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * x - 5.0; }
fn f06(x: f64) f64 { return 2.0 * x * @exp(-5.0) + 1.0 - 2.0 * @exp(-5.0 * x); }
fn f07(x: f64) f64 { return 2.0 * x * @exp(-10.0) + 1.0 - 2.0 * @exp(-10.0 * x); }
fn f08(x: f64) f64 { return 2.0 * x * @exp(-20.0) + 1.0 - 2.0 * @exp(-20.0 * x); }
fn f09(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 5.0, 2)) * x * x - std.math.pow(f64, 1.0 - 5.0 * x, 2); }
fn f10(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 10.0, 2)) * x * x - std.math.pow(f64, 1.0 - 10.0 * x, 2); }
fn f11(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 20.0, 2)) * x * x - std.math.pow(f64, 1.0 - 20.0 * x, 2); }
fn f12(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 5); }
fn f13(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 10); }
fn f14(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 20); }
fn f15(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 5.0, 4)) * x - std.math.pow(f64, 1.0 - 5.0 * x, 4); }
fn f16(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 10.0, 4)) * x - std.math.pow(f64, 1.0 - 10.0 * x, 4); }
fn f17(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 20.0, 4)) * x - std.math.pow(f64, 1.0 - 20.0 * x, 4); }
fn f18(x: f64) f64 { return @exp(-5.0 * x) * (x - 1.0) + std.math.pow(f64, x, 5); }
fn f19(x: f64) f64 { return @exp(-10.0 * x) * (x - 1.0) + std.math.pow(f64, x, 10); }
fn f20(x: f64) f64 { return @exp(-20.0 * x) * (x - 1.0) + std.math.pow(f64, x, 20); }
fn f21(x: f64) f64 { return x * x + @sin(x / 5.0) - 0.25; }
fn f22(x: f64) f64 { return x * x + @sin(x / 10.0) - 0.25; }
fn f23(x: f64) f64 { return x * x + @sin(x / 20.0) - 0.25; }
fn f24(x: f64) f64 { return (x + 2.0) * (x + 1.0) * std.math.pow(f64, x - 3.0, 3); }
fn f25(x: f64) f64 { return std.math.pow(f64, x - 4.0, 5) * @log(x); }
fn f26(x: f64) f64 { return std.math.pow(f64, @sin(x) - x / 4.0, 3); }
fn f27(x: f64) f64 {
    const q = pp(x);
    const s: f64 = if (q < 3.0) 1.0 else if (q > 3.0) -1.0 else 0.0;
    return (81.0 - q * (108.0 - q * (54.0 - q * (12.0 - q)))) * s;
}
fn f28(x: f64) f64 { return @sin(std.math.pow(f64, x - 7.143, 3)); }
fn f29(x: f64) f64 { return @exp(std.math.pow(f64, x - 3.0, 5)) - 1.0; }
fn f30(x: f64) f64 { return @exp(std.math.pow(f64, x - 3.0, 5)) - @exp(x - 1.0); }
fn f31(x: f64) f64 { return std.math.pi - 1.0 / x; }
fn f32(x: f64) f64 { return 4.0 - @tan(x); }
fn f33(x: f64) f64 { return @cos(x) - std.math.pow(f64, x, 3); }

// ============================================================================
// Test functions f34-f40: Steven A. Stage (Brent's method improvements)
// ============================================================================

fn f34(x: f64) f64 { return @cos(x) - x; }
fn f35(x: f64) f64 {
    const s: f64 = if (x <= 2.0 / 3.0) 1.0 else -1.0;
    return @sqrt(@abs(x - 2.0 / 3.0)) * s - 0.1;
}
fn f36(x: f64) f64 {
    const s: f64 = if (x <= 2.0 / 3.0) 1.0 else -1.0;
    return std.math.pow(f64, @abs(x - 2.0 / 3.0), 0.2) * s;
}
fn f37(x: f64) f64 { return std.math.pow(f64, x - 7.0 / 9.0, 3) + (x - 7.0 / 9.0) * 1e-3; }
fn f38(x: f64) f64 { return if (x <= 1.0 / 3.0) -0.5 else 0.5; }
fn f39(x: f64) f64 { return if (x <= 1.0 / 3.0) -1e-3 else 1.0 - 1e-3; }
fn f40(x: f64) f64 { return if (x == 0.0) 0.0 else 1.0 / (x - 2.0 / 3.0); }

// ============================================================================
// Test functions f41-f62: Swift & Lindfield (Continuation methods)
// ============================================================================

fn f41(x: f64) f64 { return 2.0 * x * @exp(-5.0) - 2.0 * @exp(-5.0 * x) + 1.0; }
fn f42(x: f64) f64 { return (x * x - x - 6.0) * (x * x - 3.0 * x + 2.0); }
fn f43(x: f64) f64 { return std.math.pow(f64, x, 3); }
fn f44(x: f64) f64 { return std.math.pow(f64, x, 5); }
fn f45(x: f64) f64 { return std.math.pow(f64, x, 7); }
fn f46(x: f64) f64 { return (@exp(-5.0 * x) - x - 0.5) / std.math.pow(f64, x, 5); }
fn f47(x: f64) f64 { return 1.0 / @sqrt(x) - 2.0 * @log(5e3 * @sqrt(x)) + 0.8; }
fn f48(x: f64) f64 { return 1.0 / @sqrt(x) - 2.0 * @log(5e7 * @sqrt(x)) + 0.8; }
fn f49(x: f64) f64 {
    if (x <= 0.0) {
        return -std.math.pow(f64, x, 3) - x - 1.0;
    } else {
        return std.math.cbrt(x) - x - 1.0;
    }
}
fn f50(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * x - x + 3.0; }
fn f51(x: f64) f64 { return @log(x); }
fn f52(x: f64) f64 { return (10.0 - x) * @exp(-10.0 * x) - std.math.pow(f64, x, 10) + 1.0; }
fn f53(x: f64) f64 { return @exp(@sin(x)) - x - 1.0; }
fn f54(x: f64) f64 { return 2.0 * @sin(x) - 1.0; }
fn f55(x: f64) f64 { return (x - 1.0) * @exp(-x); }
fn f56(x: f64) f64 { return std.math.pow(f64, x - 1.0, 3) - 1.0; }
fn f57(x: f64) f64 { return @exp(x * x + 7.0 * x - 30.0) - 1.0; }
fn f58(x: f64) f64 { return std.math.atan(x) - 1.0; }
fn f59(x: f64) f64 { return @exp(x) - 2.0 * x - 1.0; }
fn f60(x: f64) f64 { return @exp(-x) - x - @sin(x); }
fn f61(x: f64) f64 { return x * x - std.math.pow(f64, @sin(x), 2) - 1.0; }
fn f62(x: f64) f64 { return @sin(x) - x / 2.0; }

// ============================================================================
// Test functions f63-f83: Oliveira & Takahashi (Bisection enhancement)
// ============================================================================

fn f63(x: f64) f64 { return x * @exp(x) - 1.0; }
fn f64_fn(x: f64) f64 { return @tan(x - 0.1); }
fn f65(x: f64) f64 { return @sin(x) + 0.5; }
fn f66(x: f64) f64 { return 4.0 * std.math.pow(f64, x, 5) + x * x + 1.0; }
fn f67(x: f64) f64 { return x + std.math.pow(f64, x, 10) - 1.0; }
fn f68(x: f64) f64 { return std.math.pow(f64, std.math.pi, x) - std.math.e; }
fn f69(x: f64) f64 { return @log(@abs(x - 10.0 / 9.0)); }
fn f70(x: f64) f64 { return 1.0 / 3.0 + std.math.cbrt(x) + std.math.pow(f64, x, 3); }
fn f71(x: f64) f64 { return (x + 2.0 / 3.0) / (x + 101.0 / 100.0); }
fn f72(x: f64) f64 { return std.math.pow(f64, x * 1e6 - 1.0, 3); }
fn f73(x: f64) f64 { return @exp(x) * std.math.pow(f64, x * 1e6 - 1.0, 3); }
fn f74(x: f64) f64 { return std.math.pow(f64, x - 1.0 / 3.0, 2) * std.math.atan(x - 1.0 / 3.0); }
fn f75(x: f64) f64 {
    const t = 3.0 * x - 1.0;
    return sgn(t) * (1.0 - @sqrt(1.0 - t * t / 81.0));
}
fn f76(x: f64) f64 { return if (x > (1.0 - 1e6) / 1e6) (1.0 + 1e6) / 1e6 else -1.0; }
fn f77(x: f64) f64 { return if (x != 1.0 / 21.0) 1.0 / (21.0 * x - 1.0) else 0.0; }
fn f78(x: f64) f64 { return x * x / 4.0 + @ceil(x / 2.0) - 0.5; }
fn f79(x: f64) f64 { return @ceil(10.0 * x - 1.0) + 0.5; }
fn f80(x: f64) f64 { return x + @sin(x * 1e6) / 10.0 + 1e-3; }
fn f81(x: f64) f64 { return if (x > -1.0) 1.0 + @sin(1.0 / (x + 1.0)) else -1.0; }
fn f82(x: f64) f64 { return 202.0 * x - 2.0 * @floor((2.0 * x + 1e-2) / 2e-2) - 0.1; }
fn f83(x: f64) f64 {
    const t = 202.0 * x - 2.0 * @floor((2.0 * x + 1e-2) / 2e-2) - 0.1;
    return t * t * t;
}

// ============================================================================
// Test functions f84-f92: SciML project benchmarks
// ============================================================================

fn f84(x: f64) f64 { return (x - 1.0) * (x - 2.0) * (x - 3.0) * (x - 4.0) * (x - 5.0) - 0.05; }
fn f85(x: f64) f64 { return @sin(x) - 0.5 * x - 0.3; }
fn f86(x: f64) f64 { return @exp(x) - 1.0 - x - x * x / 2.0 - 0.005; }
fn f87(x: f64) f64 { return 1.0 / (x - 0.5) - 2.0 - 0.05; }
fn f88(x: f64) f64 { return @log(x) - x + 2.0 - 0.05; }
fn f89(x: f64) f64 { return @sin(20.0 * x) + 0.1 * x - 0.1; }
fn f90(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * std.math.pow(f64, x, 2) + x - 0.025; }
fn f91(x: f64) f64 { return x * @sin(1.0 / x) - 0.1 - 0.01; }
fn f92(x: f64) f64 { return std.math.pow(f64, x, 3) - 0.001; }

// ============================================================================
// All problems array
// ============================================================================

const all_problems = [_]Problem{
    // Galdino (f01-f33)
    .{ .name = "f01", .f = f01, .a = 0.5, .b = 1.5 },
    .{ .name = "f02", .f = f02, .a = 0.1, .b = 1.0 },
    .{ .name = "f03", .f = f03, .a = 0.1, .b = 1.0 },
    .{ .name = "f04", .f = f04, .a = -1.8, .b = 0.0 },
    .{ .name = "f05", .f = f05, .a = 2.0, .b = 3.0 },
    .{ .name = "f06", .f = f06, .a = 0.0, .b = 1.0 },
    .{ .name = "f07", .f = f07, .a = 0.0, .b = 1.0 },
    .{ .name = "f08", .f = f08, .a = 0.0, .b = 1.0 },
    .{ .name = "f09", .f = f09, .a = 0.0, .b = 1.0 },
    .{ .name = "f10", .f = f10, .a = 0.0, .b = 1.0 },
    .{ .name = "f11", .f = f11, .a = 0.0, .b = 1.0 },
    .{ .name = "f12", .f = f12, .a = 0.0, .b = 1.0 },
    .{ .name = "f13", .f = f13, .a = 0.0, .b = 1.0 },
    .{ .name = "f14", .f = f14, .a = 0.0, .b = 1.0 },
    .{ .name = "f15", .f = f15, .a = 0.0, .b = 1.0 },
    .{ .name = "f16", .f = f16, .a = 0.0, .b = 1.0 },
    .{ .name = "f17", .f = f17, .a = 0.0, .b = 1.0 },
    .{ .name = "f18", .f = f18, .a = 0.0, .b = 1.0 },
    .{ .name = "f19", .f = f19, .a = 0.0, .b = 1.0 },
    .{ .name = "f20", .f = f20, .a = 0.0, .b = 1.0 },
    .{ .name = "f21", .f = f21, .a = 0.0, .b = 1.0 },
    .{ .name = "f22", .f = f22, .a = 0.0, .b = 1.0 },
    .{ .name = "f23", .f = f23, .a = 0.0, .b = 1.0 },
    .{ .name = "f24", .f = f24, .a = 2.6, .b = 4.6 },
    .{ .name = "f25", .f = f25, .a = 3.6, .b = 5.6 },
    .{ .name = "f26", .f = f26, .a = 2.0, .b = 4.0 },
    .{ .name = "f27", .f = f27, .a = 1.0, .b = 3.0 },
    .{ .name = "f28", .f = f28, .a = 7.0, .b = 8.0 },
    .{ .name = "f29", .f = f29, .a = 2.6, .b = 4.6 },
    .{ .name = "f30", .f = f30, .a = 4.0, .b = 5.0 },
    .{ .name = "f31", .f = f31, .a = 0.05, .b = 5.0 },
    .{ .name = "f32", .f = f32, .a = 0.0, .b = 1.5 },
    .{ .name = "f33", .f = f33, .a = 0.0, .b = 4.0 },
    // Stage (f34-f40)
    .{ .name = "f34", .f = f34, .a = -11.0, .b = 9.0 },
    .{ .name = "f35", .f = f35, .a = -11.0, .b = 9.0 },
    .{ .name = "f36", .f = f36, .a = -11.0, .b = 9.0 },
    .{ .name = "f37", .f = f37, .a = -11.0, .b = 9.0 },
    .{ .name = "f38", .f = f38, .a = -11.0, .b = 9.0 },
    .{ .name = "f39", .f = f39, .a = -11.0, .b = 9.0 },
    .{ .name = "f40", .f = f40, .a = -11.0, .b = 9.0 },
    // Swift-Lindfield (f41-f62)
    .{ .name = "f41", .f = f41, .a = 0.0, .b = 10.0 },
    .{ .name = "f42", .f = f42, .a = 0.0, .b = std.math.pi },
    .{ .name = "f43", .f = f43, .a = -1.0, .b = 1.5 },
    .{ .name = "f44", .f = f44, .a = -1.0, .b = 1.5 },
    .{ .name = "f45", .f = f45, .a = -1.0, .b = 1.5 },
    .{ .name = "f46", .f = f46, .a = 0.09, .b = 0.7 },
    .{ .name = "f47", .f = f47, .a = 0.0005, .b = 0.5 },
    .{ .name = "f48", .f = f48, .a = 0.0005, .b = 0.5 },
    .{ .name = "f49", .f = f49, .a = -1.0, .b = 1.0 },
    .{ .name = "f50", .f = f50, .a = -3.0, .b = 2.0 },
    .{ .name = "f51", .f = f51, .a = 0.5, .b = 5.0 },
    .{ .name = "f52", .f = f52, .a = 0.5, .b = 8.0 },
    .{ .name = "f53", .f = f53, .a = 1.0, .b = 4.0 },
    .{ .name = "f54", .f = f54, .a = 0.1, .b = std.math.pi / 3.0 },
    .{ .name = "f55", .f = f55, .a = 0.0, .b = 1.5 },
    .{ .name = "f56", .f = f56, .a = 1.5, .b = 3.0 },
    .{ .name = "f57", .f = f57, .a = 2.6, .b = 3.5 },
    .{ .name = "f58", .f = f58, .a = 1.0, .b = 8.0 },
    .{ .name = "f59", .f = f59, .a = 0.2, .b = 3.0 },
    .{ .name = "f60", .f = f60, .a = 0.0, .b = 2.0 },
    .{ .name = "f61", .f = f61, .a = -1.0, .b = 2.0 },
    .{ .name = "f62", .f = f62, .a = std.math.pi / 2.0, .b = std.math.pi },
    // Oliveira-Takahashi (f63-f83)
    .{ .name = "f63", .f = f63, .a = -1.0, .b = 1.0 },
    .{ .name = "f64", .f = f64_fn, .a = -1.0, .b = 1.0 },
    .{ .name = "f65", .f = f65, .a = -1.0, .b = 1.0 },
    .{ .name = "f66", .f = f66, .a = -1.0, .b = 1.0 },
    .{ .name = "f67", .f = f67, .a = -1.0, .b = 1.0 },
    .{ .name = "f68", .f = f68, .a = -1.0, .b = 1.0 },
    .{ .name = "f69", .f = f69, .a = -1.0, .b = 1.0 },
    .{ .name = "f70", .f = f70, .a = -1.0, .b = 1.0 },
    .{ .name = "f71", .f = f71, .a = -1.0, .b = 1.0 },
    .{ .name = "f72", .f = f72, .a = -1.0, .b = 1.0 },
    .{ .name = "f73", .f = f73, .a = -1.0, .b = 1.0 },
    .{ .name = "f74", .f = f74, .a = -1.0, .b = 1.0 },
    .{ .name = "f75", .f = f75, .a = -1.0, .b = 1.0 },
    .{ .name = "f76", .f = f76, .a = -1.0, .b = 1.0 },
    .{ .name = "f77", .f = f77, .a = -1.0, .b = 1.0 },
    .{ .name = "f78", .f = f78, .a = -1.0, .b = 1.0 },
    .{ .name = "f79", .f = f79, .a = -1.0, .b = 1.0 },
    .{ .name = "f80", .f = f80, .a = -1.0, .b = 1.0 },
    .{ .name = "f81", .f = f81, .a = -1.0, .b = 1.0 },
    .{ .name = "f82", .f = f82, .a = -1.0, .b = 1.0 },
    .{ .name = "f83", .f = f83, .a = -1.0, .b = 1.0 },
    // SciML (f84-f92)
    .{ .name = "f84", .f = f84, .a = 0.5, .b = 5.5 },
    .{ .name = "f85", .f = f85, .a = -10.0, .b = 10.0 },
    .{ .name = "f86", .f = f86, .a = -2.0, .b = 2.0 },
    .{ .name = "f87", .f = f87, .a = 0.6, .b = 2.0 },
    .{ .name = "f88", .f = f88, .a = 0.1, .b = 3.0 },
    .{ .name = "f89", .f = f89, .a = -4.0, .b = 5.0 },
    .{ .name = "f90", .f = f90, .a = -1.0, .b = 2.0 },
    .{ .name = "f91", .f = f91, .a = 0.01, .b = 1.0 },
    .{ .name = "f92", .f = f92, .a = -10.0, .b = 10.0 },
};

pub fn main() void {
    const eps: f64 = 1e-14;
    var passed: u32 = 0;
    var failed: u32 = 0;

    std.debug.print("ModAB Root-Finding Test Results (Zig)\n", .{});
    std.debug.print("======================================================================\n", .{});
    std.debug.print("{s:>4} | {s:>22} | {s:>15} | Status\n", .{ "Func", "Root", "f(root)" });
    std.debug.print("----------------------------------------------------------------------\n", .{});

    for (all_problems) |p| {
        const root = modAB(p.f, p.a, p.b, 0.0, eps);
        const fval = p.f(root);

        var status: []const u8 = undefined;
        if (std.math.isNan(root)) {
            status = "FAIL";
            failed += 1;
        } else if (@abs(fval) < 1e-10) {
            status = "PASS";
            passed += 1;
        } else if (@abs(fval) < 1e-6) {
            status = "PASS";
            passed += 1;
        } else {
            status = "WEAK";
            passed += 1;
        }

        std.debug.print("{s:>4} | {d:>22.15} | {e:>15.6} | {s}\n", .{ p.name, root, fval, status });
    }

    std.debug.print("----------------------------------------------------------------------\n", .{});
    std.debug.print("Total: {} tests, {} passed, {} failed\n", .{ passed + failed, passed, failed });
}
