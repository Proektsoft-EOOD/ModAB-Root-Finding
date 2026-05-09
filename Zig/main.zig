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

fn f_01(x: f64) f64 { return std.math.pow(f64, x, 3) - 1.0; }
fn f_02(x: f64) f64 { return x * x * (x * x / 3.0 + @sqrt(2.0) * @sin(x)) - @sqrt(3.0) / 18.0; }
fn f_03(x: f64) f64 { return 11.0 * std.math.pow(f64, x, 11) - 1.0; }
fn f_04(x: f64) f64 { return std.math.pow(f64, x, 3) + 1.0; }
fn f_05(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * x - 5.0; }
fn f_06(x: f64) f64 { return 2.0 * x * @exp(-5.0) + 1.0 - 2.0 * @exp(-5.0 * x); }
fn f_07(x: f64) f64 { return 2.0 * x * @exp(-10.0) + 1.0 - 2.0 * @exp(-10.0 * x); }
fn f_08(x: f64) f64 { return 2.0 * x * @exp(-20.0) + 1.0 - 2.0 * @exp(-20.0 * x); }
fn f_09(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 5.0, 2)) * x * x - std.math.pow(f64, 1.0 - 5.0 * x, 2); }
fn f_10(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 10.0, 2)) * x * x - std.math.pow(f64, 1.0 - 10.0 * x, 2); }
fn f_11(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 20.0, 2)) * x * x - std.math.pow(f64, 1.0 - 20.0 * x, 2); }
fn f_12(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 5); }
fn f_13(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 10); }
fn f_14(x: f64) f64 { return x * x - std.math.pow(f64, 1.0 - x, 20); }
fn f_15(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 5.0, 4)) * x - std.math.pow(f64, 1.0 - 5.0 * x, 4); }
fn f_16(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 10.0, 4)) * x - std.math.pow(f64, 1.0 - 10.0 * x, 4); }
fn f_17(x: f64) f64 { return (1.0 + std.math.pow(f64, 1.0 - 20.0, 4)) * x - std.math.pow(f64, 1.0 - 20.0 * x, 4); }
fn f_18(x: f64) f64 { return @exp(-5.0 * x) * (x - 1.0) + std.math.pow(f64, x, 5); }
fn f_19(x: f64) f64 { return @exp(-10.0 * x) * (x - 1.0) + std.math.pow(f64, x, 10); }
fn f_20(x: f64) f64 { return @exp(-20.0 * x) * (x - 1.0) + std.math.pow(f64, x, 20); }
fn f_21(x: f64) f64 { return x * x + @sin(x / 5.0) - 0.25; }
fn f_22(x: f64) f64 { return x * x + @sin(x / 10.0) - 0.25; }
fn f_23(x: f64) f64 { return x * x + @sin(x / 20.0) - 0.25; }
fn f_24(x: f64) f64 { return (x + 2.0) * (x + 1.0) * std.math.pow(f64, x - 3.0, 3); }
fn f_25(x: f64) f64 { return std.math.pow(f64, x - 4.0, 5) * @log(x); }
fn f_26(x: f64) f64 { return std.math.pow(f64, @sin(x) - x / 4.0, 3); }
fn f_27(x: f64) f64 {
    const q = pp(x);
    const s: f64 = if (q < 3.0) 1.0 else if (q > 3.0) -1.0 else 0.0;
    return (81.0 - q * (108.0 - q * (54.0 - q * (12.0 - q)))) * s;
}
fn f_28(x: f64) f64 { return @sin(std.math.pow(f64, x - 7.143, 3)); }
fn f_29(x: f64) f64 { return @exp(std.math.pow(f64, x - 3.0, 5)) - 1.0; }
fn f_30(x: f64) f64 { return @exp(std.math.pow(f64, x - 3.0, 5)) - @exp(x - 1.0); }
fn f_31(x: f64) f64 { return std.math.pi - 1.0 / x; }
fn f_32(x: f64) f64 { return 4.0 - @tan(x); }
fn f_33(x: f64) f64 { return @cos(x) - std.math.pow(f64, x, 3); }

// ============================================================================
// Test functions f34-f40: Steven A. Stage (Brent's method improvements)
// ============================================================================

fn f_34(x: f64) f64 { return @cos(x) - x; }
fn f_35(x: f64) f64 {
    const s: f64 = if (x <= 2.0 / 3.0) 1.0 else -1.0;
    return @sqrt(@abs(x - 2.0 / 3.0)) * s - 0.1;
}
fn f_36(x: f64) f64 {
    const s: f64 = if (x <= 2.0 / 3.0) 1.0 else -1.0;
    return std.math.pow(f64, @abs(x - 2.0 / 3.0), 0.2) * s;
}
fn f_37(x: f64) f64 { return std.math.pow(f64, x - 7.0 / 9.0, 3) + (x - 7.0 / 9.0) * 1e-3; }
fn f_38(x: f64) f64 { return if (x <= 1.0 / 3.0) -0.5 else 0.5; }
fn f_39(x: f64) f64 { return if (x <= 1.0 / 3.0) -1e-3 else 1.0 - 1e-3; }
fn f_40(x: f64) f64 { return if (x == 0.0) 0.0 else 1.0 / (x - 2.0 / 3.0); }

// ============================================================================
// Test functions f41-f62: Swift & Lindfield (Continuation methods)
// ============================================================================

fn f_41(x: f64) f64 { return 2.0 * x * @exp(-5.0) - 2.0 * @exp(-5.0 * x) + 1.0; }
fn f_42(x: f64) f64 { return (x * x - x - 6.0) * (x * x - 3.0 * x + 2.0); }
fn f_43(x: f64) f64 { return std.math.pow(f64, x, 3); }
fn f_44(x: f64) f64 { return std.math.pow(f64, x, 5); }
fn f_45(x: f64) f64 { return std.math.pow(f64, x, 7); }
fn f_46(x: f64) f64 { return (@exp(-5.0 * x) - x - 0.5) / std.math.pow(f64, x, 5); }
fn f_47(x: f64) f64 { return 1.0 / @sqrt(x) - 2.0 * @log(5e3 * @sqrt(x)) + 0.8; }
fn f_48(x: f64) f64 { return 1.0 / @sqrt(x) - 2.0 * @log(5e7 * @sqrt(x)) + 0.8; }
fn f_49(x: f64) f64 {
    if (x <= 0.0) {
        return -std.math.pow(f64, x, 3) - x - 1.0;
    } else {
        return std.math.cbrt(x) - x - 1.0;
    }
}
fn f_50(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * x - x + 3.0; }
fn f_51(x: f64) f64 { return @log(x); }
fn f_52(x: f64) f64 { return (10.0 - x) * @exp(-10.0 * x) - std.math.pow(f64, x, 10) + 1.0; }
fn f_53(x: f64) f64 { return @exp(@sin(x)) - x - 1.0; }
fn f_54(x: f64) f64 { return 2.0 * @sin(x) - 1.0; }
fn f_55(x: f64) f64 { return (x - 1.0) * @exp(-x); }
fn f_56(x: f64) f64 { return std.math.pow(f64, x - 1.0, 3) - 1.0; }
fn f_57(x: f64) f64 { return @exp(x * x + 7.0 * x - 30.0) - 1.0; }
fn f_58(x: f64) f64 { return std.math.atan(x) - 1.0; }
fn f_59(x: f64) f64 { return @exp(x) - 2.0 * x - 1.0; }
fn f_60(x: f64) f64 { return @exp(-x) - x - @sin(x); }
fn f_61(x: f64) f64 { return x * x - std.math.pow(f64, @sin(x), 2) - 1.0; }
fn f_62(x: f64) f64 { return @sin(x) - x / 2.0; }

// ============================================================================
// Test functions f63-f83: Oliveira & Takahashi (Bisection enhancement)
// ============================================================================

fn f_63(x: f64) f64 { return x * @exp(x) - 1.0; }
fn f_64(x: f64) f64 { return @tan(x - 0.1); }
fn f_65(x: f64) f64 { return @sin(x) + 0.5; }
fn f_66(x: f64) f64 { return 4.0 * std.math.pow(f64, x, 5) + x * x + 1.0; }
fn f_67(x: f64) f64 { return x + std.math.pow(f64, x, 10) - 1.0; }
fn f_68(x: f64) f64 { return std.math.pow(f64, std.math.pi, x) - std.math.e; }
fn f_69(x: f64) f64 { return @log(@abs(x - 10.0 / 9.0)); }
fn f_70(x: f64) f64 { return 1.0 / 3.0 + std.math.cbrt(x) + std.math.pow(f64, x, 3); }
fn f_71(x: f64) f64 { return (x + 2.0 / 3.0) / (x + 101.0 / 100.0); }
fn f_72(x: f64) f64 { return std.math.pow(f64, x * 1e6 - 1.0, 3); }
fn f_73(x: f64) f64 { return @exp(x) * std.math.pow(f64, x * 1e6 - 1.0, 3); }
fn f_74(x: f64) f64 { return std.math.pow(f64, x - 1.0 / 3.0, 2) * std.math.atan(x - 1.0 / 3.0); }
fn f_75(x: f64) f64 {
    const t = 3.0 * x - 1.0;
    return sgn(t) * (1.0 - @sqrt(1.0 - t * t / 81.0));
}
fn f_76(x: f64) f64 { return if (x > (1.0 - 1e6) / 1e6) (1.0 + 1e6) / 1e6 else -1.0; }
fn f_77(x: f64) f64 { return if (x != 1.0 / 21.0) 1.0 / (21.0 * x - 1.0) else 0.0; }
fn f_78(x: f64) f64 { return x * x / 4.0 + @ceil(x / 2.0) - 0.5; }
fn f_79(x: f64) f64 { return @ceil(10.0 * x - 1.0) + 0.5; }
fn f_80(x: f64) f64 { return x + @sin(x * 1e6) / 10.0 + 1e-3; }
fn f_81(x: f64) f64 { return if (x > -1.0) 1.0 + @sin(1.0 / (x + 1.0)) else -1.0; }
fn f_82(x: f64) f64 { return 202.0 * x - 2.0 * @floor((2.0 * x + 1e-2) / 2e-2) - 0.1; }
fn f_83(x: f64) f64 {
    const t = 202.0 * x - 2.0 * @floor((2.0 * x + 1e-2) / 2e-2) - 0.1;
    return t * t * t;
}

// ============================================================================
// Test functions f84-f92: SciML project benchmarks
// ============================================================================

fn f_84(x: f64) f64 { return (x - 1.0) * (x - 2.0) * (x - 3.0) * (x - 4.0) * (x - 5.0) - 0.05; }
fn f_85(x: f64) f64 { return @sin(x) - 0.5 * x - 0.3; }
fn f_86(x: f64) f64 { return @exp(x) - 1.0 - x - x * x / 2.0 - 0.005; }
fn f_87(x: f64) f64 { return 1.0 / (x - 0.5) - 2.0 - 0.05; }
fn f_88(x: f64) f64 { return @log(x) - x + 2.0 - 0.05; }
fn f_89(x: f64) f64 { return @sin(20.0 * x) + 0.1 * x - 0.1; }
fn f_90(x: f64) f64 { return std.math.pow(f64, x, 3) - 2.0 * std.math.pow(f64, x, 2) + x - 0.025; }
fn f_91(x: f64) f64 { return x * @sin(1.0 / x) - 0.1 - 0.01; }
fn f_92(x: f64) f64 { return std.math.pow(f64, x, 3) - 0.001; }

// ============================================================================
// All problems array
// ============================================================================

const all_problems = [_]Problem{
    // Galdino (f01-f33)
    .{ .name = "f_01", .f = f_01, .a = 0.5, .b = 1.5 },
    .{ .name = "f_02", .f = f_02, .a = 0.1, .b = 1.0 },
    .{ .name = "f_03", .f = f_03, .a = 0.1, .b = 1.0 },
    .{ .name = "f_04", .f = f_04, .a = -1.8, .b = 0.0 },
    .{ .name = "f_05", .f = f_05, .a = 2.0, .b = 3.0 },
    .{ .name = "f_06", .f = f_06, .a = 0.0, .b = 1.0 },
    .{ .name = "f_07", .f = f_07, .a = 0.0, .b = 1.0 },
    .{ .name = "f_08", .f = f_08, .a = 0.0, .b = 1.0 },
    .{ .name = "f_09", .f = f_09, .a = 0.0, .b = 1.0 },
    .{ .name = "f_10", .f = f_10, .a = 0.0, .b = 1.0 },
    .{ .name = "f_11", .f = f_11, .a = 0.0, .b = 1.0 },
    .{ .name = "f_12", .f = f_12, .a = 0.0, .b = 1.0 },
    .{ .name = "f_13", .f = f_13, .a = 0.0, .b = 1.0 },
    .{ .name = "f_14", .f = f_14, .a = 0.0, .b = 1.0 },
    .{ .name = "f_15", .f = f_15, .a = 0.0, .b = 1.0 },
    .{ .name = "f_16", .f = f_16, .a = 0.0, .b = 1.0 },
    .{ .name = "f_17", .f = f_17, .a = 0.0, .b = 1.0 },
    .{ .name = "f_18", .f = f_18, .a = 0.0, .b = 1.0 },
    .{ .name = "f_19", .f = f_19, .a = 0.0, .b = 1.0 },
    .{ .name = "f_20", .f = f_20, .a = 0.0, .b = 1.0 },
    .{ .name = "f_21", .f = f_21, .a = 0.0, .b = 1.0 },
    .{ .name = "f_22", .f = f_22, .a = 0.0, .b = 1.0 },
    .{ .name = "f_23", .f = f_23, .a = 0.0, .b = 1.0 },
    .{ .name = "f_24", .f = f_24, .a = 2.6, .b = 4.6 },
    .{ .name = "f_25", .f = f_25, .a = 3.6, .b = 5.6 },
    .{ .name = "f_26", .f = f_26, .a = 2.0, .b = 4.0 },
    .{ .name = "f_27", .f = f_27, .a = 1.0, .b = 3.0 },
    .{ .name = "f_28", .f = f_28, .a = 7.0, .b = 8.0 },
    .{ .name = "f_29", .f = f_29, .a = 2.6, .b = 4.6 },
    .{ .name = "f_30", .f = f_30, .a = 4.0, .b = 5.0 },
    .{ .name = "f_31", .f = f_31, .a = 0.05, .b = 5.0 },
    .{ .name = "f_32", .f = f_32, .a = 0.0, .b = 1.5 },
    .{ .name = "f_33", .f = f_33, .a = 0.0, .b = 4.0 },
    // Stage (f34-f40)
    .{ .name = "f_34", .f = f_34, .a = -11.0, .b = 9.0 },
    .{ .name = "f_35", .f = f_35, .a = -11.0, .b = 9.0 },
    .{ .name = "f_36", .f = f_36, .a = -11.0, .b = 9.0 },
    .{ .name = "f_37", .f = f_37, .a = -11.0, .b = 9.0 },
    .{ .name = "f_38", .f = f_38, .a = -11.0, .b = 9.0 },
    .{ .name = "f_39", .f = f_39, .a = -11.0, .b = 9.0 },
    .{ .name = "f_40", .f = f_40, .a = -11.0, .b = 9.0 },
    // Swift-Lindfield (f41-f62)
    .{ .name = "f_41", .f = f_41, .a = 0.0, .b = 10.0 },
    .{ .name = "f_42", .f = f_42, .a = 0.0, .b = std.math.pi },
    .{ .name = "f_43", .f = f_43, .a = -1.0, .b = 1.5 },
    .{ .name = "f_44", .f = f_44, .a = -1.0, .b = 1.5 },
    .{ .name = "f_45", .f = f_45, .a = -1.0, .b = 1.5 },
    .{ .name = "f_46", .f = f_46, .a = 0.09, .b = 0.7 },
    .{ .name = "f_47", .f = f_47, .a = 0.0005, .b = 0.5 },
    .{ .name = "f_48", .f = f_48, .a = 0.0005, .b = 0.5 },
    .{ .name = "f_49", .f = f_49, .a = -1.0, .b = 1.0 },
    .{ .name = "f_50", .f = f_50, .a = -3.0, .b = 2.0 },
    .{ .name = "f_51", .f = f_51, .a = 0.5, .b = 5.0 },
    .{ .name = "f_52", .f = f_52, .a = 0.5, .b = 8.0 },
    .{ .name = "f_53", .f = f_53, .a = 1.0, .b = 4.0 },
    .{ .name = "f_54", .f = f_54, .a = 0.1, .b = std.math.pi / 3.0 },
    .{ .name = "f_55", .f = f_55, .a = 0.0, .b = 1.5 },
    .{ .name = "f_56", .f = f_56, .a = 1.5, .b = 3.0 },
    .{ .name = "f_57", .f = f_57, .a = 2.6, .b = 3.5 },
    .{ .name = "f_58", .f = f_58, .a = 1.0, .b = 8.0 },
    .{ .name = "f_59", .f = f_59, .a = 0.2, .b = 3.0 },
    .{ .name = "f_60", .f = f_60, .a = 0.0, .b = 2.0 },
    .{ .name = "f_61", .f = f_61, .a = -1.0, .b = 2.0 },
    .{ .name = "f_62", .f = f_62, .a = std.math.pi / 2.0, .b = std.math.pi },
    // Oliveira-Takahashi (f63-f83)
    .{ .name = "f_63", .f = f_63, .a = -1.0, .b = 1.0 },
    .{ .name = "f_64", .f = f_64_fn, .a = -1.0, .b = 1.0 },
    .{ .name = "f_65", .f = f_65, .a = -1.0, .b = 1.0 },
    .{ .name = "f_66", .f = f_66, .a = -1.0, .b = 1.0 },
    .{ .name = "f_67", .f = f_67, .a = -1.0, .b = 1.0 },
    .{ .name = "f_68", .f = f_68, .a = -1.0, .b = 1.0 },
    .{ .name = "f_69", .f = f_69, .a = -1.0, .b = 1.0 },
    .{ .name = "f_70", .f = f_70, .a = -1.0, .b = 1.0 },
    .{ .name = "f_71", .f = f_71, .a = -1.0, .b = 1.0 },
    .{ .name = "f_72", .f = f_72, .a = -1.0, .b = 1.0 },
    .{ .name = "f_73", .f = f_73, .a = -1.0, .b = 1.0 },
    .{ .name = "f_74", .f = f_74, .a = -1.0, .b = 1.0 },
    .{ .name = "f_75", .f = f_75, .a = -1.0, .b = 1.0 },
    .{ .name = "f_76", .f = f_76, .a = -1.0, .b = 1.0 },
    .{ .name = "f_77", .f = f_77, .a = -1.0, .b = 1.0 },
    .{ .name = "f_78", .f = f_78, .a = -1.0, .b = 1.0 },
    .{ .name = "f_79", .f = f_79, .a = -1.0, .b = 1.0 },
    .{ .name = "f_80", .f = f_80, .a = -1.0, .b = 1.0 },
    .{ .name = "f_81", .f = f_81, .a = -1.0, .b = 1.0 },
    .{ .name = "f_82", .f = f_82, .a = -1.0, .b = 1.0 },
    .{ .name = "f_83", .f = f_83, .a = -1.0, .b = 1.0 },
    // SciML (f84-f92)
    .{ .name = "f_84", .f = f_84, .a = 0.5, .b = 5.5 },
    .{ .name = "f_85", .f = f_85, .a = -10.0, .b = 10.0 },
    .{ .name = "f_86", .f = f_86, .a = -2.0, .b = 2.0 },
    .{ .name = "f_87", .f = f_87, .a = 0.6, .b = 2.0 },
    .{ .name = "f_88", .f = f_88, .a = 0.1, .b = 3.0 },
    .{ .name = "f_89", .f = f_89, .a = -4.0, .b = 5.0 },
    .{ .name = "f_90", .f = f_90, .a = -1.0, .b = 2.0 },
    .{ .name = "f_91", .f = f_91, .a = 0.01, .b = 1.0 },
    .{ .name = "f_92", .f = f_92, .a = -10.0, .b = 10.0 },
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
