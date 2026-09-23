/*
 * Edge-case tests for the shared modAB safeguards.
 *
 * The 100-problem benchmark suite never produces a non-finite or overflowing
 * residual, so the overflow/NaN branches of same_sign, safe_midpoint and
 * safe_secant, and the overflow cases of the switching test, are covered here.
 *
 * Expected values come from the C# reference in C#/Root/Node.cs and
 * C#/Root/Solvers/SgModAB.cs; every language port is held to the same table.
 *
 * Build (from the C directory):
 *     gcc -O2 -Iinclude -Isrc -o edge_test test/ModAB_edge_test.c -lm
 */
#include <stdio.h>
#include <math.h>
#include <float.h>
/* The safeguards are static inline, so the solver is pulled in directly. */
#include "../src/ModAB.c"

static int passed = 0, total = 0;

static void ck(const char *n, int i, double got, double want) {
    ++total;
    if ((isnan(got) && isnan(want)) || got == want) { ++passed; return; }
    printf("FAIL %s[%d]: got %.17g want %.17g\n", n, i, got, want);
}

static void ckb(const char *n, int i, int got, int want) {
    ++total;
    if (!got == !want) { ++passed; return; }
    printf("FAIL %s[%d]: got %d want %d\n", n, i, !!got, !!want);
}

/* f(x) = x*1e-308 - 1.2; the root is 1.2e308 and x1 + x2 overflows. */
static double f_overflow(double x) { return x * 1e-308 - 1.2; }

/* Solver-level cases for overflow and infinite residuals in the switching test. */
/* |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite. */
static double f_hump(double x) {
    return x <= 0 ? DBL_MAX * (-0.2 + 1.2 * (x + 1)) : DBL_MAX * (1 - 0.2 * x);
}
/* f2 - f1 overflows: switching is disabled until the residuals shrink. */
static double f_diff_overflow(double x) { return 1.7e308 * tanh(10 * (x - 0.3)); }
/* Infinite residuals at one or both ends of the bracket. */
static double f_inf_left(double x) { return x < -0.5 ? -INFINITY : x - 0.1; }
static double f_inf_both(double x) {
    return x < -0.5 ? -INFINITY : (x > 0.9 ? INFINITY : x - 0.1);
}
/* Subnormal residuals: AB corrections underflow towards zero. */
static double f_subnormal(double x) { return 1e-300 * (x * x * x - 0.2); }

static void ck_root(const char *n, double (*f)(double), double a, double b, double want) {
    double got = modAB_find_root(f, a, b, 1e-14, 1e-14, 200);
    ++total;
    if (fabs(got - want) <= 1e-13 * fmax(fabs(want), 1.0)) { ++passed; return; }
    printf("FAIL %s: got %.17g want %.17g\n", n, got, want);
}

int main(void) {
    double root;
    ckb("same_sign", 0, same_sign(1.0, 2.0), 1);
    ckb("same_sign", 1, same_sign(-1.0, -2.0), 1);
    ckb("same_sign", 2, same_sign(1.0, -2.0), 0);
    ckb("same_sign", 3, same_sign(-1.0, 2.0), 0);
    ckb("same_sign", 4, same_sign(0.0, 1.0), 0);
    ckb("same_sign", 5, same_sign(1.0, 0.0), 0);
    ckb("same_sign", 6, same_sign(0.0, 0.0), 0);
    ckb("same_sign", 7, same_sign(-0.0, -1.0), 0);
    ckb("same_sign", 8, same_sign(NAN, 1.0), 0);
    ckb("same_sign", 9, same_sign(1.0, NAN), 0);
    ckb("same_sign", 10, same_sign(NAN, NAN), 0);
    ckb("same_sign", 11, same_sign(INFINITY, 1.0), 1);
    ckb("same_sign", 12, same_sign(-INFINITY, -1.0), 1);
    ckb("same_sign", 13, same_sign(INFINITY, -INFINITY), 0);
    ck("safe_midpoint", 0, safe_midpoint(2.0, 4.0), 3.0);
    ck("safe_midpoint", 1, safe_midpoint(1e+308, 1e+308), 1e+308);
    ck("safe_midpoint", 2, safe_midpoint(-1e+308, 1e+308), 0.0);
    ck("safe_midpoint", 3, safe_midpoint(1e+308, 1.7e+308), 1.35e+308);
    ck("safe_midpoint", 4, safe_midpoint(-1.7e+308, -1e+308), -1.35e+308);
    ck("safe_midpoint", 5, safe_midpoint(0.0, 1.0), 0.5);
    ck("safe_midpoint", 6, safe_midpoint(-1.0, 1.0), 0.0);
    ck("safe_secant", 0, safe_secant(0.0, -1.0, 1.0, 1.0), 0.5);
    ck("safe_secant", 1, safe_secant(0.0, -1.0, 1.0, 3.0), 0.25);
    ck("safe_secant", 2, safe_secant(0.0, -1e+308, 1.0, 1e+308), 0.5);
    ck("safe_secant", 3, safe_secant(0.0, -1.7e+308, 1.0, 1.7e+308), 0.5);
    ck("safe_secant", 4, safe_secant(0.0, INFINITY, 1.0, -1.0), 0.5);
    ck("safe_secant", 5, safe_secant(0.0, -1.0, 1.0, INFINITY), 0.5);
    ck("safe_secant", 6, safe_secant(0.0, 0.0, 1.0, 0.0), 0.5);
    ck("safe_secant", 7, safe_secant(0.0, NAN, 1.0, 1.0), 0.5);
    ck("safe_secant", 8, safe_secant(0.0, -1e-300, 1.0, 1e+300), 0.0);
    ck("safe_secant", 9, safe_secant(0.0, -1e+300, 1.0, 1e-300), 1.0);
    ck("safe_secant", 10, safe_secant(1e+308, -1.0, 1.7e+308, 1.0), 1.35e+308);

    /* Solver-level regression: this returned +Inf before safe_midpoint. */
    root = modAB_find_root(f_overflow, 1e308, 1.7e308, 1e-14, 1e-14, 200);
    ++total;
    if (fabs(root - 1.2e308) <= 1e-13 * 1.2e308) ++passed;
    else printf("FAIL overflowing bracket: got %.17g want %.17g\n", root, 1.2e308);

    ck_root("hump", f_hump, -1.0, 1.0, -5.0 / 6.0);
    ck_root("f2-f1 overflow", f_diff_overflow, -1.0, 1.0, 0.3);
    ck_root("inf left", f_inf_left, -1.0, 1.0, 0.1);
    ck_root("inf both", f_inf_both, -1.0, 1.0, 0.1);
    ck_root("subnormal", f_subnormal, -1.0, 1.0, cbrt(0.2));

    printf("C edge-cases: %d/%d %s\n", passed, total, passed == total ? "PASS" : "FAIL");
    return passed == total ? 0 : 1;
}
