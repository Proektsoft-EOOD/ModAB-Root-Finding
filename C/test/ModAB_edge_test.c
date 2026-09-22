/*
 * Edge-case tests for the shared modAB safeguards.
 *
 * The 100-problem benchmark suite never produces a non-finite or overflowing
 * residual, so the overflow/NaN branches of same_nonzero_sign, safe_midpoint,
 * safe_secant, symmetry_factor and passes_switching_test are covered here.
 *
 * Expected values come from the C# reference in C#/Root/Node.cs and
 * C#/Root/Solvers/ModABCorr.cs; every language port is held to the same table.
 *
 * Build (from the C directory):
 *     gcc -O2 -Iinclude -Isrc -o edge_test test/ModAB_edge_test.c -lm
 */
#include <stdio.h>
#include <math.h>
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

int main(void) {
    double root;
    ckb("same_nonzero_sign", 0, same_nonzero_sign(1.0, 2.0), 1);
    ckb("same_nonzero_sign", 1, same_nonzero_sign(-1.0, -2.0), 1);
    ckb("same_nonzero_sign", 2, same_nonzero_sign(1.0, -2.0), 0);
    ckb("same_nonzero_sign", 3, same_nonzero_sign(-1.0, 2.0), 0);
    ckb("same_nonzero_sign", 4, same_nonzero_sign(0.0, 1.0), 0);
    ckb("same_nonzero_sign", 5, same_nonzero_sign(1.0, 0.0), 0);
    ckb("same_nonzero_sign", 6, same_nonzero_sign(0.0, 0.0), 0);
    ckb("same_nonzero_sign", 7, same_nonzero_sign(-0.0, -1.0), 0);
    ckb("same_nonzero_sign", 8, same_nonzero_sign(NAN, 1.0), 0);
    ckb("same_nonzero_sign", 9, same_nonzero_sign(1.0, NAN), 0);
    ckb("same_nonzero_sign", 10, same_nonzero_sign(NAN, NAN), 0);
    ckb("same_nonzero_sign", 11, same_nonzero_sign(INFINITY, 1.0), 1);
    ckb("same_nonzero_sign", 12, same_nonzero_sign(-INFINITY, -1.0), 1);
    ckb("same_nonzero_sign", 13, same_nonzero_sign(INFINITY, -INFINITY), 0);
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
    ck("symmetry_factor", 0, symmetry_factor(-1.0, 1.0), 1.0);
    ck("symmetry_factor", 1, symmetry_factor(-1.0, 3.0), 0.5625);
    ck("symmetry_factor", 2, symmetry_factor(-3.0, 1.0), 0.5625);
    ck("symmetry_factor", 3, symmetry_factor(-1e+308, 1e+308), 1.0);
    ck("symmetry_factor", 4, symmetry_factor(-1.7e+308, 1.0), 0.25);
    ck("symmetry_factor", 5, symmetry_factor(INFINITY, -1.0), NAN);
    ck("symmetry_factor", 6, symmetry_factor(-1.0, INFINITY), NAN);
    ck("symmetry_factor", 7, symmetry_factor(NAN, 1.0), NAN);
    ck("symmetry_factor", 8, symmetry_factor(-1e-300, 1e+300), 0.25);
    ckb("passes_switching_test", 0, passes_switching_test(1.0, 1.0, 0.5), 1);
    ckb("passes_switching_test", 1, passes_switching_test(1.0, -1.0, 0.5), 0);
    ckb("passes_switching_test", 2, passes_switching_test(1.0, 0.5, 0.5), 1);
    ckb("passes_switching_test", 3, passes_switching_test(1.0, 0.5, 0.1), 0);
    ckb("passes_switching_test", 4, passes_switching_test(1e+308, 1e+308, 0.5), 1);
    ckb("passes_switching_test", 5, passes_switching_test(1.7e+308, -1.7e+308, 0.5), 0);
    ckb("passes_switching_test", 6, passes_switching_test(1e+308, 1.6e+308, 0.5), 1);
    ckb("passes_switching_test", 7, passes_switching_test(INFINITY, 1.0, 0.5), 0);
    ckb("passes_switching_test", 8, passes_switching_test(1.0, INFINITY, 0.5), 0);
    ckb("passes_switching_test", 9, passes_switching_test(NAN, 1.0, 0.5), 0);
    ckb("passes_switching_test", 10, passes_switching_test(1.0, 1.0, NAN), 0);
    ckb("passes_switching_test", 11, passes_switching_test(0.0, 0.0, 0.5), 0);

    /* Solver-level regression: this returned +Inf before safe_midpoint. */
    root = modAB_find_root(f_overflow, 1e308, 1.7e308, 1e-14, 1e-14, 200);
    ++total;
    if (fabs(root - 1.2e308) <= 1e-13 * 1.2e308) ++passed;
    else printf("FAIL overflowing bracket: got %.17g want %.17g\n", root, 1.2e308);

    printf("C edge-cases: %d/%d %s\n", passed, total, passed == total ? "PASS" : "FAIL");
    return passed == total ? 0 : 1;
}
