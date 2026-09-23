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
 * Build and run (from the Java directory):
 *     javac ModAB.java ModABEdgeTest.java && java ModABEdgeTest
 */
public class ModABEdgeTest {
    static int passed = 0, total = 0;

    static void ck(String name, int i, double got, double want) {
        total++;
        if ((Double.isNaN(got) && Double.isNaN(want)) || got == want) { passed++; return; }
        System.out.println("FAIL " + name + "[" + i + "]: got " + got + " want " + want);
    }

    static void ckb(String name, int i, boolean got, boolean want) {
        total++;
        if (got == want) { passed++; return; }
        System.out.println("FAIL " + name + "[" + i + "]: got " + got + " want " + want);
    }

    static void ckRoot(String name, java.util.function.DoubleUnaryOperator fn, double a, double b, double want) {
        total++;
        double got = ModAB.modABRoot(fn, a, b, 0.0, 1e-14, 0.0, 200);
        if (Math.abs(got - want) <= 1e-13 * Math.max(Math.abs(want), 1.0)) { passed++; return; }
        System.out.println("FAIL " + name + ": got " + got + " want " + want);
    }

    public static void main(String[] args) {
        // --- same_sign ---
        ckb("same_sign", 0, ModAB.sameSign(1.0, 2.0), true);
        ckb("same_sign", 1, ModAB.sameSign(-1.0, -2.0), true);
        ckb("same_sign", 2, ModAB.sameSign(1.0, -2.0), false);
        ckb("same_sign", 3, ModAB.sameSign(-1.0, 2.0), false);
        ckb("same_sign", 4, ModAB.sameSign(0.0, 1.0), false);
        ckb("same_sign", 5, ModAB.sameSign(1.0, 0.0), false);
        ckb("same_sign", 6, ModAB.sameSign(0.0, 0.0), false);
        ckb("same_sign", 7, ModAB.sameSign(-0.0, -1.0), false);
        ckb("same_sign", 8, ModAB.sameSign(Double.NaN, 1.0), false);
        ckb("same_sign", 9, ModAB.sameSign(1.0, Double.NaN), false);
        ckb("same_sign", 10, ModAB.sameSign(Double.NaN, Double.NaN), false);
        ckb("same_sign", 11, ModAB.sameSign(Double.POSITIVE_INFINITY, 1.0), true);
        ckb("same_sign", 12, ModAB.sameSign(Double.NEGATIVE_INFINITY, -1.0), true);
        ckb("same_sign", 13, ModAB.sameSign(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), false);

        // --- safe_midpoint ---
        ck("safe_midpoint", 0, ModAB.safeMidpoint(2.0, 4.0), 3.0);
        ck("safe_midpoint", 1, ModAB.safeMidpoint(1e+308, 1e+308), 1e+308);
        ck("safe_midpoint", 2, ModAB.safeMidpoint(-1e+308, 1e+308), 0.0);
        ck("safe_midpoint", 3, ModAB.safeMidpoint(1e+308, 1.7e+308), 1.35e+308);
        ck("safe_midpoint", 4, ModAB.safeMidpoint(-1.7e+308, -1e+308), -1.35e+308);
        ck("safe_midpoint", 5, ModAB.safeMidpoint(0.0, 1.0), 0.5);
        ck("safe_midpoint", 6, ModAB.safeMidpoint(-1.0, 1.0), 0.0);

        // --- safe_secant ---
        ck("safe_secant", 0, ModAB.safeSecant(0.0, -1.0, 1.0, 1.0), 0.5);
        ck("safe_secant", 1, ModAB.safeSecant(0.0, -1.0, 1.0, 3.0), 0.25);
        ck("safe_secant", 2, ModAB.safeSecant(0.0, -1e+308, 1.0, 1e+308), 0.5);
        ck("safe_secant", 3, ModAB.safeSecant(0.0, -1.7e+308, 1.0, 1.7e+308), 0.5);
        ck("safe_secant", 4, ModAB.safeSecant(0.0, Double.POSITIVE_INFINITY, 1.0, -1.0), 0.5);
        ck("safe_secant", 5, ModAB.safeSecant(0.0, -1.0, 1.0, Double.POSITIVE_INFINITY), 0.5);
        ck("safe_secant", 6, ModAB.safeSecant(0.0, 0.0, 1.0, 0.0), 0.5);
        ck("safe_secant", 7, ModAB.safeSecant(0.0, Double.NaN, 1.0, 1.0), 0.5);
        ck("safe_secant", 8, ModAB.safeSecant(0.0, -1e-300, 1.0, 1e+300), 0.0);
        ck("safe_secant", 9, ModAB.safeSecant(0.0, -1e+300, 1.0, 1e-300), 1.0);
        ck("safe_secant", 10, ModAB.safeSecant(1e+308, -1.0, 1.7e+308, 1.0), 1.35e+308);

        // Solver-level regression: this returned Infinity before safeMidpoint.
        total++;
        double root = ModAB.modABRoot(x -> x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
        if (Math.abs(root - 1.2e308) <= 1e-13 * 1.2e308) passed++;
        else System.out.println("FAIL overflowing bracket: got " + root + " want 1.2e308");

        // Solver-level cases for overflow and infinite residuals in the switching test.
        final double M = Double.MAX_VALUE;
        // |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite
        ckRoot("hump", x -> x <= 0 ? M * (-0.2 + 1.2 * (x + 1)) : M * (1 - 0.2 * x), -1.0, 1.0, -5.0 / 6.0);
        // f2 - f1 overflows: switching is disabled until the residuals shrink
        ckRoot("f2-f1 overflow", x -> 1.7e308 * Math.tanh(10 * (x - 0.3)), -1.0, 1.0, 0.3);
        // infinite residuals at one or both ends of the bracket
        ckRoot("inf left", x -> x < -0.5 ? Double.NEGATIVE_INFINITY : x - 0.1, -1.0, 1.0, 0.1);
        ckRoot("inf both", x -> x < -0.5 ? Double.NEGATIVE_INFINITY : (x > 0.9 ? Double.POSITIVE_INFINITY : x - 0.1),
                -1.0, 1.0, 0.1);
        // subnormal residuals: AB corrections underflow towards zero
        ckRoot("subnormal", x -> 1e-300 * (x * x * x - 0.2), -1.0, 1.0, Math.cbrt(0.2));

        System.out.println("Java edge-cases: " + passed + "/" + total
                + " " + (passed == total ? "PASS" : "FAIL"));
        if (passed != total) System.exit(1);
    }
}
