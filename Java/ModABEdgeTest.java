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

    public static void main(String[] args) {
        // --- same_nonzero_sign ---
        ckb("same_nonzero_sign", 0, ModAB.sameNonzeroSign(1.0, 2.0), true);
        ckb("same_nonzero_sign", 1, ModAB.sameNonzeroSign(-1.0, -2.0), true);
        ckb("same_nonzero_sign", 2, ModAB.sameNonzeroSign(1.0, -2.0), false);
        ckb("same_nonzero_sign", 3, ModAB.sameNonzeroSign(-1.0, 2.0), false);
        ckb("same_nonzero_sign", 4, ModAB.sameNonzeroSign(0.0, 1.0), false);
        ckb("same_nonzero_sign", 5, ModAB.sameNonzeroSign(1.0, 0.0), false);
        ckb("same_nonzero_sign", 6, ModAB.sameNonzeroSign(0.0, 0.0), false);
        ckb("same_nonzero_sign", 7, ModAB.sameNonzeroSign(-0.0, -1.0), false);
        ckb("same_nonzero_sign", 8, ModAB.sameNonzeroSign(Double.NaN, 1.0), false);
        ckb("same_nonzero_sign", 9, ModAB.sameNonzeroSign(1.0, Double.NaN), false);
        ckb("same_nonzero_sign", 10, ModAB.sameNonzeroSign(Double.NaN, Double.NaN), false);
        ckb("same_nonzero_sign", 11, ModAB.sameNonzeroSign(Double.POSITIVE_INFINITY, 1.0), true);
        ckb("same_nonzero_sign", 12, ModAB.sameNonzeroSign(Double.NEGATIVE_INFINITY, -1.0), true);
        ckb("same_nonzero_sign", 13, ModAB.sameNonzeroSign(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY), false);

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

        // --- symmetry_factor ---
        ck("symmetry_factor", 0, ModAB.symmetryFactor(-1.0, 1.0), 1.0);
        ck("symmetry_factor", 1, ModAB.symmetryFactor(-1.0, 3.0), 0.5625);
        ck("symmetry_factor", 2, ModAB.symmetryFactor(-3.0, 1.0), 0.5625);
        ck("symmetry_factor", 3, ModAB.symmetryFactor(-1e+308, 1e+308), 1.0);
        ck("symmetry_factor", 4, ModAB.symmetryFactor(-1.7e+308, 1.0), 0.25);
        ck("symmetry_factor", 5, ModAB.symmetryFactor(Double.POSITIVE_INFINITY, -1.0), Double.NaN);
        ck("symmetry_factor", 6, ModAB.symmetryFactor(-1.0, Double.POSITIVE_INFINITY), Double.NaN);
        ck("symmetry_factor", 7, ModAB.symmetryFactor(Double.NaN, 1.0), Double.NaN);
        ck("symmetry_factor", 8, ModAB.symmetryFactor(-1e-300, 1e+300), 0.25);

        // --- passes_switching_test ---
        ckb("passes_switching_test", 0, ModAB.passesSwitchingTest(1.0, 1.0, 0.5), true);
        ckb("passes_switching_test", 1, ModAB.passesSwitchingTest(1.0, -1.0, 0.5), false);
        ckb("passes_switching_test", 2, ModAB.passesSwitchingTest(1.0, 0.5, 0.5), true);
        ckb("passes_switching_test", 3, ModAB.passesSwitchingTest(1.0, 0.5, 0.1), false);
        ckb("passes_switching_test", 4, ModAB.passesSwitchingTest(1e+308, 1e+308, 0.5), true);
        ckb("passes_switching_test", 5, ModAB.passesSwitchingTest(1.7e+308, -1.7e+308, 0.5), false);
        ckb("passes_switching_test", 6, ModAB.passesSwitchingTest(1e+308, 1.6e+308, 0.5), true);
        ckb("passes_switching_test", 7, ModAB.passesSwitchingTest(Double.POSITIVE_INFINITY, 1.0, 0.5), false);
        ckb("passes_switching_test", 8, ModAB.passesSwitchingTest(1.0, Double.POSITIVE_INFINITY, 0.5), false);
        ckb("passes_switching_test", 9, ModAB.passesSwitchingTest(Double.NaN, 1.0, 0.5), false);
        ckb("passes_switching_test", 10, ModAB.passesSwitchingTest(1.0, 1.0, Double.NaN), false);
        ckb("passes_switching_test", 11, ModAB.passesSwitchingTest(0.0, 0.0, 0.5), false);

        // Solver-level regression: this returned Infinity before safeMidpoint.
        total++;
        double root = ModAB.modABRoot(x -> x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
        if (Math.abs(root - 1.2e308) <= 1e-13 * 1.2e308) passed++;
        else System.out.println("FAIL overflowing bracket: got " + root + " want 1.2e308");

        System.out.println("Java edge-cases: " + passed + "/" + total
                + " " + (passed == total ? "PASS" : "FAIL"));
        if (passed != total) System.exit(1);
    }
}
