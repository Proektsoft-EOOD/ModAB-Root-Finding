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
 * Build and run (from the TypeScript directory):
 *     npx tsc && node dist/ModAB_edge_test.js
 */
import {
    modABRoot, sameNonzeroSign, safeMidpoint, safeSecant,
    symmetryFactor, passesSwitchingTest,
} from "./ModAB";

let passed = 0;
let total = 0;

function ck(name: string, i: number, got: number, want: number): void {
    total++;
    if ((Number.isNaN(got) && Number.isNaN(want)) || got === want) { passed++; return; }
    console.log(`FAIL ${name}[${i}]: got ${got} want ${want}`);
}

function ckb(name: string, i: number, got: boolean, want: boolean): void {
    total++;
    if (got === want) { passed++; return; }
    console.log(`FAIL ${name}[${i}]: got ${got} want ${want}`);
}

// --- same_nonzero_sign ---
ckb("same_nonzero_sign", 0, sameNonzeroSign(1.0, 2.0), true);
ckb("same_nonzero_sign", 1, sameNonzeroSign(-1.0, -2.0), true);
ckb("same_nonzero_sign", 2, sameNonzeroSign(1.0, -2.0), false);
ckb("same_nonzero_sign", 3, sameNonzeroSign(-1.0, 2.0), false);
ckb("same_nonzero_sign", 4, sameNonzeroSign(0.0, 1.0), false);
ckb("same_nonzero_sign", 5, sameNonzeroSign(1.0, 0.0), false);
ckb("same_nonzero_sign", 6, sameNonzeroSign(0.0, 0.0), false);
ckb("same_nonzero_sign", 7, sameNonzeroSign(-0.0, -1.0), false);
ckb("same_nonzero_sign", 8, sameNonzeroSign(NaN, 1.0), false);
ckb("same_nonzero_sign", 9, sameNonzeroSign(1.0, NaN), false);
ckb("same_nonzero_sign", 10, sameNonzeroSign(NaN, NaN), false);
ckb("same_nonzero_sign", 11, sameNonzeroSign(Infinity, 1.0), true);
ckb("same_nonzero_sign", 12, sameNonzeroSign(-Infinity, -1.0), true);
ckb("same_nonzero_sign", 13, sameNonzeroSign(Infinity, -Infinity), false);

// --- safe_midpoint ---
ck("safe_midpoint", 0, safeMidpoint(2.0, 4.0), 3.0);
ck("safe_midpoint", 1, safeMidpoint(1e+308, 1e+308), 1e+308);
ck("safe_midpoint", 2, safeMidpoint(-1e+308, 1e+308), 0.0);
ck("safe_midpoint", 3, safeMidpoint(1e+308, 1.7e+308), 1.35e+308);
ck("safe_midpoint", 4, safeMidpoint(-1.7e+308, -1e+308), -1.35e+308);
ck("safe_midpoint", 5, safeMidpoint(0.0, 1.0), 0.5);
ck("safe_midpoint", 6, safeMidpoint(-1.0, 1.0), 0.0);

// --- safe_secant ---
ck("safe_secant", 0, safeSecant(0.0, -1.0, 1.0, 1.0), 0.5);
ck("safe_secant", 1, safeSecant(0.0, -1.0, 1.0, 3.0), 0.25);
ck("safe_secant", 2, safeSecant(0.0, -1e+308, 1.0, 1e+308), 0.5);
ck("safe_secant", 3, safeSecant(0.0, -1.7e+308, 1.0, 1.7e+308), 0.5);
ck("safe_secant", 4, safeSecant(0.0, Infinity, 1.0, -1.0), 0.5);
ck("safe_secant", 5, safeSecant(0.0, -1.0, 1.0, Infinity), 0.5);
ck("safe_secant", 6, safeSecant(0.0, 0.0, 1.0, 0.0), 0.5);
ck("safe_secant", 7, safeSecant(0.0, NaN, 1.0, 1.0), 0.5);
ck("safe_secant", 8, safeSecant(0.0, -1e-300, 1.0, 1e+300), 0.0);
ck("safe_secant", 9, safeSecant(0.0, -1e+300, 1.0, 1e-300), 1.0);
ck("safe_secant", 10, safeSecant(1e+308, -1.0, 1.7e+308, 1.0), 1.35e+308);

// --- symmetry_factor ---
ck("symmetry_factor", 0, symmetryFactor(-1.0, 1.0), 1.0);
ck("symmetry_factor", 1, symmetryFactor(-1.0, 3.0), 0.5625);
ck("symmetry_factor", 2, symmetryFactor(-3.0, 1.0), 0.5625);
ck("symmetry_factor", 3, symmetryFactor(-1e+308, 1e+308), 1.0);
ck("symmetry_factor", 4, symmetryFactor(-1.7e+308, 1.0), 0.25);
ck("symmetry_factor", 5, symmetryFactor(Infinity, -1.0), NaN);
ck("symmetry_factor", 6, symmetryFactor(-1.0, Infinity), NaN);
ck("symmetry_factor", 7, symmetryFactor(NaN, 1.0), NaN);
ck("symmetry_factor", 8, symmetryFactor(-1e-300, 1e+300), 0.25);

// --- passes_switching_test ---
ckb("passes_switching_test", 0, passesSwitchingTest(1.0, 1.0, 0.5), true);
ckb("passes_switching_test", 1, passesSwitchingTest(1.0, -1.0, 0.5), false);
ckb("passes_switching_test", 2, passesSwitchingTest(1.0, 0.5, 0.5), true);
ckb("passes_switching_test", 3, passesSwitchingTest(1.0, 0.5, 0.1), false);
ckb("passes_switching_test", 4, passesSwitchingTest(1e+308, 1e+308, 0.5), true);
ckb("passes_switching_test", 5, passesSwitchingTest(1.7e+308, -1.7e+308, 0.5), false);
ckb("passes_switching_test", 6, passesSwitchingTest(1e+308, 1.6e+308, 0.5), true);
ckb("passes_switching_test", 7, passesSwitchingTest(Infinity, 1.0, 0.5), false);
ckb("passes_switching_test", 8, passesSwitchingTest(1.0, Infinity, 0.5), false);
ckb("passes_switching_test", 9, passesSwitchingTest(NaN, 1.0, 0.5), false);
ckb("passes_switching_test", 10, passesSwitchingTest(1.0, 1.0, NaN), false);
ckb("passes_switching_test", 11, passesSwitchingTest(0.0, 0.0, 0.5), false);

// Solver-level regression: this returned Infinity before safeMidpoint.
total++;
const root = modABRoot((x) => x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
if (Math.abs(root - 1.2e308) <= 1e-13 * 1.2e308) { passed++; }
else { console.log(`FAIL overflowing bracket: got ${root} want 1.2e308`); }

console.log(`TypeScript edge-cases: ${passed}/${total} ${passed === total ? "PASS" : "FAIL"}`);
if (passed !== total) { throw new Error("edge-case tests failed"); }
