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
 * Build and run (from the TypeScript directory):
 *     npx tsc && node dist/ModAB_edge_test.js
 */
import { modABRoot, sameSign, safeMidpoint, safeSecant } from "./ModAB";

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

// --- same_sign ---
ckb("same_sign", 0, sameSign(1.0, 2.0), true);
ckb("same_sign", 1, sameSign(-1.0, -2.0), true);
ckb("same_sign", 2, sameSign(1.0, -2.0), false);
ckb("same_sign", 3, sameSign(-1.0, 2.0), false);
ckb("same_sign", 4, sameSign(0.0, 1.0), false);
ckb("same_sign", 5, sameSign(1.0, 0.0), false);
ckb("same_sign", 6, sameSign(0.0, 0.0), false);
ckb("same_sign", 7, sameSign(-0.0, -1.0), false);
ckb("same_sign", 8, sameSign(NaN, 1.0), false);
ckb("same_sign", 9, sameSign(1.0, NaN), false);
ckb("same_sign", 10, sameSign(NaN, NaN), false);
ckb("same_sign", 11, sameSign(Infinity, 1.0), true);
ckb("same_sign", 12, sameSign(-Infinity, -1.0), true);
ckb("same_sign", 13, sameSign(Infinity, -Infinity), false);

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

// Solver-level regression: this returned Infinity before safeMidpoint.
total++;
const root = modABRoot((x) => x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
if (Math.abs(root - 1.2e308) <= 1e-13 * 1.2e308) { passed++; }
else { console.log(`FAIL overflowing bracket: got ${root} want 1.2e308`); }

// --- solver-level: overflow and infinite residuals in the switching test ---
function ckRoot(name: string, fn: (x: number) => number, a: number, b: number, want: number): void {
    total++;
    const got = modABRoot(fn, a, b, 0.0, 1e-14, 0.0, 200);
    if (Math.abs(got - want) <= 1e-13 * Math.max(Math.abs(want), 1)) { passed++; return; }
    console.log(`FAIL ${name}: got ${got} want ${want}`);
}

const M = Number.MAX_VALUE;
// |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite
ckRoot("hump", (x) => x <= 0 ? M * (-0.2 + 1.2 * (x + 1)) : M * (1 - 0.2 * x), -1.0, 1.0, -5.0 / 6.0);
// f2 - f1 overflows: switching is disabled until the residuals shrink
ckRoot("f2-f1 overflow", (x) => 1.7e308 * Math.tanh(10 * (x - 0.3)), -1.0, 1.0, 0.3);
// infinite residuals at one or both ends of the bracket
ckRoot("inf left", (x) => x < -0.5 ? -Infinity : x - 0.1, -1.0, 1.0, 0.1);
ckRoot("inf both", (x) => x < -0.5 ? -Infinity : (x > 0.9 ? Infinity : x - 0.1), -1.0, 1.0, 0.1);
// subnormal residuals: AB corrections underflow towards zero
ckRoot("subnormal", (x) => 1e-300 * (x * x * x - 0.2), -1.0, 1.0, Math.cbrt(0.2));

console.log(`TypeScript edge-cases: ${passed}/${total} ${passed === total ? "PASS" : "FAIL"}`);
if (passed !== total) { throw new Error("edge-case tests failed"); }
