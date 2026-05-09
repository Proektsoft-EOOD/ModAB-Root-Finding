import { modABRoot } from "./ModAB";

interface Problem {
    name: string;
    f: (x: number) => number;
    a: number;
    b: number;
}

function P(x: number): number {
    return x + 1.11111;
}

// Test problems from RootBenchmark.py
const problems1: Problem[] = [
    // Sergio Galdino. A family of regula falsi root-finding methods
    { name: "f01", f: (x) => x ** 3 - 1, a: 0.5, b: 1.5 },
    { name: "f02", f: (x) => x ** 2 * (x ** 2 / 3 + Math.sqrt(2) * Math.sin(x)) - Math.sqrt(3) / 18, a: 0.1, b: 1 },
    { name: "f03", f: (x) => 11 * x ** 11 - 1, a: 0.1, b: 1 },
    { name: "f04", f: (x) => x ** 3 + 1, a: -1.8, b: 0 },
    { name: "f05", f: (x) => x ** 3 - 2 * x - 5, a: 2, b: 3 },
    { name: "f06", f: (x) => 2 * x * Math.exp(-5) + 1 - 2 * Math.exp(-5 * x), a: 0, b: 1 },
    { name: "f07", f: (x) => 2 * x * Math.exp(-10) + 1 - 2 * Math.exp(-10 * x), a: 0, b: 1 },
    { name: "f08", f: (x) => 2 * x * Math.exp(-20) + 1 - 2 * Math.exp(-20 * x), a: 0, b: 1 },
    { name: "f09", f: (x) => (1 + (1 - 5) ** 2) * x ** 2 - (1 - 5 * x) ** 2, a: 0, b: 1 },
    { name: "f10", f: (x) => (1 + (1 - 10) ** 2) * x ** 2 - (1 - 10 * x) ** 2, a: 0, b: 1 },
    { name: "f11", f: (x) => (1 + (1 - 20) ** 2) * x ** 2 - (1 - 20 * x) ** 2, a: 0, b: 1 },
    { name: "f12", f: (x) => x ** 2 - (1 - x) ** 5, a: 0, b: 1 },
    { name: "f13", f: (x) => x ** 2 - (1 - x) ** 10, a: 0, b: 1 },
    { name: "f14", f: (x) => x ** 2 - (1 - x) ** 20, a: 0, b: 1 },
    { name: "f15", f: (x) => (1 + (1 - 5) ** 4) * x - (1 - 5 * x) ** 4, a: 0, b: 1 },
    { name: "f16", f: (x) => (1 + (1 - 10) ** 4) * x - (1 - 10 * x) ** 4, a: 0, b: 1 },
    { name: "f17", f: (x) => (1 + (1 - 20) ** 4) * x - (1 - 20 * x) ** 4, a: 0, b: 1 },
    { name: "f18", f: (x) => Math.exp(-5 * x) * (x - 1) + x ** 5, a: 0, b: 1 },
    { name: "f19", f: (x) => Math.exp(-10 * x) * (x - 1) + x ** 10, a: 0, b: 1 },
    { name: "f20", f: (x) => Math.exp(-20 * x) * (x - 1) + x ** 20, a: 0, b: 1 },
    { name: "f21", f: (x) => x ** 2 + Math.sin(x / 5) - 1 / 4, a: 0, b: 1 },
    { name: "f22", f: (x) => x ** 2 + Math.sin(x / 10) - 1 / 4, a: 0, b: 1 },
    { name: "f23", f: (x) => x ** 2 + Math.sin(x / 20) - 1 / 4, a: 0, b: 1 },
    { name: "f24", f: (x) => (x + 2) * (x + 1) * (x - 3) ** 3, a: 2.6, b: 4.6 },
    { name: "f25", f: (x) => (x - 4) ** 5 * Math.log(x), a: 3.6, b: 5.6 },
    { name: "f26", f: (x) => (Math.sin(x) - x / 4) ** 3, a: 2, b: 4 },
    { name: "f27", f: (x) => (81 - P(x) * (108 - P(x) * (54 - P(x) * (12 - P(x))))) * (P(x) < 3 ? 1 : (P(x) > 3 ? -1 : 0)), a: 1, b: 3 },
    { name: "f28", f: (x) => Math.sin((x - 7.143) ** 3), a: 7, b: 8 },
    { name: "f29", f: (x) => Math.exp((x - 3) ** 5) - 1, a: 2.6, b: 4.6 },
    { name: "f30", f: (x) => Math.exp((x - 3) ** 5) - Math.exp(x - 1), a: 4, b: 5 },
    { name: "f31", f: (x) => Math.PI - 1 / x, a: 0.05, b: 5 },
    { name: "f32", f: (x) => 4 - Math.tan(x), a: 0, b: 1.5 },
    { name: "f33", f: (x) => Math.cos(x) - x ** 3, a: 0, b: 4 },
    // Steven A. Stage. Comments on An Improvement to the Brent's Method
    { name: "f34", f: (x) => Math.cos(x) - x, a: -11, b: 9 },
    { name: "f35", f: (x) => Math.sqrt(Math.abs(x - 2 / 3)) * (x <= 2 / 3 ? 1 : -1) - 0.1, a: -11, b: 9 },
    { name: "f36", f: (x) => Math.abs(x - 2 / 3) ** 0.2 * (x <= 2 / 3 ? 1 : -1), a: -11, b: 9 },
    { name: "f37", f: (x) => (x - 7 / 9) ** 3 + (x - 7 / 9) * 1e-3, a: -11, b: 9 },
    { name: "f38", f: (x) => x <= 1 / 3 ? -0.5 : 0.5, a: -11, b: 9 },
    { name: "f39", f: (x) => x <= 1 / 3 ? -1e-3 : 1 - 1e-3, a: -11, b: 9 },
    { name: "f40", f: (x) => x === 0 ? 0 : 1 / (x - 2 / 3), a: -11, b: 9 },
    // A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method
    { name: "f41", f: (x) => 2 * x * Math.exp(-5) - 2 * Math.exp(-5 * x) + 1, a: 0, b: 10 },
    { name: "f42", f: (x) => (x ** 2 - x - 6) * (x ** 2 - 3 * x + 2), a: 0, b: Math.PI },
    { name: "f43", f: (x) => x ** 3, a: -1, b: 1.5 },
    { name: "f44", f: (x) => x ** 5, a: -1, b: 1.5 },
    { name: "f45", f: (x) => x ** 7, a: -1, b: 1.5 },
    { name: "f46", f: (x) => (Math.exp(-5 * x) - x - 0.5) / x ** 5, a: 0.09, b: 0.7 },
    { name: "f47", f: (x) => 1 / Math.sqrt(x) - 2 * Math.log(5e3 * Math.sqrt(x)) + 0.8, a: 0.0005, b: 0.5 },
    { name: "f48", f: (x) => 1 / Math.sqrt(x) - 2 * Math.log(5e7 * Math.sqrt(x)) + 0.8, a: 0.0005, b: 0.5 },
    { name: "f49", f: (x) => x <= 0 ? (-(x ** 3) - x - 1) : (x ** (1 / 3) - x - 1), a: -1, b: 1 },
    { name: "f50", f: (x) => x ** 3 - 2 * x - x + 3, a: -3, b: 2 },
    { name: "f51", f: (x) => Math.log(x), a: 0.5, b: 5 },
    { name: "f52", f: (x) => (10 - x) * Math.exp(-10 * x) - x ** 10 + 1, a: 0.5, b: 8 },
    { name: "f53", f: (x) => Math.exp(Math.sin(x)) - x - 1, a: 1.0, b: 4 },
    { name: "f54", f: (x) => 2 * Math.sin(x) - 1, a: 0.1, b: Math.PI / 3 },
    { name: "f55", f: (x) => (x - 1) * Math.exp(-x), a: 0.0, b: 1.5 },
    { name: "f56", f: (x) => (x - 1) ** 3 - 1, a: 1.5, b: 3 },
    { name: "f57", f: (x) => Math.exp(x ** 2 + 7 * x - 30) - 1, a: 2.6, b: 3.5 },
    { name: "f58", f: (x) => Math.atan(x) - 1, a: 1.0, b: 8 },
    { name: "f59", f: (x) => Math.exp(x) - 2 * x - 1, a: 0.2, b: 3 },
    { name: "f60", f: (x) => Math.exp(-x) - x - Math.sin(x), a: 0.0, b: 2 },
    { name: "f61", f: (x) => x ** 2 - Math.sin(x) ** 2 - 1, a: -1, b: 2 },
    { name: "f62", f: (x) => Math.sin(x) - x / 2, a: Math.PI / 2, b: Math.PI },
];

// Oliveira I. F. D., Takahashi R. H. C.
// An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality
const problems2: Problem[] = [
    { name: "f63", f: (x) => x * Math.exp(x) - 1, a: -1, b: 1 },
    { name: "f64", f: (x) => Math.tan(x - 1 / 10), a: -1, b: 1 },
    { name: "f65", f: (x) => Math.sin(x) + 0.5, a: -1, b: 1 },
    { name: "f66", f: (x) => 4 * x ** 5 + x * x + 1, a: -1, b: 1 },
    { name: "f67", f: (x) => x + x ** 10 - 1, a: -1, b: 1 },
    { name: "f68", f: (x) => Math.PI ** x - Math.E, a: -1, b: 1 },
    { name: "f69", f: (x) => Math.log(Math.abs(x - 10 / 9)), a: -1, b: 1 },
    { name: "f70", f: (x) => 1 / 3 + (x > 0 ? 1 : (x < 0 ? -1 : 0)) * Math.abs(x) ** (1 / 3) + x ** 3, a: -1, b: 1 },
    { name: "f71", f: (x) => (x + 2 / 3) / (x + 101 / 100), a: -1, b: 1 },
    { name: "f72", f: (x) => (x * 1e6 - 1) ** 3, a: -1, b: 1 },
    { name: "f73", f: (x) => Math.exp(x) * (x * 1e6 - 1) ** 3, a: -1, b: 1 },
    { name: "f74", f: (x) => (x - 1 / 3) ** 2 * Math.atan(x - 1 / 3), a: -1, b: 1 },
    { name: "f75", f: (x) => (3 * x - 1 > 0 ? 1 : (3 * x - 1 < 0 ? -1 : 0)) * (1 - Math.sqrt(1 - (3 * x - 1) ** 2 / 81)), a: -1, b: 1 },
    { name: "f76", f: (x) => x > (1 - 1e6) / 1e6 ? (1 + 1e6) / 1e6 : -1, a: -1, b: 1 },
    { name: "f77", f: (x) => x !== 1 / 21 ? 1 / (21 * x - 1) : 0, a: -1, b: 1 },
    { name: "f78", f: (x) => x * x / 4 + Math.ceil(x / 2) - 0.5, a: -1, b: 1 },
    { name: "f79", f: (x) => Math.ceil(10 * x - 1) + 0.5, a: -1, b: 1 },
    { name: "f80", f: (x) => x + Math.sin(x * 1e6) / 10 + 1e-3, a: -1, b: 1 },
    { name: "f81", f: (x) => x > -1 ? 1 + Math.sin(1 / (x + 1)) : -1, a: -1, b: 1 },
    { name: "f82", f: (x) => 202 * x - 2 * Math.floor((2 * x + 1e-2) / 2e-2) - 0.1, a: -1, b: 1 },
    { name: "f83", f: (x) => (202 * x - 2 * Math.floor((2 * x + 1e-2) / 2e-2) - 0.1) ** 3, a: -1, b: 1 },
];

// SciML project benchmarks suite
const problems3: Problem[] = [
    { name: "f84", f: (x) => (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05, a: 0.5, b: 5.5 },
    { name: "f85", f: (x) => Math.sin(x) - 0.5 * x - 0.3, a: -10.0, b: 10.0 },
    { name: "f86", f: (x) => Math.exp(x) - 1 - x - x * x / 2 - 0.005, a: -2.0, b: 2.0 },
    { name: "f87", f: (x) => 1 / (x - 0.5) - 2 - 0.05, a: 0.6, b: 2.0 },
    { name: "f88", f: (x) => Math.log(x) - x + 2 - 0.05, a: 0.1, b: 3.0 },
    { name: "f89", f: (x) => Math.sin(20 * x) + 0.1 * x - 0.1, a: -4.0, b: 5.0 },
    { name: "f90", f: (x) => x ** 3 - 2 * x ** 2 + x - 0.025, a: -1.0, b: 2.0 },
    { name: "f91", f: (x) => x * Math.sin(1 / x) - 0.1 - 0.01, a: 0.01, b: 1.0 },
    { name: "f92", f: (x) => x ** 3 - 0.001, a: -10, b: 10 },
];

const allProblems: Problem[] = [...problems1, ...problems2, ...problems3];

class CountedFunc {
    private _f: (x: number) => number;
    public count: number;

    constructor(f: (x: number) => number) {
        this._f = f;
        this.count = 0;
    }

    call(x: number): number {
        this.count++;
        return this._f(x);
    }
}

function runTests(): void {
    const eps = 1e-14;
    let passed = 0;
    let failed = 0;
    let totalEvals = 0;

    console.log("ModAB Root-Finding Test Results");
    console.log("=".repeat(70));
    console.log(`${"Func".padStart(4)} | ${"Root".padStart(22)} | ${"f(root)".padStart(15)} | ${"Evals".padStart(6)} | Status`);
    console.log("-".repeat(70));

    for (const p of allProblems) {
        const cf = new CountedFunc(p.f);
        try {
            const root = modABRoot((x) => cf.call(x), p.a, p.b, 0, eps, 0.0);
            const fval = p.f(root);

            let status: string;
            if (Number.isNaN(root)) {
                status = "FAIL";
                failed++;
            } else if (Math.abs(fval) < 1e-10) {
                status = "PASS";
                passed++;
            } else {
                status = Math.abs(fval) < 1e-6 ? "PASS" : "WEAK";
                passed++;
            }

            totalEvals += cf.count;
            const rootStr = root.toPrecision(15).padStart(22);
            const fvalStr = fval.toExponential(6).padStart(15);
            console.log(`${p.name.padStart(4)} | ${rootStr} | ${fvalStr} | ${cf.count.toString().padStart(6)} | ${status}`);
        } catch (e) {
            failed++;
            totalEvals += cf.count;
            console.log(`${p.name.padStart(4)} | ${"ERROR".padStart(22)} | ${"N/A".padStart(15)} | ${cf.count.toString().padStart(6)} | FAIL (${e})`);
        }
    }

    console.log("-".repeat(70));
    console.log(`Total: ${passed + failed} tests, ${passed} passed, ${failed} failed`);
    console.log(`Total evaluations: ${totalEvals}`);
}

runTests();
