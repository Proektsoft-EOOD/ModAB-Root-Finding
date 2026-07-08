/*
* Finds the root of "F(x) = 0" within the interval [x1, x2]
* with the specified precisions - absolute: aTol and relative: rTol,
* using an improved version of the modified Anderson Bjork's method:
*     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
*     Improvements to the Modified Anderson–Björck(modAB) Root-Finding Algorithm.
*     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
* Additional fixes proposed by L. Tomov are applied in this version:
*     1. The secant point is clamped to the interval [p1.X, p2.X] before the X-convergence exit
*     2. The original function values y1 and y2 (without A&B corrections) 
*        are stored for later use in bisection fallback
* F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))
 */
export function modABRoot(
    f: (x: number) => number,
    x1: number,
    x2: number,
    y: number,
    xtol: number = 1e-14,
    ytol: number = 0.0,
    maxiter: number = 200
): number {
    if (x2 < x1) {
        [x1, x2] = [x2, x1];
    }
    const epsy = ytol * Math.max(Math.abs(y), 1);
    let y1 = f(x1) - y;
    if (Math.abs(y1) <= epsy) {
        return x1;
    }
    let y2 = f(x2) - y;
    if (Math.abs(y2) <= epsy) {
        return x2;
    }
    if ((y1 > 0) === (y2 > 0)) {
        return NaN;
    }
    let side = 0;
    let f1 = y1, f2 = y2;
    let bisection = true;
    let threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    for (let i = 0; i < maxiter; i++) {
        let x3 = bisection ? (x1 + x2) * 0.5 : (x1 * y2 - y1 * x2) / (y2 - y1);
        const epsx = xtol * Math.max(Math.abs(x3), 1);
        if (x2 - x1 <= epsx) { // x-convergence check
            return bisection ? x3 : Math.max(x1, Math.min(x3, x2));
        }
        let y3: number;
        if (bisection) {
            y3 = f(x3) - y; // Function value at midpoint
            const ym = (f1 + f2) * 0.5; // Ordinate of chord at midpoint
            const dy = f2 - f1;
            const r = 1 - Math.abs(ym / dy); // Symmetry factor
            const k = r * r; // Deviation factor
            if (Math.abs(ym - y3) < k * (Math.abs(y3) + Math.abs(ym))) {
                bisection = false;
                threshold = (x2 - x1) * 16; // Safety factor: 4 bisection iterations = 2^4
            }
        } else {
            if (x3 <= x1) {
                x3 = x1;
                y3 = f1;
            } else if (x3 >= x2) {
                x3 = x2;
                y3 = f2;
            } else {
                y3 = f(x3) - y;
            }
            threshold *= 0.5;
        }
        if (Math.abs(y3) <= epsy) { // y-convergence check
            return x3;
        }
        if ((y1 > 0) === (y3 > 0)) { // Same sign check
            if (side === 1) {
                const m = 1 - y3 / y1;
                y2 = m > 0 ? y2 * m : y2 * 0.5;
            } else if (!bisection) {
                side = 1;
            }
            x1 = x3; f1 = y1 = y3;
        } else {
            if (side === -1) {
                const m = 1 - y3 / y2;
                y1 = m > 0 ? y1 * m : y1 * 0.5;
            } else if (!bisection) {
                side = -1;
            }
            x2 = x3; f2 = y2 = y3;
        }
        if (x2 - x1 > threshold) { // AB failed to shrink the interval enough
            bisection = true;
            side = 0;
        }
    }
    return NaN;
}
