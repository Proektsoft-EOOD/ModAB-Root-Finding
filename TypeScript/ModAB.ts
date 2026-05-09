/**
 * Finds the root of f(x) = y within [x1, x2] using
 * modified Anderson-Bjork method (Ganchovski, Traykov, 2023).
 * f(x) must be continuous and sign(f(x1) - y) != sign(f(x2) - y).
 * It also includes the recent changes from the following paper:
 * Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
 * Improvements to the Modified Anderson-Bjork (modAB) Root-Finding Algorithm. Algorithms 2026, 19, 332.
 * https://doi.org/10.3390/a19050332
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
    let bisection = true;
    let threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    for (let i = 0; i < maxiter; i++) {
        let x3 = bisection ? (x1 + x2) * 0.5 : (x1 * y2 - y1 * x2) / (y2 - y1);
        const epsx = xtol * Math.max(Math.abs(x3), 1);
        if (x2 - x1 <= epsx) { // x-convergence check
            return x3;
        }
        let y3: number;
        if (bisection) {
            y3 = f(x3) - y; // Function value at midpoint
            const ym = (y1 + y2) * 0.5; // Ordinate of chord at midpoint
            const dy = y2 - y1;
            const r = 1 - Math.abs(ym / dy); // Symmetry factor
            const k = r * r; // Deviation factor
            if (Math.abs(ym - y3) < k * (Math.abs(y3) + Math.abs(ym))) {
                bisection = false;
                threshold = (x2 - x1) * 8; // Safety factor: 4 bisection iterations = 2^4
            }
        } else {
            if (x3 <= x1) {
                x3 = x1;
                y3 = y1;
            } else if (x3 >= x2) {
                x3 = x2;
                y3 = y2;
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
            x1 = x3;
            y1 = y3;
        } else {
            if (side === -1) {
                const m = 1 - y3 / y2;
                y1 = m > 0 ? y1 * m : y1 * 0.5;
            } else if (!bisection) {
                side = -1;
            }
            x2 = x3;
            y2 = y3;
        }
        if (x2 - x1 > threshold) { // AB failed to shrink the interval enough
            bisection = true;
            side = 0;
        }
    }
    return NaN;
}
