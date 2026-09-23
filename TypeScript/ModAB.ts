/*
 * Shared safeguards
 *
 * These helpers carry the overflow/NaN postconditions that the bracketing
 * solver relies on. They are the TypeScript form of the reference
 * implementation in C#/Root/Node.cs and C#/Root/Solvers/{Solver,SgModAB}.cs,
 * so a change to the numerics belongs in one place rather than inline in the
 * solver.
 */

/**
 * Returns true only when both values have the same non-zero sign.
 *
 * Comparisons with NaN are false, so a caller must reject NaN before using this
 * predicate to update a bracket.
 */
export function sameSign(x: number, y: number): boolean {
    return (x < 0 && y < 0) || (x > 0 && y > 0);
}

/**
 * Midpoint without first forming x1 + x2, which can overflow for finite
 * endpoints of the same sign. For finite ordered endpoints the result lies in
 * the closed bracket, so callers need no extra clamp.
 */
export function safeMidpoint(x1: number, x2: number): number {
    return 0.5 * x1 + 0.5 * x2;
}

/**
 * Safeguarded false-position/secant point for an ordered bracket [x1, x2].
 *
 * The textbook formula (x1*y2 - x2*y1) / (y2 - y1) may overflow although the
 * mathematical intersection is finite. With opposite-sign ordinates, writing
 * a = |y1| and b = |y2| gives the equivalent convex combination
 *
 *     x = (b/(a+b))*x1 + (a/(a+b))*x2,
 *
 * whose weights lie in [0, 1] and sum to 1. This helper owns the complete
 * postcondition every caller needs: the returned point is finite and lies in
 * [x1, x2]. Where the secant geometry cannot deliver that, the safe midpoint is
 * returned instead.
 */
export function safeSecant(x1: number, y1: number, x2: number, y2: number): number {
    let a = Math.abs(y1);
    let b = Math.abs(y2);
    let den = a + b;

    // One test on the denominator covers every unusable case: a NaN ordinate
    // propagates into it, two zero ordinates make it zero, and an infinite
    // ordinate or an overflowing sum makes it infinite. A zero magnitude does
    // NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    // auxiliary values, so bisection is the safe and neutral fallback.
    if (!(den > 0)) {
        return safeMidpoint(x1, x2);
    }

    if (!Number.isFinite(den)) {
        // An infinite ordinate carries no usable slope. Otherwise a + b merely
        // overflowed, and halving both restores it without changing the ratio
        // that defines the weights.
        if (!Number.isFinite(a) || !Number.isFinite(b)) {
            return safeMidpoint(x1, x2);
        }
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    const x = (b / den) * x1 + (a / den) * x2;

    // In exact arithmetic the convex combination is strictly inside the
    // bracket; the projection only corrects a possible last-ulp excursion.
    return x < x1 ? x1 : (x > x2 ? x2 : x);
}

/** The Anderson-Bjorck contraction factor for the ordinate that did not move. */
function abFactor(y3: number, yMoved: number): number {
    const m = 1 - y3 / yMoved;
    return m > 0 ? m : 0.5;
}

/*
* Finds the root of "F(x) = 0" within the interval [x1, x2]
* with the specified precisions - absolute: aTol and relative: rTol,
* using an improved version of the modified Anderson Bjork's method:
*     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
*     Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
*     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
* Additional fixes proposed by L. Tomov are applied in this version:
*     1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
*     2. The original function values y1 and y2 (without A&B corrections)
*        are stored for later use in bisection fallback
* The overflow- and NaN-safe form of the interpolation lives in the shared
* safeguards above; the switching test is written so that it cannot overflow.
* F(x) must be continuous and sign(F(x1)) != sign(F(x2))
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
    // NaN has no usable sign, and sameSign is false for it, so it must
    // be rejected before the predicate is used to update a bracket.
    if (Number.isNaN(y1) || Number.isNaN(y2) || sameSign(y1, y2)) {
        return NaN;
    }
    let side = 0;
    let f1 = y1, f2 = y2; // True residuals, kept unmodified by A&B corrections
    let ymin = 0.0; // Best true residual of the bracket
    let bisection = true;
    let threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    for (let i = 0; i < maxiter; i++) {
        const x3 = bisection ? safeMidpoint(x1, x2) : safeSecant(x1, y1, x2, y2);
        const epsx = xtol * Math.max(Math.abs(x3), 1);
        if (x2 - x1 <= epsx) { // x-convergence check
            return x3;
        }
        let y3: number;
        if (bisection) {
            y3 = f(x3) - y; // Function value at midpoint
            if (Number.isFinite(f2 - f1)) { // Avoids overflow in the calculations below
                const ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                const r = 1 - Math.abs(ym / (f2 - f1)); // Symmetry factor
                const k = r * r; // Deviation factor
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (Math.abs(ym - y3) < k * Math.abs(ym) + k * Math.abs(y3)) {
                    bisection = false;
                    threshold = 2 * (x2 - x1); // Safety factor
                    y1 = f1; y2 = f2; // A&B starts from the true residuals
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            y3 = x3 === x1 ? f1 :
                 x3 === x2 ? f2 :
                 f(x3) - y;

            threshold *= 0.5;
            ymin = Math.min(Math.abs(f1), Math.abs(f2));
        }
        if (Math.abs(y3) <= epsy) { // y-convergence check
            return x3;
        }
        // A NaN residual has no usable sign, so the bracket cannot be updated.
        if (Number.isNaN(y3)) {
            return NaN;
        }
        if (bisection) {
        if (sameSign(f1, y3)) { // Same sign check
            x1 = x3; f1 = y3;
        } else {
            x2 = x3; f2 = y3;
        }
        } else {
            if (sameSign(f1, y3)) { // Same sign check
                if (side === 1) {
                    y2 *= abFactor(y3, y1);
                } else if (!bisection) {
                    side = 1;
                }
                x1 = x3; f1 = y1 = y3;
            } else {
                if (side === -1) {
                    y1 *= abFactor(y3, y2);
                } else if (!bisection) {
                    side = -1;
                }
                x2 = x3; f2 = y2 = y3;
            }        
            // Fallback if AB fails to reduce the bracket width, unless it still halves the residual
            if (x2 - x1 > threshold && Math.abs(y3) > 0.5 * ymin) {
                bisection = true;
                side = 0;
            }
        }
    }
    return NaN;
}
