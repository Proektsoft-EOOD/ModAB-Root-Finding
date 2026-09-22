/*
 * Shared safeguards
 *
 * These helpers carry the overflow/NaN postconditions that the bracketing
 * solver relies on. They are the TypeScript form of the reference
 * implementation in C#/Root/Node.cs and C#/Root/Solvers/{Solver,ModABCorr}.cs,
 * so a change to the numerics belongs in one place rather than inline in the
 * solver.
 */

/** Smallest positive subnormal double (2^-1074), the analogue of C#'s double.Epsilon. */
const MIN_SUBNORMAL = 5e-324;

/** Largest finite double, the analogue of C#'s double.MaxValue. */
const MAX_DOUBLE = Number.MAX_VALUE;

/** Math.sign-free copysign: returns |magnitude| with the sign of `sign`. */
function copySign(magnitude: number, sign: number): number {
    const m = Math.abs(magnitude);
    // Object.is distinguishes -0 from +0, which Math.sign and `< 0` do not.
    return sign < 0 || Object.is(sign, -0) ? -m : m;
}

/**
 * Returns true only when both values have the same non-zero sign.
 *
 * Comparisons with NaN are false, so a caller must reject NaN before using this
 * predicate to update a bracket.
 */
export function sameNonzeroSign(x: number, y: number): boolean {
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

/**
 * Returns k = r^2 for the symmetry-sensitive switching criterion. The
 * calculation is homogeneous in the true endpoint residuals.
 */
export function symmetryFactor(y1: number, y2: number): number {
    let a = Math.abs(y1);
    let b = Math.abs(y2);
    let den = a + b;

    if (!Number.isFinite(den)) {
        // Infinite true residuals deliberately disable switching and keep the
        // controller in bisection mode. NaN is returned rather than an infinity
        // because every exit of passesSwitchingTest is a "<" comparison, which
        // is false against NaN; an infinity would instead satisfy it and
        // switch. Residuals are never zero here, so only an overflowing sum
        // remains, and halving both restores it without changing the ratio.
        if (!Number.isFinite(a) || !Number.isFinite(b)) {
            return NaN;
        }
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    // |b-a| <= den, so the quotient lies in [0,1]; halving after the division
    // avoids forming 2*den, which could overflow.
    const r = 1 - Math.abs(b - a) / den / 2;
    return r * r;
}

/**
 * Tests whether the true midpoint value yf is close enough to the midpoint
 * value ym of the chord through the true endpoint residuals.
 */
export function passesSwitchingTest(ym: number, yf: number, symmetry: number): boolean {
    const absYm = Math.abs(ym);
    const absYf = Math.abs(yf);
    const sum = absYf + absYm;

    // Fast path. The exact-root case is handled before this is called, and a
    // non-finite ordinate or a NaN symmetry factor fails the comparison, which
    // disables switching as intended.
    if (Number.isFinite(sum)) {
        return Math.abs(ym - yf) < symmetry * sum;
    }

    // Only reached when the sum overflows. Non-finite values are unsuitable for
    // the linearity comparison.
    if (!Number.isFinite(ym) || !Number.isFinite(yf)) {
        return false;
    }

    // Normalize both sides of the homogeneous inequality to avoid overflow.
    const scale = Math.max(absYf, absYm);
    const normYm = ym / scale;
    const normYf = yf / scale;
    return Math.abs(normYm - normYf) < symmetry * (Math.abs(normYf) + Math.abs(normYm));
}

/** The Anderson-Bjorck contraction factor for the ordinate that did not move. */
function abFactor(y3: number, yMoved: number): number {
    const m = 1 - y3 / yMoved;
    return m > 0 ? m : 0.5;
}

/**
 * Multiplies an auxiliary Anderson-Bjorck ordinate by a positive factor while
 * preserving a finite non-zero sign in binary64 arithmetic. This keeps
 * sameNonzeroSign sound: an auxiliary ordinate that underflowed to zero would
 * otherwise silently change which branch of the bracket update is taken. It
 * acts only on auxiliary ordinates; an underflowed working value is never
 * accepted as a root of f.
 */
function scalePreservingNonzeroSign(value: number, positiveFactor: number): number {
    const scaled = value * positiveFactor;

    if (scaled === 0 && value !== 0) {
        return copySign(MIN_SUBNORMAL, value);
    }

    if (scaled === Infinity || scaled === -Infinity) {
        return copySign(MAX_DOUBLE, value);
    }

    return scaled;
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
* The overflow- and NaN-safe forms of the interpolation and switching
* arithmetic live in the shared safeguards above.
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
    // NaN has no usable sign, and sameNonzeroSign is false for it, so it must
    // be rejected before the predicate is used to update a bracket.
    if (Number.isNaN(y1) || Number.isNaN(y2) || sameNonzeroSign(y1, y2)) {
        return NaN;
    }
    let side = 0;
    let f1 = y1, f2 = y2, ymin = 0.0;
    let bisection = true;
    let threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    for (let i = 0; i < maxiter; i++) {
        // safeSecant already returns a point inside [x1, x2], so the separate
        // clamp on the convergence exit is no longer needed.
        const x3 = bisection ? safeMidpoint(x1, x2) : safeSecant(x1, y1, x2, y2);
        const epsx = xtol * Math.max(Math.abs(x3), 1);
        if (x2 - x1 <= epsx) { // x-convergence check
            return x3;
        }
        let y3: number;
        if (bisection) {
            y3 = f(x3) - y; // Function value at midpoint
            const ym = safeMidpoint(f1, f2); // Ordinate of chord at midpoint
            if (passesSwitchingTest(ym, y3, symmetryFactor(f1, f2))) {
                bisection = false;
                threshold = (x2 - x1) * 2; // Safety factor
            }
        } else {
            // If rounding makes the proposal coincide with an endpoint, reuse
            // the true residual already stored there.
            if (x3 === x1) {
                y3 = f1;
            } else if (x3 === x2) {
                y3 = f2;
            } else {
                y3 = f(x3) - y;
            }
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
        if (sameNonzeroSign(y1, y3)) { // Same sign check
            if (side === 1) {
                y2 = scalePreservingNonzeroSign(y2, abFactor(y3, y1));
            } else if (!bisection) {
                side = 1;
            }
            x1 = x3; f1 = y1 = y3;
        } else {
            if (side === -1) {
                y1 = scalePreservingNonzeroSign(y1, abFactor(y3, y2));
            } else if (!bisection) {
                side = -1;
            }
            x2 = x3; f2 = y2 = y3;
        }
        // Fallback if AB fails to reduce the bracket width, unless it still halves the residual
        if (!bisection && x2 - x1 > threshold && Math.abs(y3) > 0.5 * ymin) {
            bisection = true;
            side = 0;
        }
    }
    return NaN;
}
