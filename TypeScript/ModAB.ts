/*
 * Shared safeguards
 *
 * These helpers carry the overflow/NaN postconditions that the bracketing
 * solver relies on. They are the TypeScript form of the reference
 * implementation in C#/Root/Node.cs and C#/Root/Solvers/{Solver,RModAB}.cs,
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

/**
 * Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
 * (x1, y1), (x2, y2) and a third point (x3, f3) outside it. The bracket ordinates are the
 * A&B-corrected residuals and f3 is the true residual of the dropped endpoint. The model fits
 * both nearly-quadratic shapes and functions with a pole; for a straight line it gives the
 * secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative to y2,
 * q1 = y1/y2 and q3 = f3/y2 (so the step does not depend on the scale of f), and
 * c = 1 - q1 > 1, the root is
 *     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
 * two divisions per step, 1/y2 and this one. The A&B secant point is xs = x1 - w*q1/c and
 * its distance to the nearer endpoint is d = w*min(|q1|, 1)/c. The step is kept within
 * 0.7*d of xs (on a plateau the model degenerates towards a step function), so it stays
 * at least 0.3*d inside the bracket. Multiplied by c, the test reads
 *     |t*c + w*q1| <= 0.7*w*min(|q1|, 1),
 * with no division. A step outside it is clamped to xs +- 0.7*d on the side the model
 * points to; this did better than the A&B secant in randomized-bracket tests and on
 * 200,000 random functions, where factors 0.6 to 0.9 differed little.
 * A step that rounds onto x2, or a NaN step (from overflow in t, or a NaN third point
 * when no endpoint was dropped on that side yet), returns x2. Its residual is already
 * stored, so the caller reuses it instead of evaluating f; the closing step may first
 * move the point to x2 - epsx/2.
 */
export function rationalStep(x1: number, y1: number, x2: number, y2: number,
                             x3: number, f3: number): number {
    const k = 1 / y2;
    const q1 = y1 * k, q3 = f3 * k;
    const c = 1 - q1;
    const w = x2 - x1, u = x3 - x1;
    let t = w * u * q1 * (q3 - 1) / (w * (q3 - q1) - u * q3 * c);
    const a = Math.abs(q1);
    const lim = 0.7 * (a < 1 ? w * a : w); // c * 0.7 * d
    const dev = t * c + w * q1;            // c * (x - xs)
    // Clamp: t = ts +- 0.7*d. False for a NaN dev, which only comes from a NaN t (or q1, and
    // then t is NaN too), so x is NaN and the test below returns x2.
    if (Math.abs(dev) > lim) {
        t = ((dev < 0 ? -lim : lim) - w * q1) / c;
    }
    const x = x1 + t;
    // Rounding can put x on x2 when the bracket is a few ulps wide; NaN fails the test too.
    return x < x2 ? x : x2;
}

/**
 * True if the second divided differences on the left (xl, x1, x2) and on the right
 * (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
 * one-signed curvature does worse than the secant. Scaled by positive factors, they are
 *     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
 *     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
 * so only their signs are compared, with no division or product of the two.
 * False when xl or xr is NaN (no endpoint dropped on that side yet).
 */
export function curvatureChangesSign(xl: number, fl: number, x1: number, f1: number,
                                     x2: number, f2: number, xr: number, fr: number): boolean {
    const df = f2 - f1, w = x2 - x1;
    const pl = df * (x1 - xl) - (f1 - fl) * w;
    const pr = (fr - f2) * w - df * (xr - x2);
    return (pl < 0 && pr > 0) || (pl > 0 && pr < 0);
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
* The A&B phase adds two changes:
*     1. Rational correction: the A&B secant is corrected for curvature by a 3-point rational
*        model through the A&B-corrected bracket and the endpoint it most recently dropped,
*        with the step kept within 0.7 of the distance from the A&B secant point to the
*        nearer endpoint. The plain A&B secant is used when the curvature changes sign
*        within the bracket.
*     2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
*        just evaluated, so that one-sided convergence still closes the bracket.
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
    let bisection = true;
    let threshold = x2 - x1; // Threshold to fall back to bisection if AB fails to shrink the interval enough
    // A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    // than 1 % lies on a plateau; the first A&B step is then not rescaled.
    const PLATEAU_RATIO = 0.99;
    // Each return from A&B to bisection halves the margin of the switching test, so that
    // borderline switches are rejected after unsuccessful A&B phases.
    let margin = 1;
    let last = 0; // Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none
    let xl = NaN, fl = NaN, xr = NaN, fr = NaN; // Last endpoints dropped on each side
    for (let i = 0; i < maxiter; i++) {
        const dx = x2 - x1; // Bracket width
        let x3: number, epsx: number;
        if (bisection) {
            x3 = safeMidpoint(x1, x2);
            epsx = xtol * Math.max(Math.abs(x3), 1);
        } else {
            // Rational correction through the endpoint dropped by the last update, or the A&B
            // secant when the curvature changes sign within the bracket. A missing endpoint is
            // NaN, which fails both tests, so it needs no check of its own.
            if (curvatureChangesSign(xl, fl, x1, f1, x2, f2, xr, fr)) {
                x3 = safeSecant(x1, y1, x2, y2);
            } else {
                x3 = last === 1 ? rationalStep(x1, y1, x2, y2, xl, fl)
                                : rationalStep(x1, y1, x2, y2, xr, fr);
            }
            epsx = xtol * Math.max(Math.abs(x3), 1);
            // Closing step: when A&B converges from one side, its next point lands within epsx/2
            // of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
            // towards the far side then crosses the root and closes the bracket.
            const h = 0.5 * epsx;
            if (last === 1) {
                if (x3 - x1 < h && x1 + h < x2) x3 = x1 + h;
            } else if (x2 - x3 < h && x2 - h > x1) {
                x3 = x2 - h;
            }
        }
        if (dx <= epsx) { // x-convergence check
            return x3;
        }
        let y3: number;
        let switching = false; // This bisection switches to A&B
        if (bisection) {
            y3 = f(x3) - y; // Function value at midpoint
            const dy = f2 - f1;
            if (Number.isFinite(dy)) { // Avoids overflow in the calculations below
                const ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                const r = 1 - Math.abs(ym / dy); // Symmetry factor
                const k = r * r; // Deviation factor
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (Math.abs(ym - y3) < (k * Math.abs(ym) + k * Math.abs(y3)) * margin) {
                    bisection = false;
                    switching = true;
                    threshold = 2 * dx; // Safety factor
                    y1 = f1; y2 = f2; // A&B starts from the true residuals
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            y3 = x3 === x1 ? f1 :
                 x3 === x2 ? f2 :
                 f(x3) - y;

            threshold *= 0.5;
        }
        const ay3 = Math.abs(y3);
        if (!(ay3 > epsy)) { // Exit on y-convergence or NaN with a single branch on the hot path
            if (ay3 <= epsy) {
                return x3;
            }
            // A NaN residual has no usable sign, so the bracket cannot be updated.
            return NaN;
        }
        if (bisection) {
            if (sameSign(f1, y3)) { // Same sign check
                xl = x1; fl = f1;
                x1 = x3; f1 = y3; last = 1;
            } else {
                xr = x2; fr = f2;
                x2 = x3; f2 = y3; last = -1;
            }
        } else {
            const yl = f1, yr = f2; // True residuals of the bracket before the update
            if (sameSign(f1, y3)) { // Same sign check
                if (side === 1) {
                    y2 *= abFactor(y3, y1);
                } else if (!switching || ay3 <= PLATEAU_RATIO * Math.abs(yl)) {
                    side = 1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 !== x1) { xl = x1; fl = f1; }
                x1 = x3; f1 = y1 = y3; last = 1;
            } else {
                if (side === -1) {
                    y1 *= abFactor(y3, y2);
                } else if (!switching || ay3 <= PLATEAU_RATIO * Math.abs(yr)) {
                    side = -1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if (x3 !== x2) { xr = x2; fr = f2; }
                x2 = x3; f2 = y2 = y3; last = -1;
            }        
            // Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            // (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if (x2 - x1 > threshold && 42 * ay3 >= Math.min(Math.abs(yl), Math.abs(yr))) {
                bisection = true;
                side = 0;
                margin *= 0.5;
            }
        }
    }
    return NaN;
}
