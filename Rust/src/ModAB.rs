//! The modified Anderson-Bjorck (modAB) bracketing root finder.
//!
//! The `safe_*` helpers below carry the overflow/NaN postconditions the solver
//! relies on. They are the Rust form of the reference implementation in
//! `C#/Root/Node.cs` and `C#/Root/Solvers/{Solver,RModAB}.cs`, so a change to
//! the numerics belongs in one place rather than inline in the solver.

/// Returns true only when both values have the same non-zero sign.
///
/// Comparisons with NaN are false, so a caller must reject NaN before using
/// this predicate to update a bracket.
#[inline]
pub fn same_sign(x: f64, y: f64) -> bool {
    (x < 0.0 && y < 0.0) || (x > 0.0 && y > 0.0)
}

/// Midpoint without first forming `x1 + x2`, which can overflow for finite
/// endpoints of the same sign.
///
/// For finite ordered endpoints the result lies in the closed bracket, so
/// callers need no extra clamp.
#[inline]
pub fn safe_midpoint(x1: f64, x2: f64) -> f64 {
    0.5 * x1 + 0.5 * x2
}

/// Safeguarded false-position/secant point for an ordered bracket `[x1, x2]`.
///
/// The textbook formula `(x1*y2 - x2*y1) / (y2 - y1)` may overflow although the
/// mathematical intersection is finite. With opposite-sign ordinates, writing
/// `a = |y1|` and `b = |y2|` gives the equivalent convex combination
///
/// ```text
///     x = (b/(a+b))*x1 + (a/(a+b))*x2,
/// ```
///
/// whose weights lie in `[0, 1]` and sum to 1. This helper owns the complete
/// postcondition every caller needs: the returned point is finite and lies in
/// `[x1, x2]`. Where the secant geometry cannot deliver that, the safe midpoint
/// is returned instead.
#[inline]
pub fn safe_secant(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let mut a = y1.abs();
    let mut b = y2.abs();
    let mut den = a + b;

    // One test on the denominator covers every unusable case: a NaN ordinate
    // propagates into it, two zero ordinates make it zero, and an infinite
    // ordinate or an overflowing sum makes it infinite. A zero magnitude does
    // NOT indicate a root here, because the ordinates may be Anderson-Bjorck
    // auxiliary values, so bisection is the safe and neutral fallback.
    if !(den > 0.0) {
        return safe_midpoint(x1, x2);
    }

    if den.is_infinite() {
        // An infinite ordinate carries no usable slope. Otherwise a + b merely
        // overflowed, and halving both restores it without changing the ratio
        // that defines the weights.
        if a.is_infinite() || b.is_infinite() {
            return safe_midpoint(x1, x2);
        }
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    let x = (b / den) * x1 + (a / den) * x2;

    // In exact arithmetic the convex combination is strictly inside the
    // bracket; the projection only corrects a possible last-ulp excursion.
    x.clamp(x1, x2)
}

/// The Anderson-Bjorck contraction factor for the ordinate that did not move.
#[inline]
fn ab_factor(y3: f64, y_moved: f64) -> f64 {
    let m = 1.0 - y3 / y_moved;
    if m > 0.0 { m } else { 0.5 }
}

// Root of the linear-fractional model  f(x) = (x - r) / (p*x + q)  through the bracket
// (x1, y1), (x2, y2) and a third point (x3, f3) outside it. The bracket ordinates are the
// A&B-corrected residuals and f3 is the true residual of the dropped endpoint. The model fits
// both nearly-quadratic shapes and functions with a pole; for a straight line it gives the
// secant point. With w = x2 - x1, u = x3 - x1 and the residuals relative to y2,
// q1 = y1/y2 and q3 = f3/y2 (so the step does not depend on the scale of f), and
// c = 1 - q1 > 1, the root is
//     r = x1 + t,   t = w*u*q1*(q3 - 1) / (w*(q3 - q1) - u*q3*c),
// two divisions per step, 1/y2 and this one. The A&B secant point is xs = x1 - w*q1/c and
// its distance to the nearer endpoint is d = w*min(|q1|, 1)/c. The step is kept within
// 0.7*d of xs (on a plateau the model degenerates towards a step function), so it stays
// at least 0.3*d inside the bracket. Multiplied by c, the test reads
//     |t*c + w*q1| <= 0.7*w*min(|q1|, 1),
// with no division. A step outside it is clamped to xs +- 0.7*d on the side the model
// points to; this did better than the A&B secant in randomized-bracket tests and on
// 200,000 random functions, where factors 0.6 to 0.9 differed little.
// A step that rounds onto x2, or a NaN step (from overflow in t, or a NaN third point
// when no endpoint was dropped on that side yet), returns x2. Its residual is already
// stored, so the caller reuses it instead of evaluating f; the closing step may first
// move the point to x2 - epsx/2.
#[inline]
pub fn rational_step(x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, f3: f64) -> f64 {
    let k = 1.0 / y2;
    let (q1, q3) = (y1 * k, f3 * k);
    let c = 1.0 - q1;
    let (w, u) = (x2 - x1, x3 - x1);
    let mut t = w * u * q1 * (q3 - 1.0) / (w * (q3 - q1) - u * q3 * c);
    let a = q1.abs();
    let lim = 0.7 * (if a < 1.0 { w * a } else { w }); // c * 0.7 * d
    let dev = t * c + w * q1; // c * (x - xs)
    // Clamp: t = ts +- 0.7*d. False for a NaN dev, which only comes from a NaN t (or q1, and
    // then t is NaN too), so x is NaN and the test below returns x2.
    if dev.abs() > lim {
        t = (lim.copysign(dev) - w * q1) / c;
    }
    let x = x1 + t;
    // Rounding can put x on x2 when the bracket is a few ulps wide; NaN fails the test too.
    if x < x2 { x } else { x2 }
}

// True if the second divided differences on the left (xl, x1, x2) and on the right
// (x1, x2, xr) have opposite signs: an inflection within the bracket, where a model with
// one-signed curvature does worse than the secant. Scaled by positive factors, they are
//     pl = (f2 - f1)*(x1 - xl) - (f1 - fl)*(x2 - x1),
//     pr = (fr - f2)*(x2 - x1) - (f2 - f1)*(xr - x2),
// so only their signs are compared, with no division or product of the two.
// False when xl or xr is NaN (no endpoint dropped on that side yet).
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn curvature_changes_sign(
    xl: f64, fl: f64, x1: f64, f1: f64, x2: f64, f2: f64, xr: f64, fr: f64,
) -> bool {
    let (df, w) = (f2 - f1, x2 - x1);
    let pl = df * (x1 - xl) - (f1 - fl) * w;
    let pr = (fr - f2) * w - df * (xr - x2);
    (pl < 0.0 && pr > 0.0) || (pl > 0.0 && pr < 0.0)
}

// Finds the root of "F(x) = 0" within the interval [x1, x2]
// with the specified precisions - absolute: aTol and relative: rTol,
// using an improved version of the modified Anderson Bjork's method:
//     Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//     Improvements to the Modified Anderson-Bjorck (modAB) Root-Finding Algorithm.
//     Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332
// Additional fixes proposed by L. Tomov are applied in this version:
//     1. The secant point is clamped to the interval [x1, x2] before the X-convergence exit
//     2. The original function values y1 and y2 (without A&B corrections)
//        are stored for later use in bisection fallback
// The overflow- and NaN-safe form of the interpolation lives in the shared
// safeguards above; the switching test is written so that it cannot overflow.
// The A&B phase adds two changes:
//     1. Rational correction: the A&B secant is corrected for curvature by a 3-point rational
//        model through the A&B-corrected bracket and the endpoint it most recently dropped,
//        with the step kept within 0.7 of the distance from the A&B secant point to the
//        nearer endpoint. The plain A&B secant is used when the curvature changes sign
//        within the bracket.
//     2. Closing step: an A&B point is kept at least epsx/2 away from the endpoint that was
//        just evaluated, so that one-sided convergence still closes the bracket.
// F(x) must be continuous and sign(F(x1)) != sign(F(x2))
pub fn mod_ab_root<F>(
    f: F,
    mut x1: f64,
    mut x2: f64,
    y: f64,
    xtol: f64,
    ytol: f64,
    maxiter: usize,
) -> f64
where
    F: Fn(f64) -> f64,
{
    if x2 < x1 {
        std::mem::swap(&mut x1, &mut x2);
    }
    let epsy = ytol * y.abs().max(1.0);
    let mut y1 = f(x1) - y;
    if y1.abs() <= epsy {
        return x1;
    }
    let mut y2 = f(x2) - y;
    if y2.abs() <= epsy {
        return x2;
    }
    // NaN has no usable sign, and `same_sign` is false for it, so it
    // must be rejected before the predicate is used to update a bracket.
    if y1.is_nan() || y2.is_nan() || same_sign(y1, y2) {
        return f64::NAN;
    }
    // True residuals, kept unmodified by A&B corrections
    let mut f1 = y1;
    let mut f2 = y2;
    let mut side: i32 = 0;
    let mut bisection = true;
    let mut threshold = x2 - x1;
    // A switching bisection whose midpoint reduces the residual of the replaced endpoint by less
    // than 1 % lies on a plateau; the first A&B step is then not rescaled.
    const PLATEAU_RATIO: f64 = 0.99;
    // Each return from A&B to bisection halves the margin of the switching test, so that
    // borderline switches are rejected after unsuccessful A&B phases.
    let mut margin = 1.0;
    let mut last: i32 = 0; // Endpoint set by the last evaluation: 1 = x1, -1 = x2, 0 = none
    // Last endpoints dropped on each side
    let (mut xl, mut fl, mut xr, mut fr) = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    for _ in 0..maxiter {
        let dx = x2 - x1; // Bracket width
        let x3;
        let epsx;
        if bisection {
            x3 = safe_midpoint(x1, x2);
            epsx = xtol * x3.abs().max(1.0);
        } else {
            // Rational correction through the endpoint dropped by the last update, or the A&B
            // secant when the curvature changes sign within the bracket. A missing endpoint is
            // NaN, which fails both tests, so it needs no check of its own.
            let mut xc = if curvature_changes_sign(xl, fl, x1, f1, x2, f2, xr, fr) {
                safe_secant(x1, y1, x2, y2)
            } else if last == 1 {
                rational_step(x1, y1, x2, y2, xl, fl)
            } else {
                rational_step(x1, y1, x2, y2, xr, fr)
            };
            epsx = xtol * xc.abs().max(1.0);
            // Closing step: when A&B converges from one side, its next point lands within epsx/2
            // of the endpoint just evaluated and the bracket does not shrink. A step of epsx/2
            // towards the far side then crosses the root and closes the bracket.
            let h = 0.5 * epsx;
            if last == 1 {
                if xc - x1 < h && x1 + h < x2 {
                    xc = x1 + h;
                }
            } else if x2 - xc < h && x2 - h > x1 {
                xc = x2 - h;
            }
            x3 = xc;
        }
        if dx <= epsx {
            return x3;
        }
        let y3: f64;
        let mut switching = false; // This bisection switches to A&B
        if bisection {
            y3 = f(x3) - y;
            let dy = f2 - f1;
            if dy.is_finite() {
                // Avoids overflow in the calculations below
                let ym = (f1 + f2) * 0.5; // Chord ordinate at midpoint; f1, f2 have opposite signs
                let r = 1.0 - (ym / dy).abs(); // Symmetry factor
                let k = r * r; // Deviation factor
                // k*|ym| + k*|y3| cannot overflow; an infinite y3 fails the test.
                if (ym - y3).abs() < (k * ym.abs() + k * y3.abs()) * margin {
                    bisection = false;
                    switching = true;
                    threshold = 2.0 * dx;
                    (y1, y2) = (f1, f2); // A&B starts from the true residuals
                }
            }
        } else {
            // If x3 got clamped, reuse the true residual stored at the endpoint.
            if x3 == x1 {
                y3 = f1;
            } else if x3 == x2 {
                y3 = f2;
            } else {
                y3 = f(x3) - y;
            }
            threshold *= 0.5;
        }
        let ay3 = y3.abs();
        if !(ay3 > epsy) {
            // Exit on y-convergence or NaN with a single branch on the hot path
            if ay3 <= epsy {
                return x3;
            }
            // A NaN residual has no usable sign, so the bracket cannot be updated.
            return f64::NAN;
        }
        if bisection {
            if same_sign(f1, y3) {
                (xl, fl) = (x1, f1);
                (x1, f1, last) = (x3, y3, 1);
            } else {
                (xr, fr) = (x2, f2);
                (x2, f2, last) = (x3, y3, -1);
            }
        } else {
            let (yl, yr) = (f1, f2); // True residuals of the bracket before the update
            if same_sign(f1, y3) {
                if side == 1 {
                    y2 *= ab_factor(y3, y1);
                } else if !switching || ay3 <= PLATEAU_RATIO * yl.abs() {
                    side = 1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if x3 != x1 {
                    (xl, fl) = (x1, f1);
                }
                (x1, y1, f1, last) = (x3, y3, y3, 1);
            } else {
                if side == -1 {
                    y1 *= ab_factor(y3, y2);
                } else if !switching || ay3 <= PLATEAU_RATIO * yr.abs() {
                    side = -1; // A switching bisection on a plateau starts the A&B phase afresh
                }
                if x3 != x2 {
                    (xr, fr) = (x2, f2);
                }
                (x2, y2, f2, last) = (x3, y3, y3, -1);
            }
            // Fallback if AB fails to reduce the bracket width, unless it still shrinks the residual more than 42 times
            // (2^(1074/200) = 41.4: at this rate the residual underflows within 200 steps from |y| = 1)
            if x2 - x1 > threshold && 42.0 * ay3 >= yl.abs().min(yr.abs()) {
                bisection = true;
                side = 0;
                margin *= 0.5;
            }
        }
    }
    f64::NAN
}

#[cfg(test)]
mod safeguard_tests {
    //! Edge-case tests for the shared modAB safeguards.
    //!
    //! The 100-problem benchmark suite never produces a non-finite or overflowing
    //! residual, so the overflow/NaN branches of same_sign, safe_midpoint and
    //! safe_secant, and the overflow cases of the switching test, are covered here.
    //!
    //! Expected values come from the C# reference in C#/Root/Node.cs and
    //! C#/Root/Solvers/RModAB.cs; every language port is held to the same table.

    use super::*;

    fn eq(got: f64, want: f64) -> bool {
        (got.is_nan() && want.is_nan()) || got == want
    }

    #[test]
    fn same_sign_edge_cases() {
        assert_eq!(same_sign(1.0f64, 2.0f64), true, "case 0");
        assert_eq!(same_sign(-1.0f64, -2.0f64), true, "case 1");
        assert_eq!(same_sign(1.0f64, -2.0f64), false, "case 2");
        assert_eq!(same_sign(-1.0f64, 2.0f64), false, "case 3");
        assert_eq!(same_sign(0.0f64, 1.0f64), false, "case 4");
        assert_eq!(same_sign(1.0f64, 0.0f64), false, "case 5");
        assert_eq!(same_sign(0.0f64, 0.0f64), false, "case 6");
        assert_eq!(same_sign(-0.0f64, -1.0f64), false, "case 7");
        assert_eq!(same_sign(f64::NAN, 1.0f64), false, "case 8");
        assert_eq!(same_sign(1.0f64, f64::NAN), false, "case 9");
        assert_eq!(same_sign(f64::NAN, f64::NAN), false, "case 10");
        assert_eq!(same_sign(f64::INFINITY, 1.0f64), true, "case 11");
        assert_eq!(same_sign(-f64::INFINITY, -1.0f64), true, "case 12");
        assert_eq!(same_sign(f64::INFINITY, -f64::INFINITY), false, "case 13");
    }

    #[test]
    fn safe_midpoint_edge_cases() {
        assert!(eq(safe_midpoint(2.0f64, 4.0f64), 3.0f64), "case 0");
        assert!(eq(safe_midpoint(1e+308f64, 1e+308f64), 1e+308f64), "case 1");
        assert!(eq(safe_midpoint(-1e+308f64, 1e+308f64), 0.0f64), "case 2");
        assert!(eq(safe_midpoint(1e+308f64, 1.7e+308f64), 1.35e+308f64), "case 3");
        assert!(eq(safe_midpoint(-1.7e+308f64, -1e+308f64), -1.35e+308f64), "case 4");
        assert!(eq(safe_midpoint(0.0f64, 1.0f64), 0.5f64), "case 5");
        assert!(eq(safe_midpoint(-1.0f64, 1.0f64), 0.0f64), "case 6");
    }

    #[test]
    fn safe_secant_edge_cases() {
        assert!(eq(safe_secant(0.0f64, -1.0f64, 1.0f64, 1.0f64), 0.5f64), "case 0");
        assert!(eq(safe_secant(0.0f64, -1.0f64, 1.0f64, 3.0f64), 0.25f64), "case 1");
        assert!(eq(safe_secant(0.0f64, -1e+308f64, 1.0f64, 1e+308f64), 0.5f64), "case 2");
        assert!(eq(safe_secant(0.0f64, -1.7e+308f64, 1.0f64, 1.7e+308f64), 0.5f64), "case 3");
        assert!(eq(safe_secant(0.0f64, f64::INFINITY, 1.0f64, -1.0f64), 0.5f64), "case 4");
        assert!(eq(safe_secant(0.0f64, -1.0f64, 1.0f64, f64::INFINITY), 0.5f64), "case 5");
        assert!(eq(safe_secant(0.0f64, 0.0f64, 1.0f64, 0.0f64), 0.5f64), "case 6");
        assert!(eq(safe_secant(0.0f64, f64::NAN, 1.0f64, 1.0f64), 0.5f64), "case 7");
        assert!(eq(safe_secant(0.0f64, -1e-300f64, 1.0f64, 1e+300f64), 0.0f64), "case 8");
        assert!(eq(safe_secant(0.0f64, -1e+300f64, 1.0f64, 1e-300f64), 1.0f64), "case 9");
        assert!(eq(safe_secant(1e+308f64, -1.0f64, 1.7e+308f64, 1.0f64), 1.35e+308f64), "case 10");
    }

    /// Solver-level regression: this returned +Inf before `safe_midpoint`.
    #[test]
    fn overflowing_bracket() {
        let root = mod_ab_root(|x| x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
        assert!((root - 1.2e308).abs() <= 1e-13 * 1.2e308, "got {root}");
    }

    fn ck_root(name: &str, f: impl Fn(f64) -> f64, a: f64, b: f64, want: f64) {
        let got = mod_ab_root(f, a, b, 0.0, 1e-14, 0.0, 200);
        assert!((got - want).abs() <= 1e-13 * want.abs().max(1.0), "{name}: got {got}, want {want}");
    }

    /// Solver-level cases for overflow and infinite residuals in the switching test.
    #[test]
    fn switching_test_overflow() {
        let m = f64::MAX;
        // |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite
        ck_root("hump", |x| if x <= 0.0 { m * (-0.2 + 1.2 * (x + 1.0)) } else { m * (1.0 - 0.2 * x) },
                -1.0, 1.0, -5.0 / 6.0);
        // f2 - f1 overflows: switching is disabled until the residuals shrink
        ck_root("f2-f1 overflow", |x| 1.7e308 * (10.0 * (x - 0.3)).tanh(), -1.0, 1.0, 0.3);
        // infinite residuals at one or both ends of the bracket
        ck_root("inf left", |x| if x < -0.5 { f64::NEG_INFINITY } else { x - 0.1 }, -1.0, 1.0, 0.1);
        ck_root("inf both", |x| if x < -0.5 { f64::NEG_INFINITY } else if x > 0.9 { f64::INFINITY } else { x - 0.1 },
                -1.0, 1.0, 0.1);
        // subnormal residuals: AB corrections underflow towards zero
        ck_root("subnormal", |x| 1e-300 * (x * x * x - 0.2), -1.0, 1.0, 0.2f64.cbrt());
    }
}
