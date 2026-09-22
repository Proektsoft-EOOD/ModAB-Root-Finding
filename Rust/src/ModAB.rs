//! The modified Anderson-Bjorck (modAB) bracketing root finder.
//!
//! The `safe_*` helpers below carry the overflow/NaN postconditions the solver
//! relies on. They are the Rust form of the reference implementation in
//! `C#/Root/Node.cs` and `C#/Root/Solvers/{Solver,ModABCorr}.cs`, so a change to
//! the numerics belongs in one place rather than inline in the solver.

/// Returns true only when both values have the same non-zero sign.
///
/// Comparisons with NaN are false, so a caller must reject NaN before using
/// this predicate to update a bracket.
#[inline]
pub fn same_nonzero_sign(x: f64, y: f64) -> bool {
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

/// Returns `k = r^2` for the symmetry-sensitive switching criterion.
///
/// The calculation is homogeneous in the true endpoint residuals.
#[inline]
pub fn symmetry_factor(y1: f64, y2: f64) -> f64 {
    let mut a = y1.abs();
    let mut b = y2.abs();
    let mut den = a + b;

    if den.is_infinite() {
        // Infinite true residuals deliberately disable switching and keep the
        // controller in bisection mode. NaN is returned rather than an infinity
        // because every exit of `passes_switching_test` is a "<" comparison,
        // which is false against NaN; an infinity would instead satisfy it and
        // switch. Residuals are never zero here, so only an overflowing sum
        // remains, and halving both restores it without changing the ratio.
        if a.is_infinite() || b.is_infinite() {
            return f64::NAN;
        }
        a *= 0.5;
        b *= 0.5;
        den = a + b;
    }

    // |b-a| <= den, so the quotient lies in [0,1]; halving after the division
    // avoids forming 2*den, which could overflow.
    let r = 1.0 - (b - a).abs() / den / 2.0;
    r * r
}

/// Tests whether the true midpoint value `yf` is close enough to the midpoint
/// value `ym` of the chord through the true endpoint residuals.
#[inline]
pub fn passes_switching_test(ym: f64, yf: f64, symmetry: f64) -> bool {
    let abs_ym = ym.abs();
    let abs_yf = yf.abs();
    let sum = abs_yf + abs_ym;

    // Fast path. The exact-root case is handled before this is called, and a
    // non-finite ordinate or a NaN symmetry factor fails the comparison, which
    // disables switching as intended.
    if sum.is_finite() {
        return (ym - yf).abs() < symmetry * sum;
    }

    // Only reached when the sum overflows. Non-finite values are unsuitable for
    // the linearity comparison.
    if !ym.is_finite() || !yf.is_finite() {
        return false;
    }

    // Normalize both sides of the homogeneous inequality to avoid overflow.
    let scale = abs_yf.max(abs_ym);
    let norm_ym = ym / scale;
    let norm_yf = yf / scale;
    (norm_ym - norm_yf).abs() < symmetry * (norm_yf.abs() + norm_ym.abs())
}

/// Smallest positive subnormal (2^-1074), the f64 analogue of C#'s `double.Epsilon`.
const MIN_SUBNORMAL: f64 = f64::from_bits(1);

/// The Anderson-Bjorck contraction factor for the ordinate that did not move.
#[inline]
fn ab_factor(y3: f64, y_moved: f64) -> f64 {
    let m = 1.0 - y3 / y_moved;
    if m > 0.0 { m } else { 0.5 }
}

/// Multiplies an auxiliary Anderson-Bjorck ordinate by a positive factor while
/// preserving a finite non-zero sign in binary64 arithmetic.
///
/// This keeps [`same_nonzero_sign`] sound: an auxiliary ordinate that
/// underflowed to zero would otherwise silently change which branch of the
/// bracket update is taken. It acts only on auxiliary ordinates; an underflowed
/// working value is never accepted as a root of f.
#[inline]
fn scale_preserving_nonzero_sign(value: f64, positive_factor: f64) -> f64 {
    let scaled = value * positive_factor;

    if scaled == 0.0 && value != 0.0 {
        return MIN_SUBNORMAL.copysign(value);
    }

    if scaled.is_infinite() {
        return f64::MAX.copysign(value);
    }

    scaled
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
// The overflow- and NaN-safe forms of the interpolation and switching
// arithmetic live in the shared safeguards above.
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
    // NaN has no usable sign, and `same_nonzero_sign` is false for it, so it
    // must be rejected before the predicate is used to update a bracket.
    if y1.is_nan() || y2.is_nan() || same_nonzero_sign(y1, y2) {
        return f64::NAN;
    }
    let mut f1 = y1;
    let mut f2 = y2;
    let mut side: i32 = 0;
    let mut bisection = true;
    let mut threshold = x2 - x1;
    let mut ymin = 0.0;
    for _ in 0..maxiter {
        // `safe_secant` already returns a point inside [x1, x2], so the separate
        // clamp on the convergence exit is no longer needed.
        let x3 = if bisection {
            safe_midpoint(x1, x2)
        } else {
            safe_secant(x1, y1, x2, y2)
        };
        let epsx = xtol * x3.abs().max(1.0);
        if x2 - x1 <= epsx {
            return x3;
        }
        let y3: f64;
        if bisection {
            y3 = f(x3) - y;
            let ym = safe_midpoint(f1, f2);
            if passes_switching_test(ym, y3, symmetry_factor(f1, f2)) {
                bisection = false;
                threshold = (x2 - x1) * 2.0;
            }
        } else {
            // If rounding makes the proposal coincide with an endpoint, reuse
            // the true residual already stored there.
            if x3 == x1 {
                y3 = f1;
            } else if x3 == x2 {
                y3 = f2;
            } else {
                y3 = f(x3) - y;
            }
            threshold *= 0.5;
            ymin = f1.abs().min(f2.abs());
        }
        if y3.abs() <= epsy {
            return x3;
        }
        // A NaN residual has no usable sign, so the bracket cannot be updated.
        if y3.is_nan() {
            return f64::NAN;
        }
        if same_nonzero_sign(y1, y3) {
            if side == 1 {
                y2 = scale_preserving_nonzero_sign(y2, ab_factor(y3, y1));
            } else if !bisection {
                side = 1;
            }
            (x1, y1, f1) = (x3, y3, y3);
        } else {
            if side == -1 {
                y1 = scale_preserving_nonzero_sign(y1, ab_factor(y3, y2));
            } else if !bisection {
                side = -1;
            }
            (x2, y2, f2) = (x3, y3, y3);
        }
        if !bisection && x2 - x1 > threshold && y3.abs() > 0.5 * ymin {
            bisection = true;
            side = 0;
        }
    }
    f64::NAN
}

#[cfg(test)]
mod safeguard_tests {
    //! Edge-case tests for the shared modAB safeguards.
    //!
    //! The 100-problem benchmark suite never produces a non-finite or overflowing
    //! residual, so the overflow/NaN branches of same_nonzero_sign, safe_midpoint,
    //! safe_secant, symmetry_factor and passes_switching_test are covered here.
    //!
    //! Expected values come from the C# reference in C#/Root/Node.cs and
    //! C#/Root/Solvers/ModABCorr.cs; every language port is held to the same table.

    use super::*;

    fn eq(got: f64, want: f64) -> bool {
        (got.is_nan() && want.is_nan()) || got == want
    }

    #[test]
    fn same_nonzero_sign_edge_cases() {
        assert_eq!(same_nonzero_sign(1.0f64, 2.0f64), true, "case 0");
        assert_eq!(same_nonzero_sign(-1.0f64, -2.0f64), true, "case 1");
        assert_eq!(same_nonzero_sign(1.0f64, -2.0f64), false, "case 2");
        assert_eq!(same_nonzero_sign(-1.0f64, 2.0f64), false, "case 3");
        assert_eq!(same_nonzero_sign(0.0f64, 1.0f64), false, "case 4");
        assert_eq!(same_nonzero_sign(1.0f64, 0.0f64), false, "case 5");
        assert_eq!(same_nonzero_sign(0.0f64, 0.0f64), false, "case 6");
        assert_eq!(same_nonzero_sign(-0.0f64, -1.0f64), false, "case 7");
        assert_eq!(same_nonzero_sign(f64::NAN, 1.0f64), false, "case 8");
        assert_eq!(same_nonzero_sign(1.0f64, f64::NAN), false, "case 9");
        assert_eq!(same_nonzero_sign(f64::NAN, f64::NAN), false, "case 10");
        assert_eq!(same_nonzero_sign(f64::INFINITY, 1.0f64), true, "case 11");
        assert_eq!(same_nonzero_sign(-f64::INFINITY, -1.0f64), true, "case 12");
        assert_eq!(same_nonzero_sign(f64::INFINITY, -f64::INFINITY), false, "case 13");
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

    #[test]
    fn symmetry_factor_edge_cases() {
        assert!(eq(symmetry_factor(-1.0f64, 1.0f64), 1.0f64), "case 0");
        assert!(eq(symmetry_factor(-1.0f64, 3.0f64), 0.5625f64), "case 1");
        assert!(eq(symmetry_factor(-3.0f64, 1.0f64), 0.5625f64), "case 2");
        assert!(eq(symmetry_factor(-1e+308f64, 1e+308f64), 1.0f64), "case 3");
        assert!(eq(symmetry_factor(-1.7e+308f64, 1.0f64), 0.25f64), "case 4");
        assert!(eq(symmetry_factor(f64::INFINITY, -1.0f64), f64::NAN), "case 5");
        assert!(eq(symmetry_factor(-1.0f64, f64::INFINITY), f64::NAN), "case 6");
        assert!(eq(symmetry_factor(f64::NAN, 1.0f64), f64::NAN), "case 7");
        assert!(eq(symmetry_factor(-1e-300f64, 1e+300f64), 0.25f64), "case 8");
    }

    #[test]
    fn passes_switching_test_edge_cases() {
        assert_eq!(passes_switching_test(1.0f64, 1.0f64, 0.5f64), true, "case 0");
        assert_eq!(passes_switching_test(1.0f64, -1.0f64, 0.5f64), false, "case 1");
        assert_eq!(passes_switching_test(1.0f64, 0.5f64, 0.5f64), true, "case 2");
        assert_eq!(passes_switching_test(1.0f64, 0.5f64, 0.1f64), false, "case 3");
        assert_eq!(passes_switching_test(1e+308f64, 1e+308f64, 0.5f64), true, "case 4");
        assert_eq!(passes_switching_test(1.7e+308f64, -1.7e+308f64, 0.5f64), false, "case 5");
        assert_eq!(passes_switching_test(1e+308f64, 1.6e+308f64, 0.5f64), true, "case 6");
        assert_eq!(passes_switching_test(f64::INFINITY, 1.0f64, 0.5f64), false, "case 7");
        assert_eq!(passes_switching_test(1.0f64, f64::INFINITY, 0.5f64), false, "case 8");
        assert_eq!(passes_switching_test(f64::NAN, 1.0f64, 0.5f64), false, "case 9");
        assert_eq!(passes_switching_test(1.0f64, 1.0f64, f64::NAN), false, "case 10");
        assert_eq!(passes_switching_test(0.0f64, 0.0f64, 0.5f64), false, "case 11");
    }

    /// Solver-level regression: this returned +Inf before `safe_midpoint`.
    #[test]
    fn overflowing_bracket() {
        let root = mod_ab_root(|x| x * 1e-308 - 1.2, 1e308, 1.7e308, 0.0, 1e-14, 0.0, 200);
        assert!((root - 1.2e308).abs() <= 1e-13 * 1.2e308, "got {root}");
    }
}
