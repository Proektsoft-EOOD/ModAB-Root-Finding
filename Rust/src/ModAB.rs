//! Modified Anderson–Björck root-finding algorithm.
//!
//! Reference: Ganchovski, N.; Smith, O.; Rackauckas, C.; Tomov, L.; Traykov, A.
//! Improvements to the Modified Anderson–Björck (modAB) Root-Finding Algorithm.
//! Algorithms 2026, 19, 332. https://doi.org/10.3390/a19050332

/// Find a root of `f(x) = y` on `[x1, x2]` using the modified Anderson–Björck
/// method. `f` must be continuous and `sign(f(x1) - y) != sign(f(x2) - y)`.
/// Returns NaN if no root is found within `maxiter` iterations.
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
    if (y1 > 0.0) == (y2 > 0.0) {
        return f64::NAN;
    }
    let mut side: i32 = 0;
    let mut bisection = true;
    let mut threshold = x2 - x1;
    for _ in 0..maxiter {
        let mut x3 = if bisection {
            (x1 + x2) * 0.5
        } else {
            (x1 * y2 - y1 * x2) / (y2 - y1)
        };
        let epsx = xtol * x3.abs().max(1.0);
        if x2 - x1 <= epsx {
            return x3;
        }
        let y3: f64;
        if bisection {
            y3 = f(x3) - y;
            let ym = (y1 + y2) * 0.5;
            let dy = y2 - y1;
            let r = 1.0 - (ym / dy).abs();
            let k = r * r;
            if (ym - y3).abs() < k * (y3.abs() + ym.abs()) {
                bisection = false;
                threshold = (x2 - x1) * 8.0;
            }
        } else {
            if x3 <= x1 {
                x3 = x1;
                y3 = y1;
            } else if x3 >= x2 {
                x3 = x2;
                y3 = y2;
            } else {
                y3 = f(x3) - y;
            }
            threshold *= 0.5;
        }
        if y3.abs() <= epsy {
            return x3;
        }
        if (y1 > 0.0) == (y3 > 0.0) {
            if side == 1 {
                let m = 1.0 - y3 / y1;
                y2 = if m > 0.0 { y2 * m } else { y2 * 0.5 };
            } else if !bisection {
                side = 1;
            }
            x1 = x3;
            y1 = y3;
        } else {
            if side == -1 {
                let m = 1.0 - y3 / y2;
                y1 = if m > 0.0 { y1 * m } else { y1 * 0.5 };
            } else if !bisection {
                side = -1;
            }
            x2 = x3;
            y2 = y3;
        }
        if x2 - x1 > threshold {
            bisection = true;
            side = 0;
        }
    }
    f64::NAN
}
