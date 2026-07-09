//! Benchmark of bracketing root-finders.
//!
//! Compares the modified Anderson–Björck (modAB) algorithm
//! (Ganchovski et al., 2026) against the bracketing methods provided by
//! the `roots` crate (https://github.com/vorot/roots):
//!
//!   * `find_root_brent`             — Brent–Dekker
//!   * `find_root_inverse_quadratic` — Inverse quadratic interpolation
//!   * `find_root_regula_falsi`      — Regula falsi (Illinois variant)
//!   * `find_root_secant`            — Secant
//!

use roots::{
    find_root_brent, find_root_inverse_quadratic, find_root_regula_falsi,
    find_root_secant, SimpleConvergency,
};
use std::cell::Cell;
use std::f64::consts::{E, PI};
use std::time::Instant;

// modAB algorithm lives in src/ModAB.rs
#[allow(non_snake_case)]
mod ModAB;
use ModAB::mod_ab_root;

// ---------------------------------------------------------------------------
// Solver wrappers — uniform signature for the benchmark dispatch table
// ---------------------------------------------------------------------------

const MAX_ITER: usize = 100;

type SolverFn = fn(&dyn Fn(f64) -> f64, f64, f64, f64) -> f64;

fn solve_brent(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let mut conv = SimpleConvergency { eps, max_iter: MAX_ITER };
    find_root_brent(a, b, |x| f(x), &mut conv).unwrap_or(f64::NAN)
}

fn solve_invq(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let mut conv = SimpleConvergency { eps, max_iter: MAX_ITER };
    find_root_inverse_quadratic(a, b, |x| f(x), &mut conv).unwrap_or(f64::NAN)
}

fn solve_regf(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let mut conv = SimpleConvergency { eps, max_iter: MAX_ITER };
    find_root_regula_falsi(a, b, |x| f(x), &mut conv).unwrap_or(f64::NAN)
}

fn solve_secant(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let mut conv = SimpleConvergency { eps, max_iter: MAX_ITER };
    find_root_secant(a, b, |x| f(x), &mut conv).unwrap_or(f64::NAN)
}

fn solve_modab(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    mod_ab_root(|x| f(x), a, b, 0.0, eps, 0.0, MAX_ITER)
}

// ---------------------------------------------------------------------------
// Problem definitions
// ---------------------------------------------------------------------------

struct Problem {
    name: &'static str,
    f: Box<dyn Fn(f64) -> f64>,
    a: f64,
    b: f64,
}

fn make<F>(name: &'static str, f: F, a: f64, b: f64) -> Problem
where
    F: Fn(f64) -> f64 + 'static,
{
    Problem { name, f: Box::new(f), a, b }
}

#[inline] fn pp(x: f64) -> f64 { x + 1.11111 }
#[inline] fn sgn(x: f64) -> f64 {
    if x > 0.0 { 1.0 } else if x < 0.0 { -1.0 } else { 0.0 }
}

fn all_problems() -> Vec<Problem> {
    let mut v: Vec<Problem> = Vec::new();

    // ----- Galdino (problems1) -----
    v.push(make("f01", |x| x.powi(3) - 1.0, 0.5, 1.5));
    v.push(make("f02", |x| x.powi(2) * (x.powi(2) / 3.0 + 2f64.sqrt() * x.sin()) - 3f64.sqrt() / 18.0, 0.1, 1.0));
    v.push(make("f03", |x| 11.0 * x.powi(11) - 1.0, 0.1, 1.0));
    v.push(make("f04", |x| x.powi(3) + 1.0, -1.8, 0.0));
    v.push(make("f05", |x| x.powi(3) - 2.0 * x - 5.0, 2.0, 3.0));
    v.push(make("f06", |x| 2.0 * x * (-5f64).exp() + 1.0 - 2.0 * (-5.0 * x).exp(), 0.0, 1.0));
    v.push(make("f07", |x| 2.0 * x * (-10f64).exp() + 1.0 - 2.0 * (-10.0 * x).exp(), 0.0, 1.0));
    v.push(make("f08", |x| 2.0 * x * (-20f64).exp() + 1.0 - 2.0 * (-20.0 * x).exp(), 0.0, 1.0));
    v.push(make("f09", |x| (1.0 + (1.0_f64 - 5.0).powi(2)) * x * x - (1.0 - 5.0 * x).powi(2), 0.0, 1.0));
    v.push(make("f10", |x| (1.0 + (1.0_f64 - 10.0).powi(2)) * x * x - (1.0 - 10.0 * x).powi(2), 0.0, 1.0));
    v.push(make("f11", |x| (1.0 + (1.0_f64 - 20.0).powi(2)) * x * x - (1.0 - 20.0 * x).powi(2), 0.0, 1.0));
    v.push(make("f12", |x| x * x - (1.0 - x).powi(5), 0.0, 1.0));
    v.push(make("f13", |x| x * x - (1.0 - x).powi(10), 0.0, 1.0));
    v.push(make("f14", |x| x * x - (1.0 - x).powi(20), 0.0, 1.0));
    v.push(make("f15", |x| (1.0 + (1.0_f64 - 5.0).powi(4)) * x - (1.0 - 5.0 * x).powi(4), 0.0, 1.0));
    v.push(make("f16", |x| (1.0 + (1.0_f64 - 10.0).powi(4)) * x - (1.0 - 10.0 * x).powi(4), 0.0, 1.0));
    v.push(make("f17", |x| (1.0 + (1.0_f64 - 20.0).powi(4)) * x - (1.0 - 20.0 * x).powi(4), 0.0, 1.0));
    v.push(make("f18", |x| (-5.0 * x).exp() * (x - 1.0) + x.powi(5), 0.0, 1.0));
    v.push(make("f19", |x| (-10.0 * x).exp() * (x - 1.0) + x.powi(10), 0.0, 1.0));
    v.push(make("f20", |x| (-20.0 * x).exp() * (x - 1.0) + x.powi(20), 0.0, 1.0));
    v.push(make("f21", |x| x * x + (x / 5.0).sin() - 0.25, 0.0, 1.0));
    v.push(make("f22", |x| x * x + (x / 10.0).sin() - 0.25, 0.0, 1.0));
    v.push(make("f23", |x| x * x + (x / 20.0).sin() - 0.25, 0.0, 1.0));
    v.push(make("f24", |x| (x + 2.0) * (x + 1.0) * (x - 3.0).powi(3), 2.6, 4.6));
    v.push(make("f25", |x| (x - 4.0).powi(5) * x.ln(), 3.6, 5.6));
    v.push(make("f26", |x| (x.sin() - x / 4.0).powi(3), 2.0, 4.0));
    v.push(make("f27", |x| {
        let q = pp(x);
        let s = if q < 3.0 { 1.0 } else if q > 3.0 { -1.0 } else { 0.0 };
        (81.0 - q * (108.0 - q * (54.0 - q * (12.0 - q)))) * s
    }, 1.0, 3.0));
    v.push(make("f28", |x| ((x - 7.143).powi(3)).sin(), 7.0, 8.0));
    v.push(make("f29", |x| ((x - 3.0).powi(5)).exp() - 1.0, 2.6, 4.6));
    v.push(make("f30", |x| ((x - 3.0).powi(5)).exp() - (x - 1.0).exp(), 4.0, 5.0));
    v.push(make("f31", |x| PI - 1.0 / x, 0.05, 5.0));
    v.push(make("f32", |x| 4.0 - x.tan(), 0.0, 1.5));
    v.push(make("f33", |x| x.cos() - x.powi(3), 0.0, 4.0));

    // ----- Stage (problems1, continued) -----
    v.push(make("f34", |x| x.cos() - x, -11.0, 9.0));
    v.push(make("f35", |x| {
        let s = if x <= 2.0 / 3.0 { 1.0 } else { -1.0 };
        (x - 2.0 / 3.0).abs().sqrt() * s - 0.1
    }, -11.0, 9.0));
    v.push(make("f36", |x| {
        let s = if x <= 2.0 / 3.0 { 1.0 } else { -1.0 };
        (x - 2.0 / 3.0).abs().powf(0.2) * s
    }, -11.0, 9.0));
    v.push(make("f37", |x| (x - 7.0 / 9.0).powi(3) + (x - 7.0 / 9.0) * 1e-3, -11.0, 9.0));
    v.push(make("f38", |x| if x <= 1.0 / 3.0 { -0.5 } else { 0.5 }, -11.0, 9.0));
    v.push(make("f39", |x| if x <= 1.0 / 3.0 { -1e-3 } else { 1.0 - 1e-3 }, -11.0, 9.0));
    v.push(make("f40", |x| if x == 0.0 { 0.0 } else { 1.0 / (x - 2.0 / 3.0) }, -11.0, 9.0));

    // ----- Swift–Lindfield (problems1, continued) -----
    v.push(make("f41", |x| 2.0 * x * (-5f64).exp() - 2.0 * (-5.0 * x).exp() + 1.0, 0.0, 10.0));
    v.push(make("f42", |x| (x * x - x - 6.0) * (x * x - 3.0 * x + 2.0), 0.0, PI));
    v.push(make("f43", |x| x.powi(3), -1.0, 1.5));
    v.push(make("f44", |x| x.powi(5), -1.0, 1.5));
    v.push(make("f45", |x| x.powi(7), -1.0, 1.5));
    v.push(make("f46", |x| ((-5.0 * x).exp() - x - 0.5) / x.powi(5), 0.09, 0.7));
    v.push(make("f47", |x| 1.0 / x.sqrt() - 2.0 * (5e3 * x.sqrt()).ln() + 0.8, 0.0005, 0.5));
    v.push(make("f48", |x| 1.0 / x.sqrt() - 2.0 * (5e7 * x.sqrt()).ln() + 0.8, 0.0005, 0.5));
    v.push(make("f49", |x| {
        if x <= 0.0 { -x.powi(3) - x - 1.0 } else { x.cbrt() - x - 1.0 }
    }, -1.0, 1.0));
    v.push(make("f50", |x| x.powi(3) - 2.0 * x - x + 3.0, -3.0, 2.0));
    v.push(make("f51", |x| x.ln(), 0.5, 5.0));
    v.push(make("f52", |x| (10.0 - x) * (-10.0 * x).exp() - x.powi(10) + 1.0, 0.5, 8.0));
    v.push(make("f53", |x| x.sin().exp() - x - 1.0, 1.0, 4.0));
    v.push(make("f54", |x| 2.0 * x.sin() - 1.0, 0.1, PI / 3.0));
    v.push(make("f55", |x| (x - 1.0) * (-x).exp(), 0.0, 1.5));
    v.push(make("f56", |x| (x - 1.0).powi(3) - 1.0, 1.5, 3.0));
    v.push(make("f57", |x| (x * x + 7.0 * x - 30.0).exp() - 1.0, 2.6, 3.5));
    v.push(make("f58", |x| x.atan() - 1.0, 1.0, 8.0));
    v.push(make("f59", |x| x.exp() - 2.0 * x - 1.0, 0.2, 3.0));
    v.push(make("f60", |x| (-x).exp() - x - x.sin(), 0.0, 2.0));
    v.push(make("f61", |x| x * x - x.sin().powi(2) - 1.0, -1.0, 2.0));
    v.push(make("f62", |x| x.sin() - x / 2.0, PI / 2.0, PI));

    // ----- Oliveira–Takahashi (problems2) -----
    v.push(make("f63", |x| x * x.exp() - 1.0, -1.0, 1.0));
    v.push(make("f64", |x| (x - 0.1).tan(), -1.0, 1.0));
    v.push(make("f65", |x| x.sin() + 0.5, -1.0, 1.0));
    v.push(make("f66", |x| 4.0 * x.powi(5) + x * x + 1.0, -1.0, 1.0));
    v.push(make("f67", |x| x + x.powi(10) - 1.0, -1.0, 1.0));
    v.push(make("f68", |x| PI.powf(x) - E, -1.0, 1.0));
    v.push(make("f69", |x| (x - 10.0 / 9.0).abs().ln(), -1.0, 1.0));
    v.push(make("f70", |x| 1.0 / 3.0 + x.cbrt() + x.powi(3), -1.0, 1.0));
    v.push(make("f71", |x| (x + 2.0 / 3.0) / (x + 101.0 / 100.0), -1.0, 1.0));
    v.push(make("f72", |x| (x * 1e6 - 1.0).powi(3), -1.0, 1.0));
    v.push(make("f73", |x| x.exp() * (x * 1e6 - 1.0).powi(3), -1.0, 1.0));
    v.push(make("f74", |x| (x - 1.0 / 3.0).powi(2) * (x - 1.0 / 3.0).atan(), -1.0, 1.0));
    v.push(make("f75", |x| {
        let t = 3.0 * x - 1.0;
        sgn(t) * (1.0 - (1.0 - t * t / 81.0).sqrt())
    }, -1.0, 1.0));
    v.push(make("f76", |x| {
        if x > (1.0 - 1e6) / 1e6 { (1.0 + 1e6) / 1e6 } else { -1.0 }
    }, -1.0, 1.0));
    v.push(make("f77", |x| {
        if x != 1.0 / 21.0 { 1.0 / (21.0 * x - 1.0) } else { 0.0 }
    }, -1.0, 1.0));
    v.push(make("f78", |x| x * x / 4.0 + (x / 2.0).ceil() - 0.5, -1.0, 1.0));
    v.push(make("f79", |x| (10.0 * x - 1.0).ceil() + 0.5, -1.0, 1.0));
    v.push(make("f80", |x| x + (x * 1e6).sin() / 10.0 + 1e-3, -1.0, 1.0));
    v.push(make("f81", |x| {
        if x > -1.0 { 1.0 + (1.0 / (x + 1.0)).sin() } else { -1.0 }
    }, -1.0, 1.0));
    v.push(make("f82", |x| 202.0 * x - 2.0 * ((2.0 * x + 1e-2) / 2e-2).floor() - 0.1, -1.0, 1.0));
    v.push(make("f83", |x| {
        let t = 202.0 * x - 2.0 * ((2.0 * x + 1e-2) / 2e-2).floor() - 0.1;
        t * t * t
    }, -1.0, 1.0));

    // ----- SciML (problems3) -----
    v.push(make("f84", |x| (x - 1.0) * (x - 2.0) * (x - 3.0) * (x - 4.0) * (x - 5.0) - 0.05, 0.5, 5.5));
    v.push(make("f85", |x| x.sin() - 0.5 * x - 0.3, -10.0, 10.0));
    v.push(make("f86", |x| x.exp() - 1.0 - x - x * x / 2.0 - 0.005, -2.0, 2.0));
    v.push(make("f87", |x| 1.0 / (x - 0.5) - 2.0 - 0.05, 0.6, 2.0));
    v.push(make("f88", |x| x.ln() - x + 2.0 - 0.05, 0.1, 3.0));
    v.push(make("f89", |x| (20.0 * x).sin() + 0.1 * x - 0.1, -4.0, 5.0));
    v.push(make("f90", |x| x.powi(3) - 2.0 * x.powi(2) + x - 0.025, -1.0, 2.0));
    v.push(make("f91", |x| x * (1.0 / x).sin() - 0.1 - 0.01, 0.01, 1.0));
    v.push(make("f92", |x| x.powi(3) - 0.001, -10.0, 10.0));
    v.push(make("f93", |x| x.powi(7) - 0.001, -10.0, 10.0));
    v
}

// ---------------------------------------------------------------------------
// Counting wrapper — interior mutability lets the closure stay `Fn`
// ---------------------------------------------------------------------------

fn run_counting<S>(solver: S, f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64)
    -> (f64, u64)
where
    S: Fn(&dyn Fn(f64) -> f64, f64, f64, f64) -> f64,
{
    let counter = Cell::new(0u64);
    let counted = |x: f64| {
        counter.set(counter.get() + 1);
        f(x)
    };
    let r = solver(&counted, a, b, eps);
    (r, counter.get())
}

// ---------------------------------------------------------------------------
// Pretty-printing helpers
// ---------------------------------------------------------------------------

const COL_W: usize = 22; // result / f(root) columns
const CNT_W: usize = 8;  // count / time columns

fn fmt_root(v: f64) -> String {
    if v.is_nan() { format!("{:>w$}", "NaN", w = COL_W) }
    else { format!("{:>w$.15}", v, w = COL_W) }
}

fn fmt_fval(v: f64) -> String {
    if v.is_nan() { format!("{:>w$}", "NaN", w = COL_W) }
    else { format!("{:>w$.6e}", v, w = COL_W) }
}

fn median(mut xs: Vec<f64>) -> f64 {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = xs.len();
    if n == 0 { return f64::NAN; }
    if n % 2 == 0 { (xs[n / 2 - 1] + xs[n / 2]) / 2.0 } else { xs[n / 2] }
}

/// Print avg/median/min/max + factor-vs-modAB rows for one metric.
/// `data[i][j]` is the value for problem i, solver j (NaN for errors).
/// modAB is assumed to be the LAST solver column (the baseline).
fn print_stats_block(data: &[Vec<f64>], n_solvers: usize, col_w: usize,
                     integer_vals: bool) {
    let modab_idx = n_solvers - 1;
    let mut avg = vec![f64::NAN; n_solvers];
    let mut med = vec![f64::NAN; n_solvers];
    let mut mn  = vec![f64::NAN; n_solvers];
    let mut mx  = vec![f64::NAN; n_solvers];

    for j in 0..n_solvers {
        let vals: Vec<f64> = data.iter()
            .map(|row| row[j]).filter(|v| !v.is_nan()).collect();
        if vals.is_empty() { continue; }
        avg[j] = vals.iter().sum::<f64>() / vals.len() as f64;
        med[j] = median(vals.clone());
        mn[j]  = vals.iter().cloned().fold(f64::INFINITY,    f64::min);
        mx[j]  = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    }

    let cell = |v: f64| -> String {
        if v.is_nan() {
            format!("{:>w$}| ", "N/A", w = col_w)
        } else if integer_vals {
            format!("{:>w$}| ", v.round() as i64, w = col_w)
        } else {
            format!("{:>w$.4}| ", v, w = col_w)
        }
    };

    for (label, row) in [("AVG", &avg), ("MEDIAN", &med),
                         ("MIN", &mn), ("MAX", &mx)] {
        let mut line = format!("{:>6}| ", label);
        for &v in row.iter() { line.push_str(&cell(v)); }
        println!("{line}");
    }

    // factor relative to modAB's average
    let base = avg[modab_idx];
    let mut line = format!("{:>6}| ", "FACTOR");
    for j in 0..n_solvers {
        if avg[j].is_nan() || base.is_nan() || base == 0.0 {
            line.push_str(&format!("{:>w$}| ", "N/A", w = col_w));
        } else {
            let f = avg[j] / base;
            line.push_str(&format!("{:>w$.3}x| ", f, w = col_w - 1));
        }
    }
    println!("{line}");
    println!();
}

// ---------------------------------------------------------------------------
// Main benchmark
// ---------------------------------------------------------------------------

fn main() {
    let problems = all_problems();
    // modAB MUST be last — the FACTOR row uses it as baseline.
    let solvers: Vec<(&str, SolverFn)> = vec![
        ("brent",   solve_brent),
        ("invq",    solve_invq),
        ("regf",    solve_regf),
        ("secant",  solve_secant),
        (" modAB",  solve_modab),
    ];
    let n_solv = solvers.len();
    let eps = 1e-14;
    let iterations = 200;

    // ----- Roots -----
    println!("Results (root values)");
    println!("===");
    print!("{:>4}| ", "Func");
    for (name, _) in &solvers { print!("{:>w$}| ", name, w = COL_W); }
    println!();
    println!("{}", "-".repeat(6 + (COL_W + 2) * n_solv));

    for p in &problems {
        print!("{:>4}| ", p.name);
        for (_, s) in &solvers {
            let (root, _) = run_counting(s, &*p.f, p.a, p.b, eps);
            print!("{}| ", fmt_root(root));
        }
        println!();
    }
    println!();

    // ----- f(root) -----
    println!("Function values  f(root)");
    println!("===");
    print!("{:>4}| ", "Func");
    for (name, _) in &solvers { print!("{:>w$}| ", name, w = COL_W); }
    println!();
    println!("{}", "-".repeat(6 + (COL_W + 2) * n_solv));

    let mut fval_data: Vec<Vec<f64>> = Vec::with_capacity(problems.len());
    for p in &problems {
        print!("{:>4}| ", p.name);
        let mut row = Vec::with_capacity(n_solv);
        for (_, s) in &solvers {
            let (root, _) = run_counting(s, &*p.f, p.a, p.b, eps);
            let fv = if root.is_nan() { f64::NAN } else { (p.f)(root).abs() };
            row.push(fv);
            print!("{}| ", fmt_fval(fv));
        }
        fval_data.push(row);
        println!();
    }
    println!();

    // ----- Function evaluations -----
    println!("Function evaluations");
    println!("===");
    print!("{:>4}| ", "Func");
    for (name, _) in &solvers { print!("{:>w$}| ", name, w = CNT_W); }
    println!();
    println!("{}", "-".repeat(6 + (CNT_W + 2) * n_solv));

    let mut count_data: Vec<Vec<f64>> = Vec::with_capacity(problems.len());
    let mut totals = vec![0u64; n_solv];
    for p in &problems {
        print!("{:>4}| ", p.name);
        let mut row = Vec::with_capacity(n_solv);
        for (j, (_, s)) in solvers.iter().enumerate() {
            let (_, c) = run_counting(s, &*p.f, p.a, p.b, eps);
            totals[j] += c;
            row.push(c as f64);
            print!("{:>w$}| ", c, w = CNT_W);
        }
        count_data.push(row);
        println!();
    }
    print!("{:>4}| ", "SUM");
    for t in &totals { print!("{:>w$}| ", t, w = CNT_W); }
    println!();
    print_stats_block(&count_data, n_solv, CNT_W, true);

    // ----- Timings -----
    println!("Execution times  (ms per problem, {iterations} iterations)");
    println!("===");
    print!("{:>4}| ", "Func");
    for (name, _) in &solvers { print!("{:>w$}| ", name, w = CNT_W); }
    println!();
    println!("{}", "-".repeat(6 + (CNT_W + 2) * n_solv));

    let mut time_data: Vec<Vec<f64>> = Vec::with_capacity(problems.len());
    let mut total_t = vec![0.0f64; n_solv];
    for p in &problems {
        print!("{:>4}| ", p.name);
        let mut row = Vec::with_capacity(n_solv);
        for (j, (_, s)) in solvers.iter().enumerate() {
            let f_ref: &dyn Fn(f64) -> f64 = &*p.f;
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = s(f_ref, p.a, p.b, eps);
            }
            let ms = start.elapsed().as_secs_f64() * 1000.0;
            total_t[j] += ms;
            row.push(ms);
            print!("{:>w$.2}| ", ms, w = CNT_W);
        }
        time_data.push(row);
        println!();
    }
    print!("{:>4}| ", "SUM");
    for t in &total_t { print!("{:>w$.2}| ", t, w = CNT_W); }
    println!();
    print_stats_block(&time_data, n_solv, CNT_W, false);
}