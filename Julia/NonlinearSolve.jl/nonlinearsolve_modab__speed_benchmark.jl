using Printf
using NonlinearSolve.BracketingNonlinearSolve
using BenchmarkTools

# Function-call counting wrapper
mutable struct CountedFunc{F} <: Function
    f::F
    count::Int
end
CountedFunc(f) = CountedFunc(f, 0)
(cf::CountedFunc)(x) = (cf.count += 1; cf.f(x))
reset!(cf::CountedFunc) = (cf.count = 0; cf)

# NonlinearSolve.jl - benchmark solve() directly to avoid wrapper overhead
# Problem definition
struct Problem
    name::String
    f::Function
    a::Float64
    b::Float64
    value::Float64
end
Problem(name, f, a, b) = Problem(name, f, Float64(a), Float64(b), 0.0)

P(x) = x + 1.11111

# Test problems
const problems1 = [
    Problem("f01", x -> x^3 - 1, 0.5, 1.5),
    Problem("f02", x -> x^2 * (x^2 / 3 + sqrt(2) * sin(x)) - sqrt(3) / 18, 0.1, 1),
    Problem("f03", x -> 11 * x^11 - 1, 0.1, 1),
    Problem("f04", x -> x^3 + 1, -1.8, 0),
    Problem("f05", x -> x^3 - 2x - 5, 2, 3),
    Problem("f06", x -> 2x * exp(-5) + 1 - 2exp(-5x), 0, 1),
    Problem("f07", x -> 2x * exp(-10) + 1 - 2exp(-10x), 0, 1),
    Problem("f08", x -> 2x * exp(-20) + 1 - 2exp(-20x), 0, 1),
    Problem("f09", x -> (1 + (1 - 5)^2) * x^2 - (1 - 5x)^2, 0, 1),
    Problem("f10", x -> (1 + (1 - 10)^2) * x^2 - (1 - 10x)^2, 0, 1),
    Problem("f11", x -> (1 + (1 - 20)^2) * x^2 - (1 - 20x)^2, 0, 1),
    Problem("f12", x -> x^2 - (1 - x)^5, 0, 1),
    Problem("f13", x -> x^2 - (1 - x)^10, 0, 1),
    Problem("f14", x -> x^2 - (1 - x)^20, 0, 1),
    Problem("f15", x -> (1 + (1 - 5)^4) * x - (1 - 5x)^4, 0, 1),
    Problem("f16", x -> (1 + (1 - 10)^4) * x - (1 - 10x)^4, 0, 1),
    Problem("f17", x -> (1 + (1 - 20)^4) * x - (1 - 20x)^4, 0, 1),
    Problem("f18", x -> exp(-5x) * (x - 1) + x^5, 0, 1),
    Problem("f19", x -> exp(-10x) * (x - 1) + x^10, 0, 1),
    Problem("f20", x -> exp(-20x) * (x - 1) + x^20, 0, 1),
    Problem("f21", x -> x^2 + sin(x / 5) - 1 / 4, 0, 1),
    Problem("f22", x -> x^2 + sin(x / 10) - 1 / 4, 0, 1),
    Problem("f23", x -> x^2 + sin(x / 20) - 1 / 4, 0, 1),
    Problem("f24", x -> (x + 2) * (x + 1) * (x - 3)^3, 2.6, 4.6),
    Problem("f25", x -> (x - 4)^5 * log(x), 3.6, 5.6),
    Problem("f26", x -> (sin(x) - x / 4)^3, 2, 4),
    Problem("f27", x -> (81 - P(x) * (108 - P(x) * (54 - P(x) * (12 - P(x))))) * sign(P(x) - 3), 1, 3),
    Problem("f28", x -> sin((x - 7.143)^3), 7, 8),
    Problem("f29", x -> exp((x - 3)^5) - 1, 2.6, 4.6),
    Problem("f30", x -> exp((x - 3)^5) - exp(x - 1), 4, 5),
    Problem("f31", x -> π - 1 / x, 0.05, 5),
    Problem("f32", x -> 4 - tan(x), 0, 1.5),
    Problem("f33", x -> cos(x) - x^3, 0, 4),
    Problem("f34", x -> cos(x) - x, -11, 9),
    Problem("f35", x -> sqrt(abs(x - 2 / 3)) * (x <= 2 / 3 ? 1 : -1) - 0.1, -11, 9),
    Problem("f36", x -> abs(x - 2 / 3)^0.2 * (x <= 2 / 3 ? 1 : -1), -11, 9),
    Problem("f37", x -> (x - 7 / 9)^3 + (x - 7 / 9) * 1e-3, -11, 9),
    Problem("f38", x -> x <= 1 / 3 ? -0.5 : 0.5, -11, 9),
    Problem("f39", x -> x <= 1 / 3 ? -1e-3 : 1 - 1e-3, -11, 9),
    Problem("f40", x -> x == 0 ? 0.0 : 1 / (x - 2 / 3), -11, 9),
    Problem("f41", x -> 2x * exp(-5) - 2exp(-5x) + 1, 0, 10),
    Problem("f42", x -> (x^2 - x - 6) * (x^2 - 3x + 2), 0, π),
    Problem("f43", x -> x^3, -1, 1.5),
    Problem("f44", x -> x^5, -1, 1.5),
    Problem("f45", x -> x^7, -1, 1.5),
    Problem("f46", x -> (exp(-5x) - x - 0.5) / x^5, 0.09, 0.7),
    Problem("f47", x -> 1 / sqrt(x) - 2log(5e3 * sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f48", x -> 1 / sqrt(x) - 2log(5e7 * sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f49", x -> x <= 0 ? (-x^3 - x - 1) : (x^(1 / 3) - x - 1), -1, 1),
    Problem("f50", x -> x^3 - 2x - x + 3, -3, 2),
    Problem("f51", x -> log(x), 0.5, 5),
    Problem("f52", x -> (10 - x) * exp(-10x) - x^10 + 1, 0.5, 8),
    Problem("f53", x -> exp(sin(x)) - x - 1, 1.0, 4),
    Problem("f54", x -> 2sin(x) - 1, 0.1, π / 3),
    Problem("f55", x -> (x - 1) * exp(-x), 0.0, 1.5),
    Problem("f56", x -> (x - 1)^3 - 1, 1.5, 3),
    Problem("f57", x -> exp(x^2 + 7x - 30) - 1, 2.6, 3.5),
    Problem("f58", x -> atan(x) - 1, 1.0, 8),
    Problem("f59", x -> exp(x) - 2x - 1, 0.2, 3),
    Problem("f60", x -> exp(-x) - x - sin(x), 0.0, 2),
    Problem("f61", x -> x^2 - sin(x)^2 - 1, -1, 2),
    Problem("f62", x -> sin(x) - x / 2, π / 2, π),
]

const problems2 = [
    Problem("f63", x -> x * exp(x) - 1, -1, 1),
    Problem("f64", x -> tan(x - 1 / 10), -1, 1),
    Problem("f65", x -> sin(x) + 0.5, -1, 1),
    Problem("f66", x -> 4x^5 + x * x + 1, -1, 1),
    Problem("f67", x -> x + x^10 - 1, -1, 1),
    Problem("f68", x -> π^x - ℯ, -1, 1),
    Problem("f69", x -> log(abs(x - 10 / 9)), -1, 1),
    Problem("f70", x -> 1 / 3 + sign(x) * abs(x)^(1 / 3) + x^3, -1, 1),
    Problem("f71", x -> (x + 2 / 3) / (x + 101 / 100), -1, 1),
    Problem("f72", x -> (x * 1e6 - 1)^3, -1, 1),
    Problem("f73", x -> exp(x) * (x * 1e6 - 1)^3, -1, 1),
    Problem("f74", x -> (x - 1 / 3)^2 * atan(x - 1 / 3), -1, 1),
    Problem("f75", x -> sign(3x - 1) * (1 - sqrt(1 - (3x - 1)^2 / 81)), -1, 1),
    Problem("f76", x -> x > (1 - 1e6) / 1e6 ? (1 + 1e6) / 1e6 : -1.0, -1, 1),
    Problem("f77", x -> x != 1 / 21 ? 1 / (21x - 1) : 0.0, -1, 1),
    Problem("f78", x -> x * x / 4 + ceil(x / 2) - 0.5, -1, 1),
    Problem("f79", x -> ceil(10x - 1) + 0.5, -1, 1),
    Problem("f80", x -> x + sin(x * 1e6) / 10 + 1e-3, -1, 1),
    Problem("f81", x -> x > -1 ? 1 + sin(1 / (x + 1)) : -1.0, -1, 1),
    Problem("f82", x -> 202x - 2floor((2x + 1e-2) / 2e-2) - 0.1, -1, 1),
    Problem("f83", x -> (202x - 2floor((2x + 1e-2) / 2e-2) - 0.1)^3, -1, 1),
]

const problems3 = [
    Problem("f84", x -> (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05, 0.5, 5.5),
    Problem("f85", x -> sin(x) - 0.5x - 0.3, -10.0, 10.0),
    Problem("f86", x -> exp(x) - 1 - x - x * x / 2 - 0.005, -2.0, 2.0),
    Problem("f87", x -> 1 / (x - 0.5) - 2 - 0.05, 0.6, 2.0),
    Problem("f88", x -> log(x) - x + 2 - 0.05, 0.1, 3.0),
    Problem("f89", x -> sin(20x) + 0.1x - 0.1, -4.0, 5.0),
    Problem("f90", x -> x^3 - 2x^2 + x - 0.025, -1.0, 2.0),
    Problem("f91", x -> x * sin(1 / x) - 0.1 - 0.01, 0.01, 1.0),
    Problem("f92", x -> x^3 - 0.001, -10, 10),
    Problem("f93", x -> x^7 - 0.001, -10, 10),
]

const all_problems = vcat(problems1, problems2, problems3)

# Algorithms and formatting

# Format nanoseconds to a readable string
function fmt_time(ns::Float64)
    if isnan(ns)
        return "ERR"
    elseif ns < 1_000
        return @sprintf("%.0f ns", ns)
    elseif ns < 1_000_000
        return @sprintf("%.1f μs", ns / 1_000)
    else
        return @sprintf("%.2f ms", ns / 1_000_000)
    end
end

# Algorithms to benchmark (solve() directly, no wrapper overhead)
const algorithms = [
    (" bisect", Bisection()),
    ("  brent", Brent()),
    (" ridder", Ridder()),
    ("alefeld", Alefeld()),
    ("    ITP", ITP()),
    ("  modAB", ModAB()),
]

# Benchmark runner
function run_benchmark()
    eps = 1e-14
    col_w = 10

    # Warmup: run each algorithm once to trigger compilation
    println("Warming up...")
    p1 = all_problems[1]
    g1 = p1.value != 0 ? x -> p1.f(x) - p1.value : p1.f
    prob1 = IntervalNonlinearProblem((x, p) -> g1(x), (p1.a, p1.b))
    for (name, alg) in algorithms
        try
            solve(prob1, alg; abstol=eps, maxiters=200)
        catch
        end
    end
    println()

    # --- Timing: benchmark solve() only, problem constructed outside the loop ---
    println("Timing (median of @belapsed, solve() only)")
    header = lpad("Func", 4) * "; " * join([lpad(name, col_w) for (name, _) in algorithms], "; ")
    println(header)
    total_time = zeros(Float64, length(algorithms))
    for p in all_problems
        g = p.value != 0 ? x -> p.f(x) - p.value : p.f
        prob = IntervalNonlinearProblem((x, pp) -> g(x), (min(p.a, p.b), max(p.a, p.b)))
        line = lpad(p.name, 4) * "; "
        for (j, (name, alg)) in enumerate(algorithms)
            try
                t = @belapsed solve($prob, $alg; abstol=$eps, maxiters=200)
                t_ns = t * 1e9
                total_time[j] += t_ns
                line *= lpad(fmt_time(t_ns), col_w) * "; "
            catch
                total_time[j] += NaN
                line *= lpad("ERR", col_w) * "; "
            end
        end
        println(line)
    end

    # Print totals
    line = lpad("SUM", 4) * "; "
    for t in total_time
        line *= lpad(fmt_time(t), col_w) * "; "
    end
    println(line)
    println()
end

run_benchmark()