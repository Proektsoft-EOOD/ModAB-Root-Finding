import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleUnaryOperator;

public class ModABTest {
    static class Problem {
        final String name;
        final DoubleUnaryOperator f;
        final double a;
        final double b;

        Problem(String name, DoubleUnaryOperator f, double a, double b) {
            this.name = name;
            this.f = f;
            this.a = a;
            this.b = b;
        }
    }

    static class CountedFunc implements DoubleUnaryOperator {
        private final DoubleUnaryOperator f;
        int count = 0;

        CountedFunc(DoubleUnaryOperator f) {
            this.f = f;
        }

        @Override
        public double applyAsDouble(double x) {
            count++;
            return f.applyAsDouble(x);
        }
    }

    private static double P(double x) {
        return x + 1.11111;
    }

    private static List<Problem> createProblems() {
        List<Problem> problems = new ArrayList<>();

        // Sergio Galdino. A family of regula falsi root-finding methods
        problems.add(new Problem("f01", x -> Math.pow(x, 3) - 1, 0.5, 1.5));
        problems.add(new Problem("f02", x -> x * x * (x * x / 3 + Math.sqrt(2) * Math.sin(x)) - Math.sqrt(3) / 18, 0.1, 1));
        problems.add(new Problem("f03", x -> 11 * Math.pow(x, 11) - 1, 0.1, 1));
        problems.add(new Problem("f04", x -> Math.pow(x, 3) + 1, -1.8, 0));
        problems.add(new Problem("f05", x -> Math.pow(x, 3) - 2 * x - 5, 2, 3));
        problems.add(new Problem("f06", x -> 2 * x * Math.exp(-5) + 1 - 2 * Math.exp(-5 * x), 0, 1));
        problems.add(new Problem("f07", x -> 2 * x * Math.exp(-10) + 1 - 2 * Math.exp(-10 * x), 0, 1));
        problems.add(new Problem("f08", x -> 2 * x * Math.exp(-20) + 1 - 2 * Math.exp(-20 * x), 0, 1));
        problems.add(new Problem("f09", x -> (1 + Math.pow(1 - 5, 2)) * x * x - Math.pow(1 - 5 * x, 2), 0, 1));
        problems.add(new Problem("f10", x -> (1 + Math.pow(1 - 10, 2)) * x * x - Math.pow(1 - 10 * x, 2), 0, 1));
        problems.add(new Problem("f11", x -> (1 + Math.pow(1 - 20, 2)) * x * x - Math.pow(1 - 20 * x, 2), 0, 1));
        problems.add(new Problem("f12", x -> x * x - Math.pow(1 - x, 5), 0, 1));
        problems.add(new Problem("f13", x -> x * x - Math.pow(1 - x, 10), 0, 1));
        problems.add(new Problem("f14", x -> x * x - Math.pow(1 - x, 20), 0, 1));
        problems.add(new Problem("f15", x -> (1 + Math.pow(1 - 5, 4)) * x - Math.pow(1 - 5 * x, 4), 0, 1));
        problems.add(new Problem("f16", x -> (1 + Math.pow(1 - 10, 4)) * x - Math.pow(1 - 10 * x, 4), 0, 1));
        problems.add(new Problem("f17", x -> (1 + Math.pow(1 - 20, 4)) * x - Math.pow(1 - 20 * x, 4), 0, 1));
        problems.add(new Problem("f18", x -> Math.exp(-5 * x) * (x - 1) + Math.pow(x, 5), 0, 1));
        problems.add(new Problem("f19", x -> Math.exp(-10 * x) * (x - 1) + Math.pow(x, 10), 0, 1));
        problems.add(new Problem("f20", x -> Math.exp(-20 * x) * (x - 1) + Math.pow(x, 20), 0, 1));
        problems.add(new Problem("f21", x -> x * x + Math.sin(x / 5) - 0.25, 0, 1));
        problems.add(new Problem("f22", x -> x * x + Math.sin(x / 10) - 0.25, 0, 1));
        problems.add(new Problem("f23", x -> x * x + Math.sin(x / 20) - 0.25, 0, 1));
        problems.add(new Problem("f24", x -> (x + 2) * (x + 1) * Math.pow(x - 3, 3), 2.6, 4.6));
        problems.add(new Problem("f25", x -> Math.pow(x - 4, 5) * Math.log(x), 3.6, 5.6));
        problems.add(new Problem("f26", x -> Math.pow(Math.sin(x) - x / 4, 3), 2, 4));
        problems.add(new Problem("f27", x -> {
            double p = P(x);
            double val = 81 - p * (108 - p * (54 - p * (12 - p)));
            return val * (p < 3 ? 1 : (p > 3 ? -1 : 0));
        }, 1, 3));
        problems.add(new Problem("f28", x -> Math.sin(Math.pow(x - 7.143, 3)), 7, 8));
        problems.add(new Problem("f29", x -> Math.exp(Math.pow(x - 3, 5)) - 1, 2.6, 4.6));
        problems.add(new Problem("f30", x -> Math.exp(Math.pow(x - 3, 5)) - Math.exp(x - 1), 4, 5));
        problems.add(new Problem("f31", x -> Math.PI - 1 / x, 0.05, 5));
        problems.add(new Problem("f32", x -> 4 - Math.tan(x), 0, 1.5));
        problems.add(new Problem("f33", x -> Math.cos(x) - Math.pow(x, 3), 0, 4));

        // Steven A. Stage. Comments on An Improvement to the Brent's Method
        problems.add(new Problem("f34", x -> Math.cos(x) - x, -11, 9));
        problems.add(new Problem("f35", x -> Math.sqrt(Math.abs(x - 2.0 / 3)) * (x <= 2.0 / 3 ? 1 : -1) - 0.1, -11, 9));
        problems.add(new Problem("f36", x -> Math.pow(Math.abs(x - 2.0 / 3), 0.2) * (x <= 2.0 / 3 ? 1 : -1), -11, 9));
        problems.add(new Problem("f37", x -> Math.pow(x - 7.0 / 9, 3) + (x - 7.0 / 9) * 1e-3, -11, 9));
        problems.add(new Problem("f38", x -> x <= 1.0 / 3 ? -0.5 : 0.5, -11, 9));
        problems.add(new Problem("f39", x -> x <= 1.0 / 3 ? -1e-3 : 1 - 1e-3, -11, 9));
        problems.add(new Problem("f40", x -> x == 0 ? 0 : 1 / (x - 2.0 / 3), -11, 9));

        // A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method
        problems.add(new Problem("f41", x -> 2 * x * Math.exp(-5) - 2 * Math.exp(-5 * x) + 1, 0, 10));
        problems.add(new Problem("f42", x -> (x * x - x - 6) * (x * x - 3 * x + 2), 0, Math.PI));
        problems.add(new Problem("f43", x -> Math.pow(x, 3), -1, 1.5));
        problems.add(new Problem("f44", x -> Math.pow(x, 5), -1, 1.5));
        problems.add(new Problem("f45", x -> Math.pow(x, 7), -1, 1.5));
        problems.add(new Problem("f46", x -> (Math.exp(-5 * x) - x - 0.5) / Math.pow(x, 5), 0.09, 0.7));
        problems.add(new Problem("f47", x -> 1 / Math.sqrt(x) - 2 * Math.log(5e3 * Math.sqrt(x)) + 0.8, 0.0005, 0.5));
        problems.add(new Problem("f48", x -> 1 / Math.sqrt(x) - 2 * Math.log(5e7 * Math.sqrt(x)) + 0.8, 0.0005, 0.5));
        problems.add(new Problem("f49", x -> x <= 0 ? (-Math.pow(x, 3) - x - 1) : (Math.pow(x, 1.0 / 3) - x - 1), -1, 1));
        problems.add(new Problem("f50", x -> Math.pow(x, 3) - 2 * x - x + 3, -3, 2));
        problems.add(new Problem("f51", x -> Math.log(x), 0.5, 5));
        problems.add(new Problem("f52", x -> (10 - x) * Math.exp(-10 * x) - Math.pow(x, 10) + 1, 0.5, 8));
        problems.add(new Problem("f53", x -> Math.exp(Math.sin(x)) - x - 1, 1.0, 4));
        problems.add(new Problem("f54", x -> 2 * Math.sin(x) - 1, 0.1, Math.PI / 3));
        problems.add(new Problem("f55", x -> (x - 1) * Math.exp(-x), 0.0, 1.5));
        problems.add(new Problem("f56", x -> Math.pow(x - 1, 3) - 1, 1.5, 3));
        problems.add(new Problem("f57", x -> Math.exp(x * x + 7 * x - 30) - 1, 2.6, 3.5));
        problems.add(new Problem("f58", x -> Math.atan(x) - 1, 1.0, 8));
        problems.add(new Problem("f59", x -> Math.exp(x) - 2 * x - 1, 0.2, 3));
        problems.add(new Problem("f60", x -> Math.exp(-x) - x - Math.sin(x), 0.0, 2));
        problems.add(new Problem("f61", x -> x * x - Math.pow(Math.sin(x), 2) - 1, -1, 2));
        problems.add(new Problem("f62", x -> Math.sin(x) - x / 2, Math.PI / 2, Math.PI));

        // Oliveira I. F. D., Takahashi R. H. C.
        // An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality
        problems.add(new Problem("f63", x -> x * Math.exp(x) - 1, -1, 1));
        problems.add(new Problem("f64", x -> Math.tan(x - 0.1), -1, 1));
        problems.add(new Problem("f65", x -> Math.sin(x) + 0.5, -1, 1));
        problems.add(new Problem("f66", x -> 4 * Math.pow(x, 5) + x * x + 1, -1, 1));
        problems.add(new Problem("f67", x -> x + Math.pow(x, 10) - 1, -1, 1));
        problems.add(new Problem("f68", x -> Math.pow(Math.PI, x) - Math.E, -1, 1));
        problems.add(new Problem("f69", x -> Math.log(Math.abs(x - 10.0 / 9)), -1, 1));
        problems.add(new Problem("f70", x -> 1.0 / 3 + Math.signum(x) * Math.pow(Math.abs(x), 1.0 / 3) + Math.pow(x, 3), -1, 1));
        problems.add(new Problem("f71", x -> (x + 2.0 / 3) / (x + 101.0 / 100), -1, 1));
        problems.add(new Problem("f72", x -> Math.pow(x * 1e6 - 1, 3), -1, 1));
        problems.add(new Problem("f73", x -> Math.exp(x) * Math.pow(x * 1e6 - 1, 3), -1, 1));
        problems.add(new Problem("f74", x -> Math.pow(x - 1.0 / 3, 2) * Math.atan(x - 1.0 / 3), -1, 1));
        problems.add(new Problem("f75", x -> {
            double v = 3 * x - 1;
            return Math.signum(v) * (1 - Math.sqrt(1 - v * v / 81));
        }, -1, 1));
        problems.add(new Problem("f76", x -> x > (1 - 1e6) / 1e6 ? (1 + 1e6) / 1e6 : -1, -1, 1));
        problems.add(new Problem("f77", x -> x != 1.0 / 21 ? 1 / (21 * x - 1) : 0, -1, 1));
        problems.add(new Problem("f78", x -> x * x / 4 + Math.ceil(x / 2) - 0.5, -1, 1));
        problems.add(new Problem("f79", x -> Math.ceil(10 * x - 1) + 0.5, -1, 1));
        problems.add(new Problem("f80", x -> x + Math.sin(x * 1e6) / 10 + 1e-3, -1, 1));
        problems.add(new Problem("f81", x -> x > -1 ? 1 + Math.sin(1 / (x + 1)) : -1, -1, 1));
        problems.add(new Problem("f82", x -> 202 * x - 2 * Math.floor((2 * x + 1e-2) / 2e-2) - 0.1, -1, 1));
        problems.add(new Problem("f83", x -> Math.pow(202 * x - 2 * Math.floor((2 * x + 1e-2) / 2e-2) - 0.1, 3), -1, 1));

        // SciML project benchmarks suite
        problems.add(new Problem("f84", x -> (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05, 0.5, 5.5));
        problems.add(new Problem("f85", x -> Math.sin(x) - 0.5 * x - 0.3, -10.0, 10.0));
        problems.add(new Problem("f86", x -> Math.exp(x) - 1 - x - x * x / 2 - 0.005, -2.0, 2.0));
        problems.add(new Problem("f87", x -> 1 / (x - 0.5) - 2 - 0.05, 0.6, 2.0));
        problems.add(new Problem("f88", x -> Math.log(x) - x + 2 - 0.05, 0.1, 3.0));
        problems.add(new Problem("f89", x -> Math.sin(20 * x) + 0.1 * x - 0.1, -4.0, 5.0));
        problems.add(new Problem("f90", x -> Math.pow(x, 3) - 2 * x * x + x - 0.025, -1.0, 2.0));
        problems.add(new Problem("f91", x -> x * Math.sin(1 / x) - 0.1 - 0.01, 0.01, 1.0));
        problems.add(new Problem("f92", x -> Math.pow(x, 3) - 0.001, -10, 10));

        return problems;
    }

    public static void runTests() {
        double eps = 1e-14;
        int passed = 0;
        int failed = 0;
        int totalEvals = 0;

        List<Problem> allProblems = createProblems();

        System.out.println("ModAB Root-Finding Test Results");
        System.out.println("======================================================================");
        System.out.printf("%4s | %22s | %15s | %6s | Status%n", "Func", "Root", "f(root)", "Evals");
        System.out.println("----------------------------------------------------------------------");

        for (Problem p : allProblems) {
            CountedFunc cf = new CountedFunc(p.f);
            try {
                double root = ModAB.modABRoot(cf, p.a, p.b, 0, eps, 0.0, 200);
                double fval = p.f.applyAsDouble(root);
                String status;

                if (Double.isNaN(root)) {
                    status = "FAIL";
                    failed++;
                } else if (Math.abs(fval) < 1e-10) {
                    status = "PASS";
                    passed++;
                } else {
                    status = Math.abs(fval) < 1e-6 ? "PASS" : "WEAK";
                    passed++;
                }

                totalEvals += cf.count;
                System.out.printf("%4s | %22.15g | %15.6e | %6d | %s%n",
                        p.name, root, fval, cf.count, status);
            } catch (Exception e) {
                failed++;
                totalEvals += cf.count;
                System.out.printf("%4s | %22s | %15s | %6d | FAIL (%s)%n",
                        p.name, "ERROR", "N/A", cf.count, e.getMessage());
            }
        }

        System.out.println("----------------------------------------------------------------------");
        System.out.printf("Total: %d tests, %d passed, %d failed%n", passed + failed, passed, failed);
        System.out.printf("Total evaluations: %d%n", totalEvals);
    }

    public static void main(String[] args) {
        runTests();
    }
}