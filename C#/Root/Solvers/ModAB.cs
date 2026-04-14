namespace Proektsoft.Root
{
    public static partial class Solver
    {
        // Finds the root of "F(x) = 0" within the interval [x1, x2]
        // with the specified precisions - absolute: aTol and relative: rTol,
        // using an improved version of the modified Anderson Bjork's method (Ganchovski, Traykov).
        // F(x) must be continuous and sign(F(x1)) ≠ sign(F(x2))
        public static double ModAB(Func<double, double> F, double x1, double x2,
            double aTol = 1e-14, double rTol = 1e-14)
        {
            if (!Initialize(x1, x2, F, out Node p1, out Node p2))
                return double.NaN;

            var bisection = true; // Initialize the method to bisection
            var side = 0; // Store the side that moved last: -1 for left, 1 for right, 0 for none
            var threshold = p2.X - p1.X; // Threshold to reset to bisection
            const double C = 16; // Safety factor of 4 iterations behind the threshold
            for (int i = 1; i <= MaxIterations; ++i)
            {
                var x3 = bisection ? Node.Mid(p1, p2) : Node.Sec(p1, p2);
                // Check for X-convergence and return the result
                var eps2 = aTol + rTol * Math.Abs(x3);
                if (p2.X - p1.X <= eps2)
                {
                    EvaluationCount = i + 1;
                    return x3;
                }
                Node p3;
                if (bisection)
                {
                    p3 = new Node(x3, F);
                    double y1 = p1.Y, y2 = p2.Y;
                    var ym = 0.5 * (y1 + y2);
                    var r = 1 - Math.Abs(ym / (y2 - y1)); // Symmetry factor
                    var k = r * r; // Deviation factor - quadratic
                    // Check if function is close enough to straight line and switch to false-position
                    if (Math.Abs(ym - p3.Y) < k * (Math.Abs(p3.Y) + Math.Abs(ym)))
                    {
                        bisection = false;
                        threshold = (p2.X - p1.X) * C;
                    }
                }
                else
                {
                    // Clamp the secant point to the interval [p1.X, p2.X]
                    // Otherwise, for very flat functions where p1.Y and p2.Y are close
                    // floating-point round-off errors can shoot the point outside the interval
                    if (x3 <= p1.X)
                        p3 = p1;
                    else if (x3 >= p2.X)
                        p3 = p2;
                    else
                        p3 = new Node(x3, F); // Evaluate only if not clamped

                    threshold *= 0.5;
                }
                // Check for Y-convergence and return the result
                if (p3.Y == 0)
                {
                    EvaluationCount = i + 2;
                    return x3;
                }
                if (Math.Sign(p1.Y) == Math.Sign(p3.Y))
                {
                    if (side == 1) // Apply Anderson-Bjork correction to the right side
                    {
                        var m = 1 - p3.Y / p1.Y;
                        p2.Y *= m <= 0 ? 0.5 : m;
                    }
                    else if (!bisection)
                        side = 1;

                    p1 = p3;
                }
                else
                {
                    if (side == -1) // Apply Anderson-Bjork correction to the left side
                    {
                        var m = 1 - p3.Y / p2.Y;
                        p1.Y *= m <= 0 ? 0.5 : m;
                    }
                    else if (!bisection)
                        side = -1;

                    p2 = p3;
                }
                if (p2.X - p1.X > threshold) // If AB fails to shrink the interval enough
                {
                    bisection = true;        // reset to bisection
                    side = 0;
                }
            }
            EvaluationCount = MaxIterations + 2;
            return double.NaN; // When failed to converge within maxIterations
        }
    }
}