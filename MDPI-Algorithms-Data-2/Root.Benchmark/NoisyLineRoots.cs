namespace Root.Benchmark
{
    internal static class NoisyLineRoots
    {
        private const double NoisyLineFrequency = 1e6;
        private static double NoisyLine(double x) =>
            x + Math.Sin(NoisyLineFrequency * x) / 10.0 + 1e-3;

        internal static double[] Build()
        {
            // A root is possible only when |x + 0.001| <= 0.1, hence all roots
            // lie in [-0.101, 0.099].  The derivative is
            //
            //     f'(x) = 1 + 100000 cos(10^6 x).
            //
            // Its critical points satisfy cos(10^6 x) = -10^-5.  Between two
            // consecutive critical points f is monotone, so each interval
            // contains at most one root.  Testing every such interval therefore
            // produces the complete root set, independently of the algorithms
            // being compared.  For the present constants the set has 63,661
            // roots.
            const double left = -0.101;
            const double right = 0.099;
            var omega = NoisyLineFrequency;
            var twoPi = 2.0 * Math.PI;
            var alpha = Math.Acos(-1.0 / (0.1 * omega));

            var criticalPoints = new List<double> { left, right };
            var kMin = (int)Math.Floor((omega * left - alpha) / twoPi) - 2;
            var kMax = (int)Math.Ceiling((omega * right + alpha) / twoPi) + 2;

            for (var k = kMin; k <= kMax; ++k)
            {
                var xFirst = (alpha + twoPi * k) / omega;
                var xSecond = (twoPi - alpha + twoPi * k) / omega;

                if (xFirst > left && xFirst < right)
                    criticalPoints.Add(xFirst);
                if (xSecond > left && xSecond < right)
                    criticalPoints.Add(xSecond);
            }

            criticalPoints.Sort();
            var roots = new List<double>(63661);
            var x0 = criticalPoints[0];
            var y0 = NoisyLine(x0);

            for (var i = 1; i < criticalPoints.Count; ++i)
            {
                var x1 = criticalPoints[i];
                var y1 = NoisyLine(x1);

                if (y0 == 0.0)
                    AddDistinctRoot(roots, x0);

                if ((y0 < 0.0 && y1 > 0.0) || (y0 > 0.0 && y1 < 0.0))
                    AddDistinctRoot(roots, Bisect(x0, x1, y0));

                x0 = x1;
                y0 = y1;
            }

            if (y0 == 0.0)
                AddDistinctRoot(roots, x0);

            return roots.ToArray();
        }

        private static double Bisect(double left, double right, double fLeft)
        {
            // The caller supplies an interval with opposite endpoint signs.
            // Stop only when binary64 can no longer represent a distinct
            // midpoint; this makes the generated reference root as accurate as
            // the arithmetic permits without introducing a solver-dependent
            // stopping tolerance.
            while (true)
            {
                var midpoint = 0.5 * left + 0.5 * right;
                if (midpoint == left || midpoint == right)
                    return Math.Abs(NoisyLine(left)) <= Math.Abs(NoisyLine(right))
                        ? left
                        : right;

                var fMidpoint = NoisyLine(midpoint);
                if (fMidpoint == 0.0)
                    return midpoint;

                if ((fLeft < 0.0 && fMidpoint < 0.0)
                    || (fLeft > 0.0 && fMidpoint > 0.0))
                {
                    left = midpoint;
                    fLeft = fMidpoint;
                }
                else
                {
                    right = midpoint;
                }
            }
        }

        private static void AddDistinctRoot(List<double> roots, double root)
        {
            // A root at a shared interval endpoint can be encountered twice.
            // Exact comparison is sufficient because both occurrences then use
            // the same stored critical-point value.
            if (roots.Count == 0 || roots[^1] != root)
                roots.Add(root);
        }
    }
}
