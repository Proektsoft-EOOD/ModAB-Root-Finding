using System.Runtime.CompilerServices;

namespace Proektsoft.Root
{
    internal struct Node
    {
        public double X;
        public double Y;

        public Node(double x, double y)
        {
            X = x;
            Y = y;
        }

        public Node(double x, Func<double, double> F)
        {
            X = x;
            Y = F(x);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double SafeSecant(Node p1, Node p2) =>
            SafeSecant(p1.X, p1.Y, p2.X, p2.Y);

        // Secant through (x1, y1) and (x2, y2), where the ordinates need not be
        // the stored residuals (e.g. Anderson-Björck corrected values).
        public static double SafeSecant(double x1, double y1, double x2, double y2)
        {
            var a = Math.Abs(y1);
            var b = Math.Abs(y2);
            var den = a + b;
            if (!(den > 0.0))
                return 0.5 * x1 + 0.5 * x2;

            if (double.IsInfinity(den))
            {
                if (double.IsInfinity(a) || double.IsInfinity(b))
                    return 0.5 * x1 + 0.5 * x2;

                a *= 0.5;
                b *= 0.5;
                den = a + b;
            }
            var w1 = b / den;
            var w2 = a / den;
            var x = w1 * x1 + w2 * x2;
            return x < x1 ? x1 :
                   x > x2 ? x2 :
                   x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double SafeMidpoint(Node p1, Node p2) =>
            0.5 * p1.X + 0.5 * p2.X;
    }
}
