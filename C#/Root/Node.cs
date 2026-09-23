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

        public static double SafeSecant(Node p1, Node p2)
        {
            var a = Math.Abs(p1.Y);
            var b = Math.Abs(p2.Y);
            var den = a + b;
            if (!(den > 0.0))
                return SafeMidpoint(p1, p2);

            if (double.IsInfinity(den))
            {
                if (double.IsInfinity(a) || double.IsInfinity(b))
                    return SafeMidpoint(p1, p2);

                a *= 0.5;
                b *= 0.5;
                den = a + b;
            }
            var w1 = b / den;
            var w2 = a / den;
            var x = w1 * p1.X + w2 * p2.X;
            return x < p1.X ? p1.X :
                   x > p2.X ? p2.X :
                   x;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double SafeMidpoint(Node p1, Node p2) =>
            0.5 * p1.X + 0.5 * p2.X;
    }
}
