using Proektsoft.Root;
using System.Diagnostics;
namespace Root.Benchmark
{
    internal class BenchmarkProblems
    {
        private static double P(double x) => x + 1.11111;

        // Roots of the sawtooth functions f82 and f83.
        // With n = floor(100*x + 0.5) the function is 202*x - 2*n - 0.1, so x = (20*n + 1)/2020.
        // That value is consistent with its own n only for n = -45...55, which gives 101 roots.
        // Note that the sawtooth also drops through zero at the 100 cell boundaries in between,
        // x = (2*n + 1)/200, n = -45...54; those are discontinuities, not roots.
        private static readonly double[] sawtoothRoots =
            Enumerable.Range(-45, 101).Select(n => (20*n + 1)/2020d).ToArray();

        private static readonly double[] noisyLineRoots = NoisyLineRoots.Build();

        // Test examples from various publications as specified bellow
        internal static readonly Problem[] Set1 =
        {
            //Sérgio Galdino. A family of regula falsi root-finding methods
            new() {
                Name = "f01",
                F = (x) => Math.Pow(x, 3) - 1,
                a = 0.5, b = 1.5,
                Roots = [1]
            },
            new() {
                Name = "f02",
                F = (x) => Math.Pow(x, 2)*(Math.Pow(x, 2)/3 + Math.Sqrt(2)*Math.Sin(x)) - Math.Sqrt(3)/18,
                a = 0.1, b = 1,
                Roots = [0.39942229171096819]
            },
            new() {
                Name = "f03",
                F = (x) => 11*Math.Pow(x, 11) - 1,
                a = 0.1, b = 1,
                Roots = [Math.Pow(11, -1.0 / 11.0)]
            },
            new() {
                Name = "f04",
                F = (x) => Math.Pow(x, 3) + 1,
                a = -1.8, b = 0,
                Roots = [-1]
            },
            new() {       
                Name = "f05",
                F = (x) => Math.Pow(x, 3) - 2*x - 5,
                a = 2, b = 3,
                Roots = [2.09455148154232659]
            },
            new() {
                Name = "f06",
                F = (x) => 2*x*Math.Exp(-5) + 1 - 2*Math.Exp(-5*x),
                a = 0, b = 1,
                Roots = [0.138257155056824076]
            },
            new() {
                Name = "f07",
                F = (x) => 2*x*Math.Exp(-10) + 1 - 2*Math.Exp(-10*x),
                a = 0, b = 1,
                Roots = [0.069314088687023473]
            },
            new() {
                Name = "f08",
                F = (x) => 2*x*Math.Exp(-20) + 1 - 2*Math.Exp(-20*x),
                a = 0, b = 1,
                Roots = [0.034657359020853851]
            },
            new() {
                Name = "f09",
                F = (x) => (1 + Math.Pow(1 - 5, 2))*Math.Pow(x, 2) - Math.Pow(1 - 5*x, 2),
                a = 0, b = 1,
                Roots = [(5 - Math.Sqrt(17))/8]
            },
            new() {
                Name = "f10",
                F = (x) => (1 + Math.Pow(1 - 10, 2))*Math.Pow(x, 2) - Math.Pow(1 - 10*x, 2),
                a = 0, b = 1,
                Roots = [(10 - Math.Sqrt(82)) / 18]
            },
            new() {
                Name = "f11",
                F = (x) => (1 + Math.Pow(1 - 20, 2))*Math.Pow(x, 2) - Math.Pow(1 - 20*x, 2),
                a = 0, b = 1,
                Roots = [ (20 - Math.Sqrt(362))/38]
            },
            new() {
                Name = "f12",
                F = (x) => Math.Pow(x, 2) - Math.Pow(1 - x, 5),
                a = 0, b = 1,
                Roots = [0.34595481584824202]
            },
            new(){
                Name = "f13",
                F = (x) => Math.Pow(x, 2) - Math.Pow(1 - x, 10),
                a = 0, b = 1,
                Roots = [0.24512233375330725]
            },
            new() {
                Name = "f14",
                F = (x) => Math.Pow(x, 2) - Math.Pow(1 - x, 20),
                a = 0, b = 1,
                Roots = [0.16492095727644096]
            },
            new() {
                Name = "f15",
                F = (x) => (1 + Math.Pow(1 - 5, 4))*x -Math.Pow(1 - 5*x, 4),
                a = 0, b = 1,
                Roots = [0.0036171081789040634]
            },
            new() {
                Name = "f16",
                F = (x) => (1 + Math.Pow(1 - 10, 4))*x - Math.Pow(1 - 10*x, 4),
                a = 0,b = 1,
                Roots = [0.0001514713347838914]
            },
            new() {
                Name = "f17",
                F = (x) => (1 + Math.Pow(1 - 20, 4))*x - Math.Pow(1 - 20*x, 4),
                a = 0, b = 1,
                Roots = [7.668595122185337e-6]
            },
            new() {
                Name = "f18",
                F = (x) => Math.Exp(-5*x)*(x - 1) + Math.Pow(x, 5),
                a = 0, b = 1,
                Roots = [0.5161535187579336]
            },
            new() {
                Name = "f19",
                F = (x) => Math.Exp(-10*x)*(x - 1) + Math.Pow(x, 10),
                a = 0, b = 1,
                Roots = [0.5395222269084159]
            },
            new() {
                Name = "f20",
                F = (x) => Math.Exp(-20*x)*(x - 1) + Math.Pow(x, 20),
                a = 0, b = 1,
                Roots = [0.5527046666784878]
            },
            new() {
                Name = "f21",
                F = (x) => Math.Pow(x, 2) + Math.Sin(x/5) - 1d/4d,
                a = 0, b = 1,
                Roots = [0.40999201798913715]
            },
            new() {
                Name = "f22",
                F = (x) => Math.Pow(x, 2) + Math.Sin(x/10) - 1d/4d,
                a = 0, b = 1,
                Roots = [0.4525091455776412]
            },
            new() {
                Name = "f23",
                F = (x) => Math.Pow(x, 2) + Math.Sin(x/20) - 1d/4d,
                a = 0, b = 1,
                Roots = [0.4756268485960624]
            },
            new() {
                Name = "f24",
                F = (x) => (x + 2)*(x + 1)*Math.Pow(x - 3, 3),
                a = 2.6, b = 4.6,
                Roots = [3]
            },
            new() {
                Name = "f25",
                F = (x) => Math.Pow(x - 4, 5) * Math.Log(x),
                a = 3.6, b = 5.6,
                Roots = [4]
            },
            new() {
                Name = "f26",
                F = (x) => Math.Pow(Math.Sin(x) - x/4, 3),
                a = 2, b = 4,
                Roots = [2.4745767873698292]
            },
            new() {
                Name = "f27",
                F = (x) => (81 - P(x)*(108 - P(x)*(54 - P(x)*(12 - P(x)))))*Math.Sign(P(x) - 3),
                a = 1, b = 3,
                Roots = [1.88889] // F = (P - 3)^4*sign(P - 3), so P(x) = 3
            },
            new() {
                Name = "f28",
                F = (x) => Math.Sin(Math.Pow(x - 7.143, 3)),
                a = 7, b = 8,
                Roots = [7.143]
            },
            new() {
                Name = "f29",
                F = (x) => Math.Exp(Math.Pow(x - 3, 5)) - 1,
                a = 2.6, b = 4.6,
                Roots = [3]
            },
            new() {
                Name = "f30",
                F = (x) => Math.Exp(Math.Pow(x - 3, 5)) - Math.Exp(x - 1),
                a = 4, b = 5,
                Roots = [4.267168304542125]
            },
            //My functions
            new() {
                Name = "f31",
                F = (x) => Math.PI - 1/x,
                a = 0.05, b = 5,
                Roots = [1/Math.PI]
            },
            new() {
                Name = "f32",
                F = (x) => 4 - Math.Tan(x),
                a = 0, b = 1.5,
                Roots = [Math.Atan(4)]
            },
            //Steven A. Stage. Comments on An Improvement to the Brent’s Method 
            new() {
                Name = "f33",
                F = (x) => Math.Cos(x) - Math.Pow(x, 3),
                a = 0, b = 4,
                Roots = [0.8654740331016144]
            },
            new() {
                Name = "f34",
                F = (x) => Math.Cos(x) - x,
                a = -11, b = 9,
                Roots = [0.7390851332151607] // Dottie number
            },
            new() {
                Name = "f35",
                F = (x) => Math.Sqrt(Math.Abs(x - 2d/3d))*(x <= 2d/3d ? 1 : -1) - 0.1,
                a = -11, b = 9,
                Roots = [2d/3d - 0.01]
            },
            new() {
                Name = "f36",
                F = (x) => Math.Pow(Math.Abs(x - 2d/3d), 0.2)*(x <= 2d/3d ? 1 : -1),
                a = -11, b = 9,
                Roots = [2d/3d]
            },
            new() {
                Name = "f37",
                F = (x) => Math.Pow(x - 7d/9d, 3) + (x - 7d/9d) * 1e-3,
                a = -11, b = 9,
                Roots = [7d/9d]
            },
            new() {
                Name = "f38",
                F = (x) => x <= 1d/3d ? -0.5 : 0.5,
                a = -11, b = 9,
                Roots = [] // no root: jump discontinuity from -0.5 to +0.5 at x = 1/3
            },
            new() {
                Name = "f39",
                F = (x) => x <= 1d/3d ? -1e-3 : 1 - 1e-3,
                a = -11, b = 9,
                Roots = [] // no root: jump discontinuity from -1e-3 to 0.999 at x = 1/3
            },
            new() {
                Name = "f40",
                F = (x) => x == 0 ? 0 : 1 / (x - 2d/3d),
                a = -11, b = 9,
                // no root: sign change is caused by the pole at x = 2/3.
                // The coded special case F(0) = 0 is a spurious zero of the C# lambda only
                Roots = []
            },
            //A. Swift and G.R. Lindfield. Comparison of a Continuation Method with Brents Method for the Numerical Solution of a Single Nonlinear Equation
            new() {
                Name = "f41",
                F = (x) => 2*x*Math.Exp(-5) - 2*Math.Exp(-5*x) + 1,
                a = 0, b = 10,
                Roots = [0.13825715505682407] // same as f06; F is strictly increasing on [0, 10]
            },
            new() {
                Name = "f42",
                F = (x) => (Math.Pow(x, 2) - x - 6)*(Math.Pow(x, 2) - 3*x + 2),
                a = 0, b = Math.PI,
                Roots = [1, 2, 3] // -2 is outside [a, b]
            },
            new() {
                Name = "f43",
                F = (x) => Math.Pow(x, 3),
                a = -1, b = 1.5,
                Roots = [0]
            },
            new() {
                Name = "f44",
                F = (x) => Math.Pow(x, 5),
                a = -1, b = 1.5,
                Roots = [0]
            },
            new() {
                Name = "f45",
                F = (x) => Math.Pow(x, 7),
                a = -1, b = 1.5,
                Roots = [0]
            },
            new() {
                Name = "f46",
                F = (x) => (Math.Exp(-5*x) - x - 0.5)/Math.Pow(x, 5),
                a = 0.09, b = 0.7,
                Roots = [0.10162439229354954]
            },
            new() {
                Name = "f47",
                F = (x) => 1/Math.Sqrt(x) - 2*Math.Log(5e3*Math.Sqrt(x)) + 0.8,
                a = 0.0005, b = 0.5,
                Roots = [0.007732523200612798]
            },
            new() {
                Name = "f48",
                F = (x) => 1/Math.Sqrt(x) - 2*Math.Log(5e7*Math.Sqrt(x)) + 0.8,
                a = 0.0005, b = 0.5,
                Roots = [0.0012763049457355643]
            },
            new() {
                Name = "f49",
                F = (x) => x <= 0 ? -Math.Pow(x, 3) - x - 1 : Math.Pow(x, 1.0/3.0) - x - 1,
                a = -1, b = 1,
                // real root of x^3 + x + 1 = 0; the branch x > 0 stays strictly negative
                Roots = [-0.6823278038280193]
            },
            new() {
                Name = "f50",
                F = (x) => Math.Pow(x, 3) - 2*x - x + 3,
                a = -3, b = 2,
                Roots = [-2.1038034027355366] // x^3 - 3x + 3; the other two roots are imaginary
            },
            new() {
                Name = "f51",
                F = (x) => Math.Log(x),
                a = 0.5, b = 5,
                Roots = [1]
            },
            new() {
                Name = "f52",
                F = (x) => (10 - x)*Math.Exp(-10*x) - Math.Pow(x, 10) + 1,
                a = 0.5, b = 8,
                Roots = [1.0000408355647268]
            },
            new() {
                Name = "f53",
                F = (x) => Math.Exp(Math.Sin(x)) - x - 1,
                a = 1.0, b = 4,
                Roots = [1.6968123868097515]
            },
            new() {
                Name = "f54",
                F = (x) => 2*Math.Sin(x) - 1,
                a = 0.1, b = Math.PI/3,
                Roots = [0.5235987755982989] // pi/6; 5*pi/6 is outside [a, b]
            },
            new() {
                Name = "f55",
                F = (x) => (x - 1)*Math.Exp(-x),
                a = 0.0, b = 1.5,
                Roots = [1]
            },
            new() {
                Name = "f56",
                F = (x) => Math.Pow(x - 1, 3) - 1,
                a = 1.5, b = 3,
                Roots = [2]
            },
            new() {
                Name = "f57",
                F = (x) => Math.Exp(Math.Pow(x, 2) + 7*x - 30) - 1,
                a = 2.6, b = 3.5,
                Roots = [3] // (x + 10)*(x - 3) = 0; -10 is outside [a, b]
            },
            new() {
                Name = "f58",
                F = (x) => Math.Atan(x) - 1,
                a = 1.0, b = 8,
                Roots = [Math.Tan(1)]
            },
            new() {
                Name = "f59",
                F = (x) => Math.Exp(x) - 2*x - 1,
                a = 0.2, b = 3,
                Roots = [1.2564312086261697] // the second root, x = 0, is outside [a, b]
            },
            new() {
                Name = "f60",
                F = (x) => Math.Exp(-x) - x - Math.Sin(x),
                a = 0.0, b = 2,
                Roots = [0.3544631043750253]
            },
            new() {
                Name = "f61",
                F = (x) => Math.Pow(x, 2) - Math.Pow(Math.Sin(x),2)  - 1,
                a = -1, b = 2,
                Roots = [1.4044916482153411] // F is even; -1.4044916482153411 is outside [a, b]
            },
            new() {
                Name = "f62",
                F = (x) => Math.Sin(x) - x/2,
                a = Math.PI/2, b = Math.PI,
                Roots = [1.895494267033981]
            }
        };

        // Test examples from the publication
        // Oliveira I. F. D., Takahashi R. H. C.
        // An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality

        internal static readonly Problem[] Set2 =
        {
            new() { //Lambert
                Name = "f63",
                F = (x) => x * Math.Exp(x) - 1d,
                a = -1d, b = 1d,
                Roots = [0.5671432904097838] // omega constant W(1)
            },
            new() { //Trigonometric 1
                Name = "f64",
                F = (x) => Math.Tan(x - 1d/10d),
                a = -1d, b = 1d,
                Roots = [0.1] // x - 0.1 = k*pi, k != 0, falls outside [a, b]
            },
            new() { //Trigonometric 2
                Name = "f65",
                F = (x) => Math.Sin(x) + 0.5,
                a = -1d, b = 1d,
                Roots = [-Math.PI/6d]
            },
            new() { //Polynomial 1
                Name = "f66",
                F = (x) => 4 * Math.Pow(x, 5d) + x * x + 1d,
                a = -1d, b = 1d,
                Roots = [-0.8439145686492663]
            },
            new() { //Polynomial 2
                Name = "f67",
                F = (x) => x + Math.Pow(x, 10d) - 1d,
                a = -1d, b = 1d,
                Roots = [0.835079042723559]
            },
            new() { //Exponential
                Name = "f68",
                F = (x) => Math.Pow(Math.PI, x) - Math.E,
                a = -1d, b = 1d,
                Roots = [1d/Math.Log(Math.PI)]
            },
            new() { //Logarithmic
                Name = "f69",
                F = (x) => Math.Log(Math.Abs(x - 10d/9d)),
                a = -1d, b = 1d,
                Roots = [1d/9d]
            },
            new() { //Posynomial
                Name = "f70",
                F = (x) => 1d/3d + Math.Sign(x) * Math.Cbrt(Math.Abs(x)) + Math.Pow(x, 3d),
                a = -1d, b = 1d,
                Roots = [-0.03702012770786093]
            },
            new() { //Poly.Frac.
                Name = "f71",
                F = (x) => (x + 2d/3d)/(x + 101d/100d),
                a = -1d, b = 1d,
                Roots = [-2d/3d]
            },
            new() { //Polynomial 3
                Name = "f72",
                F = (x) => Math.Pow(x * 1e6 - 1d, 3d),
                a = -1d, b = 1d,
                Roots = [1e-6]
            },
            new() { //Exp. Poly.
                Name = "f73",
                F = (x) => Math.Exp(x) * Math.Pow(x * 1e6 - 1d, 3d),
                a = -1d, b = 1d,
                Roots = [1e-6]
            },
            new() { //Tan. Poly.
                Name = "f74",
                F = (x) => Math.Pow(x - 1d/3d, 2d) * Math.Atan(x - 1d/3d),
                a = -1d, b = 1d,
                Roots = [1d/3d]
            },
            new() { //Circles
                Name = "f75",
                F = (x) => Math.Sign(3d*x - 1d) * (1d - Math.Sqrt(1d - Math.Pow(3d*x - 1d, 2d)/81d)),
                a = -1d, b = 1d,
                Roots = [1d/3d]
            },
            new() { //Step Function
                Name = "f76",
                F = (x)  => x > (1d - 1e6) / 1e6 ? (1d + 1e6) / 1e6 : 0d - 1d,
                a = -1d, b = 1d,
                // no root: jump discontinuity from -1 to 1.000001 at x = -0.999999
                Roots = []
            },
            new() { //Geometric
                Name = "f77",
                F = (x) => x != 1d/21d ? 1/(21d*x - 1d) : 0d,
                a = -1d, b = 1d,
                // no root: sign change is caused by the pole at x = 1/21.
                // The coded special case F(1/21) = 0 is a spurious zero of the C# lambda only
                Roots = []
            },
            new() { //Trunc.Poly.
                Name = "f78",
                F = (x) => x * x / 4d + Math.Ceiling(x/2d) - 0.5,
                a = -1d, b = 1d,
                // no root: on [-1, 0] F = x^2/4 - 0.5 <= -0.25, on (0, 1] F = x^2/4 + 0.5 > 0,
                // so F only jumps through zero at x = 0
                Roots = []
            },
            new() { //Staircase
                Name = "f79",
                F = (x) => Math.Ceiling(10d*x - 1d) + 0.5,
                a = -1d, b = 1d,
                // no root: F is always an integer + 0.5; it jumps through zero at x = 0
                Roots = []
            },
            new() { //Noisy Line
                Name = "f80",
                F = (x) => x + Math.Sin(x*1e6)/10d + 1e-3,
                a = -1d, b = 1d,
                // 63661 roots in total, all inside (-0.101, 0.099). Listed below are the roots
                // that the benchmarked methods converge to, plus the two nearest neighbours of each
                Roots = noisyLineRoots
            },
            new() { //Warsaw
                Name = "f81",
                F = (x) => x > -1 ? 1 + Math.Sin(1d/(x + 1d)) : 0d - 1d,
                a = -1d, b = 1d,
                // Infinitely many roots x_k = 1/(3*pi/2 + 2*k*pi) - 1, k = 0, 1, 2, ..., accumulating
                // at x = -1. All are touching (double) roots - F >= 0 for x > -1, so no root produces
                // a sign change and every bracketing method converges to the discontinuity at x = -1,
                // which is not a root. Listed below are k = 23...0
                Roots = []
            },
            new() { //Sawtooth
                Name = "f82",
                F = (x) => 202d*x - 2*Math.Floor((2d*x + 1e-2)/2e-2) - 0.1,
                a = -1d, b = 1d,
                Roots = sawtoothRoots
            },
            new() { //Sawtooth Cube
                Name = "f83",
                F = (x) => Math.Pow(202d*x - 2d*Math.Floor((2d*x + 1e-2)/2e-2) - 0.1, 3d),
                a = -1d, b = 1d,
                Roots = sawtoothRoots
            },
        };

        //SciML Benchmarks test suite
        internal static readonly Problem[] Set3 =
        {
            new() { // Polynomial with multiple roots  
                Name = "f84",
                F = (x) => (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05,
                a = 0.5, b = 5.5,
                Roots = [
                    1.0020924414628258, 1.9917232770437472, 3.0125024427614666,
                    3.9916074829067334, 5.0020743558252265
                ]
            },
            new() { // Function 2: Trigonometric with multiple roots
                Name = "f85",
                F = (x) => Math.Sin(x) - 0.5*x - 0.3,
                a = -10.0, b = 10.0,
                Roots = [-2.207783144886932, 0.72245153809491, 1.352562882201154]
            },
            new() { // Function 3: Exponential function (sensitive near zero)
                Name = "f86",
                F = (x) =>  Math.Exp(x) - 1 - x - x*x/2 - 0.005,
                a = -2.0, b = 2.0,
                Roots = [0.3028041572706959] // F is strictly increasing
            },
            new() { // Function 4: Rational function with pole
                Name = "f87",
                F = (x) =>  1/(x - 0.5) - 2 - 0.05,
                a = 0.6, b = 2.0,
                Roots = [0.5 + 1d/2.05]
            },
            new() { // Function 5: Logarithmic function
                Name = "f88",
                F = (x) => Math.Log(x) - x + 2 - 0.05,
                a = 0.1, b = 3.0,
                Roots = [0.16836217326253688] // the second root, 3.0509..., is outside [a, b]
            },
            new() { // Function 6: High oscillation function
                Name = "f89",
                F = (x) =>  Math.Sin(20*x) + 0.1*x - 0.1,
                a = -4.0, b = 5.0,
                Roots = [
                    -3.952899243704171, -3.745190306935815, -3.636939548751743, -3.4327892525880244,
                    -3.3210133649214963, -3.120358291368688, -3.005116780945275, -2.8079008579419362,
                    -2.6892463442983034, -2.4954200116274747, -2.3733989763275805, -2.182918503433951,
                    -2.057571905059434, -1.870398829861522, -1.7417626112164155, -1.5578632766897285,
                    -1.4259687842239264, -1.2453139551072887, -1.1101882858457803, -0.9327528319398148,
                    -0.7944191196876997, -0.6201817553057012, -0.47865940523246825, -0.3076024767250692,
                    -0.16290735537457432, 0.004983329514230639, 0.15283874335850106, 0.3175740490899948,
                    0.4685805485527247, 0.6301681119510958, 0.784319678162274, 0.9427639777507327,
                    1.1000577257332431, 1.2553601220207036, 1.4157962750034587, 1.5679550227365961,
                    1.7315369142539843, 1.8805471469428237, 2.0472812507754727, 2.1931349371093125,
                    2.3630309258182582, 2.5057167968786542, 2.678787630419438, 2.8182910758325614,
                    2.9945531225451907, 3.130856052857107, 3.3103292460561895, 3.4434099176132276,
                    3.626117952104279, 3.7559507495159505, 3.9419213237091664, 4.068476493483119,
                    4.25774160445862, 4.380984931517933, 4.573581232545281, 4.693473648918655,
                    4.889442881728638
                ]
            },
            new() { // Function 7: Function with very flat region
                Name = "f90",
                F = (x) =>  x*x*x - 2*x*x + x - 0.025,
                a = -1.0, b = 2.0,
                Roots = [0.02637269541443003, 0.826031016050984, 1.147596288534586]
            },
            new() { // Function 8: Bessel-like function
                Name = "f91",
                F = (x) =>  x*Math.Sin(1/x) - 0.1 - 0.01,
                a = 0.01, b = 1.0,
                // |x*sin(1/x)| <= x < 0.11 for x < 0.11, so there is no root below 0.11
                Roots = [0.12077801181192625, 0.1389533072336297, 0.353913733550918]
            },
            new() { // Near symmetric function
                Name = "f92",
                F = (x) => x*x*x - 0.001,
                a = -10, b = 10,
                Roots = [0.1] // 0.001^(1/3)
            },
            new() { // Near symmetric function
                Name = "f93",
                F = (x) => Math.Pow(x, 5) - 0.001,
                a = -10, b = 10,
                Roots = [0.251188643150958] // 0.001^(1/5)
            },
            new() { // Near symmetric function
                Name = "f94",
                F = (x) => Math.Pow(x, 7) - 0.001,
                a = -10, b = 10,
                Roots = [0.372759372031494] // 0.001^(1/7)
            },
            new() { // Near symmetric function
                Name = "f95",
                F = (x) => Math.Pow(x, 9) - 0.001,
                a = -10, b = 10,
                Roots = [0.46415888336127786] // 0.001^(1/9)
            },
            new() { // Near symmetric function
                Name = "f96",
                F = (x) => Math.Pow(x, 11) - 0.001,
                a = -10, b = 10,
                Roots = [0.533669923120631] // 0.001^(1/11)
            },
            new() { // Near symmetric function
                Name = "f97",
                F = (x) => Math.Pow(x, 13) - 0.001,
                a = -10, b = 10,
                Roots = [0.5878016072274913] // 0.001^(1/13)
            },
            new() { // Near symmetric function
                Name = "f98",
                F = (x) => Math.Pow(x, 15) - 0.001,
                a = -10, b = 10,
                Roots = [0.6309573444801932] // 0.001^(1/15)
            },
            new() { // Near symmetric function
                Name = "f99",
                F = (x) => Math.Pow(x, 17) - 0.001,
                a = -10, b = 10,
                Roots = [0.6660846290809158] // 0.001^(1/17)
            },
            new() { // Verti.cal tangent at root
                Name = "f100",
                F = (x) => Math.Cbrt((4*x - 3)/x),
                a = 0, b = Math.E,
                Roots = [0.75]
            },        
        };
    }
}
