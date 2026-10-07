using System.Globalization;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text.Json;

namespace ModAB.Switching.Verification;

internal static class Program
{
    private const string DataHash = "CEDE3121D4B54F765A2CB7DC1992A9F0C68A9710EFD6AD95B5FD4BFBE461E4BB";
    private static readonly Rational Epsilon = Rational.FromDouble(Binary64.UnitRoundoff);
    private static readonly string[] Names = ["k", "h", "left", "right", "gap"];
    private static readonly int[] Bounds = [8, 3, 5, 19, 25];

    private static int Main(string[] args)
    {
        try
        {
            string? output = null;
            for (int i = 0; i < args.Length; ++i)
            {
                if (args[i] == "--help")
                {
                    Console.WriteLine("ModAB.Switching.Verification [--output results.json]");
                    return 0;
                }
                if (args[i] == "--output" && i + 1 < args.Length) output = args[++i];
                else throw new ArgumentException("Expected --output <file> or --help.");
            }

            Require(Binary64.HasRequiredArithmetic(), "Round-to-nearest and gradual-underflow checks failed.");
            CheckRationalArithmetic();
            RoundingCheck.CheckBoundaries();
            object regressions = CheckRegressionExamples();
            string dataPath = Path.Combine(AppContext.BaseDirectory, "Data", "cases.csv");
            using (FileStream stream = File.OpenRead(dataPath))
                Require(Convert.ToHexString(SHA256.HashData(stream)) == DataHash, "Test data checksum differs.");

            int total = 0, accepted = 0, rejected = 0, ambiguous = 0;
            int rawDisagreements = 0, rawRegularChecked = 0;
            Rational[] maxima = Enumerable.Repeat(Rational.Zero, Names.Length).ToArray();

            foreach (string line in File.ReadLines(dataPath).Skip(1))
            {
                string[] cells = line.Split(',');
                Require(cells.Length == 3, "Invalid test row.");
                double a = FromBits(cells[0]), b = FromBits(cells[1]), z = FromBits(cells[2]);
                RoundingCheck.Verify(a, b, z);
                SwitchingResult computed = SwitchingCriterion.Evaluate(a, b, z);
                ExactValues exact = ExactValues.Evaluate(a, b, z);
                string context = $"Case {total + 1}: {line}";
                Require(computed.Decision != SwitchingDecision.InvalidInput, context);
                double[] all = [computed.Scale, computed.Rho, computed.K, computed.P1,
                    computed.P2, computed.P3, computed.H, computed.Left, computed.Right, computed.Gap];
                Require(all.All(double.IsFinite), $"{context}: non-finite intermediate.");
                double[] actual = [computed.K, computed.H, computed.Left, computed.Right, computed.Gap];
                Rational[] expected = [exact.K, exact.H, exact.Left, exact.Right, exact.Gap];
                for (int j = 0; j < Names.Length; ++j)
                {
                    Rational error = Rational.Abs(Rational.FromDouble(actual[j]) - expected[j]);
                    Require(error <= Bounds[j] * Epsilon, $"{context}: {Names[j]} bound failed.");
                    maxima[j] = Rational.Max(maxima[j], error / Epsilon);
                }
                switch (computed.Decision)
                {
                    case SwitchingDecision.CertifiedSwitch:
                        ++accepted;
                        Require(exact.Gap > 0, $"{context}: false certified switch.");
                        break;
                    case SwitchingDecision.CertifiedNoSwitch:
                        ++rejected;
                        Require(exact.Gap < 0, $"{context}: false certified rejection.");
                        break;
                    case SwitchingDecision.Indeterminate:
                        ++ambiguous;
                        break;
                    default: throw new InvalidOperationException(context);
                }
                Require(computed.ShouldSwitch == (computed.Gap > SwitchingCriterion.Guard), context);
                if (exact.Gap > 57 * Epsilon)
                    Require(computed.Decision == SwitchingDecision.CertifiedSwitch, $"{context}: positive margin.");
                if (exact.Gap < -57 * Epsilon)
                    Require(computed.Decision == SwitchingDecision.CertifiedNoSwitch, $"{context}: negative margin.");

                RawValues? raw = RawValues.Evaluate(a, b, z);
                if ((raw?.ShouldSwitch ?? false) != (exact.Gap > 0)) ++rawDisagreements;
                if (raw is RawValues r && r.MeetsNormalRangeHypotheses(a, b, z))
                {
                    Rational af = Rational.FromDouble(a), bf = Rational.FromDouble(b), zf = Rational.FromDouble(z);
                    Rational hScale = Rational.Max(Rational.Abs(bf - af), Rational.Abs(zf));
                    Rational rawGap = Rational.FromDouble(r.Right) - Rational.FromDouble(r.Left);
                    Require(Rational.Abs(rawGap - exact.Gap * exact.Scale) <= 16 * Epsilon * hScale,
                        $"{context}: unscaled bound failed.");
                    ++rawRegularChecked;
                }
                ++total;
            }

            Require(total == 23816, "Unexpected case count.");
            Require(accepted == 5040 && rejected == 7780 && ambiguous == 10996,
                "Decision counts differ from the fixed reference dataset.");
            Require(rawDisagreements == 4163 && rawRegularChecked == 22615,
                "Unscaled counts differ from the fixed reference dataset.");

            var report = new
            {
                Status = "passed",
                Runtime = RuntimeInformation.FrameworkDescription,
                Architecture = RuntimeInformation.ProcessArchitecture.ToString(),
                OperatingSystem = OperatingSystem.IsLinux() ? "Linux" : OperatingSystem.IsWindows() ? "Windows" : "Other",
                Counts = new { Total = total, CertifiedAccept = accepted, CertifiedReject = rejected,
                    Ambiguous = ambiguous, RawDisagreements = rawDisagreements, RawRegularChecked = rawRegularChecked },
                MaximumObservedErrorsInUnitsOfEpsilon = Names.Select((name, i) => (name, value: maxima[i].ToApproximateDouble()))
                    .ToDictionary(x => x.name, x => x.value),
                ExactMaximumErrorsInUnitsOfEpsilon = Names.Select((name, i) => (name, value: maxima[i].ToString()))
                    .ToDictionary(x => x.name, x => x.value),
                TheoremBoundsInUnitsOfEpsilon = Names.Select((name, i) => (name, value: Bounds[i]))
                    .ToDictionary(x => x.name, x => x.value),
                RoundingOracle = new { CheckedOperations = RoundingCheck.CheckedOperations,
                    BoundaryOperations = RoundingCheck.BoundaryOperations,
                    AlgorithmOperations = total * 17, BitwiseMismatches = 0 },
                Guard = SwitchingCriterion.Guard,
                UnitRoundoff = Binary64.UnitRoundoff,
                DataSha256 = DataHash.ToLowerInvariant(),
                RegressionExamples = regressions
            };
            string json = JsonSerializer.Serialize(report, new JsonSerializerOptions
                { WriteIndented = true, PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower }) + Environment.NewLine;
            if (output is not null)
            {
                string fullPath = Path.GetFullPath(output);
                Directory.CreateDirectory(Path.GetDirectoryName(fullPath)!);
                File.WriteAllText(fullPath, json);
            }
            Console.Write(json);
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine($"Verification failed: {ex.Message}");
            return 1;
        }
    }

    private static object CheckRegressionExamples()
    {
        double tiny = double.Epsilon;
        var subnormal = SwitchingCriterion.Evaluate(-tiny, 2 * tiny, tiny);
        var subnormalExact = ExactValues.Evaluate(-tiny, 2 * tiny, tiny);
        Require(RawValues.Evaluate(-tiny, 2 * tiny, tiny)?.ShouldSwitch == false, "Subnormal raw example.");
        Require(subnormal.ShouldSwitch && subnormalExact.Gap == new Rational(13, 48), "Subnormal scaled example.");

        double a = FromBits("C00D674B0EF52CA9"), b = FromBits("400BAF7E83D81B4B"), z = FromBits("BF5AAB2B125E02F1");
        var boundary = SwitchingCriterion.Evaluate(a, b, z);
        var boundaryExact = ExactValues.Evaluate(a, b, z);
        var boundaryRaw = RawValues.Evaluate(a, b, z);
        Require(boundaryRaw is RawValues r && r.ShouldSwitch && r.MeetsNormalRangeHypotheses(a, b, z),
            "Normal-range raw false positive.");
        Require(boundaryExact.Gap < 0 && boundary.Decision == SwitchingDecision.Indeterminate,
            "Guard must decline the normal-range false positive.");
        Require(SwitchingCriterion.Evaluate(-1, 1, 0).Decision == SwitchingDecision.Indeterminate, "Exact equality.");

        (double, double, double)[] invalid = [(0, 1, 1), (-0.0, -1, 0), (-1, 0, 0),
            (1, 2, 1), (-1, -2, 1), (double.NaN, 1, 0), (-1, double.NaN, 0), (-1, 1, double.NaN),
            (double.NegativeInfinity, 1, 0), (-1, double.PositiveInfinity, 0),
            (-1, 1, double.PositiveInfinity), (-1, 1, double.NegativeInfinity)];
        foreach (var (f1, f2, f3) in invalid)
        {
            var result = SwitchingCriterion.Evaluate(f1, f2, f3);
            Require(result.Decision == SwitchingDecision.InvalidInput && !result.ShouldSwitch, "Invalid input accepted.");
        }

        double q = Binary64.Divide(double.MaxValue, 4);
        double endpoint = Binary64.Multiply(3, q);
        RawValues overflow = RawValues.Evaluate(-q, endpoint, double.MaxValue)!.Value;
        double grouped = Binary64.Multiply(overflow.K, Binary64.Add(Math.Abs(overflow.Y), double.MaxValue));
        Require(double.IsPositiveInfinity(grouped) && double.IsFinite(overflow.Right), "Grouped overflow example.");
        Require(double.IsFinite(SwitchingCriterion.Evaluate(-q, endpoint, double.MaxValue).Gap), "Scaled overflow example.");

        return new
        {
            Subnormal = new { ExactScaledGap = "13/48", ComputedScaledGap = subnormal.Gap,
                RawAccept = false, GuardedAccept = subnormal.ShouldSwitch },
            NormalRangeFalsePositive = new { InputBits = new[] { "C00D674B0EF52CA9", "400BAF7E83D81B4B", "BF5AAB2B125E02F1" },
                ExactScaledGap = boundaryExact.Gap.ToString(), ExactScaledGapApproximate = boundaryExact.Gap.ToApproximateDouble(),
                ComputedScaledGap = boundary.Gap, RawAccept = true, GuardedAccept = boundary.ShouldSwitch },
            GroupedOverflow = new { GroupedIsInfinite = true, DistributedIsFinite = true },
            InvalidInputCases = invalid.Length
        };
    }

    private static void CheckRationalArithmetic()
    {
        Require(Rational.FromDouble(0.1) == new Rational(3602879701896397, 36028797018963968), "Exact conversion of 0.1.");
        Require(Rational.FromDouble(double.Epsilon) == new Rational(1, BigInteger.One << 1074), "Subnormal conversion.");
        Require(Rational.FromDouble(Binary64.MinimumNormal) == new Rational(1, BigInteger.One << 1022), "Normal conversion.");
        Require(Rational.FromDouble(double.MaxValue) == new Rational((BigInteger.One << 1024) - (BigInteger.One << 971)), "Maximum finite conversion.");
        Require(Rational.FromDouble(-0.0) == Rational.Zero, "Signed-zero conversion.");
        Rational third = new(1, 3), sixth = new(1, 6);
        Require(third + sixth == new Rational(1, 2) && third / sixth == 2 && third * 3 == 1,
            "Rational arithmetic.");
        Require(new Rational(-2, -6) == third && third - sixth == sixth, "Rational normalization.");
    }

    private static double FromBits(string hex) => BitConverter.UInt64BitsToDouble(ulong.Parse(hex, NumberStyles.HexNumber, CultureInfo.InvariantCulture));
    private static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidOperationException(message);
    }
}
