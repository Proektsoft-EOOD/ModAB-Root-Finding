Fixed binary64 verification data

cases.csv contains 23,816 ordered triples. Each field is the 16-digit
hexadecimal IEEE 754 binary64 bit pattern of one input; the order is f1,f2,f3.
There is one header row. The file uses ASCII and LF line endings.
Negative zero is preserved. Hexadecimal bit patterns avoid decimal parsing
and locale-dependent changes of the inputs.

SHA-256:
cede3121d4b54f765a2cb7dc1992a9f0c68a9710efd6ad95b5fd4bfbe461e4bb

Rows 1-2016: a Cartesian family of 12 positive endpoint magnitudes,
7 midpoint choices, and both endpoint sign orientations. Magnitudes include
the least subnormal, twice that value, the least normal, 2^-500, 2^-53,
1/2, 1, the predecessor of 1, 2, 2^500, half the maximum finite value,
and the maximum finite value.

Rows 2017-8016: random finite nonzero binary64 magnitudes, with endpoint
signs negative and positive and with both signs for the midpoint residual.
The distribution is over bit patterns, not uniform over real numbers.

Rows 8017-13016: ordinary endpoint magnitudes between 0.1 and 4, with
midpoint residuals inside the exact switching region.

Rows 13017-23816: both exact switching boundaries for 1,800 endpoint
pairs, rounded to binary64, together with their predecessor and successor.

The fixed dataset was constructed with Python's random.Random(20261004)
and exact fractions for boundary values. The C# verification program reads
the recorded bit patterns and computes the mathematical reference afresh
using BigInteger rational arithmetic; Python is not required to run it.
The sample deliberately concentrates on difficult cases. Its counts do
not estimate the frequency of rounding disagreements in a root solver.
