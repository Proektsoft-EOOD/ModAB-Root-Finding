# ModAB Lean Project

Version 0.4.0. This is the standard multi-file Lean/Lake project.

## Contents

The `ModAB/` modules contain the switching, binary64 arithmetic, unscaled error,
cubic bound, exact-arithmetic convergence, width bounds, residual-to-distance
estimate, and exact Anderson-Bjoerck examples. `ModAB.lean` is the main import.
`Audit.lean` checks the transitive axioms of all 287 project lemmas and theorems.

## Build and check

Install Lean through elan and install Git. Open `ModAB.Switching.Lean` as the
project directory, or run the following commands from that directory:

```sh
lake exe cache get
lake build
lake env lean Audit.lean
```

The `lean-toolchain` file selects Lean 4.34.0. `lakefile.toml` and
`lake-manifest.json` pin FloatLib, mathlib, and the transitive dependencies.
An internet connection is needed for the initial dependency installation and
cache download. Compiler binaries and downloaded dependency/build directories
are not included in the source archive.

The audit accepts only `propext`, `Classical.choice`, and `Quot.sound`.
An unexpected axiom, including `sorryAx`, fails the audit.
The final successful message is:

```text
PASS: 287 declarations; only standard logical axioms.
```

The formalization proves its stated mathematical results under their explicit
hypotheses. Global convergence and the example trajectories concern exact real
arithmetic. The binary64 switching results do not establish convergence of the
entire floating-point root-finding loop or compiler/hardware correctness.

All Lean proof modules, `Audit.lean`, and dependency configuration files are
unchanged from the checked version 0.4.0. The FloatLib MIT license is included
in `THIRD_PARTY_NOTICES/FloatLib-LICENSE.txt`.
