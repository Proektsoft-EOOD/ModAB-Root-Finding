# ModAB Lean Project

Version 0.5.0. This is the standard multi-file Lean/Lake project.

## Contents

The `ModAB/` modules contain the switching, binary64 arithmetic, unscaled error,
cubic bound, exact-arithmetic convergence, width bounds, residual-to-distance
estimate, and exact Anderson-Bjoerck examples. `ModAB.lean` is the main import.
`Audit.lean` checks the transitive axioms of all 380 project lemmas and theorems.

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
PASS: 380 declarations; only standard logical axioms.
```

The formalization proves its stated mathematical results under their explicit
hypotheses. Global convergence and the example trajectories concern exact real
arithmetic. The binary64 switching results do not establish convergence of the
entire floating-point root-finding loop or compiler/hardware correctness.

Version 0.5.0 adds `TopEntry.lean`, `TopOneSided.lean`, `TopUncapped.lean`, and
`Eta42.lean` (93 lemmas/theorems). All 21 earlier proof modules and dependency
pins are unchanged. The proof map follows the revised manuscript and S1.
See `AXIOM_CHECKLIST.md` for the full check procedure and exact proof boundary. The FloatLib MIT license is included
in `THIRD_PARTY_NOTICES/FloatLib-LICENSE.txt`.

## New exact-arithmetic coverage

- `TopEntry.lean`: the switching midpoint records its actual moved side;
  either next-side branch preserves working signs, and the width policy
  agrees with the existing `Step.bisectStart` model.
- `TopOneSided.lean`: the top hybrid variant of Proposition 2, including
  errors 1/3, 1/15, 1/255, 1/65535 and quadratic convergence.
- `TopUncapped.lean`: the full top trajectory of f101 (Proposition 3), using
  the original invariant-region proof and proving uncapped continuation.
- `Eta42.lean`: C1 regularity, strict monotonicity, the unique simple root,
  actual secants, exact closed form, positive corrections, quadratic point
  convergence, width tending to 1, and the six-chord capped schedule.

`eta42Error 0` is the first chord, whereas `eta42FullError 0` is the endpoint
following the switching bisection. The manuscript numbers true chords from 1.
The cap schedules fallback after chord 6; the following visit is bisection.
The theorem does not claim that chord 6 itself is a bisection.

For the bundled single-file edition:

```sh
lake env lean ModAB_Live_Complete.lean
```

It is compiled locally with the pinned toolchain; an interactive Lean Live
browser session is not part of this check.
