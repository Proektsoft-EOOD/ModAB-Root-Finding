# Rational ModAB — standalone Lean project

Formal proofs accompanying **Rational ModAB: safeguarded interpolation, global convergence and local rates**, by **Latchezar Tomov and Nedelcho Ganchovski**.

This is an independent Lake project for the new RModAB repository. It contains proof release **0.8.0**, repackaged on **14 November 2026**. Upload this directory's contents to the new repository root, or retain the entire directory as a separate Lake project. The old ModAB repository is not an installation target.

The destination GitHub address is pending. All current-project links are relative. The article and its public links will be updated separately once the address is known.

## Start here

- [Interactive paper–Lean map](RModAB_Proof_Map.html)
- [English correspondence table: 12 paper results](docs/Proof_Correspondence.csv)
- [Focused single-file edition](RModAB_Live.lean)
- [Complete single-file edition](ModAB_Live_Complete.lean)
- [Axiom checklist](AXIOM_CHECKLIST.md)
- [Standalone validation record](validation/standalone/README.md)
- [Bulgarian instructions](README_BG.txt)

Extract the entire project, then open the HTML map in a browser so its relative source links resolve. GitHub's normal file view displays its source; optional GitHub Pages hosting can be configured later. [ModAB_Proof_Map.html](ModAB_Proof_Map.html) is a compatibility copy retaining the same anchors.

## Coverage

All **12 numbered results** in the accompanying paper have formal correspondence entries. The proof release audits **800 lemmas/theorems**: **380 inherited statements** and **420 rational-extension statements**. The compiled inventory contains 1027 named authored declarations when definitions and structures are counted separately.

| Paper result | Connected local theorem |
|---|---|
| Theorem 9 | `ModAB.Rational.LocalComplete.local_alternating_complete` |
| Proposition 10 and Corollary 11 | `ModAB.Rational.CycleComplete.local_pair_complete` |
| Proposition 12 | `ModAB.Rational.OneSidedComplete.smooth_rawPoint_one_sided_hasQOrder` |

The map records every statement and its hypotheses. The mathematical results concern the exact-real model. Complete C# binary64 equivalence, automatic entry into each local regime from every initial bracket, and a universal 200-iteration guarantee are outside these statements. The one-sided result retains the paper's convergent-tail and nonzero-coefficient assumptions. See [IMPLEMENTATION_SCOPE.md](IMPLEMENTATION_SCOPE.md) for the two experimental C# guards.

## Reproduce

Install Git, Python 3.9 or later, and Lean through elan. Use **Lean 4.34.0**, as pinned in `lean-toolchain`, and the exact Mathlib/FloatLib revisions in `lake-manifest.json`. Initial dependency installation requires network access. From this directory run (on Windows, `python` may replace `python3`):

```sh
lake exe cache get
lake build ModAB
lake env lean Audit.lean
lake env lean RModAB_Live.lean
lake env lean ModAB_Live_Complete.lean
python3 tools/check_audit_coverage.py
python3 tools/validate_standalone.py
```

`Audit.lean` checks transitive axioms and rejects `sorryAx`. The allowed logical axioms are `propext`, `Classical.choice` and `Quot.sound`. Structural and source-integrity checks are distinct from Lean kernel compilation; the validation record identifies each check.

For Lean Live, use `RModAB_Live.lean` in a matching Lean/Mathlib environment. This file includes the rational proofs and their exact-real dependencies. The complete edition also includes the inherited binary64/FloatLib material. Browser execution is not claimed by the local validation records.

## Independent project and inherited modules

The internal package name `modab_switching`, library name `ModAB`, and namespaces `ModAB` / `ModAB.Rational` are retained to preserve imports and theorem names. They do not determine the GitHub repository name. All inherited authored modules needed by the new proofs are included locally.

Build this directory as its own Lake root; no sibling checkout of the old project is required. Third-party Lean dependencies still need installation. Authored Lean sources and dependency pins are compared byte-for-byte with proof release 0.8.0.

The predecessor proof remains a historical reference at [its immutable revision](https://github.com/Proektsoft-EOOD/ModAB-Root-Finding/tree/ea64a126c9a06294733086824da95a0eec808051/MDPI-Algorithms-Data-2/ModAB.Switching.Lean). The inspected algorithm baseline is recorded separately in [RModAB_Source_Reference.json](RModAB_Source_Reference.json).

## Regenerate

To rebuild the map from the included compiled inventory:

```sh
python3 tools/build_rational_map.py
python3 tools/check_audit_coverage.py
python3 tools/validate_standalone.py
```

After changing proof declarations, regenerate the audit and single-file sources, run Lean, and export a fresh inventory before rebuilding the map:

```sh
python3 tools/prepare_rational_release.py
lake build ModAB
lake env lean Audit.lean
lake env lean ExportInventory.lean > validation/inventory_compiled.jsonl
python3 tools/build_rational_map.py
```

Historical validation records retain their original dates and paths. Current installation instructions and source links use this standalone layout. Third-party notices remain in `THIRD_PARTY_NOTICES/`.
