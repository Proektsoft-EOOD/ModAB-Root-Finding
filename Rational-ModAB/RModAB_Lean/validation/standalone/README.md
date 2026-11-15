# Standalone package verification — 14 November 2026

This directory contains the fresh checks of the independent RModAB Lean project. The supplied 11 November proof-release logs remain in the parent directory as historical evidence.

- [Kernel/build verification and environment](verification.json) records the pinned Lean version, dependency revisions, commands and exit codes.
- [Packaging/source-integrity check](packaging_check.json) verifies the preserved Lean files and dependency pins, local links and correspondence table. It is a structural check, not a Lean kernel run.
- [Source manifest](../../docs/standalone-source-manifest.json) records the exact input archive and the 72 preserved Lean/configuration files.

The project uses Lean 4.34.0 and the original Mathlib/FloatLib revisions. No theorem body, module name, dependency pin or mathematical hypothesis was changed for this packaging. Existing linter/deprecation warnings may occur; the proof files are retained unchanged.

The checks are local executions. They do not claim a browser run at Lean Live, GitHub publication, or a complete formal verification of the C# binary64 implementation. The standalone package does not include the local toolchain or downloaded dependency caches; normal installation retrieves them using the pinned configuration.

## Fresh execution results

| Check | Result | Seconds | Log |
|---|---|---:|---|
| modular build | PASS | 207.6 | [build.log](build.log) |
| axiom audit | PASS | 4.5 | [audit.log](audit.log) |
| focused live | PASS | 38.3 | [live_focused.log](live_focused.log) |
| complete live | PASS | 174.7 | [live_complete.log](live_complete.log) |
| compiled inventory | PASS | 12.6 | [inventory_compiled.jsonl](inventory_compiled.jsonl) |

The modular build and both single-file checks completed successfully with the original pins. The 800-statement full audit and the 420-statement focused audit passed.
