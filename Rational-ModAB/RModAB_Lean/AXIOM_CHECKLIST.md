# Axiom and correspondence audit — release 0.8.0

Validated locally on 11 November 2026 with Lean 4.34.0 and the pinned dependencies.

- Modular build and full audit: 800 proved statements.
- Focused single-file audit: 420 rational-extension statements.
- Complete single-file audit: 800 statements.
- Added in this release: 133 statements in 18 modules.
- Compiled inventory: 1027 authored declarations, including definitions.
- Paper map: 12 results; all preceding numbers and aliases retained.
- Previous 46 proof modules and dependency pins: byte-identical.
- Allowed axioms: propext, Classical.choice, Quot.sound only.
- No admitted proofs, new axioms or native_decide in authored proofs.

| Rational-extension module | Proved statements |
|---|---:|
| Rational | 38 |
| RationalBridge | 7 |
| RationalSpacing | 3 |
| RationalResidual | 3 |
| RationalUncapped | 29 |
| RationalAssembly | 12 |
| RationalLocalAlgebra | 33 |
| RationalOrderStability | 7 |
| RationalOrder | 16 |
| RationalGoldenOrder | 6 |
| RationalLocalSafeguards | 15 |
| RationalCycle | 14 |
| RationalCycleEntry | 18 |
| RationalDSlope | 15 |
| RationalLocalAnalysis | 18 |
| RationalRegularity | 24 |
| RationalAnalyticBridge | 5 |
| RationalSmoothBridge | 5 |
| RationalOneSided | 13 |
| RationalOneSidedAssembly | 4 |
| RationalLocalAssembly | 2 |
| RationalLocalNeighborhood | 5 |
| RationalLocalState | 7 |
| RationalAlternating | 15 |
| RationalCurvature | 4 |
| RationalCycleUniform | 15 |
| RationalCycleState | 6 |
| RationalCycleUniformAlgebra | 15 |
| RationalCycleUniformBridge | 8 |
| RationalSecantLocal | 11 |
| RationalCycleSafeguards | 6 |
| RationalCycleNeighborhood | 1 |
| RationalCycleIntermediate | 6 |
| RationalLocalComplete | 2 |
| RationalCycleStateGuards | 5 |
| RationalCycleStateUpdate | 3 |
| RationalCycleTransport | 5 |
| RationalCycleComplete | 7 |
| RationalOneSidedComplete | 12 |
| Total rational extension | 420 |

Audit.lean checks transitive axioms and rejects sorryAx. The connected local
theorems derive their analytic estimates and operational state invariants
from the paper's stated hypotheses. The previous Theorem9, Proposition10 and
Proposition12 correspondence boundaries are closed. The scope remains local
and exact-real; full IEEE implementation equivalence, automatic entry from
every bracket and the C=0 classification are not asserted. Browser execution
in Lean Live is not part of the recorded local checks.

See README.md, the proof map, and validation/local_completion_validation.json.

Standalone repackaging was recompiled on 14 November 2026 with the same pins. See [fresh verification](validation/standalone/README.md); proof sources are unchanged.
