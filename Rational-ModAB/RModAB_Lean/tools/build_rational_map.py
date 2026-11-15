#!/usr/bin/env python3
"""Map the standalone RModAB paper to proof release 0.8.0.
The earlier ModAB article, all its anchors and the draft R1-R6 aliases survive.
This documentation generator does not alter Lean sources or proof counts.
"""
from pathlib import Path
import csv,json,re,html,hashlib,shutil,subprocess
root=Path(__file__).resolve().parents[1]
esc=html.escape
# A standalone repository has no assigned public URL yet. Current-project
# links are local so this map works in a checkout and after relocation.
MAP_FILE='RModAB_Proof_Map.html'
PREDECESSOR_URL='https://github.com/Proektsoft-EOOD/ModAB-Root-Finding/tree/ea64a126c9a06294733086824da95a0eec808051/MDPI-Algorithms-Data-2/ModAB.Switching.Lean'
INSPECTED_ALGORITHM_URL='https://github.com/Proektsoft-EOOD/ModAB-Root-Finding/tree/64a95b1dce05551feaa2934e45d37fee23a9b2e8'
legacy=json.loads((root/'docs/legacy-proof-inventory.json').read_text(encoding='utf-8'))
new=json.loads((root/'docs/rational-source-index.json').read_text(encoding='utf-8'))
compiled={}
for line in (root/'validation/inventory_compiled.jsonl').read_text(encoding='utf-8').splitlines():
    if line.startswith('{'):
        row=json.loads(line);compiled[row['name']]=row
rows=legacy+new
assert set(compiled)=={r['name'] for r in rows}
for row in rows:row.update(compiled[row['name']])
groups=[
 {'number':'1','kind':'Proposition','aliases':['R1'],'label':'prop:rational_algebra_limit','title':'Rational algebra and bounded displacement','status':'Proved algebra and interval bound',
  'hypotheses':'a < b; opposite nonzero working ordinates; stated nonzero divisions for interpolation equations. Function-value interpolation additionally needs nonzero model denominators.',
  'proof':'Normalize the three data points, verify the cross-multiplied equations algebraically, then use the existing secant interiority theorem and d = min(s-a,b-s). Projection gives the 0.3d margins.',
  'names':['ModAB.Rational.normalized_interpolation_equations','ModAB.Rational.secant_normalized','ModAB.Rational.limited_bounds','ModAB.Rational.limited_margins','ModAB.Rational.limited_interior','ModAB.Rational.candidate_interior','ModAB.Rational.numerator_ne_zero_of_denominator_zero']},
 {'number':'2','kind':'Proposition','aliases':['R2'],'label':'prop:rational_closed_update','title':'Closed proposals and guarded closing','status':'Proved; source branches assembled',
  'hypotheses':'Sign-changing bracket; candidate in the closed interval; nonzero true evaluation; sign update. Positive supplied tolerance for closing.',
  'proof':'A collapsed retained interval would contradict the strict sign product. SourceProposal implies closed membership; guarded closing preserves it. The operational assembly derives candidate membership from these facts.',
  'names':['ModAB.RationalBridge.updated_invariants_closed','ModAB.Rational.sourceProposal_mem_closed','ModAB.Rational.closeCandidate_mem_closed','ModAB.Rational.closeCandidate_interior','ModAB.RationalAssembly.source_closed_update','ModAB.RationalAssembly.source_candidate_mem']},
 {'number':'3','kind':'Theorem','aliases':['R3'],'label':'thm:rational_capped','title':'Inheritance of the capped width theorem','status':'Proved by reuse of the original phase/schedule proof',
  'hypotheses':'Continuous f on initial sign bracket; closed proposals and true updates; original capped Step B M with initial bisection and B >= 1/2.',
  'proof':'Reuse the existing width Run, global width bound, finite tolerance theorem and nested-interval theorem. Endpoint reuse consumes a step; it does not reset the clock. B=2, M=3 yields Wn <= W0 / 2^floor(n/7), an amortized bound from initialization.',
  'names':['ModAB.RationalBridge.ClosedExactRun.widthRun','ModAB.RationalBridge.closed_invariants','ModAB.RationalBridge.closed_zeros','ModAB.RationalBridge.capped_global_convergence','ModAB.RationalBridge.capped_global_width','ModAB.RationalBridge.capped_width_seven','ModAB.RationalBridge.capped_finite_tolerance']},
 {'number':'4','kind':'Theorem','aliases':['R4'],'label':'thm:rational_uncapped_alternative','title':'Uncapped width–residual alternative','status':'Proved from transitions; not a model assumption',
  'hypotheses':'Nonnegative widths/residuals; nonincreasing widths; real bisections; advancing phase clock; continuation only through width success or residual decrease; minimum <= previous residual; 0 <= eta < 1. No continuity needed.',
  'proof':'If widths are arbitrarily small, monotonicity implies their convergence to zero. Otherwise their positive lower bound forces finitely many bisections. The threshold in the last phase tends to zero, so eventually every step satisfies R(n+1) < eta R(n).',
  'names':['ModAB.RationalUncapped.decision_step','ModAB.RationalUncapped.continuing_suffix','ModAB.RationalUncapped.bisection_count_bound','ModAB.RationalUncapped.width_or_residual','ModAB.RationalUncapped.finite_combined_tolerance','ModAB.RationalAssembly.source_policy']},
 {'number':'5','kind':'Theorem','aliases':['R5'],'label':'thm:rational_uncapped_roots','title':'Residual, zero-set and unique-root convergence','status':'Proved for the continuous sign-bracketed model',
  'hypotheses':'Theorem 4 (legacy alias R4) plus continuity on the compact initial interval and true sign updates. Uniqueness is needed only for convergence to one specified root.',
  'proof':'In the width branch reuse the old nested-interval theorem and continuity. In the residual branch no further estimate is needed. A compact-subsequence contradiction proves eventual proximity to a true zero. Uniqueness turns that proximity into point convergence.',
  'names':['ModAB.RationalUncapped.global_residual_convergence','ModAB.RationalUncapped.eta42_global_residual','ModAB.RationalUncapped.finite_residual_tolerance','ModAB.RationalResidual.limit_is_zero','ModAB.RationalResidual.eventually_near_zero_set','ModAB.RationalResidual.tendsto_unique_zero','ModAB.RationalAssembly.source_global_residual','ModAB.RationalAssembly.source_unique_root_convergence']},
 {'number':'6','kind':'Proposition','aliases':['R6'],'label':'prop:rational_closing_termination','title':'Width stopping from the actual closing operation','status':'Proved with dynamic tolerances; spacing is derived',
  'hypotheses':'Operational SourceRun; unique zero on the initial interval; 0 <= eta < 1; tolerance(n) >= tau0 > 0. Exact real arithmetic, no claimed 200-visit bound.',
  'proof':'SourceRun.exactRun derives membership and schedule from the step equations. Infinitely many bisections force a width stop. Otherwise the interpolation tail converges to the unique root, while failure of the width test makes closing enforce steps of length at least tau0/2. This contradicts point convergence.',
  'names':['ModAB.RationalAssembly.SourceRun.exactRun','ModAB.RationalAssembly.source_closing_separation','ModAB.RationalUncapped.eventually_short_steps','ModAB.RationalUncapped.no_convergent_separated_tail','ModAB.RationalAssembly.source_finite_width_tolerance','ModAB.RationalAssembly.eta42_source_finite_width']}
]
groups.extend(json.loads((root/'docs/local-map-groups.json').read_text(encoding='utf-8')))
index={r['name']:r for r in rows}
for g in groups:
    for name in g['names']:
        assert name in index,name
        if g['label'] not in index[name]['paper_labels']:index[name]['paper_labels'].append(g['label'])
theorems=[r for r in rows if r['kind'] in ('theorem','lemma')]
allowed={'propext','Classical.choice','Quot.sound'}
for r in rows:assert set(r['axioms'])<=allowed,(r['name'],r['axioms'])
(root/'docs/proof-inventory.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False), encoding='utf-8')
(root/'docs/rational-map.json').write_text(json.dumps(groups,indent=2,ensure_ascii=False), encoding='utf-8')
with (root/'docs/Proof_Correspondence.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['Paper result','Preserved alias','LaTeX label','Title','Hypotheses','Proof correspondence','Lean declarations'])
    for g in groups:
        writer.writerow([g['kind']+' '+g['number'],','.join(g['aliases']),g['label'],g['title'],g['hypotheses'],g['proof'],'; '.join(g['names'])])

dot='''digraph G {
rankdir=TB; graph [bgcolor="transparent",pad="0.2",nodesep="0.3",ranksep="0.55"];
node [shape=box,style="rounded,filled",fillcolor="#eef5f6",color="#678a95",fontname="DejaVu Sans",fontsize=13,margin="0.15,0.12"];
edge [color="#567581",arrowsize=0.7];
old [label="Existing secant and sign lemmas"];
rational [label="Propositions 1–2: rational limit / closing\nclosed-interval admissibility"];
phase [label="Existing capped phase / schedule"];
capped [label="Theorem 3: capped width bound\nWn <= W0 / 2^floor(n/7)"];
run [label="SourceRun.exactRun\nmembership and schedule derived"];
alternative [label="Theorem 4: width OR residual tends to zero"];
roots [label="Theorem 5: residual / zero-set convergence\nunique root => point convergence"];
stop [label="Proposition 6: dynamic width tolerance reached\nunique root + positive tolerance floor"];
old -> rational; rational -> capped; phase -> capped; rational -> run;
run -> alternative; alternative -> roots; roots -> stop; rational -> stop;
localalgebra [label="Proposition 7: mixed-ordinate error identity"];
localanalysis [label="C3 regularity / divided-difference limits\\nChecked true-data error coefficient"];
localorder [label="Lemma 8: tribonacci Q-order\\nChecked stable logarithmic recurrence"];
localbounds [label="Uniform C3 local estimates\\nClamp / closing / residual42 decision"];
localassembly [label="Theorem 9: actual guarded state orbit\\nHistory invariance / tribonacci Q-order"];
cycle [label="Proposition 10: uniform interpolation cone\\nContraction / separation / pair recurrence"];
cycleorder [label="Corollary 11: alternating orders 9/5, 5/3\\nCycle-normalized sqrt(3)"];
onesided [label="Proposition 12: actual AB tail\\nPositive factor limit / golden-ratio order"];
rational -> localalgebra; localalgebra -> localanalysis;
localanalysis -> localbounds [label="proved C3 neighbourhood"];
localbounds -> localassembly; localorder -> localassembly;
localalgebra -> cycle [label="proved uniform bounds"];
localanalysis -> cycle; localbounds -> cycle;
cycle -> cycleorder; localalgebra -> onesided;
{rank=same;rational;phase;} {rank=same;capped;run;}
}'''
(root/'docs/rational-dependencies.dot').write_text(dot, encoding='utf-8')
svgpath=root/'docs/rational-dependencies.svg'
if shutil.which('dot'):
    svg=subprocess.check_output(['dot','-Tsvg'],input=dot.encode()).decode()
    svg=svg[svg.index('<svg'):]
    svgpath.write_text(svg, encoding='utf-8')
elif svgpath.is_file():
    # The bundled overview is sufficient for a normal documentation refresh.
    # Install Graphviz when changing the graph specification itself.
    svg=svgpath.read_text(encoding='utf-8')
else:
    raise SystemExit('Graphviz dot is required to create the missing dependency overview.')
parts=['<section id="rational-extension"><span id="standalone-rmodab-paper"></span><p class="meta">STANDALONE RMODAB REPOSITORY · PROOF RELEASE 0.8.0 · PACKAGED 14 NOVEMBER 2026</p><h1>RModAB paper ↔ Lean proof map</h1>',
 f'<p class="lead">{len(new)} new named declarations; {len(theorems)-380} new proved statements. The complete project audit covers {len(theorems)} statements.</p>',
 '<p class="note"><strong>Separate follow-up paper.</strong> The standalone RModAB manuscript preserves results 1–6 and adds local-analysis results 7–12. The earlier ModAB revision is a separate submitted article and is not renumbered or replaced. Draft aliases R1–R6, their links and all historical article anchors remain valid. Aliases L1–L6 identify the new local results. Release 0.8.0 completes the local C3-to-state and interpolation-to-uniform-cone assemblies, and derives the one-sided scaling bound. Each entry states its exact mathematical hypotheses and scope.</p>',
 '<p><a href="README.md">Project instructions</a> · <a href="ModAB.lean">Modular proof entry point</a> · <a href="RModAB_Live.lean">Focused single-file proof</a> · <a href="ModAB_Live_Complete.lean">Complete single-file proof</a> · <a href="docs/Proof_Correspondence.csv">Paper correspondence table</a>. These links are relative to this independent project. Its future public repository address has not yet been assigned.</p>',
 f'<p>Historical references: <a href="{PREDECESSOR_URL}">separate predecessor ModAB proof (immutable revision)</a> · <a href="{INSPECTED_ALGORITHM_URL}">original RModAB source inspected for this proof release</a>. The predecessor is not replaced by this package.</p>',
 '<p>The algorithm was extracted from the supplied description and checked against the public source. The proofs develop the earlier phase, bracket and compactness arguments. The local arguments are independently developed from the algorithm and the coauthor’s analysis; their exact Lean coverage is recorded below. Modular and single-file proofs were compiled with Lean 4.34.0; the original compiler records are preserved.</p>',
 '<p><strong>Scope:</strong> the global results retain their exact-real SourceRun hypotheses. The connected local theorems derive uniform estimates and actual state transitions from C3, the simple root and the stated initial geometry. Both discarded-side histories, AB updates, true-ordinate curvature, ordered candidate, zero closing and uncapped continuation are included. Uniform negative pair coefficients and the positive limiting one-sided weight are derived. The local assembly boundaries of release 0.7.0 are closed. Universal entry from arbitrary brackets, the degenerate C=0 boundary, full binary64 program equivalence and a universal 200-visit bound are not claimed.</p>',
 '<div class="diagram">'+svg+'</div><p class="note">All arrows in this overview connect checked results. The local C3-to-state and interpolation-to-cone links are proved in this release; individual declaration dependencies are listed below.</p>',
 '<div class="tablewrap"><table><tr><th>Standalone paper result</th><th>Preserved alias</th><th>Statement</th><th>Formal route</th></tr>']
for g in groups:
    parts.append(f'<tr><td><a href="#standalone-paper-{g["number"]}">{g["kind"]} {g["number"]}</a></td><td><a href="#rational-paper-{g["aliases"][0]}">{g["aliases"][0]}</a></td><td>{esc(g["title"])}</td><td>{esc(g["status"])}</td></tr>')
parts.append('</table></div>')
for g in groups:
    parts.append(f'<article class="entry" id="rational-paper-{g["aliases"][0]}"><h3 id="standalone-paper-{g["number"]}">{g["kind"]} {g["number"]}. {esc(g["title"])}</h3><p>Map alias: <code>{g["aliases"][0]}</code> · <code>{g["label"]}</code> · <span class="tag">{esc(g["status"])}</span></p><p><strong>Hypotheses.</strong> {esc(g["hypotheses"])}</p><p><strong>Proof correspondence.</strong> {esc(g["proof"])}</p><p>')
    parts.append(', '.join(f'<a href="#decl-{n}"><code>{n}</code></a>' for n in g['names']))
    parts.append('</p></article>')
parts.append('<h3>New declaration inventory and compiled dependencies</h3><p>Dependencies below are extracted from the elaborated Lean environment, not inferred from the prose. Model structures and definitions are not counted as proved theorems.</p>')
for r in new:
    rr=index[r['name']];deps=', '.join(f'<a href="#decl-{n}"><code>{n}</code></a>' for n in rr['dependencies']) or 'No authored declaration dependencies.'
    parts.append(f'<article class="entry" id="decl-{r["name"]}"><header><h3><code>{r["name"]}</code></h3><span class="tag">{r["kind"]}</span></header><p><a href="{esc(r["file"])}#L{r["line"]}"><code>{r["file"]}</code></a>, line {r["line"]}</p><p>{deps}</p><details><summary>Formal declaration and proof</summary><pre>{esc(r["source"])}</pre></details><p class="meta">Axioms: {esc(", ".join(rr["axioms"]) or "none")}.</p></article>')
parts.append('<h3>Boundaries of the source correspondence</h3><p>The inspected source commit is <code>64a95b1dce05551feaa2934e45d37fee23a9b2e8</code>. The explicit SourceRun equations model midpoint proposals, post-limit closing, true sign updates and the uncapped nextMode decision. Positive working signs and correct last-endpoint history are stated premises, supported by separate state lemmas. This is not a machine-checked equivalence of the full C# binary64 program. Its 200-visit limit and function evaluation accuracy are not proved by these theorems.</p></section><hr><h2>Preserved map of the separate original ModAB article</h2><p class="note">The following entries retain the original 0.5.0 statements, numbers, labels, anchors and proof excerpts. They belong to the previously submitted ModAB revision, not to results 1–12 of the standalone RModAB paper above. Their source locations and historical manifest identify that separate article. Current source hashes and the combined inventory are in <code>docs/proof-manifest.json</code> and <code>docs/proof-inventory.json</code>.</p>')
legacyhtml=(root/'docs/legacy-proof-map.html').read_text(encoding='utf-8')
legacyhtml=legacyhtml.replace('<title>ModAB — manuscript and Lean proof map</title>','<title>Standalone RModAB paper — proof map 0.8.0 with preserved ModAB map</title>')
legacyhtml=legacyhtml.replace('<main>','<main>'+''.join(parts),1)
(root/MAP_FILE).write_text(legacyhtml, encoding='utf-8')
# Preserve all old anchors for existing private bookmarks without redirecting
# the browser to a different repository or dropping its URL fragment.
(root/'ModAB_Proof_Map.html').write_text(legacyhtml, encoding='utf-8')

manifest={'project_version':'0.8.0','date':'2026-11-11','lean':'4.34.0',
    'audited_theorems':len(theorems),'authored_declarations':len(rows),
    'new_theorems':len(theorems)-380,'new_declarations':len(new),
    'previous_release':'0.7.0','added_local_theorems':len(theorems)-472,
    'preserved_previous_proof_modules':46,'local_modules':len({v['module'] for v in new})-6,
    'added_this_release':len(theorems)-667,'local_correspondence_complete':True,
    'legacy_proof_step_records':91,'legacy_manuscript_entries':23,'rational_manuscript_entries':len(groups),
    'standalone_manuscript_entries':len(groups),'manuscript_role':'Standalone RModAB follow-up paper; previous ModAB article remains separate',
    'standalone_results':[g['kind']+' '+g['number'] for g in groups],
    'preserved_aliases':{g['aliases'][0]:g['kind']+' '+g['number'] for g in groups},
    'documentation_refresh':'Completed local state assemblies, uniform negative coefficients and one-sided scaling; all 12 statement numbers preserved',
    'repository_url':None,'project_url':None,'proof_map_url':None,
    'proof_map_file':MAP_FILE,'packaging_date':'2026-11-14',
    'packaging_role':'Independent repository; no overwrite or deployment of the predecessor',
    'predecessor_proof_url':PREDECESSOR_URL,'inspected_algorithm_url':INSPECTED_ALGORITHM_URL,
    'scope':'Inherited global proofs; exact-real rational/uncapped model, stopping, local error algebra, divided-difference limits, recurrence Q-orders, pair-cone dynamics and safeguard estimates. Local C3-to-state, interpolation-to-cone and one-sided scaling assemblies are proved; no complete C# IEEE proof.',
    'source_commit':'64a95b1dce05551feaa2934e45d37fee23a9b2e8','files':{}}
files=list((root/'ModAB').glob('*.lean'))+[root/n for n in ['ModAB.lean','Audit.lean','ExportInventory.lean','ModAB_Live_Complete.lean','RModAB_Live.lean','lakefile.toml','lake-manifest.json','lean-toolchain']]
for f in files:manifest['files'][str(f.relative_to(root))]=hashlib.sha256(f.read_bytes()).hexdigest()
(root/'docs/proof-manifest.json').write_text(json.dumps(manifest,indent=2), encoding='utf-8')
print(json.dumps({k:v for k,v in manifest.items() if k!='files'},indent=2))
