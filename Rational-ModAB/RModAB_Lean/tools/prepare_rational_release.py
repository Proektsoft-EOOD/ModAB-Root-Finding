#!/usr/bin/env python3
"""Prepare proof release 0.8.0 for the standalone RModAB follow-up paper.
Regenerate explicit audit/export lists and dependency-ordered live files.
This script prepares sources; kernel checks are separate commands. A paper-only
numbering refresh uses build_rational_map.py instead and does not run this file.
Run from any directory: all inputs belong to this independent project root.
The Lake package identity and dependency pins are deliberately left unchanged.
"""
from pathlib import Path
import json,re

root=Path(__file__).resolve().parents[1]
modules=['Rational','RationalBridge','RationalSpacing','RationalResidual',
         'RationalUncapped','RationalAssembly','RationalLocalAlgebra',
         'RationalOrderStability','RationalOrder','RationalGoldenOrder',
         'RationalLocalSafeguards','RationalCycle','RationalCycleEntry',
         'RationalDSlope','RationalLocalAnalysis','RationalRegularity','RationalAnalyticBridge',
         'RationalSmoothBridge',
         'RationalOneSided','RationalOneSidedAssembly',
         'RationalLocalAssembly']

# Extend the existing proof sequence with every new rational module, ordered
# by actual imports. This preserves the old proof bodies and supports the
# completed analytic/state assemblies without a manually duplicated list.
module_candidates = modules + sorted(p.stem for p in (root/'ModAB').glob('Rational*.lean')
    if p.stem not in modules)
module_seen = set()
module_order = []
def order_module(mod):
    if mod in module_seen:
        return
    module_seen.add(mod)
    src = (root/'ModAB'/f'{mod}.lean').read_text(encoding='utf-8')
    for dep in re.findall(r'^import ModAB\.(Rational[\w]+|Rational)$', src, re.M):
        order_module(dep)
    module_order.append(mod)
for mod in module_candidates:
    order_module(mod)
modules = module_order
rows=[]
for mod in modules:
    path=root/'ModAB'/f'{mod}.lean';s=path.read_text(encoding='utf-8')
    ns=re.search(r'^namespace ([\w.]+)',s,re.M).group(1)
    matches=list(re.finditer(r'^(?:@\[[^\n]*\][ \t]*)*(?:noncomputable )?(def|structure|inductive|theorem|lemma)\s+([\w.]+)',s,re.M))
    for i,m in enumerate(matches):
        end=matches[i+1].start() if i+1<len(matches) else s.rfind('\nend')
        rows.append(dict(name=ns+'.'+m.group(2),kind=m.group(1),module=mod,
            source=s[m.start():end].strip(),file=f'ModAB/{mod}.lean',
            line=s[:m.start()].count('\n')+1,paper_labels=[]))
(root/'docs/rational-source-index.json').write_text(json.dumps(rows,indent=2), encoding='utf-8')

def replace_names(filename,variable,newnames):
    path=root/filename;s=path.read_text(encoding='utf-8')
    pat=r'(let '+variable+r' : Array Name := #\[)(.*?)(\n  \])'
    m=re.search(pat,s,re.S);assert m,filename
    old=[n for n in re.findall(r'``([\w.]+)',m.group(2)) if not n.startswith('ModAB.Rational')]
    allnames=old+newnames
    body='\n'+',\n'.join('    ``'+n for n in allnames)
    s=s[:m.start(2)]+body+s[m.end(2):];path.write_text(s, encoding='utf-8')
    return allnames

new_theorems=[r['name'] for r in rows if r['kind'] in ('theorem','lemma')]
audit_names=replace_names('Audit.lean','declarations',new_theorems)
replace_names('ExportInventory.lean','names',[r['name'] for r in rows])
p=root/'ModAB.lean';s=p.read_text(encoding='utf-8')
s='\n'.join(line for line in s.splitlines() if not line.startswith('import ModAB.Rational')
    and not line.startswith('-- Standalone RModAB paper: rational proposal and uncapped convergence'))
s=s.rstrip()+'\n'
s+='\n-- Standalone RModAB paper: rational proposal and uncapped convergence, 11 November 2026.\n'
s+='\n'.join('import ModAB.'+m for m in modules)+'\n';p.write_text(s, encoding='utf-8')

def closure(mod,seen,ordered):
    if mod in seen:return
    seen.add(mod)
    s=(root/(mod.replace('.','/')+'.lean')).read_text(encoding='utf-8')
    for imp in re.findall(r'^import (ModAB\.[\w.]+)',s,re.M):closure(imp,seen,ordered)
    ordered.append(mod)
ordered=[];seen=set()
for mod in modules:closure('ModAB.'+mod,seen,ordered)
external=set()
for mod in ordered:
    for imp in re.findall(r'^import ([\w.]+)',(root/(mod.replace('.','/')+'.lean')).read_text(encoding='utf-8'),re.M):
        if not imp.startswith('ModAB.'):external.add(imp)
external.add('Lean.Util.CollectAxioms')
header='''/-
Standalone RModAB paper: exact-real proofs, release 0.8.0 (11 November 2026).
Self-contained source for a Lean 4.34.0 / matching Mathlib environment.
Includes inherited exact-real dependencies and all new rational modules.
It does NOT include the old binary64/FloatLib layer; use the full project
or ModAB_Live_Complete.lean for that layer. This is not an IEEE solver proof.
The file is checked locally with the pinned project; no browser execution is implied.
-/\n'''
parts=[header]+['import '+imp+'\n' for imp in sorted(external)]
for i,mod in enumerate(ordered):
    s=(root/(mod.replace('.','/')+'.lean')).read_text(encoding='utf-8')
    s=re.sub(r'^import .*\n','',s,flags=re.M)
    parts.append(f'\n/- Original module: {mod} -/\nsection RModABBundle{i}\n'+s+f'\nend RModABBundle{i}\n')
audit=(root/'Audit.lean').read_text(encoding='utf-8');auditbody=audit[audit.index('open Lean Elab Command'):]
start=auditbody.index('  let declarations : Array Name := #[')
end=auditbody.index('\n  ]',start)+len('\n  ]')
newarr='  let declarations : Array Name := #[\n'+',\n'.join('    ``'+n for n in new_theorems)+'\n  ]'
focusaudit=auditbody[:start]+newarr+auditbody[end:]
parts.append('\n'+focusaudit)
(root/'RModAB_Live.lean').write_text(''.join(parts), encoding='utf-8')

# Preserve all original proof bodies in the complete edition and append the new modules.
full=root/'ModAB_Live_Complete.lean'
s=full.read_text(encoding='utf-8')
if '/- ===== Rational extension modules ===== -/' in s:
    s=s.split('/- ===== Rational extension modules ===== -/')[0]
else:s=s.split('/- ===== Final transitive axiom audit ===== -/')[0]
s=re.sub(r'version 0\.[5678]\.0','version 0.8.0',s)
s=re.sub(r'PASS: \d+ declarations;',f'PASS: {len(audit_names)} declarations;',s)
s=s.replace('  Stable Release (v4.34.0 with mathlib, cslib)', '  a project matching Lean 4.34.0 and the pinned Mathlib revision')
s=s.replace('cubic lemma, exact-arithmetic global convergence and width bound, both switching conventions for the exact AB', 'cubic lemma, exact-arithmetic global convergence and width bound, local C3 rational error,\ntribonacci/golden-ratio Q-orders and pair-cycle consequences; both switching conventions for the exact AB')
s=s.replace('All 511 authored ModAB declarations are retained with their original names.',
    'All original authored ModAB declarations are retained; the rational extension is appended.')
s=s.replace('public import Mathlib.Topology.Sequences\n','')
s=s.replace('@[expose] public section','public import Mathlib.Topology.Sequences\n\n@[expose] public section',1)
existing_imports=set(re.findall(r'^(?:public )?import ([\w.]+)',s,re.M))
needed=sorted(external-existing_imports)
s=s.replace('@[expose] public section',''.join('public import '+imp+'\n' for imp in needed)+'\n@[expose] public section',1)
head,marker,body=s.partition('@[expose] public section')
s=re.sub(r'\n(?:[ \t]*\n)+', '\n\n', head)+marker+body
s=s.rstrip()+'\n'
s+='\n/- ===== Rational extension modules ===== -/\n'
for i,mod in enumerate(modules):
    t=(root/'ModAB'/f'{mod}.lean').read_text(encoding='utf-8');t=re.sub(r'^import .*\n','',t,flags=re.M)
    s+=f'\nsection RationalCompleteBundle{i}\n'+t+f'\nend RationalCompleteBundle{i}\n'
s+='\n/- ===== Final transitive axiom audit ===== -/\n'+auditbody
full.write_text(s, encoding='utf-8')
print(json.dumps({'new_theorems':len(new_theorems),'total_audited':len(audit_names),
    'new_declarations':len(rows),'focus_live_lines':len(''.join(parts).splitlines()),
    'complete_live_lines':len(s.splitlines()),'dependency_modules':ordered},indent=2))
