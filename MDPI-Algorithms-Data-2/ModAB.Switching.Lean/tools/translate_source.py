#!/usr/bin/env python3
"""Fail-closed translation of this project's finite C# arithmetic body to a Lean expression.
This small front end is reviewable infrastructure, not a verified C# compiler.
"""
from pathlib import Path
import argparse, hashlib, json, re
root=Path(__file__).resolve().parents[1]
source=root/'csharp/src/ModAB.Switching/SwitchingCriterion.cs'
primitives=root/'csharp/src/ModAB.Switching/Binary64.cs'
target=root/'ModAB/Source.lean'
FUNCTIONS={'Math.Abs':('abs',1),'Math.Min':('min',2),'Math.Max':('max',2),
 'Binary64.Add':('op .add',2),'Binary64.Subtract':('op .sub',2),
 'Binary64.Multiply':('op .mul',2),'Binary64.Divide':('op .div',2)}

def split_top(s,sep=','):
 depth=0; start=0; out=[]
 for i,c in enumerate(s):
  if c=='(': depth+=1
  if c==')': depth-=1
  if depth<0: raise ValueError('Unbalanced parentheses')
  if c==sep and depth==0: out.append(s[start:i]);start=i+1
 if depth: raise ValueError('Unbalanced parentheses')
 return out+[s[start:]]

def parse(s,known):
 s=s.strip()
 if s in ('f1','f2','f3'): return f'(.input {int(s[1])-1})'
 if s in ('1.0','2.0'): return '.one' if s=='1.0' else '.two'
 if s in known: return s+'Expr'
 m=re.fullmatch(r'([A-Za-z0-9_.]+)\((.*)\)',s)
 if not m or m[1] not in FUNCTIONS: raise ValueError('Unsupported expression: '+s)
 name,arity=FUNCTIONS[m[1]]; args=split_top(m[2])
 if len(args)!=arity: raise ValueError('Wrong arity')
 return f'(.{name} '+ ' '.join(parse(a,known) for a in args)+')'

def produce():
 text=source.read_text()
 m=re.search(r'        double e1 =.*?(?=        SwitchingDecision decision =)',text,re.S)
 if not m: raise ValueError('Arithmetic block not found')
 body=m[0]; rows=[]; known=set()
 for statement in body.strip().split(';'):
  if not statement.strip(): continue
  if not statement.strip().startswith('double '): raise ValueError('Unsupported statement')
  for binding in split_top(statement.strip()[7:]):
   name,expr=binding.split('=',1);name=name.strip()
   if not re.fullmatch(r'[a-z][A-Za-z0-9]*',name) or name in known: raise ValueError('Invalid local')
   rows.append(f'def {name}Expr : Expr := {parse(expr,known)}')
   known.add(name)
 if len(rows)!=18 or 'gap' not in known: raise ValueError('Unexpected arithmetic body shape')
 # Check the wrapper operators and the guard literal exactly, not by floating conversion.
 primitive_text=primitives.read_text()
 for name,op in [('Add','+'),('Subtract','-'),('Multiply','*'),('Divide','/')]:
  expected=f'[MethodImpl(MethodImplOptions.NoInlining)]\n    public static double {name}(double x, double y) => x {op} y;'
  if expected not in primitive_text: raise ValueError('Primitive wrapper changed: '+name)
 guard='3.552713678800500929355621337890625e-15'
 if f'public const double Guard = {guard};' not in text: raise ValueError('Guard changed')
 tail='''SwitchingDecision decision = gap > Guard
            ? SwitchingDecision.CertifiedSwitch
            : gap < -Guard
                ? SwitchingDecision.CertifiedNoSwitch
                : SwitchingDecision.Indeterminate;'''
 if tail not in text: raise ValueError('Decision tree changed')
 hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,primitives]}
 out='import ModAB.Expressions\n\n/- Arithmetic expressions extracted from the accompanying C# source.\nRegenerate with tools/translate_source.py; check with --check. -/\nnamespace ModAB.Source\nopen IEEE64\n\n'+'\n'.join(rows)+'''\n
noncomputable def decision (gap : ℝ) : Nat :=
  if 1/(2:ℝ)^48 < gap then 3 else if gap < -(1/(2:ℝ)^48) then 2 else 1

end ModAB.Source
'''
 return out,hashes

parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
out,hashes=produce();manifest=root/'source-manifest.json'
traced=root/'csharp/tests/ModAB.Switching.Verification/TracedCriterion.cs'
original=source.read_text()
traced_text='namespace ModAB.Switching.Verification;\n\n'+original[original.index('public static class SwitchingCriterion'):]
traced_text=traced_text.replace('class SwitchingCriterion','class TracedCriterion').replace('Binary64.','RoundingCheck.')
if args.check:
 if traced.read_text()!=traced_text: raise SystemExit('Instrumented C# differs; regenerate.')
 if target.read_text()!=out: raise SystemExit('Generated Lean source differs; regenerate and rebuild.')
 if json.loads(manifest.read_text())!=hashes: raise SystemExit('C# source hashes differ; regenerate and rebuild.')
 print('PASS: C# arithmetic AST, primitive wrappers, decision tree and source hashes.')
else:
 traced.write_text(traced_text)
 target.write_text(out);manifest.write_text(json.dumps(hashes,indent=2)+'\n')
 print('Wrote ModAB/Source.lean and source-manifest.json')
