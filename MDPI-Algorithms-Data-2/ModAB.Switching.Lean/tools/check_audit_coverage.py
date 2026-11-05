#!/usr/bin/env python3
"""Check the explicit audit list, compiled inventory, map links and source hashes.
Run lake build and lake env lean Audit.lean first to re-check kernel proofs.
This script checks packaging/coverage consistency; it is not a theorem prover.
"""
from pathlib import Path
from html.parser import HTMLParser
import hashlib,json,re
p=Path(__file__).resolve().parents[1]
inv=json.loads((p/'docs/proof-inventory.json').read_text())
manifest=json.loads((p/'docs/proof-manifest.json').read_text())
audit=re.findall(r'``(ModAB\.[A-Za-z0-9_.]+)',(p/'Audit.lean').read_text())
expected={v['name'] for v in inv if v['kind'] in {'theorem','lemma'}}
assert len(audit)==len(set(audit)), 'Duplicate audit entries'
assert set(audit)==expected, {'missing':sorted(expected-set(audit)), 'extra':sorted(set(audit)-expected)}
assert len(audit)==manifest['audited_theorems']==380
assert len(inv)==len({v['name'] for v in inv})==manifest['authored_declarations']==511
allowed={'propext','Classical.choice','Quot.sound'}
for row in inv:
 assert set(row['axioms'])<=allowed,(row['name'],row['axioms'])
 assert set(row['dependencies'])<={v['name'] for v in inv},row['name']
for file in (p/'ModAB').glob('*.lean'):
 # Project proof sources have no escapes/admitted results or native proof oracles.
 text=file.read_text()
 text=re.sub(r'/\-[\s\S]*?\-/','',text)
 text=re.sub(r'--[^\n]*','',text)
 assert not re.search(r'\b(sorry|admit|axiom|native_decide)\b',text),file
class MapParser(HTMLParser):
 def __init__(self):super().__init__();self.ids=[];self.links=[]
 def handle_starttag(self,tag,attrs):
  attrs=dict(attrs)
  if 'id' in attrs:self.ids.append(attrs['id'])
  if tag=='a' and attrs.get('href','').startswith('#'):self.links.append(attrs['href'][1:])
m=MapParser();m.feed((p/'ModAB_Proof_Map.html').read_text())
assert len(m.ids)==len(set(m.ids)), 'Duplicate map IDs'
assert set(m.links)<=set(m.ids),'Broken map anchors'
assert {x[5:] for x in m.ids if x.startswith('decl-')}=={v['name'] for v in inv}
assert len([x for x in m.ids if x.startswith('step-')])==91
assert len([x for x in m.ids if x.startswith('paper-')])==23
for source,expected_hash in manifest['files'].items():
 prefix='MDPI-Algorithms-Data-2/ModAB.Switching.Lean/'
 file=p/source[len(prefix):] if source.startswith(prefix) else p.parents[1]/source
 # A standalone proof archive need not include the whole manuscript.
 if not file.exists() and source.startswith('manuscript/'):continue
 assert file.is_file(),file
 assert hashlib.sha256(file.read_bytes()).hexdigest()==expected_hash,f'Source changed: {source}'
print('PASS: 380 audited statements, 511 mapped declarations, 91 proof steps, 23 manuscript entries; source hashes and anchors agree.')
