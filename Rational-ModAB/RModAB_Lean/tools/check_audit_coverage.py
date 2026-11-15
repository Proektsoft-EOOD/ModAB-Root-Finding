#!/usr/bin/env python3
"""Check the standalone paper map, preserved aliases, audit and source hashes.
Run lake build and lake env lean Audit.lean first to re-check kernel proofs.
This script checks packaging/coverage consistency; it is not a theorem prover.
"""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import csv,hashlib,json,re
p=Path(__file__).resolve().parents[1]
inv=json.loads((p/'docs/proof-inventory.json').read_text(encoding='utf-8'))
manifest=json.loads((p/'docs/proof-manifest.json').read_text(encoding='utf-8'))
audit=re.findall(r'``(ModAB\.[A-Za-z0-9_.]+)',(p/'Audit.lean').read_text(encoding='utf-8'))
expected={v['name'] for v in inv if v['kind'] in {'theorem','lemma'}}
assert len(audit)==len(set(audit)), 'Duplicate audit entries'
assert set(audit)==expected, {'missing':sorted(expected-set(audit)), 'extra':sorted(set(audit)-expected)}
assert len(audit)==manifest['audited_theorems']
assert len(inv)==len({v['name'] for v in inv})==manifest['authored_declarations']
new=[v for v in inv if v['name'].startswith('ModAB.Rational')]
assert len(new)==manifest['new_declarations']
assert sum(v['kind'] in {'theorem','lemma'} for v in new)==manifest['new_theorems']
allowed={'propext','Classical.choice','Quot.sound'}
for row in inv:
 assert set(row['axioms'])<=allowed,(row['name'],row['axioms'])
 assert set(row['dependencies'])<={v['name'] for v in inv},row['name']
for file in (p/'ModAB').glob('*.lean'):
 # Project proof sources have no escapes/admitted results or native proof oracles.
 text=file.read_text(encoding='utf-8')
 text=re.sub(r'/\-[\s\S]*?\-/','',text)
 text=re.sub(r'--[^\n]*','',text)
 assert not re.search(r'\b(sorry|admit|axiom|native_decide)\b',text),file
class MapParser(HTMLParser):
 def __init__(self):super().__init__();self.ids=[];self.links=[];self.files=[]
 def handle_starttag(self,tag,attrs):
  attrs=dict(attrs)
  if 'id' in attrs:self.ids.append(attrs['id'])
  if tag=='a':
   href=attrs.get('href','')
   if href.startswith('#'):self.links.append(href[1:])
   elif href and not urlsplit(href).scheme:self.files.append(href)
map_file=manifest.get('proof_map_file','RModAB_Proof_Map.html')
map_text=(p/map_file).read_text(encoding='utf-8')
m=MapParser();m.feed(map_text)
assert (p/'ModAB_Proof_Map.html').read_bytes()==(p/map_file).read_bytes(),'Compatibility map differs'
assert len(m.ids)==len(set(m.ids)), 'Duplicate map IDs'
assert set(m.links)<=set(m.ids),'Broken map anchors'
assert {x[5:] for x in m.ids if x.startswith('decl-')}=={v['name'] for v in inv}
assert len([x for x in m.ids if x.startswith('step-')])==manifest['legacy_proof_step_records']
assert len([x for x in m.ids if x.startswith('paper-')])==manifest['legacy_manuscript_entries']
assert len([x for x in m.ids if x.startswith('rational-paper-')])==manifest['rational_manuscript_entries']
groups=json.loads((p/'docs/rational-map.json').read_text(encoding='utf-8'))
expected_kinds=['Proposition','Proposition','Theorem','Theorem','Theorem','Proposition','Proposition','Lemma','Theorem','Proposition','Corollary','Proposition']
assert len(groups)==manifest['standalone_manuscript_entries']==12
assert [g['number'] for g in groups]==[str(n) for n in range(1,13)]
assert [g['kind'] for g in groups]==expected_kinds
assert [g['aliases'] for g in groups]==[[f'R{n}'] for n in range(1,7)]+[[f'L{n}'] for n in range(1,7)]
assert manifest['standalone_results']==[g['kind']+' '+g['number'] for g in groups]
assert manifest['preserved_aliases']=={g['aliases'][0]:g['kind']+' '+g['number'] for g in groups}
assert {x for x in m.ids if x.startswith('standalone-paper-')}=={f'standalone-paper-{n}' for n in range(1,13)}
assert {f'rational-paper-R{n}' for n in range(1,7)}<=set(m.ids)
legacy=MapParser();legacy.feed((p/'docs/legacy-proof-map.html').read_text(encoding='utf-8'))
assert set(legacy.ids)<=set(m.ids),'Original ModAB map anchors were removed'
previous=json.loads((p/'docs/previous-release-map-anchors.json').read_text(encoding='utf-8'))
assert set(previous['ids'])<=set(m.ids),'Previous RModAB map anchors were removed'
assert len(audit)>667 and len(inv)>850,'Completed local proof extension missing'
assert manifest['local_correspondence_complete'] is True
assert manifest['preserved_previous_proof_modules']==46
assert manifest['project_version']=='0.8.0','Unexpected proof release'
assert '2^floor(n/7)' in groups[2]['proof'],'Capped bound is missing the floor'
for g in groups:
 assert set(g['names'])<={v['name'] for v in inv}
 for name in g['names']:
  row=next(v for v in inv if v['name']==name)
  assert g['label'] in row['paper_labels'],(name,g['label'])
for href in m.files:
 parts=urlsplit(href)
 target=p/unquote(parts.path)
 assert target.resolve().is_relative_to(p.resolve()),f'Link leaves the independent project: {href}'
 assert target.is_file(),f'Missing local link target: {href}'
 if target.suffix=='.lean' and parts.fragment:
  assert re.fullmatch(r'L[1-9][0-9]*',parts.fragment),href
  assert int(parts.fragment[1:])<=len(target.read_text(encoding='utf-8').splitlines()),href
assert manifest['packaging_role'].startswith('Independent repository')
assert manifest['repository_url'] is None and manifest['project_url'] is None
assert manifest['proof_map_url'] is None,'A future repository URL has not been assigned'
assert manifest['predecessor_proof_url'] in map_text
assert manifest['inspected_algorithm_url'] in map_text
with (p/'docs/Proof_Correspondence.csv').open(encoding='utf-8',newline='') as f:
 table=list(csv.DictReader(f))
assert [r['Paper result'] for r in table]==manifest['standalone_results']
for source,expected_hash in manifest['files'].items():
 file=p/source
 assert not Path(source).is_absolute() and file.resolve().is_relative_to(p.resolve()),source
 assert file.is_file(),file
 assert hashlib.sha256(file.read_bytes()).hexdigest()==expected_hash,f'Source changed: {source}'
print(f"PASS: {len(audit)} audited statements, {len(inv)} mapped declarations, {manifest['legacy_proof_step_records']} inherited proof steps, {manifest['legacy_manuscript_entries']} original ModAB article entries and {manifest['standalone_manuscript_entries']} standalone RModAB paper entries (1–12; aliases R1–R6 preserved, L1–L6 added); source hashes, local links, CSV correspondence and all historical/current anchors agree. Independent packaging of proof release 0.8.0.")
