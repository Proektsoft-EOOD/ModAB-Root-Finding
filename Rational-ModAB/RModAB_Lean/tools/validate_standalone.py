#!/usr/bin/env python3
"""Check release integrity and local links; this does not run the Lean kernel."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import csv
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / 'docs/standalone-source-manifest.json').read_text(encoding='utf-8'))
failures = []

def check(condition, message):
    if not condition:
        failures.append(message)

for item in manifest['preserved_files']:
    path = ROOT / item['path']
    check(path.is_file(), 'Missing preserved file: ' + item['path'])
    if path.is_file():
        check(hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'],
              'Changed proof source or pin: ' + item['path'])

required = ['README.md', 'README_BG.txt', 'IMPLEMENTATION_SCOPE.md',
            'RModAB_Proof_Map.html', 'ModAB_Proof_Map.html',
            'docs/Proof_Correspondence.csv', 'Audit.lean',
            'RModAB_Live.lean', 'ModAB_Live_Complete.lean',
            'validation/standalone/README.md']
for name in required:
    check((ROOT / name).is_file(), 'Missing required file: ' + name)

class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.links = []
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key == 'id':
                self.ids.add(value)
            if key in ('href', 'src'):
                self.links.append(value)

def check_links(path, links, ids=None):
    count = 0
    for link in links:
        url = urlsplit(link)
        if url.scheme or url.netloc:
            continue
        target = path.parent / unquote(url.path) if url.path else path
        check(target.is_file(), f'Broken local link in {path.name}: {link}')
        if not url.path and url.fragment and ids is not None:
            check(unquote(url.fragment) in ids, f'Unknown map anchor: {link}')
        count += 1
    return count

primary = ROOT / 'RModAB_Proof_Map.html'
link_count = 0
if primary.is_file():
    page = Page()
    page.feed(primary.read_text(encoding='utf-8'))
    link_count += check_links(primary, page.links, page.ids)
    alias = ROOT / 'ModAB_Proof_Map.html'
    check(alias.is_file() and alias.read_bytes() == primary.read_bytes(),
          'Compatibility map is not identical to the primary map')
readme = ROOT / 'README.md'
if readme.is_file():
    link_count += check_links(readme, re.findall(r'\[[^\]]+\]\(([^)]+)\)', readme.read_text(encoding='utf-8')))

for path in (ROOT / 'ModAB').glob('*.lean'):
    for name in re.findall(r'^import\s+(ModAB(?:\.[\w]+)+)', path.read_text(encoding='utf-8'), re.M):
        check((ROOT / (name.replace('.', '/') + '.lean')).is_file(),
              f'Missing authored import {name} from {path.name}')

correspondence = ROOT / 'docs/Proof_Correspondence.csv'
if correspondence.is_file():
    with correspondence.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    check(len(rows) == 12, f'Expected 12 paper correspondence rows, got {len(rows)}')

report = {
    'status': 'PASS' if not failures else 'FAIL',
    'kind': 'Structural/source-integrity check, not Lean kernel compilation',
    'preserved_files_checked': len(manifest['preserved_files']),
    'local_links_checked': link_count,
    'failures': failures,
}
out = ROOT / 'validation/standalone/packaging_check.json'
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, indent=2))
raise SystemExit(bool(failures))
