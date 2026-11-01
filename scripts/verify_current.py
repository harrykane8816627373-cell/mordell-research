"""Check maintained navigation; replay the two current suites without dirtying Git."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent.parent


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def main():
    pages = [ROOT/'README.md', ROOT/'configuration-bounds/PROOF_STATUS.md',
             *sorted((ROOT/'docs').glob('*.md')), *sorted((ROOT/'eight-point').glob('*.md'))]
    links = 0
    for page in pages:
        text = page.read_text(encoding='utf-8')
        for target in re.findall(r'\[[^\]\n]+\]\(([^)\n]+)\)', text):
            parsed = urlsplit(target.strip('<>'))
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            resolved = (page.parent/unquote(parsed.path)).resolve()
            require(resolved.is_relative_to(ROOT) and resolved.exists(),
                    'missing or escaping local link: '+str(page.relative_to(ROOT))+' -> '+target)
            links += 1
    print('PASS maintained navigation: '+str(len(pages))+' pages, '+str(links)+' local links', flush=True)
    inputs = {}
    suites = ('power-moment', 'eight-point')
    manifests = {}
    for suite in suites:
        manifest_path = ROOT/suite/'source_manifest.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        manifests[suite+'/source_manifest.json'] = manifest_path.read_bytes()
        for row in manifest['files']:
            name = row['path']
            resolved = (ROOT/name).resolve()
            require(resolved.is_relative_to(ROOT), 'manifest input in repository')
            data = resolved.read_bytes()
            require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                    'fixed bytes '+name)
            require(name not in inputs or inputs[name] == data, 'consistent shared dependency')
            inputs[name] = data
    inputs.update(manifests)
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-current-', dir=ROOT) as temp:
        stage = Path(temp)
        for name, data in inputs.items():
            dest = stage/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
        for suite in suites:
            proc = subprocess.run([sys.executable, '-X', 'int_max_str_digits=0', '-B', '-S', '-O',
                str(stage/suite/'verify_all.py')], cwd=stage, env=env,
                capture_output=True, text=True, encoding='utf-8', timeout=180)
            require(proc.returncode == 0, suite+': '+proc.stdout[-2000:]+proc.stderr[-3000:])
            print('PASS '+suite+'/verify_all.py', flush=True)
    print('PASS current certificates. Complete eight-point exclusion remains unproved.')


if __name__ == '__main__':
    main()
