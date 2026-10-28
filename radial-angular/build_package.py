"""Build a deterministic ZIP containing this update and its pinned inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def build(output):
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    members = {}
    for row in manifest['files']:
        data = (ROOT/row['path']).read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('manifest mismatch '+row['path'])
        members[row['path']] = data
    for name in ('source_manifest.json', 'verification_report.json'):
        members['radial-angular/'+name] = (HERE/name).read_bytes()
    members['START_HERE.txt'] = (
        'Mordell radial/angular research, 2026-10-28.\n'
        'From this extracted directory run:\n'
        'python -B -S -O radial-angular/verify_all.py\n'
        'Read radial-angular/README.md and its four proof documents.\n'
        'Only Python standard library is required. No network or installer is used.\n'
        'This is this update with fixed runtime/proof inputs, not the entire historical archive.\n'
        'Repository: https://github.com/harrykane8816627373-cell/mordell-research\n'
    ).encode('utf-8')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            item = zipfile.ZipInfo(name, (2026, 10, 28, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, data)
    data = output.read_bytes()
    return {'path': str(output), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'members': len(members)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'mordell-radial-angular-2026-10-28.zip')
    print(json.dumps(build(parser.parse_args().output), indent=2))
