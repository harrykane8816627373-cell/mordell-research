"""Build the pinned-input offline replay ZIP without changing old packages."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def build(output):
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    members = {row['path']: (ROOT/row['path']).read_bytes() for row in manifest['files']}
    for row in manifest['files']:
        data = members[row['path']]
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('manifest mismatch '+row['path'])
    for name in ('source_manifest.json', 'verification_report.json'):
        members['joint-bounds/'+name] = (HERE/name).read_bytes()
    members['START_HERE.txt'] = (
        'Mordell joint bounds, 2026-10-28.\n'
        'From this extracted directory run:\n'
        'python -B -S -O joint-bounds/verify_all.py\n'
        'Read joint-bounds/README.md and its three proof documents.\n'
        'Only Python standard library is required; no network is used.\n'
        'This contains this update and its fixed runtime inputs, not the full historical archive.\n'
        'Historical context: https://github.com/harrykane8816627373-cell/mordell-research\n'
    ).encode('utf-8')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 28, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    data = output.read_bytes()
    return {'path': str(output), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'members': len(members)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'mordell-joint-bounds-2026-10-28.zip')
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2))
