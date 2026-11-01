"""Build a deterministic offline package including the unchanged supplied ZIP."""
from pathlib import Path, PurePosixPath
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
        name = row['path']
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
            raise ValueError('unsafe fixed path')
        if name in members:
            raise ValueError('duplicate fixed path')
        data = (ROOT/name).read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('fixed bytes '+name)
        members[name] = data
    members['radial-gap/source_manifest.json'] = (HERE/'source_manifest.json').read_bytes()
    members['START_HERE.txt'] = (
        'Mordell eight-point radial/angular exclusion, 2026-11-01.\n'
        'Run: python -B -S -O radial-gap/verify_all.py\n'
        'Only the Python standard library is needed.\n'
        'The new theorem gives log(Delta_q/7^7)<=-E/32 for E<=1/1200.\n'
        'Read radial-gap/README.md for the analytic proof and the exact seven-point dependency.\n'
        'The unchanged supplied unified ZIP is preserved under radial-gap/inherited/.\n'
        'Default replay: 14 supplied entries, 2 inherited objective entries, 1 new entry.\n'
        'Optional SymPy root-pairing tests are not run.\n'
        'The full eight-point sharp bound is still unproved.\n'
        'Links to other repository topics require the complete repository.\n'
    ).encode('utf-8')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            info = zipfile.ZipInfo(name, (2026, 11, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    data = output.read_bytes()
    return {'path': str(output), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest(), 'members': len(members)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(build(parser.parse_args().output), indent=2))
