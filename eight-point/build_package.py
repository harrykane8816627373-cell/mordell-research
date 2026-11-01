"""Build an offline ZIP from the already checked, fixed input manifest."""
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
        rel = PurePosixPath(name)
        if rel.is_absolute() or '..' in rel.parts or '\\' in name or ':' in name or name != rel.as_posix():
            raise ValueError('unsafe fixed path')
        if name in members:
            raise ValueError('duplicate fixed path')
        data = (ROOT/name).read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('fixed bytes '+name)
        members[name] = data
    members['eight-point/source_manifest.json'] = (HERE/'source_manifest.json').read_bytes()
    members['START_HERE.txt'] = (
        'Mordell eight-point research, 2026-11-01.\n'
        'From this extracted directory run:\n'
        'python -B -S -O eight-point/verify_all.py\n'
        'Read eight-point/MOMENT_TUBE.md for the analytic proof.\n'
        'The p7-free tube is certified; the complete eight-point sharp bound remains open.\n'
        'Only the current proof and its required legacy dependencies are included.\n'
    ).encode('utf-8')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            item = zipfile.ZipInfo(name, (2026, 11, 1, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, data)
    data = output.read_bytes()
    return {'path': str(output), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest(), 'members': len(members)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'mordell-eight-point-2026-11-01.zip')
    print(json.dumps(build(parser.parse_args().output), indent=2))
