"""Build a deterministic, self-contained algebraic-gap ZIP with pinned sources."""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import zipfile

from verify_all import read_sources

HERE = Path(__file__).resolve().parent


def build(output, source_archive):
    config = json.loads((HERE/'certificate.json').read_text(encoding='utf-8'))
    manifest = json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    members = {}
    for row in manifest['files']:
        name = row['path']
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
            raise ValueError('unsafe fixed path')
        data = (HERE/name).read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('fixed current input '+name)
        members['algebraic-gap/'+name] = data
    members['algebraic-gap/source_manifest.json'] = (HERE/'source_manifest.json').read_bytes()
    for name, data in read_sources(config, archive=source_archive).items():
        members['inherited-first3/'+name] = data
    members['START_HERE.txt'] = (
        'Mordell eight-point algebraic moment gap, 2026-11-01.\n'
        'python -B -S -O algebraic-gap/verify_all.py\n'
        'This rebuilds the inherited resultants and checks all 14179 integer coefficients.\n'
        'Then it recomputes seven multiplier norms and the quantitative moment bound.\n'
        'Read algebraic-gap/README.md for all hypotheses and the analytic proof.\n'
        'The certified lower bound is |p2| >= 10^-36 |p3| for stationary challengers in the old ellipse.\n'
        'The complete eight-point sharp bound remains unproved.\n'
        'Only this proof and three required inherited source files are included.\n'
        'Links to other repository topics require the full repository.\n'
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
    parser.add_argument('--source-archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.source_archive), indent=2))
