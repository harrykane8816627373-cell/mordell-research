"""Build the self-contained PARTIAL joint-cover package from pinned bytes."""
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
        if rel.is_absolute() or '..' in rel.parts or ':' in name or '\\' in name or name in members:
            raise ValueError('unsafe or duplicate fixed path')
        data = (ROOT/name).read_bytes()
        if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
            raise ArithmeticError('fixed bytes '+name)
        members[name] = data
    members['joint-closure/source_manifest.json'] = (HERE/'source_manifest.json').read_bytes()
    members['START_HERE.txt'] = (
        'Eight-point Mordell joint system, 2026-11-01: PARTIAL COVER ONLY.\n'
        'Run: python -B -S -O joint-closure/verify_all.py\n'
        'The degree-14 polynomial PSD model is equivalent to true stationarity plus full stability.\n'
        'The stored tree has 60021 nodes, 29990 excluded leaves and 21 UNRESOLVED regions.\n'
        'The complete eight-point sharp bound is NOT proved.\n'
        'A separate exact physical stationary orbit passes four old power modes but is fully unstable.\n'
        'Replay uses the Python standard library; numerical proposal tools optionally use NumPy.\n'
        'Read joint-closure/README.md for the analytic proof, global-only cuts and dependencies.\n'
        'The supplied older unified ZIP and required radial proof inputs are preserved.\n'
        'Other repository navigation links require the complete repository.\n'
    ).encode('utf-8')
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            item = zipfile.ZipInfo(name, (2026, 11, 1, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, data)
    data = output.read_bytes()
    return {'path': str(output), 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'members': len(members), 'full_eight_point_exclusion_proved': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(build(parser.parse_args().output), indent=2))
