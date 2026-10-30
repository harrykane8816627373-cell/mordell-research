"""Deterministic offline overlap-budget package with pinned inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def build(output):
    manifest=json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    members={}
    for row in manifest['files']:
        rel=Path(row['path'])
        if rel.is_absolute() or rel.drive or '..' in rel.parts:
            raise ValueError('unsafe manifest path')
        data=(ROOT/rel).read_bytes()
        if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise ArithmeticError('manifest mismatch '+row['path'])
        members[row['path']]=data
    for name in ('source_manifest.json','verification_report.json'):
        members['overlap-budget/'+name]=(HERE/name).read_bytes()
    members['START_HERE.txt']=(
        'Mordell configurations: parabolic overlap, vertex covers and joint budgets, 2026-10-30.\n'
        'From this extracted directory run:\n'
        'python -B -S -O overlap-budget/verify_all.py\n'
        'Read overlap-budget/README.md and its proof documents.\n'
        'Acceptance uses only Python standard library; no optimizer or network is needed.\n'
        'Four new and eight prior entries use exact arithmetic and outward intervals.\n'
        'This package pins this update and needed prior inputs, not the complete historical archive.\n'
        'https://github.com/harrykane8816627373-cell/mordell-research\n').encode('utf-8')
    output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name,data in sorted(members.items()):
            item=zipfile.ZipInfo(name,(2026,10,30,0,0,0))
            item.compress_type=zipfile.ZIP_DEFLATED
            item.external_attr=0o100644<<16
            archive.writestr(item,data)
    data=output.read_bytes()
    return {'path':str(output),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'members':len(members)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'mordell-overlap-budget-2026-10-30.zip')
    print(json.dumps(build(parser.parse_args().output),indent=2))
