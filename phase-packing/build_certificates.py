"""Regenerate the finite examples; verify_all.py accepts fixed input files."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

from packing import build_certificate, check_certificate


def examples():
    square_radii = [F(9, 10), F(1), F(11, 10), F(6, 5)]
    mixed_radii = [F(4, 5), F(9, 10), F(1), F(11, 10)]
    width = F(1, 10000)
    return {
        'four_squares_triangles.json': build_certificate([4]*4, [(r, r) for r in square_radii]),
        'four_squares_radius_box.json': build_certificate([4]*4, [(r-width, r+width) for r in square_radii]),
        'mixed_twentyfive_triangles.json': build_certificate([8, 6, 3, 8], [(r, r) for r in mixed_radii], 3),
        'mixed_twentyfive_cycles.json': build_certificate([8, 6, 3, 8], [(r, r) for r in mixed_radii], 4),
        'compatible_fourteen.json': build_certificate([2, 4, 8], [(F(1), F(1)), (F(2), F(2)), (F(3), F(3))]),
    }


def input_paths(base):
    return sorted(p for p in base.rglob('*') if p.is_file()
                  and (p.suffix in ('.py', '.md') or p.name == '.gitattributes'
                       or (p.parent.name == 'certificates' and p.suffix == '.json')))


if __name__ == '__main__':
    base = Path(__file__).resolve().parent
    (base/'certificates').mkdir(exist_ok=True)
    for name, cert in examples().items():
        check_certificate(cert)
        (base/'certificates'/name).write_text(json.dumps(cert, indent=2)+'\n', encoding='utf-8')
    hashes = {p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths(base)}
    (base/'source_manifest.json').write_text(json.dumps({'algorithm': 'sha256', 'files': hashes}, indent=2)+'\n', encoding='utf-8')
    print('Generated five independently checkable LP certificates and source manifest.')
