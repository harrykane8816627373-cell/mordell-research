"""Floating search proposes supports; integer Bernstein covers accept them.

No numerical optimum, residual or sampled positivity is part of acceptance.
"""
from fractions import Fraction as F
from pathlib import Path
from math import ceil
import hashlib
import json

from core import (check_certificate, cycle_geometry, edge_minorant,
                  jensen_reinforcement, old_gamma, pair_box_upper,
                  positive_cover, support_polynomial, validate_geometry)


def horner(poly, x):
    value = 0.0
    for a in reversed(poly):
        value = value*x+a
    return value


def golden_maximum(fun, lo, hi, count=55):
    ratio = (5**0.5-1)/2
    a, b = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
    fa, fb = fun(a), fun(b)
    for _ in range(count):
        if fa > fb:
            hi, b, fb = b, a, fa
            a = hi-ratio*(hi-lo)
            fa = fun(a)
        else:
            lo, a, fa = a, b, fb
            b = lo+ratio*(hi-lo)
            fb = fun(b)
    x = (lo+hi)/2
    return x, fun(x)


def conjugate(poly, coefficient):
    def fun(x): return coefficient*x-horner(poly, x)
    grid = [fun(i/128) for i in range(129)]
    values = [(0.0, grid[0]), (1.0, grid[-1])]
    for i in range(1, 128):
        if grid[i] >= grid[i-1] and grid[i] >= grid[i+1]:
            values.append(golden_maximum(fun, (i-1)/128, (i+1)/128, 42))
    return max(values, key=lambda pair: pair[1])


def propose_cycle(edges):
    polys = [[float(c) for c in edge_minorant(e['g'], e['eta'])] for e in edges]
    def dual(p):
        return p-sum(conjugate(poly, p*edge['a'])[1] for poly, edge in zip(polys, edges))
    pfloat, value = golden_maximum(dual, 0, 128)
    slope = F(round(pfloat*10**6), 10**6)
    supports = []
    for edge, poly in zip(edges, polys):
        approximate = conjugate(poly, float(slope)*edge['a'])[1]
        b = F(ceil(approximate*10**7)+1, 10**7)
        for trial in range(12):
            try:
                cover = positive_cover(support_polynomial(edge['g'], edge['eta'], edge['a'], slope, b))
                break
            except ValueError:
                b += F(1, 10**6)
        else:
            raise ArithmeticError('unable to certify proposed affine support')
        supports.append({'edge': list(edge['edge']), 'intercept': str(b), 'cover': cover})
    return slope, supports


def make_certificate(ms, rs, cycles):
    ms, rs = validate_geometry(ms, rs, 0)
    records = []
    old = jensen = nonlinear = F(0)
    for vertices, weight in cycles:
        edges = cycle_geometry(ms, rs, vertices)
        p, supports = propose_cycle(edges)
        records.append({'vertices': list(vertices), 'weight': str(weight), 'slope': str(p), 'supports': supports})
        old += weight*old_gamma(edges)
        jensen += weight*jensen_reinforcement(edges)
        nonlinear += weight*(p-sum((F(v['intercept']) for v in supports), F(0)))
    U = pair_box_upper(ms, rs, 0)
    cert = {'format': 'mordell-nonlinear-phase-v1', 'date': '2026-11-08',
            'populations': list(ms), 'radius_intervals': [[str(a), str(b)] for a, b in rs],
            'center': 0, 'cycles': records, 'old_packing_gap': str(old), 'jensen_gap': str(jensen),
            'nonlinear_gap': str(nonlinear), 'pair_box_upper': str(U), 'all_phase_upper': str(U/(1+nonlinear)),
            'scope': 'specified regular-ring populations and radius intervals, all continuous phases; not unrestricted M_n'}
    check_certificate(cert)
    return cert


def examples():
    sq = [F(9, 10), F(1), F(11, 10), F(6, 5)]
    box = [(r-F(1, 10000), r+F(1, 10000)) for r in sq]
    triangles = [([0, 1, 2], F(1, 2)), ([0, 1, 3], F(1, 2)),
                 ([0, 2, 3], F(1, 2)), ([1, 2, 3], F(1, 2))]
    mix = [F(4, 5), F(9, 10), F(1), F(11, 10)]
    return {
        'four_squares.json': make_certificate([4]*4, [(r, r) for r in sq], triangles),
        'four_squares_box.json': make_certificate([4]*4, box, triangles),
        'mixed_twentyfive_cycle.json': make_certificate([8, 6, 3, 8], [(r, r) for r in mix], [([0, 1, 2, 3], F(1))]),
    }


def source_paths(base):
    return sorted(p for p in base.rglob('*') if p.is_file() and
                  (p.suffix in ('.py', '.md') or p.name == '.gitattributes' or
                   (p.parent.name == 'certificates' and p.suffix == '.json')))


if __name__ == '__main__':
    base = Path(__file__).resolve().parent
    (base/'certificates').mkdir(exist_ok=True)
    for name, cert in examples().items():
        (base/'certificates'/name).write_text(json.dumps(cert, indent=2)+'\n', encoding='utf-8')
        result = check_certificate(cert)
        print(name, 'old', float(result['old_gap']), 'Jensen', float(result['jensen_gap']),
              'nonlinear', float(result['nonlinear_gap']), 'leaves', result['bernstein_leaves'], flush=True)
    hashes = {p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths(base)}
    (base/'source_manifest.json').write_text(json.dumps({'algorithm': 'sha256', 'files': hashes}, indent=2)+'\n', encoding='utf-8')
