"""Exact inner-pair localization using the sharper two-band Ward kernel."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INNER_LOWER = F(9, 200)
INNER_UPPER = F(49, 200)
DOT_LOWER = F(-7, 50)
DOT_UPPER = F(-47, 500)
MIDPOINT_SQUARED_UPPER = F(1, 25)


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def inside(box):
    (u0,u1), (v0,v1), (t0,t1) = box
    return (min(u0,v0)>INNER_LOWER and max(u1,v1)<INNER_UPPER
            and DOT_LOWER<t0<=t1<DOT_UPPER
            and u1+v1+2*t1<4*MIDPOINT_SQUARED_UPPER and 400*t1*t1>121*u1*v1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for base,key in ((HERE, 'files'), (ROOT, 'dependencies')):
        for file,digest in manifest[key].items():
            need(hashlib.sha256((base/file).read_bytes()).hexdigest()==digest, 'fixed '+file)
    banded_dir = ROOT/'ten-point-banded-kernel'
    banded = load('new_banded_kernel', banded_dir/'verify_all.py')
    counted = load('prior_counted_radii', ROOT/'ten-point-counted-radii/verify_all.py')
    row = load('prior_ward_rows', ROOT/'ten-point-ward-radii/verify_all.py')
    upper = banded.upper
    data = json.loads((banded_dir/'kernel.json').read_bytes())
    gram,force = upper.corrections(data)
    base,shape = upper.polynomial(data, gram, force)
    poly = banded.root_polynomial(base, shape, banded.parsed(data), 0, 0)
    pack = upper.Packed(shape, True)
    upper.old.selftest()
    gap = F(1, 10**6)
    need(upper.log.log_lower(upper.witness())>upper.hx(data['B'])-gap, 'certified witness log gap')
    budget = 44*banded.C+gap
    need(banded.C==F(1, 10000) and budget==F(4401, 10**6), 'fixed banded pair slack and new budget')
    other = json.loads((ROOT/'ten-point-counted-radii/verification_report.json').read_bytes())
    need(other['passed'] and F(other['inner_pair']['dot_interval'][1])<0, 'strictly negative inner dot product')
    lo,hi = banded.BANDS[0]
    box = ((lo,hi),(lo,hi),(-hi,F(0)))
    stack = [(box,upper.initial(poly,shape,box),0)]
    counts = dict.fromkeys('AGP012', 0)
    depth = 0
    for node in (HERE/'cover.txt').read_text(encoding='ascii').strip():
        need(bool(stack) and node in counts, 'valid action in nonempty tree')
        counts[node] += 1
        box,coefficients,level = stack.pop()
        depth = max(depth, level)
        if node == 'A':
            need(pack.dominates(coefficients, box, -budget), 'Bernstein exact inner-pair budget separation')
        elif node == 'G':
            need(inside(box), 'all four inner geometry bounds')
        elif node == 'P':
            need(box[0][0]>box[1][1] or counted.prunable(row, upper, box), 'exact symmetric physical prune')
        else:
            axis = int(node)
            a,b = row.split(box, axis)
            ca,cb = pack.split(coefficients, axis)
            stack.extend(((b,cb,level+1),(a,ca,level+1)))
    need(not stack and sum(counts[k] for k in 'AGP')==1+sum(counts[k] for k in '012'), 'complete inner-pair tree')
    result = {'passed': True, 'normalization_mass': 9, 'applies_to_every_global_maximizer': True,
              'inner_squared_radius_open_interval': [str(INNER_LOWER), str(INNER_UPPER)],
              'inner_dot_product_open_interval': [str(DOT_LOWER), str(DOT_UPPER)],
              'inner_midpoint_squared_strict_upper': str(MIDPOINT_SQUARED_UPPER),
              'inner_angle_cosine_strict_upper': '-11/20', 'pair_deficit_strict_upper': str(budget),
              'nodes': sum(counts.values()), 'counts': counts, 'maximum_depth': depth, 'unresolved': 0}
    output = json.dumps(result, indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(output.encode('utf-8'))
    print(output, end='')


if __name__ == '__main__':
    main()
