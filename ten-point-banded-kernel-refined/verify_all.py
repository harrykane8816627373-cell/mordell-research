"""Exact refined band certificate and inner-pair localization."""
import argparse
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
C = F(3, 40000)
C_INNER = F(9, 100000)
RATIO = F(5017, 5000)
PAIR_BUDGET = F(3301, 10**6)


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def load(path):
    spec = importlib.util.spec_from_file_location(path.replace('/', '_'), ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def geometry(box):
    (u0,u1), (v0,v1), (t0,t1) = box
    return (min(u0,v0)>F(27,500) and max(u1,v1)<F(217,1000)
            and F(-27,200)<t0<=t1<F(-1,10)
            and u1+v1+2*t1<F(1,10) and 9*t1*t1>4*u1*v1
            and u1+v1-2*t0<F(3,5))


def check_pair_tree(banded, data, tree):
    upper = banded.upper
    row = load('ten-point-ward-radii/verify_all.py')
    counted = load('ten-point-counted-radii/verify_all.py')
    previous = json.loads((ROOT/'ten-point-inner-pair/verification_report.json').read_bytes())
    need(previous['passed'] and previous['applies_to_every_global_maximizer']
         and F(previous['inner_dot_product_open_interval'][1])<0, 'inherited negative inner dot product')
    A, force = upper.corrections(data)
    poly, shape = upper.polynomial(data, A, force)
    poly = banded.root_polynomial(poly, shape, banded.parsed(data), 0, 0)
    pack = upper.Packed(shape,True)
    band = banded.BANDS[0]
    box = (band,band,(-band[1],F(0)))
    stack = [(box,upper.initial(poly,shape,box),0)]
    totals = dict.fromkeys('AGP012',0)
    depth_max = 0
    for action in tree:
        need(stack and action in totals,'valid action in inner-pair tree')
        box,poly,depth = stack.pop()
        depth_max = max(depth_max,depth)
        totals[action] += 1
        if action=='A':
            need(pack.dominates(poly,box,-PAIR_BUDGET),'exact deficit above the inner-pair budget')
        elif action=='G':
            need(geometry(box),'strict inner-pair geometry bounds')
        elif action=='P':
            need(box[0][0]>box[1][1] or counted.prunable(row,upper,box),'valid necessary-domain exclusion')
        else:
            stack.extend(banded.split(pack,poly,box,depth,int(action)))
    need(not stack and sum(totals[k] for k in 'AGP')==1+sum(totals[k] for k in '012'),
         'complete exact inner-pair tree')
    return {'visited':sum(totals.values()),'unresolved':0,'counts':totals,'maximum_depth':depth_max,
            'applies_to_every_global_maximizer':True,
            'inner_squared_radius_open_interval':['27/500','217/1000'],
            'inner_dot_product_open_interval':['-27/200','-1/10'],
            'inner_midpoint_squared_strict_upper':'1/40','angle_cosine_strict_upper':'-2/3',
            'inner_squared_distance_strict_upper':'3/5',
            'pair_deficit_strict_upper':str(PAIR_BUDGET)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    parser.add_argument('--progress',action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for directory,section in ((HERE,'files'),(ROOT,'dependencies')):
        for name,digest in manifest[section].items():
            need(sha256((directory/name).read_bytes()).hexdigest()==digest,'fixed '+name)
    banded = load('ten-point-banded-kernel/verify_all.py')
    upper = banded.upper
    banded.verify_dependencies()
    data = json.loads((HERE/'kernel.json').read_bytes())
    tree = json.loads((HERE/'cover.json').read_bytes())
    previous = json.loads((ROOT/'ten-point-banded-kernel/kernel.json').read_bytes())
    need(all(data[k]==previous[k] for k in ('radial_nodes','radial_moment_centers','transforms','Q_blocks','B')),
         'same exact basis, positive Gram blocks and total budget')
    need(len(banded.parsed(data))==19,'exact piecewise directions')
    need(F(tree['c'])==C and F(tree['c_inner'])==C_INNER,'fixed unequal slacks')
    witness = upper.witness()
    B = upper.hx(data['B'])
    need(upper.actual_root_lower()>witness,'accepted candidate exceeds the feasible witness')
    log_upper = B+44*C+C_INNER
    need(upper.log.log_lower(witness*RATIO)>log_upper+F(1,250000),'strict logarithmic ratio')
    need(B+44*C-upper.log.log_lower(witness)<PAIR_BUDGET,'strict single inner-pair budget')
    with redirect_stdout(sys.stderr if args.progress else StringIO()):
        report = banded.verify(data,tree,C,C_INNER)
    report['global_over_witness_strict_upper'] = str(RATIO)
    report['exact_inner_count'] = 2
    report['exact_outer_count'] = 8
    report['sharp_global_optimum_proved'] = False
    report['inner_pair'] = check_pair_tree(banded,data,(HERE/'inner-pair.txt').read_text(encoding='ascii').strip())
    output = json.dumps(report,indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_bytes(output.encode('utf-8'))
    print(output,end='')


if __name__=='__main__':
    main()
