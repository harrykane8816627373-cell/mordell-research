"""Certify the radial Hessian bound throughout the mass-nine half-turn domain."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
TICKS = 2048
SCALE = 65536


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def prerequisites(repo):
    manifest = read(HERE / 'source_manifest.json')
    for base, key in ((HERE, 'files'), (repo, 'dependencies')):
        for name, digest in manifest[key].items():
            require(sha256((base / name).read_bytes()).hexdigest() == digest, 'fixed input ' + name)

    def proof(folder):
        row = read(repo / folder / 'verification_report.json')
        require(row['passed'] and row['normalization_mass'] == 9
                and row['applies_to_every_global_maximizer'] and row['unresolved'] == 0,
                'accepted global bound ' + folder)
        return row

    inner = proof('ten-point-inner-edge')
    outer = proof('ten-point-angular-radii')
    gap = proof('ten-point-outer-gap')
    a, b = (F(d)/4 for d in inner['inner_squared_distance_open_interval'])
    l, u = map(F, outer['outer_squared_radius_open_interval'])
    gl, gu = map(F, gap['outer_successive_angle_cosine_open_interval'])
    cl, cu = 2*gl**2-1, 2*gu**2-1
    require((a, b, l, u) == (F(43,400), F(53,400), F(89,100), F(1271,1000))
            and (outer['exact_inner_count'],outer['exact_outer_count']) == (2,8),
            'squared-radius bands and populations')
    require((cl, cu) == (F(-1519,5000), F(271,1250)) and max(2*cl**2-1, 2*cu**2-1) < F(-81,100),
            'adjacent and opposite doubled-angle cosines')

    # The concave sum of logs takes its minimum over the mass polytope at
    # a vertex. Monotonicity reduces this to the lowest allowed outer mass.
    cut = F(9,2)-b
    products = []
    for triple in product((l,u), repeat=3):
        last = cut - sum(triple)
        if l <= last <= u:
            products.append(math.prod(triple)*last)
    require(products and min(products) > F(43,40)**4, 'log-convex mass-compatible domain')

    # The sign of the c derivative of b(d,c) is the sign of 2-d*d-c*d.
    # For negative c and the inner-outer pairs we use c >= -1.
    outer_d = (l/u+u/l)/2
    inner_d = (b/l+l/b)/2
    require(2-outer_d**2-cu*outer_d > 0 and inner_d > 3 and 2-inner_d**2+inner_d < 0,
            'monotonicity of the edge curvatures in the cosine')
    return min(products)


def atan_interval(q):
    x = F(1,q)
    series = sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(18))
    return series, series+x**37/37


def cosine_sum(x, order):
    term = total = F(1)
    for j in range(1,order+1):
        term *= -x*x / (2*j*(2*j-1))
        total += term
    return total


def cosine_table():
    # Machin's identity pi/4 = 4 atan(1/5) - atan(1/239).
    al, au = atan_interval(5)
    bl, bu = atan_interval(239)
    pi_lo, pi_hi = 16*al-4*bu, 16*au-4*bl
    # A shorter rational bracket keeps the Taylor fractions small.
    p = 10**16
    pi_lo, pi_hi = F(math.floor(pi_lo*p),p), F(math.ceil(pi_hi*p),p)
    require(F(314,100)<pi_lo<pi_hi<F(315,100), 'pi interval')
    quarter = [(SCALE,SCALE)]
    for k in range(1,TICKS//4):
        xlo, xhi = 2*k*pi_lo/TICKS, 2*k*pi_hi/TICKS
        require(xhi < pi_lo/2, 'cosine decreases on the first quadrant')
        quarter.append((math.floor(SCALE*cosine_sum(xhi,9)),math.ceil(SCALE*cosine_sum(xlo,10))))
    quarter.append((0,0))
    table = []
    for k in range(TICKS+1):
        j = min(k,TICKS-k)
        if j <= TICKS//4:
            table.append(quarter[j])
        else:
            lo, hi = quarter[TICKS//2-j]
            table.append((-hi,-lo))
    require(F(table[64][0],SCALE)>F(271,1250), 'distinct initial bins for all outer angles')
    return table


def replay(compiler):
    table = cosine_table()
    with TemporaryDirectory(prefix='mordell-radial-hessian-') as scratch:
        scratch = Path(scratch)
        angles = scratch/'cosines.txt'
        angles.write_text(''.join(f'{lo} {hi}\n' for lo,hi in table),encoding='ascii')
        binary = scratch/'verify.exe'
        subprocess.run([compiler, '-std=c++20', '-O3', str(HERE/'verify.cpp'), '-o', str(binary)],
                       check=True,timeout=120)
        output = subprocess.run([str(binary),str(angles)],capture_output=True,encoding='ascii',
                                check=True,timeout=600).stdout
    nodes, geometry, matrices, unresolved, depth = map(int,output.split())
    roots = math.comb(32,4)
    leaves = geometry+matrices+unresolved
    require(2*leaves-roots == nodes and matrices>0 and unresolved==0, 'complete matrix covering tree')
    return {'roots':roots, 'nodes':nodes, 'geometry_leaves':geometry, 'matrix_leaves':matrices,
            'unresolved':unresolved, 'maximum_subdivision_depth':depth}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=HERE.parent)
    parser.add_argument('--compiler',default='g++')
    parser.add_argument('--record-report',action='store_true')
    args = parser.parse_args()
    floor = prerequisites(args.repo.resolve())
    result = {'passed':True, 'normalization_mass':9, 'assumption':'half-turn symmetry',
              'feasible_outer_modulus_product_lower':str(floor),
              'outer_modulus_geometric_mean_domain_lower':'43/40',
              'negative_log_radial_hessian_strict_lower':'1/4', 'radial_variables_eliminated':5,
              'remaining_angle_parameters':4, 'ceiling_squared_gradient_coefficient':'2',
              'radial_stationarity_implies_unique_radii_at_fixed_angles':True, 'cover':replay(args.compiler)}
    encoded = json.dumps(result,indent=2)+'\n'
    path = HERE/'verification_report.json'
    if args.record_report:
        path.write_text(encoded,encoding='ascii',newline='\n')
    else:
        require(result==read(path),'matches the frozen matrix cover')
    print(encoded,end='')


if __name__=='__main__':
    main()
