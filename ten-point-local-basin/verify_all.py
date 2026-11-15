"""Exact enlarged orbit basin for the pinned ten-point candidate."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def identity(n, zero, one):
    return [[one if i == j else zero for j in range(n)] for i in range(n)]


def rational_spd(M, scale, label):
    """A midpoint and a uniform diagonal reduction dominate interval errors."""
    n = len(M)
    unit = 1 << 40
    raw = [[0]*n for _ in range(n)]
    row_errors = [0]*n
    for i in range(n):
        for j in range(i+1):
            L, R = M[i][j], M[j][i]
            lo, hi = min(L.lo, R.lo), max(L.hi, R.hi)
            mid = (lo+hi)*unit//(2*scale)
            raw[i][j] = raw[j][i] = mid
            err = max(mid*scale-lo*unit, hi*unit-mid*scale)
            row_errors[i] += err
            if j != i:
                row_errors[j] += err
    shift = (max(row_errors)+scale-1)//scale
    A = [[F(raw[i][j]-(shift if i == j else 0),unit) for j in range(n)] for i in range(n)]
    T = [row.copy() for row in A]
    pivots = []
    for k in range(n):
        p = T[k][k]
        need(p > 0, label + ' strict rational LDL pivot')
        pivots.append(p)
        for i in range(k+1,n):
            for j in range(k+1,i+1):
                T[i][j] -= T[i][k]*T[j][k]/p
    record = {'dimension': n, 'dyadic_diagonal_shift': str(F(shift,unit)),
              'minimum_ldl_pivot_floor': str(F(min(pivots)*unit//1,unit))}
    payload = {'matrix': {'domain': 'QQ', 'entries': [
        [{'num': str(q.numerator), 'den': str(q.denominator)} for q in row] for row in A]}}
    return record, payload


def verify():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for path, expected in manifest['files'].items():
        need(hashlib.sha256((HERE/path).read_bytes()).hexdigest() == expected, 'fixed own input '+path)
    for path, expected in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, 'fixed source '+path)
    spec = importlib.util.spec_from_file_location('symmetry', ROOT/'candidate-global-screen/verify_symmetry.py')
    symmetry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(symmetry)
    prior = symmetry.load_prior()
    common = json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    catalog = json.loads((ROOT/'candidate-global-screen/symmetry_certificates.json').read_bytes())
    n = 10
    row = next(r for r in common['rows'] if r['n'] == n)
    entry = next(r for r in catalog['rows'] if r['n'] == n)
    sym = symmetry.verify_row(entry, prior, common)
    need(sym['reflection_exact'], 'replayed exact reflection of actual root')
    B, S = prior.Ball, prior.SCALE
    zero, one = B.integer(0), B.integer(1)
    err = F(sym['root_refinement_error_upper'])*S
    err = (err.numerator+err.denominator-1)//err.denominator
    mid = [B.dyadic(v,common['center_exponent']) for v in row['center_numerators']]
    z = [B(v.lo-err,v.hi+err) for v in mid]
    z.insert(1,zero)
    q = [v for i in range(n) for v in (-z[2*i+1],z[2*i])]
    G = [[q[i]*q[j]/n for j in range(2*n)] for i in range(2*n)]
    for i in range(2*n):
        G[i][i] += B.integer(n-1)/2
    L = {k: identity(n,zero,B.integer(bound)) for k,bound in ((3,14),(4,17))}
    for i in range(n):
        for j in range(i):
            d = [z[2*i+a]-z[2*j+a] for a in range(2)]
            s = sum(x.square() for x in d)
            need(100*s.lo > 49*S, 'all candidate distances exceed 7/10')
            D = symmetry.sqrt_ball(s,B,S)
            for k,mat in L.items():
                w = 1/(D*s if k == 3 else s.square())
                mat[i][i] -= w; mat[j][j] -= w
                mat[i][j] += w; mat[j][i] += w
            for a in range(2):
                for b in range(2):
                    v = (one if a == b else zero)/s-2*d[a]*d[b]/s.square()
                    G[2*i+a][2*i+b] -= v; G[2*j+a][2*j+b] -= v
                    G[2*i+a][2*j+b] += v; G[2*j+a][2*i+b] += v
    psi = entry['reflection_permutation']
    alpha,beta = z[2*psi[0]],z[2*psi[0]+1]
    norm = symmetry.sqrt_ball(alpha.square()+beta.square(),B,S)
    R = [[alpha/norm,beta/norm],[beta/norm,-alpha/norm]]
    gamma_minus, gamma_plus = F(3,40), F(22,25)
    mid_gamma = (gamma_minus+gamma_plus)/2
    gap_gamma = (gamma_plus-gamma_minus)/2
    for i in range(2*n):
        G[i][i] -= B.integer(mid_gamma.numerator)/mid_gamma.denominator
    for i in range(n):
        for a in range(2):
            for b in range(2):
                G[2*i+a][2*psi[i]+b] -= gap_gamma.numerator*R[a][b]/gap_gamma.denominator
    records, payloads = {}, {}
    for key,M in (('reflection_split_curvature',G),('length_laplacian_3',L[3]),('length_laplacian_4',L[4])):
        records[key],payloads[key] = rational_spd(M,S,key)
    radius, two_norm, min_distance = F(1,100), F(10,7), F(7,10)
    need(two_norm**2 > 2, 'upper bound on pair difference operator')
    r = two_norm*radius/min_distance
    c = 2*two_norm*14
    tail = two_norm**2*radius**2*17/(2*(1-r))
    a = gamma_minus-tail
    b = gamma_plus-tail-c*radius/3
    d = c*radius/2
    need(r < 1 and a > 0 and a*b-d*d > 0, 'strict upper Taylor remainder budget')
    report = {'passed': True, 'normalization_mass': 10, 'euclidean_orbit_radius': str(radius),
              'previous_euclidean_orbit_radius': '1/262144', 'candidate_value_cap_proved': True,
              'equality_only_candidate_orbit_within_this_radius': True, 'global_complement_covered': False,
              'reflection_gaps_odd_even': [str(gamma_minus),str(gamma_plus)],
              'normalized_pair_displacement_upper': str(r),
              'quartic_tail_coefficient_upper': str(tail),
              'quadratic_remainder_matrix': [[str(a),str(-d)],[str(-d),str(b)]],
              'quadratic_remainder_determinant': str(a*b-d*d), 'matrix_checks': records}
    return report,payloads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-jacobian', type=Path)
    args = parser.parse_args()
    report, payloads = verify()
    if args.export_jacobian is not None:
        args.export_jacobian.write_text(json.dumps(payloads,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
