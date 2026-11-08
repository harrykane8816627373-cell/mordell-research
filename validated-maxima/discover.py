"""Optional mpmath discovery; only verify_all.py supplies rigorous acceptance."""
from pathlib import Path
import json
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
mp.mp.dps = 160
CENTER_EXPONENT = 256
MATRIX_EXPONENT = 192


def force_jacobian(v, n):
    z = list(v[:1])+[mp.mpf(0)]+list(v[1:])
    lam = mp.mpf(n-1)/2
    f = mp.matrix([-lam*q for q in z])
    J = -lam*mp.eye(2*n)
    for i in range(n):
        for j in range(i):
            dx, dy = z[2*i]-z[2*j], z[2*i+1]-z[2*j+1]
            d2 = dx*dx+dy*dy
            for a, u in enumerate((dx,dy)):
                f[2*i+a] += u/d2
                f[2*j+a] -= u/d2
                for b, w in enumerate((dx,dy)):
                    q = (int(a == b)/d2)-2*u*w/d2**2
                    J[2*i+a,2*i+b] += q
                    J[2*j+a,2*j+b] += q
                    J[2*i+a,2*j+b] -= q
                    J[2*j+a,2*i+b] -= q
    keep = [i for i in range(2*n) if i != 1]
    return mp.matrix([f[i] for i in keep]), mp.matrix([[J[i,j] for j in keep] for i in keep])


def rounded(value, exponent):
    return int(mp.nint(value*mp.mpf(2)**exponent))


def matrix_numbers(matrix):
    return [[rounded(matrix[i,j], MATRIX_EXPONENT) for j in range(matrix.cols)] for i in range(matrix.rows)]


def discover(row, source):
    n = row['n']
    z = [mp.mpc(x,y) for x,y in row['integer_points']]
    z = [q-sum(z)/n for q in z]
    z = [q*mp.sqrt(n/sum(abs(w)**2 for w in z)) for q in z]
    first = max(range(n), key=lambda i: abs(z[i]))
    z = z[first:]+z[:first]
    rotation = mp.conj(z[0])/abs(z[0])
    z = [q*rotation for q in z]
    full = [q for p in z for q in (mp.re(p),mp.im(p))]
    v = mp.matrix(full[:1]+full[2:])
    for step in range(40):
        f,J = force_jacobian(v,n)
        if max(abs(q) for q in f) < mp.mpf(2)**-350:
            break
        v -= mp.lu_solve(J,f)
    else:
        raise ArithmeticError('discovery Newton did not converge for '+str(n))
    # Recompute all matrices at the exact dyadic midpoint written to the certificate.
    numbers = [rounded(q,CENTER_EXPONENT) for q in v]
    mid = mp.matrix([mp.mpf(q)/mp.mpf(2)**CENTER_EXPONENT for q in numbers])
    f,J = force_jacobian(mid,n)
    result = {'n':n, 'discovery_seed':source, 'discovery_integer_points':row['integer_points'], 'center_numerators':numbers,
              'inverse_numerators':matrix_numbers(mp.inverse(J)),
              'discovery_newton_steps':step,
              'floating_quantities_are_not_proof_inputs':True}
    try:
        L = mp.cholesky(-J)
        result['classification'] = 'strict_local_maximum_mod_rotation'
        result['whitener_numerators'] = matrix_numbers(mp.inverse(L.T))
    except ValueError:
        values, vectors = mp.eigsy(J)
        if values[-1] <= 0:
            # mpmath matrix negative indexing is not Python-list indexing.
            largest = values[values.rows-1]
        else:
            largest = values[-1]
        if largest <= 0:
            raise ArithmeticError('neither Cholesky nor a positive direction was found')
        vector = vectors[:,vectors.cols-1]
        vector /= max(abs(q) for q in vector)
        result['classification'] = 'unstable_stationary_orbit'
        result['positive_vector_numerators'] = [rounded(q,MATRIX_EXPONENT) for q in vector]
    print(n,result['classification'],'Newton steps',step,flush=True)
    return result


def main():
    old = json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_text(encoding='utf-8'))
    oldrows = {row['n']:row for row in old['rows']}
    data = {'date':'2026-11-05', 'arithmetic_precision':384, 'center_exponent':CENTER_EXPONENT,
            'matrix_exponent':MATRIX_EXPONENT, 'radius_exponent':80, 'rows':[]}
    for n in (10,11,12,13,14,15,19):
        data['rows'].append(discover(oldrows[n],'symmetry-breaking/integer_witnesses.json, n='+str(n)))
    new = json.loads((ROOT/'finite-search-nov5/integer_witnesses.json').read_text(encoding='utf-8'))
    newrows = {row['n']:row for row in new['rows']}
    for n in (17,18):
        data['rows'].append(discover(newrows[n],'finite-search-nov5/integer_witnesses.json, n='+str(n)))
    (HERE/'certificates.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
