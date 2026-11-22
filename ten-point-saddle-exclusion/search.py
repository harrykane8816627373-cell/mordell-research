"""Finite floating n=10 exploration with no symmetry restrictions."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize, root

N = 10
HERE = Path(__file__).resolve().parent
ROW, COL = np.triu_indices(N, 1)
KEEP = np.arange(2*N) != 1
TOL = 1e-7


def normalize(z):
    z = z-z.mean()
    return z*np.sqrt(N/np.vdot(z,z).real)


def log_product(z):
    z = normalize(z)
    return float(2*np.log(np.abs(z[ROW]-z[COL])).sum())


def objective(flat):
    z = flat[::2]+1j*flat[1::2]
    diff = z[:, None]-z[None, :]
    distances = diff.real**2 + diff.imag**2
    np.fill_diagonal(distances, 1.0)
    if not (distances > 0).all():
        return np.inf, np.zeros_like(flat)
    energy = -np.log(distances[ROW, COL]).sum()+4.5*np.vdot(z,z).real
    grad = 9*z-2*np.sum(diff/distances, axis=1)
    return float(energy), grad.view(float)


def gauge(z):
    i = np.argmax(abs(z))
    z = np.r_[z[i], z[:i], z[i+1:]]
    z = z*z[0].conjugate()/abs(z[0])
    return z.view(float)[KEEP]


def expand(part):
    flat = np.insert(part, 1, 0.0)
    return flat[::2]+1j*flat[1::2]


def jacobian(z):
    diff = z[:, None]-z[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        weights = -1/np.conj(diff)**2
    np.fill_diagonal(weights, 0)
    lap = np.diag(weights.sum(axis=1))-weights
    jac = np.empty((2*N,2*N))
    jac[::2,::2] = lap.real
    jac[::2,1::2] = jac[1::2,::2] = lap.imag
    jac[1::2,1::2] = -lap.real
    jac -= 4.5*np.eye(2*N)
    return jac


def root_solve(z):
    def force(x):
        return objective(expand(x).view(float))[1][KEEP]/-2
    def force_jac(x):
        return jacobian(expand(x))[np.ix_(KEEP,KEEP)]
    result = root(force, gauge(z), jac=force_jac, method='hybr',
                  options={'xtol':2e-11,'maxfev':350})
    output = expand(result.x)
    residual = np.linalg.norm(objective(output.view(float))[1], ord=np.inf)
    return output, float(residual)


def optimize(z):
    data = minimize(objective, normalize(z).view(float), method='L-BFGS-B', jac=True,
                    options={'ftol': 2e-15, 'gtol': 1e-10, 'maxiter': 1000})
    z = data.x[::2]+1j*data.x[1::2]
    return root_solve(z)


def align(z, anchor):
    best_err = np.inf
    best = None
    for trial in (z,z.conjugate()):
        for shift in np.linspace(0,2*np.pi,36,endpoint=False):
            v = trial*np.exp(1j*shift)
            for _ in range(3):
                _, col = linear_sum_assignment(abs(anchor[:, None]-v[None,:])**2)
                v = v[col]
                phase = np.vdot(v,anchor)
                v *= phase/abs(phase)
            error = np.linalg.norm(v-anchor)
            if error < best_err:
                best_err = error
                best = v.copy()
    return best


def read_known(repo):
    data = json.loads((repo/'candidate-global-screen/data/n10.json').read_bytes())
    out = []
    for row in data['runs']:
        points = np.array(row['terminal_integer_points_grid_1e10'])*1e-10
        z = normalize(points[:,0]+1j*points[:,1])
        if all(abs(log_product(z)-log_product(w))>TOL for w in out):
            out.append(z)
    return sorted(out,key=log_product,reverse=True)


def seeds(known, rng, scale):
    candidate = known[0]
    theta = 2*np.pi*np.arange(9)/9
    centered_ring = normalize(np.r_[0,np.exp(1j*theta)])
    decagon = np.exp(2j*np.pi*np.arange(10)/10)
    random = lambda n: rng.normal(size=n)+1j*rng.normal(size=n)
    for known_i, v in enumerate([*known,centered_ring,decagon]):
        for s in (.004,.015,.04,.08,.12,.18,.28,.40,.60):
            for _ in range(scale):
                yield f'distort_{known_i}', v+s*random(10)

    for k in range(3*scale):
        yield 'gaussian', random(10)
        yield 'disk', np.sqrt(rng.uniform(size=10))*np.exp(2j*np.pi*rng.uniform(size=10))
        radii = np.r_[np.sqrt(rng.uniform(.074,.172,2)),np.sqrt(rng.uniform(.845,1.312,8))]
        gap = rng.dirichlet(np.full(8,3))*2*np.pi
        angles = np.r_[rng.uniform(0,np.pi),rng.uniform(np.pi,2*np.pi),
                       np.cumsum(gap)+rng.uniform(0,2*np.pi)]
        yield 'banded_2_8', radii*np.exp(1j*angles)
        yield 'banded_1_9', centered_ring+rng.uniform(.01,.6)*random(10)
        yield 'slit', rng.normal(size=10)+1j*rng.uniform(.02,.3)*rng.normal(size=10)
        rep_radii=np.sqrt(np.r_[rng.uniform(.074,.172),rng.uniform(.845,1.312,4)])
        pair_a=rep_radii*np.exp(1j*(2*np.pi*rng.uniform(size=5)))
        pair_m=rng.uniform(.02,.45)*random(5)
        yield 'noisy_antipodal', np.r_[pair_m+pair_a,pair_m-pair_a]
    for k in range(scale):
        corners = np.exp(2j*np.pi*rng.uniform(size=3))
        centers=corners[rng.integers(0,3,size=10)]
        yield 'clusters', centers+.18*random(10)

    for i,w in enumerate(known[1:]+[centered_ring,decagon]):
        v=align(w,candidate)
        for t in np.linspace(0,1,11):
            for _ in range(max(1, scale//5)):
                yield f'interpolated_{i}', (1-t)*candidate+t*v+.04*random(10)


def branch_record(z, residual):
    z=expand(gauge(z))
    spectrum=np.linalg.eigvalsh(jacobian(z)[np.ix_(KEEP,KEEP)])
    sorted_radii=sorted(abs(z)**2*.9)
    return {
        'log_product_mass_10': log_product(z),
        'force_gradient_max': residual,
        'hessian_positive_count': int((spectrum>1e-7).sum()),
        'smallest_absolute_reduced_force_eigenvalue':float(np.min(abs(spectrum))),
        'positive_reduced_force_eigenvalues':spectrum[spectrum>1e-7].tolist(),
        'squared_radii_mass_9':sorted_radii,
        'strict_radial_band_population':[
            sum(1 for x in sorted_radii if .074<x<.172),
            sum(1 for x in sorted_radii if .845<x<1.312)],
        'coordinates_mass_10':z.view(float).reshape(N,2).tolist(),
        'optimization_counts':{},
        'root_counts':{},
    }


def push(records,z,residual,kind,label):
    if residual > 1e-8 or not np.isfinite(residual):
        return False
    log=log_product(z)
    for row in records:
        if abs(row['log_product_mass_10']-log)<TOL:
            break
    else:
        records.append(row:=branch_record(z,residual))
    counters=row[kind+'_counts']
    counters[label]=counters.get(label,0)+1
    return True


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',type=Path,default=HERE.parent)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--scale',type=int,default=20)
    parser.add_argument('--start-seed',type=int,default=20261116010)
    parser.add_argument('--root-only',action='store_true')
    args=parser.parse_args()
    rng=np.random.default_rng(args.start_seed)
    known=read_known(args.repo)
    opt_counts=Counter()
    root_counts=Counter()
    branches=[]
    def save():
        branches.sort(key=lambda r:-r['log_product_mass_10'])
        data={'runs_by_stratum':dict(opt_counts),'root_attempts_by_stratum':dict(root_counts),
              'random_seed':args.start_seed,'scale':args.scale,'branches':branches}
        args.output.write_text(json.dumps(data,separators=(',',':'))+'\n',encoding='utf-8')
    for i,(label,z) in enumerate(seeds(known,rng,args.scale)):
        if not args.root_only:
            v,r=optimize(z)
            opt_counts[label]+=1
            push(branches,v,r,'optimization',label)
        if args.root_only or label.startswith('interpolated') or label.startswith('distort') and i%4==0:
            root_counts[label]+=1
            v,r=root_solve(normalize(z))
            push(branches,v,r,'root',label)
        if i%200==0:
            save()
            print('attempts',i+1,'observed branches',len(branches),flush=True)
    save()
    print('optimization attempts',sum(opt_counts.values()),'root attempts',sum(root_counts.values()),
          'observed branches',len(branches),flush=True)


if __name__=='__main__':
    main()
