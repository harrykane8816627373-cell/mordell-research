"""Exploratory C2 stationary-point search with reflection unconstrained."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize, root

KEEP = np.arange(10) != 1
I, J = np.triu_indices(5, 1)


def coordinates(x):
    return np.insert(x, 1, 0.).view(complex)


def force(a):
    wm = a[:, None] - a[None, :]
    wp = a[:, None] + a[None, :]
    for v in (wm, wp):
        np.fill_diagonal(v, 1+0j)
    dminus, dplus = 1/wm, 1/wp
    np.fill_diagonal(dminus, 0)
    np.fill_diagonal(dplus, 0)
    return (0.5/a + np.sum(dminus + dplus, axis=1)).conjugate() - 5*a, dminus, dplus


def derivative(a, dm, dp):
    left, right = dm**2, dp**2
    w = left-right
    np.fill_diagonal(w, -0.5/a**2 - np.sum(left+right, axis=1))
    Q = np.empty((10, 10))
    Q[::2, ::2] = w.real
    Q[::2, 1::2] = Q[1::2, ::2] = -w.imag
    Q[1::2, 1::2] = -w.real
    Q -= 5*np.eye(10)
    return Q


def log_delta(a):
    w = a*a
    return 10*np.log(2.) + 2*np.log(abs(a)).sum() + 4*np.log(abs(w[I]-w[J])).sum()


def objective(x):
    a = coordinates(x)
    g, _, _ = force(a)
    return -log_delta(a) + 10*np.vdot(a, a).real, (-4*g).view(float)[KEEP]


def root_solve(a):
    def residual(x):
        a = coordinates(x)
        f, dm, dp = force(a)
        return f.view(float)[KEEP], derivative(a, dm, dp)[np.ix_(KEEP, KEEP)]

    result = root(residual, a.view(float)[KEEP], jac=True, method='hybr',
                  options={'xtol': 2e-11, 'maxfev': 225})
    a = coordinates(result.x)
    residual_size = float(np.linalg.norm(force(a)[0].view(float), ord=np.inf))
    return a, residual_size


def normalize(a):
    a = a*np.sqrt(4.5/np.vdot(a, a).real)
    return a*a[0].conjugate()/abs(a[0])


def canonical(a):
    j = np.argmin(abs(a))
    a = np.r_[a[j], a[:j], a[j+1:]]
    w = (a*a)*(a[0].conjugate()/abs(a[0]))**2
    order = np.argsort(np.angle(w[1:])%(2*np.pi))
    return w[np.r_[0, order+1]]


def report(a, residual):
    Q = 2*np.vdot(a, a).real
    delta = log_delta(a) + 45*np.log(9/Q)
    f, dm, dp = force(a)
    e = np.linalg.eigvalsh(derivative(a, dm, dp)[np.ix_(KEEP, KEEP)])
    w = canonical(a)
    radial = abs(w)
    costs = abs(w[1:, None]-w[None, 1:].conjugate())
    rr, cc = linear_sum_assignment(costs)
    return {
        'mass9_log_delta': float(delta), 'residual': residual,
        'positive_hessian_count_in_C2_slice': int(np.sum(e>1e-6)),
        'pair_squared_radii': radial.tolist(),
        'within_global_bands': bool(.1025 < radial[0] < .1375 and ((.845 < radial[1:]) & (radial[1:] < 1.312)).all()),
        'distance_to_conjugation': float(np.linalg.norm(costs[rr, cc])),
        'square_images': [[float(t.real), float(t.imag)] for t in w],
    }


def typical_start(rng, run):
    if run%3 == 0:
        gaps = rng.dirichlet(np.full(4, 4.))*np.pi
    elif run%3 == 1:
        gaps = (np.pi/4)*(rng.uniform(.80, 1.20, size=4))
    else:
        gaps = rng.dirichlet(np.full(4, 1.))*np.pi
    angles = (rng.uniform(0., np.pi)+np.r_[0., np.cumsum(gaps[:3])]) % np.pi
    radii = np.sqrt(rng.uniform(.78, 1.35, size=4))
    return normalize(np.r_[rng.uniform(.305, .39), radii*np.exp(1j*angles)])


def seeded_start(rng, run):
    pairs = ((.11183669980128009,0.), (1.0074790085402798,np.pi),
             (1.188144595492442,0.), (1.096269848082956,1.6374), (1.096269848082956,-1.6374))
    root = np.array([np.sqrt(r)*np.exp(0.5j*t) for r,t in pairs])
    scale = [.005,.015,.025,.05,.08,.12,.2][run%7]
    return normalize(root+scale*(rng.standard_normal(5)+1j*rng.standard_normal(5)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--starts', type=int, default=20000)
    p.add_argument('--optimize-first', type=int, default=2000)
    p.add_argument('--seed', type=int, default=1011202616)
    p.add_argument('--output', type=Path, default=Path('half-turn-search-replay.json'))
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)
    clusters = []
    success = 0
    for i in range(args.starts):
        a = typical_start(rng, i) if i%5 < 4 else seeded_start(rng, i)
        if i < args.optimize_first:
            opt = minimize(objective, a.view(float)[KEEP], jac=True, method='L-BFGS-B',
                           options={'ftol':1e-12, 'gtol':1e-8, 'maxiter':500})
            a = coordinates(opt.x)
        a, error = root_solve(a)
        if not 0<=error<1e-8:
            continue
        success += 1
        value = log_delta(a) + 45*np.log(9/(2*np.vdot(a, a).real))
        found = next((cluster for cluster in clusters if abs(cluster['mass9_log_delta']-value)<1e-7), None)
        if found is None:
            candidate = report(a, error)
            candidate['occurrences'] = 1
            clusters.append(candidate)
            if candidate['within_global_bands']:
                print('New within global bands', json.dumps(candidate), flush=True)
        else:
            found['occurrences'] += 1
        if (i+1)%1000 == 0:
            print('progress', i+1, success, len(clusters), flush=True)
    records = {'settings': {k:v for k,v in vars(args).items() if k!='output'}, 'accepted': success,
               'clusters': sorted(clusters, key=lambda c: -c['mass9_log_delta'])}
    args.output.write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
