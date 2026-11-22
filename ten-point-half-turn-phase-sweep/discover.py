"""Propose rational radial references and dyadic phase cells."""

from fractions import Fraction as F
import json
from math import atan

import numpy as np
from scipy.optimize import root

from verify_all import TOLERANCES, HERE, Sweep


def candidate_references(sweep, mode):
    box = sweep.box
    directions = [complex(*(float(box.endpoint(v)) for v in u)) for u in sweep.directions[mode]]
    initial = np.log([.11, 1.1, 1.1, 1.1, 1.1])
    theta = np.angle(directions)
    i, j = np.triu_indices(5, 1)

    def propose(left, right):
        phase = 2*atan(float((left+right)/2))
        angles = theta + np.array([0,phase,phase,phase,phase])
        cosine = np.cos(angles[i]-angles[j])

        def gradient(s):
            difference = s[i]-s[j]
            term = 2*np.sinh(difference)/(np.cosh(difference)-cosine)
            return (9-10*np.exp(s)+np.bincount(i,term,minlength=5)
                    -np.bincount(j,term,minlength=5))

        solution = root(gradient, initial)
        if not solution.success:
            raise ArithmeticError('reference root proposal')
        radii = np.exp(solution.x)
        inner = np.clip(radii[0], .10750001, .13249999)
        radii[1:] += (4.5-inner-sum(radii[1:]))/4
        values = [F(round(x*10**9), 10**9) for x in [inner,*radii[1:4]]]
        return values+[F(9,2)-sum(values)]

    return propose


def cover(sweep, mode):
    propose = candidate_references(sweep, mode)
    stack = [(F(-1), F(1), 0)]
    records = []
    while stack:
        left, right, depth = stack.pop()
        if sweep.near_orbit(mode, left, right):
            values = None
        else:
            values = propose(left, right)
        action = sweep.check(mode, left, right, values)
        if action is not None:
            records.append([str(left),str(right),action,list(map(str,values)) if values else []])
        elif depth < 18:
            middle = (left+right)/2
            stack.extend([(middle,right,depth+1),(left,middle,depth+1)])
        else:
            raise ArithmeticError(f'{mode}: unresolved cell {left} {right}')
    return records


if __name__ == '__main__':
    sweep = Sweep()
    for mode in TOLERANCES:
        records = cover(sweep, mode)
        print(mode, len(records), flush=True)
        data = ',\n'.join('  '+json.dumps(record) for record in records)
        (HERE/(mode+'.json')).write_bytes(('[\n'+data+'\n]\n').encode())
