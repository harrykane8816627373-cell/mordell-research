"""Floating choices only. All accepted facts are checked by verify_critical_box.py."""
from pathlib import Path
import json
import numpy as np
import interval_engine as E
from propose_cover import midpoint_objects, numeric_columns
from verify_orbit import check

HERE = Path(__file__).resolve().parent
E.setup(60)
row = json.loads((HERE/'exploration/scout.json').read_text())['roots'][5]
q = np.array([complex(*z) for z in row['q']])
q = np.r_[q[0], np.array(sorted(q[1:], key=lambda z: (z.real, z.imag)))]
if q[1].imag < 0:
    q = q.conjugate()
xx = [q[0].real]+[v for z in q[1:7] for v in [z.real, z.imag]]
center = [int(round(v*E.S)) for v in xx]
root_radius = 1 << 30
qm = E.positions([E.point(v) for v in center])
J = E.polynomial_jacobian(qm)
J = np.array([[(a+b)/(2*E.S) for a, b in row] for row in J])
R = [[int(round(v*E.S)) for v in row] for row in np.linalg.inv(J)]
dummy_box = [(v-root_radius, v+root_radius) for v in center]
H = midpoint_objects(dummy_box, center)[0]
values, vectors = np.linalg.eigh(H)
v = vectors[:, 0]
direction = [[int(round(v[i]*1024)), int(round(v[8+i]*1024))] for i in range(8)]
V, F, sv = numeric_columns(center)
M = 64*(V[:, 2:].conj().T@V[:, 2:])-8*np.outer(sv.conj(), sv)-F[:, 2:].conj().T@F[:, 2:]
values, vectors = np.linalg.eigh((M+M.conj().T)/2)
c = np.zeros(8, dtype=complex)
c[2:] = vectors[:, 0]
c[0] = -sv@c[2:]/8
coeff = [[int(round(z.real*1024)), int(round(z.imag*1024))] for z in c]
radii = [E.S//n for n in [100, 200, 400, 800, 1600, 3200, 6400]]
hr, pr = None, None
for r in radii:
    box = [(x-r, x+r) for x in center]
    qq = E.positions(box)
    inv, _ = E.invert_edges(qq)
    if hr is None and E.rayleigh(qq, inv, direction)[1] < 0:
        hr = r
    if pr is None and E.polynomial_rayleigh(qq, coeff)[1] < 0:
        pr = r
if hr is None or pr is None:
    raise RuntimeError('No certified broad negative-curvature box')
data = {'precision_bits': 60, 'center': center, 'root_radius': root_radius,
        'preconditioner': R, 'integer_displacement': direction,
        'integer_polynomial_coefficients': coeff, 'Hessian_box_radius': hr,
        'polynomial_box_radius': pr,
        'proposer': 'Numerical scout index 5; not a complete stationary classification.'}
path = HERE/'certificates/critical_orbit.json'
path.write_text(json.dumps(data, indent=2)+'\n')
result = check(path)
(HERE/'reports/orbit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
