"""Exact outer harmonics and rank-antipodal angular pairing at n=10."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from math import comb, isqrt
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
REFINED=ROOT/'ten-point-banded-kernel-refined'
RADII=ROOT/'ten-point-refined-radii'
BANDS=((F(37,500),F(43,250)),(F(169,200),F(164,125)))
BUDGET=F(3391,500000)
TARGET={1:F(39,100),2:F(9,10),3:F(39,40),4:F(37,25)}
ALPHA_COS=F(106,125)
BETA_COS=F(26,125)
ODD_PROJECTION_BOUND=F(377,1000)
ODD_WEIGHT=F(11,100)
DISTANCE_SQUARED_LOWER=F(2,5)

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

upper=load('upper',ROOT/'ten-point-ward-upper/verify_all.py')
radial=load('radial',ROOT/'ten-point-ward-inner-count/verify_radial.py')

def need(condition, why):
    if not condition:
        raise ArithmeticError(why)

def add(left,right,scale=F(1)):
    out=[F(0)]*max(len(left),len(right))
    for i,x in enumerate(left):out[i]+=x
    for i,x in enumerate(right):out[i]+=scale*x
    return out

def expand(c):
    return [sum(c[j]*comb(j,i)*(-1)**(j-i) for j in range(i,len(c)))
            for i in range(len(c))]

def coefficients(kernel,k,record):
    E=[[upper.hx(x) for x in row] for row in kernel['transforms'][k-1]]
    Q=[[upper.hx(x) for x in row] for row in kernel['Q_blocks'][k]]
    y=tuple(map(F,record['q_basis']))
    W=tuple(map(F,record['ward']))
    need(len(E)==6 and all(len(row)==len(y) for row in E)
         and Q==[list(q) for q in zip(*Q)] and len(Q)==len(y),'compatible exact Gram basis')
    M=[row[:]+[v] for row,v in zip(Q,y)]
    upper.gram.positive_ldl(Q)
    for j in range(len(y)):
        need(M[j][j]!=0,'nonzero exact solve pivot')
        M[j]=[x/M[j][j] for x in M[j]]
        for i in range(len(y)):
            if i!=j:
                a=M[i][j]
                M[i]=[b-a*c for b,c in zip(M[i],M[j])]
    z=[row[-1] for row in M]
    need(all(sum(a*b for a,b in zip(row,z))==v for row,v in zip(Q,y)),'exact Gram solve substitution')
    g=[sum(a*b for a,b in zip(row,y)) for row in E]
    if k==1:
        need(len(W)==2,'translation & first closure Ward weights')
        g[0]+=W[0];g[1]+=W[1]
    elif k in (2,3):
        need(len(W)==1,'mixed Ward weight')
        g[0]-=F(9-k,10)*W[0];g[1]+=W[0]
    else:
        need(not W,'no extra Ward equations')
    return expand(g),sum(a*b for a,b in zip(y,z))

def harmonics(kernel,rows):
    report={}
    for k,target in TARGET.items():
        row=rows[str(k)]
        g,norm=coefficients(kernel,k,row)
        square=radial.mul([F(0)]*k+[F(1)],radial.mul(g,g))
        e=tuple(map(F,row['e']))
        slope=F(row['slope'])
        covers=[]
        for i,(lo,hi) in enumerate(BANDS):
            affine=[e[i]-F(9,10)*slope,slope]
            edge=[affine[0]+affine[1]*u for u in (lo,hi)]
            need(min(edge)>0,'positive affine approximation radius')
            if i==0:
                bound=radial.mul(affine,affine)
                covers.append(radial.strictly_above(add(bound,square,-1),(lo,hi),F(0)))
            else:
                need(max(edge)<1,'outer lower envelope may be squared')
                covers.append(radial.strictly_above(g,(lo,hi),F(0)))
                lower=add([F(1)],affine,-1)
                upper_bound=add([F(1)],affine)
                covers.append(radial.strictly_above(add(square,radial.mul(lower,lower),-1),(lo,hi),F(0)))
                covers.append(radial.strictly_above(add(radial.mul(upper_bound,upper_bound),square,-1),(lo,hi),F(0)))
        common_error=2*e[0]+8*e[1]
        cap=target-common_error
        need(norm>0 and cap>0 and cap**2>BUDGET*norm,'harmonic error and Gram bound fit target')
        report[str(k)]={'modulus_strict_upper':str(target),'fixed_approximation_error':str(common_error),
                        'dual_gram_norm_squared_decimal_display':float(norm),
                        'one_variable_covers':covers,'one_variable_nodes':sum(c['visited'] for c in covers)}
    return report

def split_vector(values):
    left, right = [values[0]], [values[-1]]
    while len(values) > 1:
        values = [(a+b)/2 for a,b in zip(values,values[1:])]
        left.append(values[0])
        right.append(values[-1])
    return left, right[::-1]

def split_control(control,axis):
    rows = list(zip(*control)) if axis == 0 else control
    parts = [split_vector(row) for row in rows]
    left, right = list(map(list,zip(*parts)))
    return (list(zip(*left)),list(zip(*right))) if axis == 0 else (left,right)

def check_bivariate(polynomial,radius_band):
    # The first Bernstein axis is the radius and the second is the cosine.
    columns=[radial.bernstein(list(col),*radius_band) for col in zip(*polynomial)]
    tensor=[radial.bernstein(list(row),F(0),F(1)) for row in zip(*columns)]
    todo=[(tensor,0)]
    nodes=splits=leaves=depth_max=0
    while todo:
        control,depth=todo.pop()
        nodes+=1
        depth_max=max(depth_max,depth)
        if all(min(row)>0 for row in control):
            leaves+=1
        else:
            need(depth<30,'unresolved joint trigonometric approximation')
            lo,hi=split_control(control,depth%2)
            todo.extend(((lo,depth+1),(hi,depth+1)))
            splits+=1
    need(nodes==2*splits+1 and leaves==splits+1,'complete two-variable Bernstein tree')
    return {'visited':nodes,'maximum_depth':depth_max}

def joint_projection(kernel,record):
    g,norm1=coefficients(kernel,1,record['1'])
    h,norm3=coefficients(kernel,3,record['3'])
    e=tuple(map(F,record['e']))
    slope=F(record['slope'])
    radius_bands=((F(272,1000),F(415,1000)),(F(919,1000),F(1146,1000)))
    covers=[]
    for band,(rlo,rhi) in enumerate(radius_bands):
        ulo,uhi=BANDS[band]
        need(0<rlo*rlo<ulo<uhi<rhi*rhi,'rational support contains square roots of band endpoints')
        # (r*g(r虏)-I)c - lambda*(r鲁*h(r虏)-I)(4c鲁-3c).
        first=[-F(band)]+[F(0)]*13
        third=first[:]
        for j,x in enumerate(g):first[2*j+1]=x
        for j,x in enumerate(h):third[2*j+3]=x
        affine=[e[band]-F(9,10)*slope,F(0),slope]+[F(0)]*11
        for sign in (-1,1):
            pol=[[a,sign*(b+3*ODD_WEIGHT*c),F(0),-4*sign*ODD_WEIGHT*c]
                 for a,b,c in zip(affine,first,third)]
            covers.append(check_bivariate(pol,(rlo,rhi)))
    dual=norm1+ODD_WEIGHT**2*norm3
    common_error=2*e[0]+8*e[1]
    cap=ODD_PROJECTION_BOUND-common_error
    need(norm1>0 and norm3>0 and cap>0 and cap**2>BUDGET*dual,
         'combined Gram and approximation error fit the odd directional bound')
    return {'odd_directional_projection_absolute_strict_upper':str(ODD_PROJECTION_BOUND),
            'fixed_approximation_error':str(common_error),
            'joint_dual_gram_norm_squared_decimal_display':float(dual),
            'two_variable_covers':covers,
            'two_variable_nodes':sum(row['visited'] for row in covers)}

def f_at_cos(y):
    return y-ODD_WEIGHT*(4*y**3-3*y)

def pairing():
    a,b=ALPHA_COS,BETA_COS
    # The checked polynomial is separately convex in both radial variables.
    need(all((u+v-DISTANCE_SQUARED_LOWER)**2 <4*a*a*u*v
             for u in BANDS[1] for v in BANDS[1]),'strict angular separation at corners')
    # Rational square-root enclosure for cos(beta-alpha).
    val=(1-b*b)*(1-a*a)
    pow2=1<<60
    x=isqrt((val*pow2**2).numerator//(val*pow2**2).denominator)
    lo,hi=F(x,pow2),F(x+1,pow2)
    need(lo*lo<=val<hi*hi,'exact square root enclosure')
    cb=[(b,b),(a*b+lo,a*b+hi)]
    for j in range(2,5):
        cb.append((2*a*cb[j-1][0]-cb[j-2][1],2*a*cb[j-1][1]-cb[j-2][0]))
    need(all(0<c<=d<1 for c,d in cb),'each projected vertex cosine lies in the positive range')
    need(ODD_WEIGHT<F(1,9),'odd trigonometric projection is concave on the forward semicircle')
    values=[f_at_cos(pair[0]) for pair in cb]
    vertices=[sum(values[:j])+sum(values[:5-j]) for j in range(6)]
    other_upper=f_at_cos(F(1))+2*f_at_cos(a)
    need(other_upper>2*f_at_cos(F(1)),'triple upper also covers a nonpositive forward projection')
    required=ODD_PROJECTION_BOUND
    net=min(vertices)-other_upper
    need(net>required,'five ray arc would violate the first and third Fourier bounds')
    kappa=1-2*b*b
    paired_midpoint_cap=F(73,1250)
    for u in BANDS[1]:
        for v in BANDS[1]:
            w=u+v-4*paired_midpoint_cap
            need(w>0 and w*w<4*kappa*kappa*u*v,'paired midpoint upper endpoint')
    return {'angular_separation_cosine_strict_upper':str(a),
            'closed_arc_half_angle_cosine':str(b),'closed_arc_occupancy_upper':4,
            'cyclic_rank_pairing_cosine_strict_upper':str(2*b*b-1),
            'matching_outer_pairs':4,'five_point_arc_vertices_checked':6,
            'paired_squared_distance_strict_lower':str(2*BANDS[1][0]*(1+kappa)),
            'paired_midpoint_squared_strict_upper':str(paired_midpoint_cap),
            'arc_odd_directional_moment_lower_decimal_display':float(net),
            'global_odd_directional_moment_strict_upper':str(required)}

def verify():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for names, directory in (('files', HERE), ('dependencies', ROOT)):
        for name, expected in manifest[names].items():
            need(sha256((directory/name).read_bytes()).hexdigest() == expected, 'fixed input '+name)
    refined = json.loads((REFINED/'verification_report.json').read_bytes())
    radii = json.loads((RADII/'verification_report.json').read_bytes())
    need(refined['passed'] is True and refined['unresolved'] == 0
         and F(refined['c']) == F(3,40000) and F(refined['c_inner']) == F(9,100000)
         and refined['PSD_blocks'] == 10
         and refined['exact_inner_count'] == 2 and refined['exact_outer_count'] == 8,
         'refined global kernel and exact population')
    need(radii['passed'] is True and radii['applies_to_every_global_maximizer'] is True
         and radii['normalization_mass'] == 9 and radii['unresolved'] == 0
         and radii['exact_inner_count'] == 2 and radii['exact_outer_count'] == 8,
         'refined radial conclusion applies to all maximizers')
    need(tuple(map(F,radii['inner_squared_radius_open_interval'])) == BANDS[0]
         and tuple(map(F,radii['outer_squared_radius_open_interval'])) == BANDS[1],
         'fixed exact radial intervals')
    data = json.loads((REFINED/'kernel.json').read_bytes())
    need((data['n'],data['mass']) == (10,9) and upper.DMIN == DISTANCE_SQUARED_LOWER,
         'same normalization and necessary spacing')
    gap = F(1,10**6)
    need(F(refined['B']) == upper.hx(data['B'])
         and upper.hx(data['B']) - upper.log.log_lower(upper.witness()) < gap
         and 2*(gap + 44*F(refined['c']) + F(refined['c_inner'])) == BUDGET,
         'strict common ordered Gram budget')
    rows = json.loads((HERE/'harmonics.json').read_bytes())
    joint = json.loads((HERE/'odd-projection.json').read_bytes())
    return {
        'passed': True,
        'normalization_mass': 9,
        'applies_to_every_global_maximizer': True,
        'outer_count': 8,
        'common_ordered_gram_budget_strict_upper': str(BUDGET),
        'outer_normalized_harmonics': harmonics(data, rows),
        'joint_projection': joint_projection(data, joint),
        'pairing': pairing(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    args = parser.parse_args()
    report = verify()
    path = HERE/'verification_report.json'
    raw = (json.dumps(report, indent=2)+'\n').encode()
    if args.record_report:
        path.write_bytes(raw)
    else:
        need(report == json.loads(path.read_bytes()), 'default replay matches the saved report')
    print(raw.decode(),end='')


if __name__ == '__main__':
    main()
