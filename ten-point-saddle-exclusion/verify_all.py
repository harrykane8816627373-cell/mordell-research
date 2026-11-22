"""Replay exact n=10 saddle root and second-order exclusion."""
import argparse
from hashlib import sha256
import importlib.util
from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCIES = ('ten-point-contact-geometry/check_root_box.py',
                'validated-maxima/dyadic.py', 'validated-maxima/verify_all.py')
EVIDENCE = ('saddle-jacobian.json', 'saddle-inertia-jacobian.json')
NAMES = tuple('abcrx')
RADIUS = F(21, 500)
DIRECTION = [(0,-630),(0,630),(63,0),(-63,0),(0,224),(0,-224),
             (-118,-105),(118,105),(118,-105),(-118,105)]
FORCE_MULTIPLIER = [(-2195,0),(2195,0),(0,1306),(0,-1306),(-1446,0),(1446,0),
                    (641,-612),(-641,612),(641,612),(-641,-612)]


def need(condition, detail):
    if not condition:
        raise ArithmeticError(detail)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def root_image(data, h):
    a,b,c,r,x = ({tuple(int(i==j) for i in range(5)):1} for j in range(5))
    one = {(0,)*5:1}
    A,B,C,R = [one],[one],[one],[h.scaled(one,2),h.scaled(x,2)]
    for k in range(1,5):
        A.append(h.mul(a,A[-1]))
        B.append(h.mul(h.scaled(b,-1),B[-1]))
        C.append(h.mul(c,C[-1]))
        if k >= 2:
            R.append(h.plus(h.scaled(h.mul(x,R[-1]),2),h.scaled(h.mul(h.mul(r,r),R[-2]),-1)))
    S = [h.plus(A[k],B[k],C[k],R[k]) for k in range(5)]
    T = [h.plus(h.mul(a,A[k]),h.mul(b,B[k]),h.mul(c,C[k]),h.mul(r,R[k])) for k in range(5)]
    P = [h.plus(h.scaled(S[k],-(2*k+1)),h.scaled(T[k],-10),
                *(h.scaled(h.mul(S[j],S[k-j]),2) for j in range(k+1))) for k in range(5)]

    result = data['result']
    out = result['output']
    need(result['operation_id']=='polynomial.system.real_root_box.certify', 'root operation')
    need(data['input']['polynomial_map']==out['polynomial_map'] and data['input']['box']==out['box'],
         'root input matches result source')
    system = out['polynomial_map']
    need(tuple(system['input_variables'])==NAMES, 'ordered root variables')
    returned = system['output_polynomials']
    need(len(returned)==len(P), 'number of equations')
    for equation, component in zip(P,returned):
        terms = component['polynomial']['terms']
        values = {tuple(term['exponents']):h.rat(term['coefficient']) for term in terms}
        need(component['domain']=='QQ' and tuple(component['variables'])==NAMES and
             len(terms)==len(values) and equation==values, 'independently reconstructed polynomial')
    box = h.check_box(out['box'],NAMES)
    a,b,c,r,x = box
    need(0<a[0]<a[1]<b[0]<b[1]<r[0]<r[1]<c[0]<c[1] and
         r[0]>max(abs(x[0]),abs(x[1])), 'five distinct squared-plane roots')
    cert = out['conclusion']
    need(cert['status']=='CERTIFIED_UNIQUE_NONSINGULAR','recorded root result')
    ev = cert['evidence']
    center = [(l+u)/2 for l,u in box]
    need(tuple(ev['center']['variables'])==NAMES and list(map(h.rat, ev['center']['values']))==center,
         'rational center')
    point = [(v,v) for v in center]
    val = [h.ival(p,point)[0] for p in P]
    need(val==list(map(h.rat,ev['value_at_center'])), 'point value')
    jac = [[h.dif(p,k) for k in range(5)] for p in P]
    mid_jac = [[h.ival(p,point)[0] for p in row] for row in jac]
    interval_jac = [[h.ival(p,box) for p in row] for row in jac]
    need(mid_jac==[list(map(h.rat,row)) for row in ev['jacobian_at_center']['entries']] and
         interval_jac==[list(map(h.interval,row)) for row in ev['jacobian_enclosure']['entries']],
         'exact point and natural interval derivatives')
    inverse = [list(map(h.rat,row)) for row in ev['preconditioner']['entries']]
    need(len(inverse)==5 and all(len(row)==5 for row in inverse), 'inverse shape')
    need(all(sum(inverse[i][k]*mid_jac[k][j] for k in range(5))==int(i==j)
             for i in range(5) for j in range(5)), 'exact inverse')
    image = []
    norm = F(0)
    for i in range(5):
        contracted_row = []
        for j in range(5):
            l,u = h.isum(h.imul((inverse[i][k],)*2,interval_jac[k][j]) for k in range(5))
            contracted_row.append((int(i==j)-u,int(i==j)-l))
        norm = max(norm,sum(max(abs(l),abs(u)) for l,u in contracted_row))
        l,u = h.isum(h.imul(entry,(a0-c0,b0-c0))
                     for entry,(a0,b0),c0 in zip(contracted_row,box,center))
        fixed_center = center[i]-sum(w*v for w,v in zip(inverse[i],val))
        image.append((fixed_center+l,fixed_center+u))
    need(norm<F(1,1000) and image==h.check_box(ev['krawczyk_image'],NAMES) and
         all(l<a0<=b0<u for (a0,b0),(l,u) in zip(image,box)), 'strict Krawczyk contraction')
    return image, norm


def ball_range(lo, hi, Ball, scale):
    return Ball(lo.numerator*scale//lo.denominator,-((-hi.numerator*scale)//hi.denominator))


def sqrt_ball(v, Ball, scale):
    need(v.lo>=0,'nonnegative square-root input')
    return Ball(isqrt(v.lo*scale),isqrt(v.hi*scale)+1)


def coords(image, Ball, scale):
    a,b,c,r,x = [ball_range(l,u,Ball,scale) for l,u in image]
    A,B,C,p,q = [sqrt_ball(v,Ball,scale) for v in (a,b,c,(r+x)/2,(r-x)/2)]
    zero=Ball.integer(0)
    return [(A,zero),(-A,zero),(zero,B),(zero,-B),(C,zero),(-C,zero),
            (p,q),(-p,-q),(p,-q),(-p,q)]


def sub(z,w):
    return tuple(a-b for a,b in zip(z,w))


def add(z,w):
    return tuple(a+b for a,b in zip(z,w))


def mul(z,w):
    a,b=z
    c,d=w
    return a*c-b*d,a*d+b*c


def scale_complex(z,c):
    return tuple(c*a for a in z)


def norm2(z):
    return sum(a.square() for a in z)


def inv(z):
    a,b=z
    d=norm2(z)
    return a/d,-b/d


def dot(z,w):
    return sum(a*b for a,b in zip(z,w))


def inertia_signs(matrix):
    a=[row.copy() for row in matrix]
    signs=[0,0]
    for k in range(len(a)):
        pivot=a[k][k]
        need(pivot!=0,'nonzero symmetric elimination pivot')
        signs[int(pivot<0)]+=1
        for i in range(k+1,len(a)):
            for j in range(i,len(a)):
                a[j][i]=a[i][j]=a[i][j]-a[i][k]*a[j][k]/pivot
    return signs+[0]


def payload(matrix):
    return {'matrix':{'domain':'QQ','entries':[
        [{'num':str(x.numerator),'den':str(x.denominator)} for x in row] for row in matrix]}}


def check_sources(repo):
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for base, key, paths in [(HERE,'files',EVIDENCE),(repo,'dependencies',DEPENDENCIES)]:
        need(set(manifest[key])==set(paths),'all fixed proof paths')
        for name in paths:
            need(sha256((base/name).read_bytes()).hexdigest()==manifest[key][name],
                 'fixed proof input '+name)


def check_jacobian_inertia(matrices):
    saved = json.loads((HERE/'saddle-inertia-jacobian.json').read_bytes())
    need(set(saved)==set(matrices),'all four independent inertia comparisons')
    for key, inp in matrices.items():
        row = saved[key]
        out = row['output']
        want = {'remainder_laplacian_shift':(10,0,0), 'quadratic_term_shift':(20,0,0)}.get(key,(1,18,0))
        need(row['operation_id']=='matrix.inertia.compute' and
             out['matrix']['domain']==inp['matrix']['domain'] and
             out['matrix']['entries']==inp['matrix']['entries'] and
             tuple(out[q] for q in ('n_positive','n_negative','n_zero'))==want,
             'fixed Jacobian inertia result for '+key)


def directional_exclusion(z, Ball, scale):
    v = [(Ball.integer(a)/1000,Ball.integer(b)/1000) for a,b in DIRECTION]
    w = [(Ball.integer(a)/1000,Ball.integer(b)/1000) for a,b in FORCE_MULTIPLIER]
    G = -10*sum(map(norm2,v))-5*sum(dot(p,q) for p,q in zip(w,z))
    complex_gradient = [scale_complex((a,-b),-5) for a,b in w]
    hessian_lap=[[F(0) for _ in z] for _ in z]
    n=len(z)
    quadratic=[[Ball.integer(0) for _ in range(2*n)] for _ in range(2*n)]
    delta10=Ball.integer(1)
    pair_shift=Ball.integer(99)/70*Ball.integer(RADIUS.numerator)/RADIUS.denominator
    need(F(99,70)**2>2,'strict square-root-two upper bound')
    for i,p in enumerate(z):
        for j,q in enumerate(z[:i]):
            d=sub(p,q)
            k=sub(v[i],v[j])
            l=sub(w[i],w[j])
            squared_d=norm2(d)
            delta10=delta10*(squared_d*10/9)
            square_k=mul(k,k)
            reciprocal=inv(d)
            reciprocal2=mul(reciprocal,reciprocal)
            reciprocal3=mul(reciprocal2,reciprocal)
            reciprocal4=mul(reciprocal2,reciprocal2)
            value=add(scale_complex(mul(square_k,reciprocal2),-2),mul(l,reciprocal))
            G+=value[0]
            derivative=sub(scale_complex(mul(square_k,reciprocal3),4),mul(l,reciprocal2))
            complex_gradient[i]=add(complex_gradient[i],derivative)
            complex_gradient[j]=sub(complex_gradient[j],derivative)
            a,b=mul(sub(mul(l,d),scale_complex(square_k,6)),reciprocal4)
            qmat=((a,-b),(-b,-a))
            for x in range(2):
                for y in range(2):
                    qxy=qmat[x][y]
                    quadratic[2*i+x][2*i+y]+=qxy
                    quadratic[2*j+x][2*j+y]+=qxy
                    quadratic[2*i+x][2*j+y]-=qxy
                    quadratic[2*j+x][2*i+y]-=qxy
            numerator=sqrt_ball(norm2(sub(scale_complex(square_k,24),scale_complex(mul(l,d),3))),Ball,scale)
            numerator+=3*sqrt_ball(norm2(l),Ball,scale)*pair_shift
            lower_distance=sqrt_ball(squared_d,Ball,scale)-pair_shift
            need(lower_distance.lo>0,'collision-free exclusion ball')
            capacity=numerator/(lower_distance.square().square()*lower_distance)
            capacity=F(-((-capacity.hi*10**6)//scale),10**6)
            hessian_lap[i][i]+=capacity
            hessian_lap[j][j]+=capacity
            hessian_lap[i][j]-=capacity
            hessian_lap[j][i]-=capacity
    gradient_norm2=sum(map(norm2,complex_gradient))
    lower=F(1697,10**4)
    upper_gradient=F(13,2000)
    need(F(G.lo,scale)>lower,'corrected witness strictly at saddle')
    need(F(gradient_norm2.hi,scale)<upper_gradient**2,'small corrected witness gradient')
    quadratic_center, error = rational_center(quadratic,scale)
    shifted=[[a+(74-error)*int(i==j) for j,a in enumerate(row)]
             for i,row in enumerate(quadratic_center)]
    need(inertia_signs(shifted)==[20,0,0],'signed quadratic matrix strictly exceeds -74I')
    matrix=[[970*int(i==j)-value for j,value in enumerate(row)] for i,row in enumerate(hessian_lap)]
    need(inertia_signs(matrix)==[10,0,0],'dominant third derivative Laplacian norm < 970')
    comparison=lower-upper_gradient*RADIUS-74*RADIUS**2-F(970*99,3*70)*RADIUS**3
    need(comparison>F(1,250),'strict positive second variation at all stationary points in ball')
    lo,hi=F(delta10.lo,scale),F(delta10.hi,scale)
    need(53744000000<lo<=hi<53746000000,'saddle product is below established witness')
    return {
        'mass_nine_excluded_orbit_radius':str(RADIUS),
        'corrected_saddle_witness_lower':str(lower),
        'corrected_witness_gradient_upper':str(upper_gradient),
        'quadratic_term_operator_lower':-74,
        'third_derivative_laplacian_norm_upper':970,
        'positive_second_variation_lower':str(comparison),
        'saddle_delta_mass_ten_coarse_interval':[lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator)],
    }, {'remainder_laplacian_shift': payload(matrix), 'quadratic_term_shift':payload(shifted)}


def rational_center(jac,scale):
    denom=1<<30
    midpoint=[[F((v.lo+v.hi)*denom//(2*scale),denom) for v in row] for row in jac]
    err=F(1,1000000)
    for row,mid in zip(jac,midpoint):
        row_error=sum(max(abs(F(v.lo,scale)-a),abs(F(v.hi,scale)-a)) for v,a in zip(row,mid))
        need(row_error<err, 'rational midpoint matrix operator norm error')
    return midpoint,err


def sliced_inertia(z, repo, dyadic):
    sys.path.insert(0,str(repo/'validated-maxima'))
    force = load('n10_force_equations',repo/'validated-maxima/verify_all.py')
    Ball=dyadic.Ball
    scale=dyadic.SCALE
    factor=sqrt_ball(Ball.integer(10)/9,Ball,scale)
    order=[4,5,0,1,2,3,6,7,8,9]
    flat=[v*factor for i in order for v in z[i]]
    _,jac,_=force.equations(flat[:1]+flat[2:],10)
    lower_precision,err=rational_center(jac,scale)
    matrices={}
    for sign in (-1,1):
        shifted=[[v+sign*err*int(i==j) for j,v in enumerate(row)] for i,row in enumerate(lower_precision)]
        need(inertia_signs(shifted)==[1,18,0], 'same inertia for both shifts')
        matrices[f'force_jacobian_shift_{sign}']=payload(shifted)
    return matrices


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=ROOT)
    parser.add_argument('--jacobian-inputs',type=Path)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    repo=args.repo.resolve()
    check_sources(repo)
    h = load('original_root_replay_helpers',repo/'ten-point-contact-geometry/check_root_box.py')
    arithmetic=load('n10_intervals',repo/'validated-maxima/dyadic.py')
    # Use a common Ball class in the imported full-force verifier as well.
    sys.modules['dyadic']=arithmetic
    data=json.loads((HERE/'saddle-jacobian.json').read_bytes())
    image,_=root_image(data,h)
    z=coords(image,arithmetic.Ball,arithmetic.SCALE)
    report, matrices=directional_exclusion(z,arithmetic.Ball,arithmetic.SCALE)
    matrices.update(sliced_inertia(z,repo,arithmetic))
    check_jacobian_inertia(matrices)
    if args.jacobian_inputs is not None:
        args.jacobian_inputs.write_text(json.dumps(matrices,indent=2)+'\n',encoding='utf-8')
    report['reduced_force_inertia']=[1,18,0]
    report['unique_nonsingular_five_variable_root']=True
    report['contraction_norm_strict_upper']='1/1000'
    report['jacobian_inertia_crosschecks']=4
    report['passed']=True
    report_path=HERE/'verification_report.json'
    if args.record_report:
        report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    else:
        need(report==json.loads(report_path.read_bytes()),'frozen replay report')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
