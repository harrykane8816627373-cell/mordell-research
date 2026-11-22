"""Isolate the radial stationary profile in two angle-only neighborhoods."""

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SETTINGS = {
    'candidate': (F(11, 500), [4388, 619, 1049, 1049, 619]),
    'saddle': (F(3, 125), [3142, 427, 984, 1750, 984]),
}


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def run(repo, mode, delta, widths):
    spec = importlib.util.spec_from_file_location('phase',repo/'ten-point-phase-sweep/verify_all.py')
    phase = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(phase)
    box = phase.AngularBox()
    zs = box.points[mode]
    reps = [zs[i] for i in (0,2,3,4,5)]
    r = [box.saddle.norm2(z) for z in reps]
    u = [box.saddle.scale_complex(box.saddle.mul(z,z),1/rr) for z,rr in zip(reps,r)]

    def weights(shifts, angular_errors):
        radius=[]
        for rr,q in zip(r,shifts):
            hi=phase.exp_upper(q)
            radius.append(rr*box.interval(1/hi,hi))
        hessian={}
        mixed={}
        for i,j in combinations(range(5),2):
            e=angular_errors[i]+angular_errors[j]
            cosine = box.interval(1-e**2/2,1)
            sine = box.interval(-e,e)
            relative = box.saddle.mul(u[i],(u[j][0],-u[j][1]))
            c,s=box.saddle.mul(relative,(cosine,sine))
            ratio=radius[i]/radius[j]
            u_lo,u_hi=box.endpoint(ratio,False),box.endpoint(ratio)
            d_lo,d_hi=(u_lo+1/u_lo)/2,(u_hi+1/u_hi)/2
            d=box.interval(1 if u_lo<=1<=u_hi else min(d_lo,d_hi),max(d_lo,d_hi))
            sinh=box.interval((u_lo-1/u_lo)/2,(u_hi-1/u_hi)/2)
            den=(d-c).square()
            require(den.lo>0,'positive square distance')
            hessian[i,j]=2*(1-c*d)/den
            mixed[i,j]=-2*sinh*s/den
        return radius,hessian,mixed

    def apply(precond, rr, edges):
        out=[[box.zero for j in range(5)] for i in range(5)]
        for i in range(5):
            for (j,k),h in edges.items():
                excess=(precond[i][j]-precond[i][k])*h
                out[i][j]+=excess
                out[i][k]-=excess
            if rr:
                for j in range(5):
                    out[i][j]-=10*precond[i][j]*rr[j]
        return out

    zeros = [F(0)]*5
    eye=box.basin.eye(5)
    _,h0,_=weights(zeros,zeros)
    negH=box.basin.subtract([[box.zero]*5 for _ in range(5)],apply(eye,r,h0))
    P=box.basin.inverse(negH)
    # The exact inverse at the true algebraic anchor is enclosed by P.
    # It is nonsingular by positive definiteness.
    angle=[F(0)]+[delta]*4
    _,_,at_anchor=weights(zeros,angle)
    cr=apply(P,None,at_anchor)
    rr,hall,_=weights(widths,angle)
    H=apply(P,rr,hall)
    def mag(b):
        return max(abs(box.endpoint(b,False)),abs(box.endpoint(b)))
    err=[sum(mag(cr[i][j])*angle[j]+mag(eye[i][j]+H[i][j])*widths[j] for j in range(5)) for i in range(5)]
    bound=sum(box.endpoint(rad)*( (phase.exp_upper(w/2)-1)**2 + phase.exp_upper(w/2)*e**2/4)
              for rad,w,e in zip(r,widths,angle))
    bands = [(F(43,400),F(53,400))]+[(F(89,100),F(1271,1000))]*4
    for radial, (lo, hi) in zip(rr, bands):
        require(lo < box.endpoint(radial,False) <= box.endpoint(radial) < hi,
                mode+' profile box lies in radial bounds')
    outer_product = box.product(rr[1:])
    require(box.endpoint(outer_product,False) > F(43,40)**4, mode+' profile box lies in log mass domain')
    require(all(x < w for x,w in zip(err,widths)), mode+' radial fixed-point inclusion')
    distance_sq = bound if mode=='candidate' else 2*bound
    limit = F(7,250) if mode=='candidate' else F(21,500)
    require(distance_sq < limit**2, mode+' local sharp comparison')
    return {'squared_image_angle_tolerance':str(delta), 'log_radius_box_half_widths':list(map(str,widths)),
            'induced_norm_upper_millionths':-((-box.sqrt(box.exact(distance_sq)).hi*10**6)//box.scale),
            'local_comparison_radius':str(limit)}


def main():
    manifest = json.loads((HERE/'source_manifest.json').read_bytes())
    for directory, items in ((HERE, manifest['files']), (ROOT, manifest['dependencies'])):
        for name, digest in items.items():
            require(sha256((directory/name).read_bytes()).hexdigest() == digest, 'hash '+name)
    proof = json.loads((ROOT/'ten-point-half-turn-radial-concavity/verification_report.json').read_bytes())
    require(proof['passed'] and proof['radial_stationarity_implies_unique_radii_at_fixed_angles'] and
            proof['cover']['unresolved']==0, 'unique radial stationarity in global domain')
    results = {mode:run(ROOT,mode,delta,[F(x,100000) for x in widths])
               for mode,(delta,widths) in SETTINGS.items()}
    print(json.dumps({'passed':True, 'normalization_mass':9, 'assumption':'half-turn symmetry',
                      'angle_neighborhoods':results},indent=2))


if __name__=='__main__':
    main()
