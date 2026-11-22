"""Angular phase cells at fixed normalized outer radii with free inner separation."""
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path

S=F(35,4)
EPS={'candidate':F(1,20000), 'saddle':F(3,4000)}

def load(repo):
    path=repo/'ten-point-phase-sweep/verify_all.py'
    spec=importlib.util.spec_from_file_location('phase_base',path)
    phase=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(phase)
    return phase

class Tube:
    def __init__(self,repo):
        self.phase=load(repo)
        self.eng=self.phase.AngularBox()
        self.ar=self.eng.saddle

    def data(self,mode):
        e=self.eng; a=self.ar
        points=e.points[mode]
        original=[a.norm2(z) for z in points[2:]]
        factor=e.exact(S)/sum(original,e.zero)
        # The anchor has zero inner midpoint, whence its eight normalized radii sum to S.
        radii=[u*factor for u in original]
        dirs=[a.scale_complex(z,1/e.sqrt(u)) for z,u in zip(points[2:],original)]
        return points,radii,dirs,EPS[mode]

    def accept(self, mode, data, lo, hi, refs):
        e=self.eng; a=self.ar
        zs,radii,dirs,eps=data
        w,*slopes=refs
        self.phase.require(F(107,1000)<=w<=F(133,1000),'reference inner half-difference')
        u=[a.mul(x,y) for x,y in zip(dirs,e.shifts(slopes))]
        for (x,y),(v,b) in zip(u,u[1:]+u[:1]):
            dot=x*v+y*b
            self.phase.require((x*b-y*v).lo>0 and e.endpoint(dot,False)>F(59,100)
                               and e.endpoint(dot)<F(39,50),'admissible reference gaps')
        half=e.sqrt(e.exact(w)); m=e.interval(-eps,eps)
        obs=[e.sqrt(e.expand(r,eps)) for r in radii]
        middle=(lo+hi)/2
        ref=e.corner([(half+m,m),(m-half,m)],obs,u,e.interval(lo,hi))
        g=e.angle_gradient(ref,u,obs,(m,m),half)
        if g is None:return None

        def maxabs(x):return max(-e.endpoint(x,False),e.endpoint(x))
        def denom(point):
            p=a.sub(point,(m,m))
            return p,a.sub(a.mul(p,p),(e.exact(w),e.zero))

        midpoint_ref=e.corner(ref[:2],obs,u,e.exact(middle))
        gw=e.exact(1/w-90/(S+2*w))
        total_cross=e.zero
        for z,q in zip(midpoint_ref[2:],ref[2:]):
            _,dc=denom(z)
            p,d=denom(q)
            if min(a.norm2(dc).lo,a.norm2(d).lo)<=0:return None
            gw-=2*a.inv(dc)[0]
            rec=a.inv(d)
            total_cross-=4*a.mul(a.mul(q,p),a.mul(rec,rec))[1]
        gw_bound=maxabs(gw) + maxabs(total_cross)*(hi-lo)
        g2=gw_bound**2/5+2*g[0]
        massfactor=e.sqrt(e.exact(9/(S+2*w)))
        physical=[a.scale_complex(z,massfactor) for z in ref]
        if mode=='candidate':
            pairs=((0,1),(2,6),(3,7),(4,8),(5,9))
            observed=[a.scale_complex(a.sub(physical[i],physical[j]),e.exact(F(1,2))) for i,j in pairs]
            dist=e.endpoint(e.norm([a.sub(x,zs[i]) for x,(i,j) in zip(observed,pairs)]))
            if dist+F(33,25)*e.endpoint(e.sqrt(e.exact(g2)))<=F(7,250):
                return 'candidate_basin'
        base=[r*massfactor.square() for r in radii]
        exact_half=half*massfactor
        anchor=e.corner([(exact_half,e.zero),(-exact_half,e.zero)],[e.sqrt(r) for r in base],u,e.exact(middle))
        ds=e.distances(anchor)
        if any(d.lo<=0 for d in ds):return None
        qa=e.endpoint(e.product(ds))
        bump=e.energy_increase(anchor,ds,[r*massfactor for r in obs],base,
                            m*massfactor,exact_half,hi-lo)+g2/2
        if bump<2 and qa*self.phase.exp_upper(bump)<e.lower:return 'energy_excluded'
        return None
