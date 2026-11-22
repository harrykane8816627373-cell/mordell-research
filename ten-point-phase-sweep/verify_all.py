"""Replay angular elimination for two ten-point radial neighborhoods."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STIFFNESS = F(161,50)
AMPLITUDE = F(7,250)
MODES = {'candidate':F(3,20000), 'saddle':F(1,500)}
PAIRS = tuple(combinations(range(10),2))


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def read(path):
    return json.loads(path.read_bytes())


def load(folder, name='verify_all.py'):
    spec = importlib.util.spec_from_file_location(folder,ROOT/folder/name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exp_upper(q):
    require(0<=q<8, 'geometric exponential series remainder')
    return sum(q**k/factorial(k) for k in range(7)) + q**7/factorial(7)/(1-q/8)


class AngularBox:
    def __init__(self):
        self.basin=load('ten-point-weighted-pair-basin')
        self.saddle=load('ten-point-saddle-exclusion')
        self.Ball=self.basin.Ball
        self.scale=self.basin.SCALE
        self.exact=self.basin.exact
        self.zero=self.basin.ZERO
        self.sqrt=self.basin.sqrt
        rep=self.basin.algebraic_data()
        neg=lambda z: self.saddle.scale_complex(z,-1)
        a,bp,bm,cp,cm=rep
        self.candidate=[a,neg(a),cp,bp,neg(bm),neg(cm),neg(cp),neg(bp),bm,cm]
        root=read(ROOT/'ten-point-saddle-exclusion/saddle-jacobian.json')
        h=load('ten-point-contact-geometry','check_root_box.py')
        image,_=self.saddle.root_image(root,h)
        zs=self.saddle.coords(image,self.Ball,self.scale)
        self.saddle_points=[zs[i] for i in (0,1,4,6,2,9,5,7,3,8)]
        self.points={'candidate':self.candidate, 'saddle':self.saddle_points}
        self.lower=self.endpoint(self.product(self.distances(self.candidate)), False)

    def endpoint(self, a, upper=True):
        return F(a.hi if upper else a.lo,self.scale)

    def interval(self,l,r):
        return self.Ball(self.exact(l).lo,self.exact(r).hi)

    def expand(self,v,eps):
        return self.interval(self.endpoint(v,False)-eps,self.endpoint(v)+eps)

    def rotation(self,t):
        t2=t.square()
        return (1-t2)/(1+t2), (2*t)/(1+t2)

    def distance(self,z,w):
        return self.saddle.norm2(self.saddle.sub(z,w))

    def norm(self, zs):
        return self.sqrt(sum((self.saddle.norm2(z) for z in zs),self.zero))

    def distances(self,zs):
        return [self.distance(zs[i],zs[j]) for i,j in PAIRS]

    def product(self, entries):
        result=self.exact(1)
        for entry in entries:
            result*=entry
        return result

    def shifts(self, slopes):
        require(len(slopes)==3 and max(map(abs,slopes))<=F(1,12), 'small rotation corrections')
        rs=[self.rotation(self.exact(t)) for t in slopes]
        total=self.saddle.mul(self.saddle.mul(rs[0],rs[1]),rs[2])
        rs.append((total[0],-total[1]))
        return rs*2

    def setup_mode(self,mode):
        zs=self.points[mode]
        radii=[self.saddle.norm2(z) for z in zs[2:]]
        directions=[self.saddle.scale_complex(z,1/self.sqrt(u)) for z,u in zip(zs[2:],radii)]
        inner=self.saddle.norm2(zs[0]); eps=MODES[mode]
        return zs, radii, directions, self.sqrt(self.expand(inner,eps)), self.interval(-eps,eps)

    def corner(self, zs, radii, directions, phase):
        rot=self.rotation(phase)
        return [*zs[:2],*(self.saddle.mul(rot,self.saddle.scale_complex(u,r)) for r,u in zip(radii,directions))]

    def angle_gradient(self, zs, directions, radii, midpoint, half):
        outer=[self.zero]*8
        for i,j in combinations(range(8),2):
            a,b=directions[i]; c,e=directions[j]
            rho,sigma=radii[i],radii[j]
            q=1+(rho-sigma).square()/(2*rho*sigma)
            value=(b*c-a*e)/(q-a*c-b*e)
            outer[i]+=value
            outer[j]-=value
        inner=[]
        for point in zs[2:]:
            offset=self.saddle.sub(point,midpoint)
            w=(half.square(),self.zero)
            denom=self.saddle.sub(self.saddle.mul(offset,offset),w)
            if self.saddle.norm2(denom).lo<=0:
                return None
            num=self.saddle.add(w,self.saddle.mul(midpoint,offset))
            inner.append(-4*self.saddle.mul(num,self.saddle.inv(denom))[1])
        projected=[v+(7*x-sum((y for j,y in enumerate(inner) if i!=j),self.zero))/8
                   for i,(v,x) in enumerate(zip(outer,inner))]
        squared=min(self.endpoint(sum(((v+x).square() for v,x in zip(outer,inner)),self.zero)),
                    self.endpoint(sum((x.square() for x in projected),self.zero)))
        return squared,self.endpoint(self.sqrt(self.exact(squared)))

    def force(self, zs, ds):
        g=[self.saddle.scale_complex(z,-5) for z in zs]
        for (i,j),d in zip(PAIRS,ds):
            value=self.saddle.scale_complex(self.saddle.sub(zs[i],zs[j]),1/d)
            g[i]=self.saddle.add(g[i],value)
            g[j]=self.saddle.sub(g[j],value)
        return g

    def energy_increase(self, anchor, distances, rbox, radii, midpoint, half, phase_error):
        force=self.force(anchor,distances)
        phase_linear=(self.zero,self.zero)
        change=2*(midpoint.square()+midpoint.square()+(half-anchor[0][0]).square())
        inner_error=max(abs(self.endpoint(half-anchor[0][0],False)),abs(self.endpoint(half-anchor[0][0])))
        abs_upper=lambda x: max(abs(self.endpoint(x,False)),abs(self.endpoint(x)))
        linear=4*abs_upper(force[0][0])*inner_error
        for point,(a,b),observed,ideal in zip(anchor[2:],force[2:],rbox,map(self.sqrt,radii)):
            product=self.saddle.mul((a,-b),point)
            phase_linear=self.saddle.add(phase_linear,product)
            delta=abs_upper(observed-ideal)
            radial_slope=abs_upper(product[0]/ideal)
            force_size=self.endpoint(self.norm([(a,b)]))
            linear+=2*delta*(radial_slope+phase_error*force_size)
            change+=(observed-ideal).square()+observed*ideal*self.exact(phase_error**2)
        linear+=2*abs_upper(phase_linear[1])*phase_error+abs_upper(phase_linear[0])*phase_error**2
        row_sums=[F(0)]*10
        for (i,j),d in zip(PAIRS,distances):
            dl=self.endpoint(d,False)
            require(dl>0,'separated energy anchor')
            row_sums[i]+=1/dl
            row_sums[j]+=1/dl
        return linear+max(F(0),2*max(row_sums)-5)*self.endpoint(change)

    def accept(self, mode, supplied, low, high, slopes):
        zs,radii,directions,half,m= supplied
        eps=MODES[mode]
        u=[self.saddle.mul(a,b) for a,b in zip(directions,self.shifts(slopes))]
        for (a,b),(c,d) in zip(u,u[1:]+u[:1]):
            dot=a*c+b*d
            require((a*d-b*c).lo>0 and self.endpoint(dot,False)>F(59,100)
                    and self.endpoint(dot)<F(39,50), 'reference has all eight admissible consecutive gaps')
        inner=[(m+half,m),(m-half,m)]
        rbox=[self.sqrt(self.expand(r,eps)) for r in radii]
        ref=self.corner(inner,rbox,u,self.interval(low,high))
        ds=self.distances(ref)
        if any(d.lo<=0 for d in ds):
            return None
        residual=self.angle_gradient(ref,u,rbox,(m,m),half)
        if residual is None:
            return None
        g2,g_norm=residual

        if mode=='candidate':
            half_differences=[self.saddle.scale_complex(self.saddle.sub(ref[i],ref[j]),self.exact(F(1,2)))
                              for i,j in ((0,1),(2,6),(3,7),(4,8),(5,9))]
            ideal=[zs[i] for i in (0,2,3,4,5)]
            error=self.endpoint(self.norm([self.saddle.sub(z,w) for z,w in zip(half_differences,ideal)]))
            if error+F(4,5)*g_norm/STIFFNESS<=AMPLITUDE:
                return 'candidate_basin'

        middle=(low+high)/2
        anchor=self.corner(zs,[self.sqrt(r) for r in radii],u,self.exact(middle))
        ds0=self.distances(anchor)
        aq=self.endpoint(self.product(ds0))
        budget=self.energy_increase(anchor,ds0,rbox,radii,m,half,high-low)+g2/(2*STIFFNESS)
        if budget<2 and aq*exp_upper(budget)<self.lower:
            return 'energy_excluded'
        return None


def verify_inputs():
    manifest=read(HERE/'source_manifest.json')
    for base,group in ((HERE,'files'),(ROOT,'dependencies')):
        for name,expected in manifest[group].items():
            require(sha256((base/name).read_bytes()).hexdigest()==expected,'fixed input '+name)
    gaps=read(ROOT/'ten-point-outer-gap/verification_report.json')
    inner=read(ROOT/'ten-point-inner-midpoint/verification_report.json')
    radii=read(ROOT/'ten-point-angular-radii/verification_report.json')
    quartic=read(ROOT/'ten-point-quartic-midpoints/verification_report.json')
    saddle=read(ROOT/'ten-point-saddle-exclusion/verification_report.json')
    require(all(row['passed'] and row.get('normalization_mass',9)==9 for row in (gaps,inner,radii,quartic,saddle)),
            'accepted angular domain, sharp basin and saddle')
    require(F(quartic['half_difference_radius'])==AMPLITUDE
            and F(inner['outer_mean_zero_negative_angle_hessian_lower_on_successive_gap_polytope'])==STIFFNESS
            and F(radii['outer_squared_radius_open_interval'][1])/2<F(4,5)**2
            and tuple(map(F,gaps['outer_successive_angle_cosine_open_interval']))==(F(59,100),F(39,50)),
            'fixed angular residual and amplitude bounds')


def replay(engine, mode, rows):
    data=engine.setup_mode(mode)
    counts={'candidate_basin':0,'energy_excluded':0}
    cursor=F(-1)
    maximum=F(0)
    for row in rows:
        lo,hi=map(F,row['interval'])
        require(lo==cursor and lo<hi<=1, 'continuous phase coverage')
        outcome=engine.accept(mode,data,lo,hi,list(map(F,row['corrections'])))
        require(outcome==row['result'] and outcome in counts, f'accepted phase cell {mode} {lo} {hi}')
        counts[outcome]+=1
        maximum=max(maximum,hi-lo)
        cursor=hi
    require(cursor==1,'phase chart fully covered')
    return {'cell_counts':counts,'maximum_cell_width':str(maximum),'shape_parameter_tolerance':str(MODES[mode])}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    verify_inputs()
    engine=AngularBox()
    rows={mode:replay(engine,mode,read(HERE/f'{mode}.json')) for mode in MODES}
    result={'passed':True,'normalization_mass':9,'phase_tangent_interval':['-1','1'],
            'angular_hessian_gap':str(STIFFNESS),
            'phase_cells':sum(sum(record['cell_counts'].values()) for record in rows.values()),
            'modes':rows,'unresolved_within_these_radial_neighborhoods':0}
    text=json.dumps(result,indent=2)+'\n'
    output=HERE/'verification_report.json'
    if args.record_report:
        output.write_text(text,encoding='ascii',newline='\n')
    else:
        require(result==read(output),'frozen report')
    print(text,end='')


if __name__=='__main__':
    main()
