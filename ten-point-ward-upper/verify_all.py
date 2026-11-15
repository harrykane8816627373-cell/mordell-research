"""Exact replay of the ten-point PSD, Ward and radial-force upper bound."""
import argparse
from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
from math import comb, lcm, prod
from pathlib import Path
import hashlib, importlib.util, json, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
U=F(15403,10000)
U_NEW=F(1539,1000)
DMIN=F(2,5)
L=F(97,80)*9**9
N=45
C=F(3,10000)
RATIO=F(507,500)
LOG_MARGIN=F(1,3000)

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

old=load('twelve-point-kernel-upper/verify_cover.py','n10_exact_cover')
gram=load('candidate-kernel/verify_psd.py','n10_exact_gram')
base=load('candidate-kernel/verify_all.py','n10_exact_base')
log=load('twelve-point-kernel-upper/verify_all.py','n10_exact_log')
need=old.need
hx=old.hx
S,BIAS,W=old.S,old.BIAS,old.W

def corrections(data):
    w={k:hx(v) for k,v in data['ward_parameters'].items()}
    s=data['ward_feature_order']
    A=[[[F(0)]*s for i in range(s)] for k in range(4)]
    for i in (0,1):
        for j in range(i,s):A[1][i][j]=A[1][j][i]=w.pop(f'A1_{i}{j}')
    for k,key,r in ((0,'mass',F(1,10)),(2,'a2',-F(7,10)),(3,'a3',-F(3,5))):
        v=[r,F(1)]+[F(0)]*(s-2);a=[w.pop(f'{key}_{i}') for i in range(s)]
        A[k]=[[v[i]*a[j]+a[i]*v[j] for j in range(s)] for i in range(s)]
    need(not w,'exact permitted Ward correction parameters only')
    unit=[[F(i==j) for j in range(s)] for i in range(s)]
    for k,r in enumerate((F(1,10),None,-F(7,10),-F(3,5))):
        vectors=unit[2:] if r is None else [[F(1),-r]+[F(0)]*(s-2)]+unit[2:]
        need(all(base.bilinear(x,A[k],y)==0 for x in vectors for y in vectors),'exact vanishing Ward restriction')
    return A, list(map(hx,data['force_parameters']))

def polynomial(data,A,forces):
    C,shape=old.polynomial(data)
    s=data['ward_feature_order']
    T=[[comb(j,i)*(-1)**(j-i) if j>=i else 0 for j in range(s)] for i in range(s)]
    H=[{(0,0,0):1}, {(0,0,1):1}, {(0,0,2):2,(1,1,0):-1}, {(0,0,3):4,(1,1,1):-3}]
    def put(D,i,j,k,value):D[(i,j,k)]=D.get((i,j,k),F(0))+value
    for k in range(4):
        for i in range(s):
            for j in range(s):
                G=sum(T[i][r]*A[k][r][q]*T[j][q] for r in range(s) for q in range(s))
                for (a,b,c),f in H[k].items():put(C,i+a,j+b,c,-f*G)
                put(C,i+j+k,0,0,-G/18);put(C,0,i+j+k,0,-G/18)
    if forces:
        Z={}
        for (i,j,k),value in C.items():
            put(Z,i+1,j,k,value);put(Z,i,j+1,k,value);put(Z,i,j,k+1,-2*value)
        for j,w in enumerate(forces,1):
            for i in range(j+1):
                p=comb(j,i)*(-1)**(j-i)*w
                put(Z,i,0,1,-p);put(Z,0,i,1,-p)
                put(Z,i+1,0,0,p);put(Z,0,i+1,0,p)
                for a,b in ((i+1,0),(0,i+1)):
                    put(Z,a+1,b,0,-5*p/9);put(Z,a,b+1,0,-5*p/9);put(Z,a,b,1,10*p/9)
        C=Z;shape=tuple(v+1 for v in shape)
    C={index:value for index,value in C.items() if value}
    need(all(all(0<=j<dim for j,dim in zip(index,shape)) for index in C),'powers fit exact tensor')
    need(all(C.get((j,i,k),0)==p for (i,j,k),p in C.items()),'symmetric pair polynomial')
    return C,shape

def crosscheck(data,A,forces,C):
    s=data['ward_feature_order']
    def k(u,v,t):
        phi=lambda r:[(r-1)**j for j in range(s)]
        harmonics=(F(1),t,2*t*t-u*v,4*t**3-3*u*v*t)
        return base.encoded_kernel(u,v,t,data)+sum(h*base.bilinear(phi(u),Ak,phi(v)) for h,Ak in zip(harmonics,A))
    b=(hx(data['B'])-45)/10
    for u,v,t in ((F(0),F(1),F(0)),(F(1,3),F(7,8),-F(1,4)),(U,U,-U),(F(0),U,F(0)),(F(1),F(7,6),F(1,4))):
        d=u+v-2*t
        value=(5*(u+v)+2*b-(k(u,u,u)+k(v,v,v))/2)/9-k(u,v,t)
        if forces:
            value*=d
            for j,w in enumerate(forces,1):
                value+=w*((u-t-5*u*d/9)*(u-1)**j+(v-t-5*v*d/9)*(v-1)**j)
        need(value==sum(p*u**i*v**j*t**k for (i,j,k),p in C.items()),'independent exact pair expression matches polynomial')

def initial(C,shape,box):
    den=lcm(*(v.denominator for v in C.values()));size=prod(shape);A=[0]*size
    for (i,j,k),v in C.items():A[(i*shape[1]+j)*shape[2]+k]=v.numerator*(den//v.denominator)
    for axis in range(3):
        M,d=old.transform(shape[axis]-1,*box[axis]);stride=prod(shape[axis+1:]);group=shape[axis]*stride;out=[0]*size
        for base in range(0,size,group):
            for offset in range(stride):
                vals=[A[base+j*stride+offset] for j in range(shape[axis])]
                for i,row in enumerate(M):out[base+i*stride+offset]=sum(x*y for x,y in zip(row,vals))
        A=out;den*=d
    values=[S*a//den+BIAS for a in A]
    need(0<=min(values)<=max(values)<1<<62,'initial packed lane range')
    return sum(v<<(W*i) for i,v in enumerate(values))

def prune(box):
    (ul,uh),(vl,vh),(tl,th)=box
    amin=max(tl,-th,0);dhi=uh+vh-2*tl
    if ul>vh or ul>U_NEW or vl>U_NEW or amin**2>uh*vh or dhi<DMIN:return True
    q=9-F(9,8)*(ul+vl)-tl/4
    if q<F(142,25):return True
    delta=5*q-28
    if delta>=2:return False
    b=max(F(0),81*uh+vh+18*th)*max(F(0),uh+81*vh+18*th)/4096
    return (delta+2)**2*dhi*b**8*q**28<(2-delta)**2*L*7**21

class Packed(old.Packed):
    def __init__(self,shape,force):
        super().__init__(shape);self.force=force
        self.index=[sum(i*m for i,m in enumerate(masks)) for masks in self.axisones]
        if force:
            self.choose2=[sum(i*(i-1)*m for i,m in enumerate(masks)) for masks in self.axisones]
            self.cross={(i,j):sum(x*y*(a&b) for x,a in enumerate(self.axisones[i]) for y,b in enumerate(self.axisones[j]))
                        for i in range(3) for j in range(i+1,3)}

    def dominates(self,B,box,c):
        (ul,uh),(vl,vh),(tl,th)=box
        a=(max(DMIN,ul+vl-2*th)+uh+vh-2*tl)/2
        need(a>0,'positive log tangent point')
        d0=ul+vl-2*tl
        ws=(uh-ul,vh-vl,-2*(th-tl));ceil=old.ceilq
        Lhat=F(old.log_upper_scaled(a),S)-1-c
        if self.force:
            # We only require the multiplied inequality at physical d>0.
            # Hence replacing log(a) by its upper endpoint is safe there
            # even when the ambient cuboid contains nonphysical d<0.
            c0=ceil(S*(d0*d0/a+Lhat*d0))
        else:
            c0=ceil(S*(d0/a+Lhat))
        low=high=BIAS+c0;T=(BIAS+c0)*self.ones
        def add(slope,masks,maxval):
            nonlocal low,high,T
            v=ceil(slope);T+=v*masks
            low+=min(0,v*maxval);high+=max(0,v*maxval)
        for axis,w in enumerate(ws):
            n=self.shape[axis]-1
            if self.force:
                add(S*w*(2*d0/a+Lhat)/n,self.index[axis],n)
                add(S*w*w/(a*n*(n-1)),self.choose2[axis],n*(n-1))
            else:add(S*F(w,a*n),self.index[axis],n)
        if self.force:
            for (i,j),masks in self.cross.items():
                ni,nj=self.shape[i]-1,self.shape[j]-1
                slope=2*S*ws[i]*ws[j]/(a*ni*nj)
                add(slope,masks,ni*nj)
        need(0<=low<=high<1<<62,'exact upper-threshold lane range')
        return ((B+self.high-T)&self.high)==self.high

def witness():
    data=json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())
    points=next(r['integer_points'] for r in data['rows'] if r['n']==10)
    sx=sum(x for x,y in points);sy=sum(y for x,y in points)
    D=10*sum(x*x+y*y for x,y in points)-sx*sx-sy*sy
    P=prod((x-u)**2+(y-v)**2 for i,(x,y) in enumerate(points) for u,v in points[:i])
    V=F(90**N*P,D**N);need(V>L,'exact feasible witness above threshold')
    return V

def actual_root_lower():
    dyadic=load('validated-maxima/dyadic.py','n10_root_dyadic')
    prev=sys.modules.get('dyadic');sys.modules['dyadic']=dyadic
    try:prior=load('validated-maxima/verify_all.py','n10_actual_root')
    finally:
        if prev is None:del sys.modules['dyadic']
        else:sys.modules['dyadic']=prev
    common=json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    row=next(r for r in common['rows'] if r['n']==10)
    accepted=prior.verify_row(row,common)
    need(accepted['classification']=='strict_local_maximum_mod_rotation','independently replayed candidate root')
    return F(accepted['delta_interval'][0])*F(9,10)**N

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--record-report',action='store_true');p.add_argument('--progress',action='store_true')
    args=p.parse_args()
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for label,directory in (('files',HERE),('dependencies',ROOT)):
        for path,sha in manifest[label].items():
            need(hashlib.sha256((directory/path).read_bytes()).hexdigest()==sha,f'fixed {label} source {path}')
    domain=load('maximizer-domain/verify_all.py','n10_replay_necessary_domain')
    captured=StringIO()
    with redirect_stdout(captured):domain.main()
    dr=json.loads(captured.getvalue())
    need(dr['passed'] is True and dr['symmetry_assumed'] is False
         and F(dr['point_squared_radius_strict_upper'])==U
         and F(dr['pair_squared_distance_strict_lower'])==DMIN
         and F(dr['eight_point_complement_centered_mass_strict_lower'])==F(142,25)
         and domain.L==L,'accepted necessary conditions for all global maximizers')
    collapse=load('ten-point-collapse-domain/verify_all.py','n10_replay_collapse_domain')
    cr=collapse.verify()
    need(cr['passed'] is True and F(cr['radius_strict_upper'])==U_NEW
         and collapse.L==L,'accepted strict improved radius for all global maximizers')
    raw=(HERE/'kernel.json').read_bytes();data=json.loads(raw);action_raw=(HERE/'cover.txt').read_bytes()
    need((data['n'],data['mass'])==(10,9),'correct problem and normalization')
    Qs=gram.decode(data['Q_blocks'])
    need(len(Qs)==data['harmonics_through']+1==10,'full PSD block sequence')
    pivots=[gram.positive_ldl(Q) for Q in Qs]
    A,fs=corrections(data);poly,shape=polynomial(data,A,fs);crosscheck(data,A,fs,poly)
    need(len(fs)==5 and shape==(21,21,11),'fixed force and tensor degrees')
    box=((F(0),U),(F(0),U),(-U,U))
    B=initial(poly,shape,box);packer=Packed(shape,bool(fs));old.selftest();stack=[(box,B,0)]
    counts={'A':0,'P':0,'S':0};depth=0
    for i,act in enumerate(action_raw.decode('ascii').strip(),1):
        need(bool(stack),'every action has a remaining region')
        box,B,k=stack.pop();depth=max(depth,k)
        if act=='P':
            need(prune(box),f'exact domain prune {i}');counts['P']+=1
        elif act=='A':
            need(packer.dominates(B,box,C),f'exact Bernstein/log domination {i}, {box}');counts['A']+=1
        else:
            need(act in '012',f'no unresolved node {i}: {act}')
            j=int(act);l,r=packer.split(B,j);m=sum(box[j])/2;left=list(box);right=list(box)
            left[j]=(box[j][0],m);right[j]=(m,box[j][1]);stack.extend(((tuple(right),r,k+1),(tuple(left),l,k+1)))
            counts['S']+=1
        if args.progress and i%10000==0:print(f'exact cover replayed {i}',file=sys.stderr,flush=True)
    need(not stack and counts['A']+counts['P']==counts['S']+1,'exact complete partition')
    V=witness();lower=actual_root_lower();need(lower>V,'true root exceeds old integer witness')
    upper=hx(data['B'])+N*C
    gap=log.log_lower(V*RATIO)-upper
    need(gap>LOG_MARGIN,'strict rational ratio exceeds upper bound by the stated logarithmic margin')
    report={'passed':True,'kernel_sha256':hashlib.sha256(raw).hexdigest(),'tree_sha256':hashlib.sha256(action_raw).hexdigest(),
            'scope':'global nonsharp ten-point upper bound; no maximizer classification',
            'arithmetic':'Python standard library, exact rational and unbounded integer arithmetic',
            'c':str(C),'global_maximum_over_feasible_witness_strict_upper':str(RATIO),
            'actual_candidate_over_global_maximum_strict_lower':str(1/RATIO),
            'B_plus_N_c':str(upper),'strict_log_ratio_margin_lower':str(LOG_MARGIN),
            'visited':sum(counts.values()),'counts':counts,'maximum_depth':depth,'unresolved':0,'PSD_blocks':len(Qs),
            'minimum_exact_LDL_pivot_decimal_display':float(min(map(min,pivots))),'tensor_shape':shape,
            'Ward_parameters':len(data['ward_parameters']),'radial_force_parameters':len(fs),
            'necessary_domain_replayed':True,'radius_strict_upper':str(U_NEW),'root_exceeds_witness':True,
            'sharp_global_optimality_proved':False}
    text=json.dumps(report,indent=2)+'\n';print(text,end='')
    if args.record_report:(HERE/'verification_report.json').write_text(text,encoding='utf-8',newline='\n')

if __name__=='__main__':main()
