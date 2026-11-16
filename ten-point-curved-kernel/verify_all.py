"""Exact arithmetic replay for piecewise Ward kernels on the known ten-point bands."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
C=F(53,1000000)
C_INNER=F(3,50000)
C_MIXED=C
RATIO=F(1253,1250)
BANDS=((F(27,500),F(217,1000)),(F(37,50),F(139,100)))
BOXES=[(0,0,-1),(0,1,-1),(1,1,-1),(0,1,1),(1,1,1)]
WARD={'ward10':(1,(F(1),F(0))),'ward11':(1,(F(0),F(1))),'ward2':(2,(F(-17,10),F(1))),'ward3':(3,(F(-8,5),F(1)))}

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

upper=load('ten-point-ward-upper/verify_all.py','next_n10_old_upper')
collapse=load('ten-point-collapse-domain/verify_all.py','next_n10_collapse')
need=upper.need

def parsed(data):
    out=[]
    for kind,band,power,raw in data['band_parameters']:
        need(kind in ('count','mass','center','force',*WARD) and band in ('all','inner')
             and isinstance(power,int) and 0<=power<=8,'known stationary multiplication direction')
        out.append((kind,band,power,upper.hx(raw)))
    need(len(out)==len(set(row[:3] for row in out)),'distinct zero-sum directions')
    return out

def root_box(i,j,sign):
    if i==j==0:return BANDS[i],BANDS[j],(F(-27,200),F(-1,10))
    product=BANDS[i][1]*BANDS[j][1]*10000
    w=math.isqrt(product.numerator//product.denominator)
    w=F(w+(w*w<product),100)
    return BANDS[i],BANDS[j],(-w,F(0)) if sign==-1 else (F(0),w)

def weight(band,index):
    return band=='all' or band=='inner' and index==0 or band=='high' and index==2

def multiply(left,right):
    out={}
    for a,x in left.items():
        for b,y in right.items():
            key=tuple(p+q for p,q in zip(a,b));out[key]=out.get(key,F(0))+x*y
    return out

def root_polynomial(base,shape,params,i,j):
    result=base.copy()
    d={(1,0,0):F(1),(0,1,0):F(1),(0,0,1):F(-2)}
    for kind,band,power,coefficient in params:
        if kind in WARD:
            k,(a,b)=WARD[kind]
            H=(None,{(0,0,1):1},{(0,0,2):2,(1,1,0):-1},{(0,0,3):4,(1,1,1):-3})[k]
            for axis,b1 in ((0,i),(1,j)):
                if not weight(band,b1):continue
                g={(n,0,0):math.comb(power,n)*(-1)**(power-n) for n in range(power+1)}
                terms=multiply(H,{(0,0,0):a,(0,1,0):b})
                for key,value in (((k,0,0),a/9),((k+1,0,0),b/9)):terms[key]=terms.get(key,F(0))+value
                for (u,v,t),x in multiply(multiply(g,terms),d).items():
                    key=(u,v,t) if axis==0 else (v,u,t)
                    need(all(q<n for q,n in zip(key,shape)),'partial Ward polynomial fits tensor')
                    result[key]=result.get(key,F(0))-coefficient*x
            continue
        for axis,b1,b2 in ((0,i,j),(1,j,i)):
            if not weight(band,b1):continue
            if kind=='force':
                numerator={(1,0,0):F(1),(0,0,1):F(-1),(2,0,0):F(-5,9),
                           (1,1,0):F(-5,9),(1,0,1):F(10,9)}
            else:
                if kind=='count':row={(0,0,0):F(b2==0)-F(2-int(b1==0),9)}
                elif kind=='mass':row={(0,0,0):F(-1),(1,0,0):F(1,9),(0,1,0):F(1)}
                else:row={(0,0,1):F(1),(1,0,0):F(1,9)}
                numerator={}
                for a,x in row.items():
                    for b,y in d.items():
                        key=tuple(p+q for p,q in zip(a,b));numerator[key]=numerator.get(key,F(0))+x*y
            for deg in range(power+1):
                scale=coefficient*math.comb(power,deg)*(-1)**(power-deg)
                for (a,b,c),p in numerator.items():
                    key=(a+deg,b,c) if axis==0 else (b,a+deg,c)
                    need(all(k<n for k,n in zip(key,shape)),'expanded polynomial fits tensor')
                    result[key]=result.get(key,F(0))+scale*p
    return {key:value for key,value in result.items() if value}

def check_root(base,result,params,i,j):
    def at(poly,u,v,t):return sum(co*u**a*v**b*t**c for (a,b,c),co in poly.items())
    n=0
    for x,y in ((F(1,2),F(1,3)),(F(0),F(1)),(F(3,4),F(2,3))):
        u=x*BANDS[i][0]+(1-x)*BANDS[i][1]
        v=y*BANDS[j][0]+(1-y)*BANDS[j][1]
        for t in (-(u+v)/3,(u+v)/5):
            d=u+v-2*t;need(d>0,'nonzero crosscheck distance')
            change=F(0)
            for kind,band,power,q in params:
                for a,b,b1,b2 in ((u,v,i,j),(v,u,j,i)):
                    if not weight(band,b1):continue
                    if kind in WARD:
                        k,(c,e)=WARD[kind]
                        H=(0,t,2*t*t-u*v,4*t**3-3*u*v*t)[k]
                        change-=q*(a-1)**power*(H*(c+e*b)+a**k*(c+e*a)/9)
                        continue
                    g={'count':F(b2==0)-F(2-int(b1==0),9), 'mass':b-(9-a)/9,
                       'center':t+a/9,'force':(a-t)/d-5*a/9}
                    change+=q*(a-1)**power*g[kind]
            need(at(result,u,v,t)==at(base,u,v,t)+d*change,'exact independent piecewise value crosscheck');n+=1
    return n

def prune(box,inner=False):
    if upper.prune(box) or collapse.pair_box_excluded(tuple(a for a,b in box),tuple(b for a,b in box)):
        return True
    if not inner:return False
    (ul,uh),(vl,vh),(tl,th)=box
    need(th<0,'fixed negative inner dot product')
    return ul+vl-2*th>F(3,5) or ul+vl+2*tl>F(1,10) or 9*tl*tl<4*ul*vl

def fact(coeff,index):
    return math.prod(coeff-j for j in range(index))

class Packed(upper.Packed):
    def __init__(self,shape):
        super().__init__(shape,True)
        self.powers=[]
        for i in range(4):
            for j in range(4-i):
                for k in range(4-i-j):
                    power=(i,j,k);degree=i+j+k
                    slots=[(self.ones,1)]
                    for exponent,axis in zip(power,self.axisones):
                        if exponent:
                            slots=[(mask&lane,val*fact(n,exponent)) for n,lane in enumerate(axis) if n>=exponent for mask,val in slots]
                    mask=sum(mask*value for mask,value in slots)
                    denom=math.prod(fact(n-1,e) for n,e in zip(shape,power))
                    multinomial=math.factorial(degree)//math.prod(math.factorial(e) for e in power)
                    self.powers.append((power,degree,mask,denom,multinomial))

    def curved_dominates(self,B,box,c):
        threshold=self.curved_threshold(box,c)
        return ((B+self.high-threshold)&self.high)==self.high

    def curved_threshold(self,box,c):
        (ul,uh),(vl,vh),(tl,th)=box
        dhi=uh+vh-2*tl;dlo=max(upper.DMIN,ul+vl-2*th);a=(dhi+dlo)/2
        need(dhi>=dlo>0,'positive physical d interval for log curvature bound')
        S=upper.S;BIAS=upper.BIAS
        log_a=F(upper.old.log_upper_scaled(a),S)
        q=1/(2*dhi*dhi)
        # d(log(a)-1-c)+d^2/a - d(d-a)^2/(2*dhi^2).
        cd=(F(0),log_a-1-c-q*a*a,1/a+2*q*a,-q)
        d0=ul+vl-2*tl;ds=(uh-ul,vh-vl,-2*(th-tl))
        cs=[sum(cd[k]*math.comb(k,j)*d0**(k-j) for k in range(j,4)) for j in range(4)]
        T=BIAS*self.ones;low=high=BIAS
        for power,j,mask,denom,multi in self.powers:
            v=upper.old.ceilq(S*cs[j]*multi*math.prod(s**e for s,e in zip(ds,power))/denom)
            if j==0:low+=v;high+=v
            else:low+=min(0,v*denom);high+=max(0,v*denom)
            T+=v*mask
        need(0<=low<=high<1<<62,'packed cubic threshold lane bounds')
        return T

def curvature_selftest():
    """Compare the packed threshold with a separate exact tensor conversion."""
    shape=(4,5,4)
    pack=Packed(shape)
    box=((F(1,10),F(1,4)),(F(3,4),F(5,4)),(F(-1,5),F(2,5)))
    c=F(1,16000)
    (ul,uh),(vl,vh),(tl,th)=box
    dhi=uh+vh-2*tl;dlo=max(upper.DMIN,ul+vl-2*th)
    a=(dhi+dlo)/2
    d={(1,0,0):F(1),(0,1,0):F(1),(0,0,1):F(-2)}
    delta=dict(d);delta[(0,0,0)]=-a
    correction=multiply(d,multiply(delta,delta))
    d2=multiply(d,d)
    log_a=F(upper.old.log_upper_scaled(a),upper.S)
    raw={key:d.get(key,0)*(log_a-1-c)+d2.get(key,0)/a-correction.get(key,0)/(2*dhi*dhi)
         for key in set(d)|set(d2)|set(correction)}
    negative=upper.initial({key:-value for key,value in raw.items() if value},shape,box)
    threshold=pack.curved_threshold(box,c)
    mask=(1<<upper.W)-1
    rounding_bound=sum(denom for _,_,_,denom,_ in pack.powers)
    for i in range(math.prod(shape)):
        shift=i*upper.W
        rounded_up=2*upper.BIAS-((negative>>shift)&mask)
        got=(threshold>>shift)&mask
        need(0<=got-rounded_up<rounding_bound,
             f'packed curved bound equals direct convex coefficients up to outward rounding: {i}, {got-rounded_up}')

def split(packer,B,box,depth,j):
    l,r=packer.split(B,j);m=sum(box[j])/2
    left=list(box);right=list(box);left[j]=(box[j][0],m);right[j]=(m,box[j][1])
    return (tuple(right),r,depth+1),(tuple(left),l,depth+1)

def verify(data,tree,c,c_inner=None,c_mixed=None,progress=False):
    need((data['n'],data['mass'])==(10,9),'fixed problem normalization')
    n_q=data['harmonics_through']+1
    need(len(data['Q_blocks'])==n_q,'all harmonic PSD blocks')
    pivots=[upper.gram.positive_ldl(Q) for Q in upper.gram.decode(data['Q_blocks'])]
    A,force=upper.corrections(data)
    need(bool(force),'fixed denominator cleared by multiplying by d')
    C,shape=upper.polynomial(data,A,force);upper.crosscheck(data,A,force,C)
    params=parsed(data);packer=Packed(shape);upper.old.selftest();curvature_selftest()
    need([tuple(row['id']) for row in tree['roots']]==BOXES,'all band pairs in both t signs')
    counts=dict(A=0,Q=0,P=0,S=0);crosschecks=0;statistics=[];max_depth=0
    for (i,j,sgn),row in zip(BOXES,tree['roots']):
        c_region=c_inner if i==j==0 and c_inner is not None else c_mixed if i!=j and c_mixed is not None else c
        box=root_box(i,j,sgn);g=root_polynomial(C,shape,params,i,j)
        crosschecks+=check_root(C,g,params,i,j)
        B=upper.initial(g,shape,box);stack=[(box,B,0)];seq=[]
        for node in row['tree']:
            need(bool(stack),'tree node has a remaining box')
            box,B,k=stack.pop();max_depth=max(max_depth,k)
            if node in '012':
                stack.extend(split(packer,B,box,k,int(node)));seq.append(node);continue
            if (node=='P' and prune(box,i==j==0) or node=='A' and packer.dominates(B,box,c_region)
                    or node=='Q' and packer.curved_dominates(B,box,c_region)):
                seq.append(node);continue
            raise ArithmeticError(f'invalid {node} leaf at region {i,j,sgn}, box {box}')
        need(not stack,'fully expanded root and exact coverage')
        seq=''.join(seq);row['tree']=seq
        stats={k:seq.count(k) for k in 'AQP012'}
        splits=sum(stats[k] for k in '012')
        need(stats['A']+stats['Q']+stats['P']==1+splits,'binary tree conservation')
        for key in 'AQP':counts[key]+=stats[key]
        counts['S']+=splits
        statistics.append({'id':[i,j,sgn],'counts':stats,'visited':len(seq)})
        if progress:
            print('exact region',i,j,sgn,len(seq),'nodes',file=sys.stderr,flush=True)
    B=upper.hx(data['B'])
    c_inner=c if c_inner is None else c_inner
    c_mixed=c if c_mixed is None else c_mixed
    return {'passed':True,'scope':'exact ten-point curved majorant on the accepted two bands','c':str(c),'c_inner':str(c_inner),'c_mixed':str(c_mixed),'B':str(B),
            'visited':sum(counts.values()),'counts':counts,'unresolved':0,'maximum_depth':max_depth,
            'regions':statistics,'tensor_shape':shape,'PSD_blocks':n_q,'min_gram_pivot_display':float(min(map(min,pivots))),
            'polynomial_crosschecks':crosschecks,'piecewise_zero_sum_directions':len(params),
            'logarithmic_gap_upper_display':float(B+28*c+16*c_mixed+c_inner-upper.log.log_lower(upper.witness()))}

def read_json(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_dependencies():
    manifest = read_json(HERE/'source_manifest.json')
    for group, folder in (('files', HERE), ('dependencies', ROOT)):
        for path, value in manifest[group].items():
            need(digest(folder/path)==value, 'fixed proof input '+path)
    banded = load('ten-point-banded-kernel/verify_all.py', 'curved_prior_banded')
    banded.verify_dependencies()
    need(tuple(map(F, banded.BANDS[1]))==BANDS[1], 'inherited outer radius interval')
    saved = read_json(ROOT/'ten-point-banded-kernel-refined/verification_report.json')
    need(saved['passed'] is True and saved['unresolved']==0
         and (saved['exact_inner_count'],saved['exact_outer_count'])==(2,8),
         'accepted refined proof of the shell counts')
    inner = saved['inner_pair']
    need(inner['applies_to_every_global_maximizer'] is True and inner['unresolved']==0,
         'accepted complete inner pair localization')
    need(tuple(map(F,inner['inner_squared_radius_open_interval']))==BANDS[0] and
         tuple(map(F,inner['inner_dot_product_open_interval']))==(F(-27,200),F(-1,10)) and
         F(inner['inner_midpoint_squared_strict_upper'])==F(1,40) and
         F(inner['angle_cosine_strict_upper'])==F(-2,3) and
         F(inner['inner_squared_distance_strict_upper'])==F(3,5),
         'the inherited localization exactly implies the inner-pair domain')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report', action='store_true')
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args()
    verify_dependencies()
    kernel_path = HERE/'kernel.json'
    cover_path = HERE/'cover.json'
    data = read_json(kernel_path)
    tree = read_json(cover_path)
    need((F(tree['c']),F(tree['c_inner']),F(tree['c_mixed']))==(C,C_INNER,C_MIXED),
         'fixed three pair bounds')
    need(len(parsed(data))==54, 'fixed piecewise kernel directions')
    V = upper.witness()
    B = upper.hx(data['B'])
    budget = B+C_INNER+16*C_MIXED+28*C
    need(upper.actual_root_lower()>V, 'independently accepted candidate exceeds integer witness')
    margin = upper.log.log_lower(V*RATIO)-budget
    need(margin>F(1,200000), 'strict global/witness ratio with rational logarithmic margin')
    inner_pair_budget = B-upper.log.log_lower(V)+16*C_MIXED+28*C
    need(inner_pair_budget<F(2333,10**6), 'strict single inner-pair deficit bound')
    report = verify(data,tree,C,C_INNER,C_MIXED,progress=args.progress)
    report.update(kernel_sha256=digest(kernel_path), proof_sha256=digest(cover_path),
                  global_over_witness_strict_upper=str(RATIO), candidate_over_global_strict_lower=str(1/RATIO),
                  log_ratio_margin_strict_lower='1/200000', log_upper_mass_nine=str(budget),
                  inner_pair_deficit_strict_upper='2333/1000000', exact_inner_count=2, exact_outer_count=8,
                  curvature_reference_lanes_checked=80, sharp_global_optimum_proved=False)
    output = json.dumps(report,indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_text(output,encoding='utf-8',newline='\n')
    print(output,end='')


if __name__=='__main__':
    main()
