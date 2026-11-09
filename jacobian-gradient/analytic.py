"""Generate exact block gradients with Jacobian and reuse them for search.

Generation requires installed Jacobian. Residual evaluation below requires
only Python's standard library and is exploratory floating arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from importlib.metadata import version
import argparse,hashlib,importlib.util,json,math,time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'twelve-point-kernel-refinement/kernel.json'
DEFAULT_OUTPUT=ROOT/'work/jacobian-gradient-generation'
AXES=('u','v','t')


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def clean(p):return {e:c for e,c in p.items() if c}
def add(target,e,c):target[e]=target.get(e,F(0))+c
def add_poly(target,source):
    for e,c in source.items():add(target,e,c)


def derivative(p,axis):
    out={}
    for exps,c in p.items():
        if exps[axis]:
            e=list(exps);e[axis]-=1;add(out,tuple(e),c*exps[axis])
    return clean(out)


def evaluate(p,point):
    return sum((c*point[0]**e[0]*point[1]**e[1]*point[2]**e[2]
                for e,c in p.items()),F(0))


def from_native(p):
    return {tuple(t.exponents):t.coefficient.as_fraction() for t in p.polynomial.terms}


def to_plain(p):
    return [{'exponents':list(e),'coefficient':{'num':str(c.numerator),'den':str(c.denominator)}}
            for e,c in sorted(p.items(),reverse=True) if c]


def from_plain(terms):
    return {tuple(t['exponents']):F(int(t['coefficient']['num']),int(t['coefficient']['den'])) for t in terms}


def block_polynomials(data,helper):
    modes=data['harmonics_through'];n=data['n'];mass=n-1;B=helper.hx(data['B'])
    b=(B-F(n*mass,2))/n
    base={(0,0,0):2*b/mass,(1,0,0):F(n,2*mass),(0,1,0):F(n,2*mass)}
    H=[{(0,0,0):1},{(0,0,1):1}]
    for k in range(2,modes+1):
        h={}
        for (a,b,c),v in H[-1].items():add(h,(a,b,c+1),2*v)
        for (a,b,c),v in H[-2].items():add(h,(a+1,b+1,c),-v)
        H.append(clean(h))
    blocks=[]
    for k,raw in enumerate(data['Q_blocks']):
        Q=[[helper.hx(q) for q in row] for row in raw];P=helper.features(data,k);dim=len(Q)
        G=[[sum((P[i][a]*Q[a][b]*P[j][b] for a in range(dim) for b in range(dim)),F(0))
            for j in range(len(P))] for i in range(len(P))]
        C={}
        for i,row in enumerate(G):
            for j,g in enumerate(row):
                for (a,b,c),v in H[k].items():add(C,(a+i,b+j,c),-v*g)
                add(C,(i+j+k,0,0),-g/(2*mass))
                add(C,(0,i+j+k,0),-g/(2*mass))
        blocks.append(clean(C))
    return clean(base),blocks


def scaled_native_gradient(p,label,native,gradient,RationalPolynomial):
    """Compose legal bounded native gradients by exact Q-linearity.

    A direct block call exceeds Jacobian's 128-digit input guard.  Use a
    common denominator and signed base-2^250 integer coefficient limbs.
    Every native input has <=76-digit integer coefficients; every omitted
    scale is restored exactly after the native producer returns.
    """
    D=math.lcm(*(c.denominator for c in p.values()))
    integers={e:c.numerator*(D//c.denominator) for e,c in p.items()}
    count=max(abs(c).bit_length() for c in integers.values())
    count=(count+249)//250;mask=(1<<250)-1
    answer=[{}, {}, {}];reconstruction={}
    for limb in range(count):
        terms={e:((abs(c)>>(250*limb))&mask)*(1 if c>=0 else -1) for e,c in integers.items()}
        terms=clean(terms)
        if not terms:continue
        scale=F(1<<(250*limb),D)
        ordered=sorted(terms.items(),reverse=True)
        # The gradient admission guard counts total output terms across all
        # three coordinates.  Sixty-four source terms give at most 192.
        for offset in range(0,len(ordered),64):
            part=dict(ordered[offset:offset+64])
            source=RationalPolynomial.model_validate({'domain':'QQ','variables':AXES,
                'polynomial':{'terms':[{'coefficient':{'num':c,'den':1},'exponents':list(e)}
                  for e,c in part.items()]}})
            result=gradient(source)
            for e,c in part.items():add(reconstruction,e,scale*c)
            for j,component in enumerate(result.components):
                for e,c in from_native(component).items():add(answer[j],e,scale*c)
            native.append({'label':label,'integer_limb':limb,'source_term_offset':offset,
                'exact_restoring_scale':{'num':str(scale.numerator),'den':str(scale.denominator)},
                'operation':'polynomial_field.scalar.gradient.compute','output':result.model_dump(mode='json')})
    need(clean(reconstruction)==p,'signed integer limbs exactly reconstruct source block')
    return list(map(clean,answer)),count


class Residual:
    """Exploratory float residual and analytic gradient from generated data.

    No symbolic differentiation or finite differences occur in this evaluator.
    Global acceptance still belongs to the exact continuous cover.
    """
    def __init__(self,path=HERE/'gradient.json'):
        data=json.loads(Path(path).read_bytes());self.metadata=data
        self.h=self._convert(data['h'])
        self.g=[self._convert(t) for t in data['gradient']]
        self.degrees=[max(e[j] for terms in [self.h,*self.g] for e,c in terms) for j in range(3)]

    @staticmethod
    def _convert(terms):
        return [(tuple(t['exponents']),float(F(int(t['coefficient']['num']),int(t['coefficient']['den'])))) for t in terms]

    @staticmethod
    def _eval(terms,point,degrees):
        powers=[[float(x)**k for k in range(degree+1)] for x,degree in zip(point,degrees)]
        return math.fsum(c*powers[0][e[0]]*powers[1][e[1]]*powers[2][e[2]] for e,c in terms)

    def value(self,point):
        u,v,t=point;d=u+v-2*t
        if d<=0:return math.inf
        return self._eval(self.h,point,self.degrees)-math.log(d)

    def jac(self,point):
        u,v,t=point;d=u+v-2*t
        if d<=0:raise ValueError('positive squared distance required')
        return [self._eval(p,point,self.degrees)-a/d for p,a in zip(self.g,(1,1,-2))]


def build(kernel_path=SOURCE,output_dir=DEFAULT_OUTPUT):
    start=time.monotonic()
    kernel_path=Path(kernel_path).resolve();output_dir=Path(output_dir).resolve()
    need(not output_dir.is_relative_to(HERE),'generation cannot overwrite the accepted directory')
    output_dir.mkdir(parents=True,exist_ok=True)
    from jacobian.math.polynomials.values import RationalPolynomial
    from jacobian.math.polynomials.vector_calculus.operations import gradient
    need(version('jacobian')=='0.23.0','pinned Jacobian')
    helper_path=ROOT/'twelve-point-kernel-upper/verify_cover.py'
    spec=importlib.util.spec_from_file_location('gradient_strict_helper',helper_path)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    raw=kernel_path.read_bytes();data=json.loads(raw);base,blocks=block_polynomials(data,helper)
    combined=dict(base)
    for p in blocks:add_poly(combined,p)
    combined=clean(combined);whole,_=helper.polynomial(data)
    need(combined==whole,'harmonic decomposition reconstructs independently expanded complete h')
    grad=[{}, {}, {}];native=[];summary=[]
    for label,p in [('affine_base',base)]+[(str(k),p) for k,p in enumerate(blocks)]:
        need(len(p)<=256,'one block fits the bounded exact polynomial API')
        exact,limbs=scaled_native_gradient(p,label,native,gradient,RationalPolynomial)
        for j,g in enumerate(exact):
            need(g==derivative(p,j),'Jacobian block gradient exact coefficient identity')
            add_poly(grad[j],g)
        summary.append({'label':label,'source_support':len(p),'gradient_supports':list(map(len,exact)),
                        'integer_coefficient_limbs':limbs})
    grad=list(map(clean,grad))
    for j,g in enumerate(grad):need(g==derivative(whole,j),'combined native gradient equals complete monomial derivative')
    points=[(F(1,4),F(2,3),-F(1,5)),(F(7,6),F(3,2),-F(3,4)),
            (F(5264,10000),F(121,72),F(9033,10000)),tuple(map(F,data['test_pair']))]
    crosschecks=[]
    for point in points:
        d=point[0]+point[1]-2*point[2];need(d>0,'rational test squared distance')
        expected=[evaluate(derivative(whole,j),point)-a/d for j,a in enumerate((1,1,-2))]
        actual=[evaluate(g,point)-a/d for g,a in zip(grad,(1,1,-2))]
        need(actual==expected,'exact residual gradient at independent rational point')
        crosschecks.append({'point':list(map(str,point)),'exact_gradient_agrees':True})
    artifact={'variables':list(AXES),'kernel_sha256':hashlib.sha256(raw).hexdigest(),
        'jacobian_version':version('jacobian'),'operation':'polynomial_field.scalar.gradient.compute',
        'h':to_plain(whole),'gradient':[to_plain(g) for g in grad]}
    (output_dir/'gradient.json').write_text(json.dumps(artifact,indent=2)+'\n',encoding='utf-8',newline='\n')
    (output_dir/'native_blocks.json').write_text(json.dumps(native,indent=2)+'\n',encoding='utf-8',newline='\n')
    r=Residual(output_dir/'gradient.json')
    float_errors=[]
    for point in points:
        d=point[0]+point[1]-2*point[2]
        exact=[float(evaluate(g,point)-a/d) for g,a in zip(grad,(1,1,-2))]
        observed=r.jac(list(map(float,point)))
        float_errors.append(max(abs(a-b) for a,b in zip(exact,observed)))
    report={'passed':True,'kernel_source':kernel_path.relative_to(ROOT).as_posix() if kernel_path.is_relative_to(ROOT) else kernel_path.name,
        'kernel_sha256':hashlib.sha256(raw).hexdigest(),
        'helper_sha256':hashlib.sha256(helper_path.read_bytes()).hexdigest(),
        'actual_Jacobian_native_gradient_calls':len(native),'blocks':summary,
        'admission_strategy':'exact common denominator, signed 250-bit integer coefficient limbs and <=64 source terms per call; each native input has <=76-digit coefficients and <=192 gradient terms',
        'direct_block_coefficient_guard_observed':'128 digits; original block0 has coefficients exceeding it',
        'direct_total_gradient_term_guard_observed':256,
        'complete_h_support':len(whole),'complete_gradient_supports':list(map(len,grad)),
        'global_coefficient_identity_checked':True,'rational_point_crosschecks':crosschecks,
        'exploratory_float_gradient_absolute_errors':float_errors,
        'gradient_artifact_sha256':hashlib.sha256((output_dir/'gradient.json').read_bytes()).hexdigest(),
        'native_results_sha256':hashlib.sha256((output_dir/'native_blocks.json').read_bytes()).hexdigest(),
        'seconds':time.monotonic()-start,
        'scope':'exact derivative construction and exploratory numerical evaluator; no new global optimum claim'}
    (output_dir/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true')
    parser.add_argument('--kernel',type=Path,default=SOURCE)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    if args.build:print(json.dumps(build(args.kernel,args.output),indent=2))
    else:
        residual=Residual();point=[.5264,121/72,.9033]
        print(json.dumps({'point':point,'residual':residual.value(point),'gradient':residual.jac(point)},indent=2))
