"""Restricted affine pair fans and finite annular population/layer duals."""
from fractions import Fraction as F
from functools import lru_cache
import sys
from acceptance import ROOT,require
sys.path.insert(0,str(ROOT/'joint-variational'))
from affine_pair import support_line,tail_interval,pi_bounds
sys.path.insert(0,str(ROOT/'radial-angular'))
from radial_layers import smearing_penalty_upper


def floor_mesh(x,mesh=10**9):
    x=F(x);return F((x*mesh).numerator//(x*mesh).denominator,mesh)


@lru_cache(maxsize=64)
def line_at(center):
    center=F(center);eps=F(1,100000)
    if not F(143,100)<=center<=F(19,10):raise ValueError('implemented center range 1.43..1.90')
    lo,hi=tail_interval(center)
    beta=floor_mesh((lo+hi)/center)
    row=support_line(center-eps,center+eps,beta)
    return {'center':center,'left':row['left'],'right':row['right'],
        'alpha':row['alpha'],'beta':row['beta']}


def fan(h,centers):
    h=F(h)
    if not F(1)<=h<=F(11,10):raise ValueError('implemented 1<=h<=1.1')
    centers=tuple(map(F,centers))
    if not centers or len(set(centers))!=len(centers):raise ValueError('nonempty distinct centers')
    rows=sorted((line_at(c) for c in centers),key=lambda r:r['beta'])
    cmin=floor_mesh(pi_bounds()[0]/(2*h),100000)
    require(cmin*cmin>=2,'restricted domain remains in implemented lens range')
    return {'h':h,'domain_lower':cmin,'lines':rows}


def restricted_line(packet,beta):
    """Positive mixtures remain global; upper-slope tilt is domain-restricted."""
    beta=F(beta)
    if beta<0:raise ValueError('nonnegative slope')
    rows=[{'alpha':F(0),'beta':F(0)}]+packet['lines']
    for left,right in zip(rows,rows[1:]):
        if beta<=right['beta']:
            t=(beta-left['beta'])/(right['beta']-left['beta'])
            return (1-t)*left['alpha']+t*right['alpha']
    last=rows[-1]
    return last['alpha']+(beta-last['beta'])*packet['domain_lower']


def annular_certificate(packet,layer_price,boundaries):
    """Finite exact coefficients; all asymptotic quantifiers are analytic."""
    lam=F(layer_price);ends=list(map(F,boundaries))
    if lam<=0 or len(ends)<2 or ends[0]!=0 or ends[-1]!=1:
        raise ValueError('positive price and a partition from 0 to 1')
    if any(a>=b for a,b in zip(ends,ends[1:])):raise ValueError('strictly increasing partition')
    h=packet['h'];pu=pi_bounds()[1];rows=[];integral=F(0)
    for a,b in zip(ends,ends[1:]):
        beta=lam*h/(2*pu*b)
        alpha=restricted_line(packet,beta)
        require(alpha>=0 and 2*pu*beta*b/h==lam,'band price and sign')
        integral+=(b*b-a*a)*alpha
        rows.append({'a':a,'b':b,'alpha':alpha,'beta':beta})
    return {'h':h,'layer_price':lam,'bands':rows,'uniform_radial_integral_lower':integral,
        'maximum_alpha':max(row['alpha'] for row in rows),
        'restricted_domain_lower':packet['domain_lower']}


def asymptotic_lower(certificate,gap_per_n_upper):
    e=F(gap_per_n_upper)
    if e<0:raise ValueError('nonnegative proven asymptotic gap')
    return (certificate['uniform_radial_integral_lower']-e-
            smearing_penalty_upper(certificate['h']))/certificate['layer_price']


def finite_lower(certificate,n,K,populations,dirichlet_lower=0,layer_counts=None):
    """Actual band populations; the returned delta lower is one shared budget."""
    if type(n)is not int or n<2 or type(K)is not int or not 1<=K<=n:
        raise ValueError('integer n>=2 and 1<=K<=n')
    populations=list(populations);A=F(dirichlet_lower)
    if len(populations)!=len(certificate['bands']) or any(type(m)is not int or m<0 for m in populations) or sum(populations)!=n:
        raise ValueError('nonnegative integer band populations sum to n')
    if A<0:raise ValueError('nonnegative independently proved same-h Dirichlet lower')
    from affine_pair import sqrt_bounds
    measured=sum((row['alpha']*m for row,m in zip(certificate['bands'],populations)),F(0))
    if layer_counts is None:
        correction=certificate['maximum_alpha']*K
    else:
        counts=list(layer_counts)
        if len(counts)!=len(populations) or any(type(k)is not int or not 0<=k<=m or bool(k)!=bool(m) for k,m in zip(counts,populations)) or sum(counts)!=K:
            raise ValueError('actual integer band layer counts, compatible with populations, sum to K')
        correction=sum((row['alpha']*k for row,k in zip(certificate['bands'],counts)),F(0))
    return A+measured-certificate['layer_price']*sqrt_bounds(n)[1]*K-correction-n*smearing_penalty_upper(certificate['h'])
