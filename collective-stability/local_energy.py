"""Localized angular Dirichlet certificates with exact harmonic disk means.

The analytic theorem is LOCAL_ENERGY.md. Squared radii use Coulomb scale
sum|x_i|²=n-1. A supplied local moment lower bound is a hypothesis that the
caller must establish; radial bands alone do not imply angular alignment.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'joint-bounds'))
from subset_moment import log_interval, outward_decimal


def require_n(n):
    if type(n) is not int or n < 2:
        raise ValueError('integer n>=2 required')


def validate_hat(k, a, b, c):
    if type(k) is not int or k < 1:
        raise ValueError('integer angular mode k>=1 required')
    a,b,c=map(F,(a,b,c))
    if not 0 <= a < b < c:
        raise ValueError('squared radii 0<=a<b<c required')
    return k,a,b,c


def hat_norm(k, a, b, c):
    """Dirichlet norm squared of EACH real/imaginary annular hat."""
    k,a,b,c=validate_hat(k,a,b,c)
    return k*(b**k-a**k)*(c**k-a**k)/(c**k-b**k)


def center_weight(k, a, b, c, q):
    """Phi(x)=x^k*weight(|x|²). At zero the weight can be chosen as zero."""
    k,a,b,c=validate_hat(k,a,b,c)
    q=F(q)
    if q<0:
        raise ValueError('nonnegative squared radius required')
    if q<=a or q>=c:
        return F(0)
    if q<=b:
        return 1-(a/q)**k
    t=(b**k-a**k)/(c**k-b**k)
    return t*((c/q)**k-1)


def upper_support_le(q, r2, s):
    q,r2,s=map(F,(q,r2,s))
    if min(q,r2,s)<0:
        raise ValueError('nonnegative squared inputs required')
    d=s-q-r2
    return d>=0 and d*d>=4*q*r2


def lower_support_ge(q, r2, s):
    q,r2,s=map(F,(q,r2,s))
    if min(q,r2,s)<0:
        raise ValueError('nonnegative squared inputs required')
    if s==0:
        return True
    d=q-r2-s
    return d>=0 and d*d>=4*r2*s


def harmonic_region(qlo, qhi, r2, a, b, c):
    """Prove every disk in a q-band stays in one harmonic piece.

    Boundary tangency is allowed. Unknown crossings are rejected, rather
    than applying the harmonic mean formula to a nonharmonic whole disk.
    """
    qlo,qhi,r2,a,b,c=map(F,(qlo,qhi,r2,a,b,c))
    if not 0<=qlo<=qhi or r2<=0 or not 0<=a<b<c:
        raise ValueError('valid q-band, positive smear and ordered hat required')
    if upper_support_le(qhi,r2,a):
        return 'inside_hole'
    if lower_support_ge(qlo,r2,c):
        return 'outside_support'
    if lower_support_ge(qlo,r2,a) and upper_support_le(qhi,r2,b):
        return 'inner_harmonic_piece'
    if lower_support_ge(qlo,r2,b) and upper_support_le(qhi,r2,c):
        return 'outer_harmonic_piece'
    raise ValueError('a smeared disk may cross a hat interface')


def validate_orthogonal_hats(hats):
    """Different modes are orthogonal; equal modes need disjoint supports."""
    hats=[validate_hat(*hat) for hat in hats]
    if not hats:
        raise ValueError('at least one hat required')
    for i,(k,a,b,c) in enumerate(hats):
        for ell,aa,bb,cc in hats[:i]:
            if k==ell and max(a,aa)<min(c,cc):
                raise ValueError('same-mode overlapping hats require their full Gram matrix')
    return hats


def smear_penalty_upper(n,h,terms=28):
    require_n(n)
    h=F(h)
    if h<=0:
        raise ValueError('positive smearing multiplier required')
    lo,hi=log_interval(h,terms)
    return n*max(F(0),(h*h-1)/2-lo)


def sector_population_upper(n,proven_gap_upper,k,a,b,c,q_lower,gamma,
                            h=1,radial_delta_lower=0,terms=28):
    """Necessary cap on an angularly aligned local population.

    Each counted center must lie in the inner harmonic piece at q>=q_lower
    and satisfy cos(k*(theta-theta0))>=gamma>0. All OTHER centers must have
    zero hat measurement (or nonnegative projected measurement). Disk
    containment must be proved separately. See LOCAL_ENERGY.md.
    """
    require_n(n)
    k,a,b,c=validate_hat(k,a,b,c)
    gap,qlow,gamma,h,radial=map(F,(proven_gap_upper,q_lower,gamma,h,radial_delta_lower))
    if gap<0 or not a<qlow<=b or not 0<gamma<=1 or h<=0 or radial<0:
        raise ValueError('valid gap, inner radius, alignment and smear bounds required')
    budget=gap+smear_penalty_upper(n,h,terms)-radial
    if budget<0:
        raise ValueError('supplied radial lower bound already exceeds the energy budget')
    base=gamma*gamma*qlow**k*(1-(a/qlow)**k)**2
    allowed_square=hat_norm(k,a,b,c)*budget/base
    maximum=min(n,isqrt(allowed_square.numerator//allowed_square.denominator))
    return {'n':n,'maximum_population':maximum,'allowed_population_square_upper':allowed_square,
            'per_point_projected_moment_square_lower':base,
            'scope':'Necessary under stated alignment, noncancellation and harmonic-disk-containment hypotheses; not a universal disk population bound.'}


def band_local_certificate(n,bands,hats,moment_squared_lowers,h=1,
                           radial_delta_lower=0,terms=28):
    """Conditional class bound using independently PROVED local moments.

    bands=(qlo,qhi,population), hats=(k,a,b,c), and each supplied moment
    lower bounds |sum_i Phi_hat(x_i)|². The band checks establish harmonic
    disk means. Nonzero radial_delta_lower must use this same smear h.
    """
    require_n(n)
    hats=validate_orthogonal_hats(hats)
    bands=[(F(a),F(b),m) for a,b,m in bands]
    lows=list(map(F,moment_squared_lowers))
    h,radial=F(h),F(radial_delta_lower)
    if h<=0 or radial<0 or len(lows)!=len(hats) or any(v<0 for v in lows):
        raise ValueError('positive smear, matching nonnegative bounds required')
    if not bands or any(not 0<=a<=b or type(m)is not int or m<1 for a,b,m in bands):
        raise ValueError('valid bands and positive integer populations required')
    if sum(m for a,b,m in bands)!=n:
        raise ValueError('band populations must sum to n')
    if not sum(m*a for a,b,m in bands)<=n-1<=sum(m*b for a,b,m in bands):
        raise ValueError('bands incompatible with Coulomb quadratic moment')
    r2=2*h*h/n
    rows=[]
    angular=F(0)
    for (k,a,b,c),low in zip(hats,lows):
        regions=[harmonic_region(qa,qb,r2,a,b,c) for qa,qb,m in bands]
        norm=hat_norm(k,a,b,c)
        value=low/norm
        angular+=value
        rows.append({'k':k,'a':a,'b':b,'c':c,'norm':norm,
                     'band_regions':regions,'moment_squared_lower':low,
                     'angular_delta_lower':value})
    penalty=smear_penalty_upper(n,h,terms)
    return {'n':n,'smear_multiplier':h,'smear_squared':r2,
            'angular_delta_lower':angular,'radial_delta_lower':radial,
            'smear_penalty_upper':penalty,
            'delta_lower':max(F(0),angular+radial-penalty),'hats':rows,
            'scope':'Conditional necessary energy bound. Local moment inputs and nonzero radial input require independent proof; no full domain covering.'}


def orthogonal_holomorphic_bonus(n,bands,hats,moment_squared_lowers,
                                 support_squared,h=1):
    """Additional global j=0 angular terms orthogonal to compact local hats.

    moment_squared_lowers is a mapping {k: proved lower bound for |p_k(x)|²}.
    The support radius must contain EVERY smeared disk and EVERY hat. This
    is a bonus to n²D angular energy, with no second smear-penalty deduction.
    Global radial-weighted j>=1 tests are generally NOT orthogonal to hats.
    """
    require_n(n)
    hats=validate_orthogonal_hats(hats)
    S,h=F(support_squared),F(h)
    if S<=0 or h<=0:
        raise ValueError('positive support squared and smear multiplier required')
    bands=[(F(a),F(b),m) for a,b,m in bands]
    if not bands or any(not 0<=a<=b or type(m)is not int or m<1 for a,b,m in bands):
        raise ValueError('valid bands and positive integer populations required')
    if sum(m for a,b,m in bands)!=n:
        raise ValueError('band populations must sum to n')
    if not sum(m*a for a,b,m in bands)<=n-1<=sum(m*b for a,b,m in bands):
        raise ValueError('bands incompatible with Coulomb quadratic moment')
    r2=2*h*h/n
    if any(not upper_support_le(b,r2,S) for a,b,m in bands):
        raise ValueError('global support fails to contain a smeared disk')
    if any(c>S for k,a,b,c in hats):
        raise ValueError('global support fails to contain a local hat')
    if not hasattr(moment_squared_lowers,'items'):
        raise ValueError('moment lower bounds must be a mapping indexed by unique k')
    rows=[]
    total=F(0)
    for k,low in sorted(moment_squared_lowers.items()):
        if type(k)is not int or k<1 or F(low)<0:
            raise ValueError('positive integer modes and nonnegative moment bounds required')
        value=F(low)/(k*S**k)
        total+=value
        rows.append({'k':k,'moment_squared_lower':F(low),'angular_bonus':value})
    return {'n':n,'support_squared':S,'smear_multiplier':h,
            'holomorphic_angular_bonus':total,'modes':rows,
            'scope':'Global j=0 angular bonus, orthogonal to all listed compact hats. Add to their angular bound and subtract the common smearing penalty only once.'}


def cmul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def cpow(z,k):
    ans=F(1),F(0)
    for _ in range(k):
        ans=cmul(ans,z)
    return ans


def shape_local_certificate(points,hats,h=1,radial_delta_lower=0,terms=28):
    """Normalize centered rational coordinates and calculate local moments.

    Each disk must avoid every relevant hat interface. No polynomial mean
    formula is substituted when this geometric requirement fails.
    """
    pts=[(F(a),F(b)) for a,b in points]
    n=len(pts)
    require_n(n)
    if len(set(pts))!=n:
        raise ValueError('distinct input coordinates required')
    hats=validate_orthogonal_hats(hats)
    h=F(h)
    if h<=0:
        raise ValueError('positive smearing multiplier required')
    mx=sum(a for a,b in pts)/n
    my=sum(b for a,b in pts)/n
    v=[(a-mx,b-my) for a,b in pts]
    rawq=[a*a+b*b for a,b in v]
    factor=F(n-1)/sum(rawq)
    qs=[factor*q for q in rawq]
    bands=[(q,q,1) for q in qs]
    lows=[]
    moments=[]
    for k,a,b,c in hats:
        re=im=F(0)
        for point,q in zip(v,qs):
            harmonic_region(q,q,2*h*h/n,a,b,c)
            p=cpow(point,k)
            weight=center_weight(k,a,b,c,q)
            re+=p[0]*weight
            im+=p[1]*weight
        low=factor**k*(re*re+im*im)
        lows.append(low)
        moments.append({'common_factor_squared':factor**k,
                        'raw_real':re,'raw_imaginary':im})
    result=band_local_certificate(n,bands,hats,lows,h,radial_delta_lower,terms)
    result['normalization_factor_squared']=factor
    result['exact_local_moments']=moments
    return result


if __name__=='__main__':
    import json
    result=band_local_certificate(100,
        [(F(7,20),F(7,20),35),(F(347,260),F(347,260),65)],
        [(2,0,F(27,50),F(51,50))],[121])
    print(json.dumps(result,default=str,indent=2))
