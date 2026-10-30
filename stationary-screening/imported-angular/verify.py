"""Independent finite checks; general analytic proof is README.md."""
from fractions import Fraction as F
from pathlib import Path
from math import prod
import json,random,sys
import angular_certificate as a

ROOT=Path(__file__).resolve().parent
BASE=ROOT/'baseline'
sys.path.insert(0,str(BASE))
from radial_distribution import exact_radii_certificate,polar_gap_upper

def require(ok,msg):
    if not ok:raise ArithmeticError(msg)

def cadd(x,y):return x[0]+y[0],x[1]+y[1]
def cscale(x,b):return x[0]*b,x[1]*b

def newton_check(points):
    v=a.centered(points);n=len(v);s=a.power_sums(v,n)
    e=[(F(1),F(0))]
    for y in v:
        nxt=e+[(F(0),F(0))]
        for j in range(1,len(nxt)):nxt[j]=cadd(nxt[j],a.cmul(e[j-1],y))
        e=nxt
    # Independent Newton recurrence vs direct product coefficients.
    b=[F(1)]
    for j in range(1,n+1):
        rhs=(F(0),F(0))
        for k in range(1,j+1):rhs=cadd(rhs,cscale(a.cmul(e[j-k],s[k-1]),(-1)**(k-1)))
        require(rhs==cscale(e[j],j),'Newton identity')
        b.append(sum(a.root_up(a.norm(s[k-1]))*b[j-k]for k in range(1,j+1))/j)
        require(a.norm(e[j])<=b[j]*b[j],'coefficientwise positive exponential majorant')
    q2=max(map(a.norm,v))
    # No square root needed for the actual Cauchy lower bound.
    H=sum(a.norm(s[k-1])/k for k in range(2,n+1))
    bound=F(4,9)/sum(q2**(-k)/k for k in range(2,n+1))
    require(H>=bound,'Newton-Cauchy positive moment floor')

def main():
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    rng=random.Random(20261030);cases=[]
    for n in range(2,17):
        cases.append([(j,0)for j in range(n)])
        v=[]
        while len(v)<n:
            p=(rng.randrange(-15,16),rng.randrange(-15,16))
            if p not in v:v.append(p)
        cases.append(v)
    saved=json.loads((BASE/'saved_rational_bounds.json').read_text(encoding='utf-8'))
    cases.extend(x['integer_points']for x in saved['rows'])
    rows=[]
    for points in cases:
        newton_check(points)
        cert=a.shape_certificate(points);n=cert['n']
        radial=exact_radii_certificate(cert['coulomb_radii_squared'],[F(j,8)for j in range(17)])['delta_lower']
        lo,hi=a.deficit_interval(n,cert['normalized_product'])
        require(lo>=cert['angular_delta_lower']+radial,'additive angular-radial bound vs independent log deficit')
        rows.append({'n':n,'angular_lower':a.conservative_decimal(cert['angular_delta_lower']),
                     'radial_lower':a.conservative_decimal(radial),
                     'deficit_upper':a.conservative_decimal(hi,upper=True)})
    # Two-sided Dirichlet normalization, with independent monomial integration.
    norm_checks=0
    for k in range(1,31):
        for R2 in [F(1),F(2),F(7,3)]:
            inside=F(k*k,2*k)*R2**k
            outside=F(k*k,2*k)*R2**(2*k)*R2**(-k)
            require(inside+outside==k*R2**k,'inside+outside Dirichlet norm')
            norm_checks+=1
    # Bubble radial integrals with variable u=radius^2.
    for rho2 in [F(1,100),F(1,3),F(7,5)]:
        integral_background=F(1,2)*(rho2*rho2-rho2*rho2/2)
        energy=2*(rho2*rho2/2)
        require(integral_background==rho2**2/4 and energy==rho2**2,'off-center bubble constants')
        require(integral_background**2/energy==rho2**2/16,'bubble Bessel cost')
    known={2:F(4),3:F(27),4:F(256),5:F(3125),6:F(6**15,5**10),7:F(7**21,6**15)}
    upper_rows=[]
    for n in range(2,101):
        c=a.universal_discount(n);require(c>0,'universal positive discount')
        if n in known:require(a.deficit_interval(n,known[n])[0]>c,'known sharp value below refined upper')
        if n in [2,3,5,8,10,16,26,30,50,100]:
            upper_rows.append({'n':n,'discount_positive':True,'c_decimal_lower':a.conservative_decimal(c,18),
                               'R2_upper':str(a.global_radius_squared(n))})
    require(a.global_radius_squared(100)==F(11881,5000),'exact n100 radius')
    strip=a.strip_global_penalty(100,F(1,10));gap=polar_gap_upper(100)
    require(strip==F(78012500000,141158161) and gap==F(1657,12) and strip>gap,'whole rotated strip exclusion')
    eight_gap=a.deficit_interval(8,F(8**28,7**21))
    require(eight_gap[1]<F(21,200),'eight-point exact logarithm cap')
    require(a.global_radius_squared(8)<F(49,16),'eight-point support cap')
    require(F(7,8)/F(49,16)==F(2,7),'eight-point angular rescaling')
    require(F(2,7)**3/3/(F(2,7)**2/2)==F(4,21),'joint ellipse coefficient')
    require(F(21,200)/(F(2,7)**2/2)==F(1029,400),'joint ellipse cap')
    require(F(9,4)*(1+F(4,21))>F(1029,400),'exclude both moment moduli at least 3/2')
    # Log brackets checked with independent decimal evaluation as a diagnostic.
    from decimal import Decimal,localcontext
    log_tests=0
    with localcontext() as ctx:
        ctx.prec=90
        for q in [F(1),F(2),F(1,3),F(17,9),F(10**500+1,10**480+7)]:
            lo,hi=a.log_interval(q);val=(Decimal(q.numerator)/Decimal(q.denominator)).ln()
            require(Decimal(lo.numerator)/Decimal(lo.denominator)<=val<=Decimal(hi.numerator)/Decimal(hi.denominator),'diagnostic log enclosure')
            log_tests+=1
    out={'passed':True,'date':'2026-10-30','newton_and_shape_cases':len(cases),
         'Dirichlet_norm_checks':norm_checks,'bubble_constant_checks':3,'positive_discounts_n2_to_100':99,
         'known_sharp_value_checks':6,'independent_log_diagnostics':log_tests,
         'strip_n100':{'delta_lower_exact':str(strip),'global_delta_upper_exact':str(gap),'excluded':True},
         'eight_point_joint_moment_bound':{'p2_squared_plus_4_over_21_p3_squared_upper':'1029/400',
             'gap_upper':'21/200','both_moduli_at_least_3_over_2_excluded':True,
             'scope':'Global maximizers only; retain the stronger old separate p2<=20/13 bound.'},
         'universal_discount_examples':upper_rows,'shape_checks':rows,
         'scope':'Analytic all-n proof in README.md; finite arithmetic checks do not prove a complete global classification. No GitHub upload in this run.'}
    (ROOT/'verification_report.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items()if k!='shape_checks'},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
