"""Independent exact checks for the radial distribution/energy interface."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from radial_distribution import (root_upper,cdf_error_upper,polar_gap_upper,
    disk_population,annulus_population,upper_support_le,lower_support_le,
    cell_cost,band_certificate,exact_radii_certificate,inner_radius_upper)


def require(ok,name):
    if not ok:raise ArithmeticError(name)


def log_interval(q,terms=48):
    q=F(q)
    if q<=0:raise ValueError('positive logarithm input')
    if q<1:
        lo,hi=log_interval(1/q,terms)
        return -hi,-lo
    if q==1:return F(0),F(0)
    exponent=0
    while q>=2:q/=2;exponent+=1
    def core(x):
        t=(x-1)/(x+1)
        low=2*sum((t**(2*j+1)/(2*j+1) for j in range(terms)),F(0))
        tail=2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
        return low,low+tail
    lo,hi=core(q);l2,h2=core(F(2))
    return lo+exponent*l2,hi+exponent*h2


def discrete_radial_energy(atoms):
    """Independent exact log interval for radial ring probability charges.

    Each atom is (positive squared radius, rational probability).
    This exercises the CDF-to-energy lemma, not the planar optimizer.
    """
    require(all(q>0 and p>0 for q,p in atoms),'positive test atoms')
    require(sum(p for q,p in atoms)==1,'test probability mass')
    stops=sorted(set([F(0),F(2)]+[q for q,p in atoms]))
    low=high=F(0)
    sup=F(0)
    for s in stops:
        target=min(s/2,F(1))
        before=sum(p for q,p in atoms if q<s)
        after=sum(p for q,p in atoms if q<=s)
        sup=max(sup,abs(before-target),abs(after-target))
    for a,b in zip(stops,stops[1:]):
        c=sum(p for q,p in atoms if q<=a)
        if b<=2:
            constant=-c*(b-a)/2+(b*b-a*a)/16
            if c:
                require(a>0,'nonzero CDF at radius zero')
                l,h=log_interval(b/a)
                low+=c*c*l/2+constant;high+=c*c*h/2+constant
            else:
                low+=constant;high+=constant
        else:
            l,h=log_interval(b/a)
            low+=(c-1)**2*l/2;high+=(c-1)**2*h/2
    return low,high,sup


def main():
    root_checks=0
    for a in range(101):
        for b in range(1,8):
            for degree in (2,3,4):
                bound=root_upper(F(a,b),degree,10000)
                require(bound**degree>=F(a,b),'outward rational root')
                if bound:
                    require((bound-F(1,10000))**degree<F(a,b),'minimal grid root')
                root_checks+=1
    root_compare_checks=0
    for Rnum,rnum,anum in product(range(8),range(1,7),range(14)):
        R=F(Rnum,5);r=F(rnum,7);a=F(anum,5)
        require(upper_support_le(R*R,r*r,a)==((R+r)**2<=a),'upper support radical elimination')
        require(lower_support_le(R*R,r*r,a)==(max(R-r,F(0))**2<=a),'lower support radical elimination')
        root_compare_checks+=2
    energy_cases=[]
    for n in range(1,11):
        energy_cases.append([(F(2*j+1,n),F(1,n)) for j in range(n)])
        energy_cases.append([(F(j+1,3*n),F(1,n)) for j in range(n)])
        energy_cases.append([(F(2+j,n),F(j+1,n*(n+1)//2)) for j in range(n)])
    energy_cases += [[(F(1,10),F(1,2)),(F(3),F(1,2))],[(F(4),F(1))]]
    for atoms in energy_cases:
        lo,hi,sup=discrete_radial_energy(atoms)
        require(lo>=sup**3/6,'CDF cubic energy inequality')
        require(lo<=hi,'energy interval ordering')
    # Direct triangle integral, including arbitrary rational epsilon.
    for e in (F(1,100),F(1,7),F(1,2),F(1)):
        h=2*e
        integral=e*e*h-e*h*h/2+h**3/12
        require(integral/4==e**3/6,'exact triangular constant')
    cell_checks=0
    for a,b in ((F(0),F(1,4)),(F(1,7),F(3,5)),(F(1),F(2))):
        for L,U in ((F(0),F(0)),(F(0),F(1,3)),(F(1,4),F(3,4)),(F(1),F(1))):
            cost=cell_cost(a,b,L,U)
            require(cost>=0,'nonnegative radial cell')
            for c in (L,(L+U)/2,U):
                integral=(c*c*(b-a)-c*(b*b-a*a)/2+(b**3-a**3)/12)/(2*b)
                require(cost<=integral,'cell lower bound for every allowed constant CDF')
                cell_checks+=1
    wide_band=band_certificate(100,[(F(49,50),F(1),100)],[F(0),F(7,10)])
    require(wide_band['delta_lower']==F(1225,6),'wide-band exact conditional upper certificate')
    gap100=polar_gap_upper(100)
    require(gap100==F(1657,12),'existing finite polar gap at 100')
    require(wide_band['delta_lower']>gap100,'whole radial region excluded from global max')
    require(band_certificate(100,[(F(49,50),F(99),100)],[F(0),F(7,10)])['delta_lower']==F(1225,6),'no radial upper cap required')
    # Analytic monotonicity starts, not a finite sample replacement for n>=100.
    require(49*100**2-4200*100+2500>0,'all-n support polynomial start')
    require(98*100-4200>0,'support polynomial increasing beyond 100')
    require(F(49,24)>F(1,100)+F(77,60)+F(35,400),'all-n upper/lower gap start')
    examples=[]
    for n in (100,1000,10000,1000000):
        gap=polar_gap_upper(n)
        disk=disk_population(n,F(1),gap)
        annulus=annulus_population(n,F(1,2),F(3,2),gap)
        require(0<=disk['lower']<=disk['upper']<=n,'disk integer population interval')
        require(0<=annulus['lower']<=annulus['upper']<=n,'annulus integer population interval')
        examples.append({'n':n,'gap_upper':str(gap),'cdf_error_upper':str(cdf_error_upper(n,gap)),
                         'disk_s_1':[disk['lower'],disk['upper']],
                         'annulus_half_to_three_halves':[annulus['lower'],annulus['upper']],
                         'minimum_radius_upper':str(inner_radius_upper(n,gap))})
    # Original-scale correction supremum is exactly 1/n at s=2.
    scale_checks=0
    for n in range(2,101):
        a2=F(n-1,n)
        for s in [F(j,10) for j in range(101)]+[F(2)/a2]:
            require(abs(min(a2*s/2,F(1))-min(s/2,F(1)))<=F(1,n),'scale correction 1/n')
            scale_checks+=1
        require(1-a2==F(1,n),'scale correction equality')
    report={'passed':True,'arithmetic':'integer/Fraction and outward logarithm intervals',
            'outward_root_checks':root_checks,'support_radical_checks':root_compare_checks,
            'CDF_energy_distributions':len(energy_cases),'cell_comparisons':cell_checks,
            'original_scale_checks':scale_checks,
            'wide_band_delta_lower':str(wide_band['delta_lower']),'global_gap_upper_n100':str(gap100),
            'wide_band_globally_excluded':True,'examples':examples,
            'scope':'Finite arithmetic checks; radial_audit.md proves the all-n analytic claims. CDF constraints are necessary, not sufficient, and regional certificates do not cover the global domain.'}
    (Path(__file__).parent/'radial_distribution_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
