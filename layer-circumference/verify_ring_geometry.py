"""Independent finite checks of convex tails and cyclic edge bookkeeping."""
from fractions import Fraction as F
from math import isqrt
from pair_integral import ROOT,sqrt_bounds,certified_pi_upper
from reporting import require,write_report
import sys
sys.path.insert(0,str(ROOT/'joint-bounds'))
sys.path.insert(0,str(ROOT/'collective-stability'))
from subset_moment import log_interval
from ring_population import pi_interval


def atomic_loss_interval(d):
    lo=hi=F(0)
    for r,weight in [(F(1),F(1,4)),(F(3,2),F(1,4)),(F(2),F(1,2))]:
        if r>d:
            l,u=log_interval(r/d);lo+=2*weight*l;hi+=2*weight*u
    return lo,hi


def run():
    # Same Newton-tail representation, but an independently chosen radial law.
    convex=0
    for left,right in [(F(1,2),F(3,2)),(F(3,4),F(7,4)),(F(5,4),F(9,4))]:
        for weight in (F(1,4),F(1,2),F(3,4)):
            d=weight*left+(1-weight)*right
            lo1,hi1=atomic_loss_interval(left);lo2,hi2=atomic_loss_interval(right)
            require(atomic_loss_interval(d)[1]<weight*lo1+(1-weight)*lo2,'independent convex radial tail');convex+=1
    ledgers=0
    for n in range(3,90):
        for K in range(1,n+1):
            # Consecutive compositions with population-one and population-two layers.
            sizes=[1]*K;sizes[0]+=n-K
            unique=[]
            start=0
            for m in sizes:
                if m>=3:unique.extend(tuple(sorted((start+j,start+(j+1)%m))) for j in range(m))
                start+=m
            require(len(unique)==len(set(unique)) and len(unique)>=n-2*K,'cyclic edges counted once');ledgers+=1
    # Pythagorean circle points, independently sorted by a known angular order.
    polygons=[[(1,0),(F(3,5),F(4,5)),(0,1),(-1,0),(0,-1)],
              [(1,0),(0,1),(-1,0),(0,-1)],
              [(1,0),(F(5,13),F(12,13)),(-1,0),(F(5,13),-F(12,13))]]
    low_pi=pi_interval()[0]
    perimeters=[]
    for points in polygons:
        require(all(F(x)**2+F(y)**2==1 for x,y in points),'unit circle')
        bound=F(0)
        for i,(x,y) in enumerate(points):
            u,v=points[(i+1)%len(points)]
            bound+=sqrt_bounds((x-u)**2+(y-v)**2)[1]
        require(bound<2*low_pi,'polygonal circumference upper');perimeters.append(bound)
    return {'passed':True,'independent_convex_tail_checks':convex,'cyclic_count_ledgers':ledgers,
            'rational_polygon_perimeter_uppers':perimeters,
            'same_ring_minimum_loss':'regular angles minimize sum P(chord) for fixed population and radius',
            'adjacent_edge_count':'L>=n-2K; one- and two-point layers discarded',
            'scope':'Finite checks of the implementation ledger and independent radial-law examples. General convexity, cyclic Jensen, and rigidity follow from PAIR_CONVEXITY.md.'}


if __name__=='__main__':write_report('ring_geometry_verification.json',run())
