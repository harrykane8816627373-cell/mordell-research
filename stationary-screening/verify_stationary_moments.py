"""Independent pair symmetrization, stationary examples and exact Gram checks."""
from fractions import Fraction as F
from pathlib import Path
import json
import random
import sys
from stationary_moments import (cadd,csub,cmul,cscale,norm,convolution,
    compact_measurements,stationary_energy,low_mode_coefficient,
    residual_mode_lower,eight_point_certificate,eliminated_convolution_coefficient)

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'radial-angular'))
sys.path.insert(0,str(HERE.parent/'distribution-geometry'))
from mixed_energy import projection_value
from microscopic_covering import sqrt_interval


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def divide(z,w):return ((z[0]*w[0]+z[1]*w[1])/norm(w),(z[1]*w[0]-z[0]*w[1])/norm(w))


def powers(points,K):
    out=[];current=[(F(1),F(0)) for p in points]
    for k in range(K):
        current=[cmul(p,v) for p,v in zip(points,current)]
        out.append(tuple(sum((v[a] for v in current),F(0)) for a in range(2)))
    return out


def residuals(points):
    n=len(points);out=[]
    for i,z in enumerate(points):
        total=(F(0),F(0))
        for j,w in enumerate(points):
            if i!=j:total=cadd(total,divide((F(1),F(0)),csub(z,w)))
        out.append(csub(total,(F(n,2)*z[0],-F(n,2)*z[1])))
    return out


def run():
    rng=random.Random(20261030)
    identity=robust=0
    for n in range(2,9):
        points=[(F(j-2),F(rng.randrange(-9,10),7)) for j in range(n)]
        mean=tuple(sum((p[a] for p in points),F(0))/n for a in range(2))
        points=[csub(p,mean) for p in points]
        p=powers(points,9);R=residuals(points)
        for k in range(1,9):
            Q=tuple(sum((norm(z)*v[a] for z,v in zip(points,powers_each(points,k))),F(0)) for a in range(2))
            weighted=tuple(sum((cmul(v,r)[a] for v,r in zip(powers_each(points,k+1),R)),F(0)) for a in range(2))
            rhs=csub(cadd(cscale(p[k-1],F(2*n-k-1,n)),cscale(convolution(p,k),F(1,n))),cscale(weighted,F(2,n)))
            require(Q==rhs,'independent pair-symmetrized mixed closure with residual')
            identity+=1
            if k in (2,3):
                plo=sqrt_interval(norm(p[k-1]))[0]
                eup=sqrt_interval(norm(weighted))[1]
                S=F(7,2)
                A1=cadd(Q,cscale(p[k-1],F(k+1,n)))
                actual=norm(p[k-1])/(k*S**k)+(k+2)*norm(csub(A1,cscale(p[k-1],S)))/S**(k+2)
                require(residual_mode_lower(k,S,plo,eup,n)<=actual,'near-stationary reverse-triangle bound')
                robust+=1
    stationary=0
    for t in (F(0),F(1,5),F(1,2),F(2,3),F(1),F(3,2)):
        rotation=((1-t*t)/(1+t*t),2*t/(1+t*t))
        for base in ([(F(1,2),F(1,2)),(F(-1,2),F(-1,2))],
                     [(F(0),F(0)),(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1))]):
            pts=[cmul(z,rotation) for z in base];n=len(pts)
            require(sum((norm(z) for z in pts),F(0))==n-1,'known Coulomb normalization')
            require(all(r==(0,0) for r in residuals(pts)),'known rational stationary example')
            p=powers(pts,12);a=compact_measurements(n,5,p)
            for k,v in enumerate(a,1):
                Q=tuple(sum((norm(z)*v0[j] for z,v0 in zip(pts,powers_each(pts,k))),F(0)) for j in range(2))
                direct=csub(cscale(p[k-1],5),cadd(Q,cscale(p[k-1],F(k+1,n))))
                require(v==direct,'actual disk-averaged compact stationary measurement')
                stationary+=1
    gram=monotone=0
    for n in (5,8,13):
        p=[(F(0),F(0))]+[(F(k-3,5),F(2-k,7)) for k in range(2,9)]
        for S in (F(2),F(151,50),F(7,2)):
            result=stationary_energy(n,S,p)
            for k,row in enumerate(result['rows'],1):
                A1=cadd(cscale(p[k-1],2),cscale(convolution(p,k),F(1,n)))
                expected=projection_value(k,S,[p[k-1][0],A1[0]],[p[k-1][1],A1[1]])
                require(expected==row['holomorphic_energy']+row['stationary_extra'],'independent old Gram projection vs new closure')
                gram+=1
    # Exact derivative numerator is a negative square for h=1, k=2,3.
    for k in (2,3):
        for S in (F(3,2),F(2),F(151,50),F(4)):
            derivative=-S*S+(k+2)*(2*S*(S-2)-(k+2)*(S-2)**2)
            require(derivative==-((k+1)*S-2*(k+2))**2,'support coefficient monotonicity')
            monotone+=1
    eliminated=0
    for k in range(2,9):
        for S in (F(2),F(151,50),F(7,2)):
            a,b,c=1/(k*S**k),F(k+2)/S**(k+2),2-S
            conv=(F(3,7),F(-2,5));n=8
            optimum=cscale(conv,-b*c/(n*(a+b*c*c)))
            minimum=a*norm(optimum)+b*norm(cadd(cscale(optimum,c),cscale(conv,F(1,n))))
            original_scale=F(n-1,n)**k
            require(original_scale*minimum==eliminated_convolution_coefficient(n,k,S)*norm(conv),'independent quadratic elimination of higher moment')
            eliminated+=1
    eight=eight_point_certificate()
    require(F(8,7)<F(20,13),'new separate p2 condition improves the supplied old bound')
    require(F(6,5)**2<F(1029,400) and F(6,5)>F(8,7),'strict old/new moment-region comparison, not a configuration witness')
    invalid=0
    for call in (lambda:stationary_energy(1,3,[(0,0)]),lambda:stationary_energy(8,0,[(0,0)]),
                 lambda:stationary_energy(8,3,[(1,0)]),lambda:low_mode_coefficient(8,4,3),
                 lambda:residual_mode_lower(2,3,1,-1,8)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid stationary certificate accepted')
    return {'passed':True,'pair_symmetrization_with_residual_checks':identity,
        'robust_residual_checks':robust,'known_stationary_disk_moment_checks':stationary,
        'independent_old_Gram_checks':gram,'coefficient_derivative_checks':monotone,
        'eliminated_higher_moment_checks':eliminated,
        'invalid_inputs_rejected':invalid,'eight_point_certificate':eight,
        'scope':'Exact identities, implementation and eight-point necessary moment region; no stationary configuration enumeration or full sharp proof.'}


def powers_each(points,k):
    out=[]
    for z in points:
        v=(F(1),F(0))
        for j in range(k):v=cmul(v,z)
        out.append(v)
    return out


if __name__=='__main__':
    result=run()
    (HERE/'stationary_moments_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
