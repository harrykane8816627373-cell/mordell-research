"""Direct full Gram replay, polynomial integrals and a finite hole exclusion."""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial
import json
import sys
from screened_holes import conditioned_energy,one_hole_dipole
from stationary_moments import cmul,norm,compact_measurements,stationary_energy

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'imported-angular'))
import angular_certificate as imported


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def direct_solve(A,b):
    a=[list(row)+[rhs] for row,rhs in zip(A,b)]
    for j in range(len(a)):
        require(a[j][j]>0,'independent full Gram positive pivot')
        p=a[j][j];a[j]=[v/p for v in a[j]]
        for i in range(len(a)):
            if i==j:continue
            q=a[i][j];a[i]=[v-q*w for v,w in zip(a[i],a[j])]
    return [r[-1] for r in a]


def disk_monomial(a,b,radius):
    if a%2 or b%2:return F(0)
    p,q=a//2,b//2
    return F(factorial(2*p)*factorial(2*q),4**(p+q)*factorial(p)*factorial(q)*factorial(p+q+1))*radius**(a+b+2)


def run():
    integral=0
    for b in (F(1,20),F(1,5),F(2,3)):
        for c in ((F(0),F(0)),(F(1,2),F(1,4))):
            t=b**4
            # Integral/2pi of grad(S-|w|^2) dot grad(b^2-|w-c|^2).
            cross=2*(disk_monomial(2,0,b)+disk_monomial(0,2,b))
            require(cross==t,'independent radial/hole cross integral')
            # For Re[w(S-|w|^2)]: -Delta f = 8 Re w.
            bubble_integral=b*b*disk_monomial(0,0,b)-disk_monomial(2,0,b)-disk_monomial(0,2,b)
            require(4*c[0]*bubble_integral==2*t*c[0],'independent dipole/hole polynomial integral')
            integral+=2
    full=0
    holes=[(F(1,2),F(1,4),F(1,10)),(F(-1,3),F(-1,4),F(1,8))]
    for n in (8,100):
        for K in range(1,5):
            S=F(3);mom=[(F(0),F(0))]+[(F(k,9),F(1-k,11)) for k in range(2,K+1)]
            result=conditioned_energy(n,S,holes,mom)
            a=compact_measurements(n,S,mom)
            diagonal=[S*S]
            rhs=[F(0)]
            for k,z in enumerate(a,1):diagonal.extend([S**(k+2)/F(k+2)]*2);rhs.extend(z)
            M=len(diagonal);J=len(holes)
            G=[[F(0)]*(M+J) for j in range(M+J)]
            for i,d in enumerate(diagonal):G[i][i]=d
            for j,(x,y,b) in enumerate(holes):
                t=b**4;v=[F(1)];p=(F(1),F(0))
                for k in range(1,K+1):p=cmul(p,(x,y));v.extend([(k+1)*p[0],(k+1)*p[1]])
                for i,value in enumerate(v):G[i][M+j]=G[M+j][i]=t*value
                G[M+j][M+j]=t;rhs.append(-F(n,4)*t)
            sol=direct_solve(G,rhs)
            direct=sum((u*v for u,v in zip(rhs,sol)),F(0))+result['harmonic_energy']
            require(direct==result['joint_energy_lower'],'full Gram vs implemented Schur complement')
            no_holes=conditioned_energy(n,S,[],mom)['joint_energy_lower']
            require(no_holes==stationary_energy(n,S,mom)['angular_energy_lower'],'no-hole specialization is stationary mixed energy')
            require(direct>=no_holes and direct>=result['basic_hole_energy'],'joint projection nesting')
            full+=1
    # Real stationary five-point shape and genuinely empty smoothed disks.
    points=[(F(0),F(0)),(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1))]
    b=F(1,20);q=F(2,5)
    for x,y in points:
        d=(x-F(1,2))**2+(y-F(1,2))**2
        # d >= (sqrt(q)+b)^2, using exact squaring.
        margin=d-q-b*b
        require(margin>0 and margin*margin>4*q*b*b,'actual physical smoothed hole clearance')
    actual=conditioned_energy(5,3,[(F(1,2),F(1,2),b)],[(0,0),(0,0),(0,0),(4,0)])
    old=imported.shape_certificate(points)
    require(imported.deficit_interval(5,old['normalized_product'])[0]>actual['joint_energy_lower'],'joint energy below independent true shape deficit')
    n=100;S=F(11881,5000);b=F(19073,100000)
    screened=one_hole_dipole(n,S,(1,0),b)
    radial=conditioned_energy(n,S,[(1,0,b)])
    stored=json.loads((HERE.parent/'symmetry-breaking/witness_verification.json').read_text(encoding='utf-8'))
    gap=F(next(r for r in stored['rows'] if r['n']==100)['global_optimizer_gap_upper'])
    require(screened['basic_hole_energy']<radial['joint_energy_lower']<gap<screened['joint_energy_lower'],'dipole-screened class exceeds actual certified n100 gap')
    # Global support bound is exact at n=100: sqrt(S)=sqrt(2)*(99/100+1/10).
    require(S==2*(F(99,100)+F(1,10))**2,'n100 proven global smeared support')
    invalid=0
    for call in (lambda:conditioned_energy(1,3,[]),lambda:conditioned_energy(8,1,[]),
                 lambda:conditioned_energy(8,3,[(2,0,F(1,10))]),
                 lambda:conditioned_energy(8,3,[(0,0,1),(F(1,2),0,1)]),
                 lambda:conditioned_energy(8,3,[(0,0,0)])):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid joint geometry accepted')
    return {'passed':True,'independent_polynomial_integrals':integral,
        'direct_full_Gram_replays':full,'physical_stationary_shape_joint_energy':actual['joint_energy_lower'],
        'n100_energy_hole_center':['1','0'],'n100_energy_hole_radius':b,
        'n100_point_hole_radius_description':'19073/100000 + sqrt(2/100), Coulomb x scale',
        'n100_basic_hole_energy':screened['basic_hole_energy'],
        'n100_second_moment_screened_energy':radial['joint_energy_lower'],
        'n100_stationary_dipole_screened_energy':screened['joint_energy_lower'],
        'n100_certified_global_gap_upper':gap,'n100_hole_class_excluded':True,
        'invalid_inputs_rejected':invalid,
        'scope':'Joint moment-conditioned energy is stronger than the listed individual projections. The example is not claimed to pass all previous geometric criteria or classify optimizers.'}


if __name__=='__main__':
    result=run()
    (HERE/'screened_holes_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
