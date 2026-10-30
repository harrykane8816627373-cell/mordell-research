"""Independent radial integration, Schur formula and nested exact projections."""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial
import json
from radial_hierarchy import radial_gram,radial_data,eight_radial_coefficient

HERE=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def monomial(a,b,r2):
    if a%2 or b%2:return F(0)
    p,q=a//2,b//2
    return F(factorial(2*p)*factorial(2*q),4**(p+q)*factorial(p)*factorial(q)*factorial(p+q+1))*r2**(p+q)


def run():
    gram=disks=nested=0
    for degree in range(1,8):
        for S in (F(2),F(151,50),F(4)):
            H=radial_gram(degree,S)
            for j in range(1,degree+1):
                for k in range(1,degree+1):
                    # Expand derivatives in powers of q, integrate 2q dq.
                    direct=F(0)
                    for a in range(j):
                        for b in range(k):
                            direct+=2*j*k*(-1)**(a+b)*comb(j-1,a)*comb(k-1,b)*S**(j+k-2-a-b)*S**(a+b+2)/F(a+b+2)
                    require(H[j-1][k-1]==direct,'independent expanded radial Gram integral')
                    gram+=1
    for m in range(1,7):
        for x in (F(0),F(1,3),F(2)):
            r2=F(1,5)
            # Expand ((x+u)^2+v^2)^m in Cartesian monomials.
            direct=sum((comb(m,j)*comb(2*(m-j),a)*x**(2*(m-j)-a)*monomial(a,2*j,r2)
                        for j in range(m+1) for a in range(2*(m-j)+1)),F(0))
            closed=sum((F(comb(m,l)**2,l+1)*r2**l*(x*x)**(m-l) for l in range(m+1)),F(0))
            require(direct==closed,'independent shifted-disk radial moment')
            disks+=1
    examples=[]
    for n,q in ((2,[F(1,2)]*2),(5,[F(0)]+[F(1)]*4),(8,[F(0)]+[F(1)]*7),
                (10,[F(9*j,45) for j in range(10)])):
        S=F(4);prev=F(0)
        for d in range(1,7):
            result=radial_data(n,q,d,S)
            require(result['radial_energy_lower']>=prev,'radial projection nesting')
            if d==2:
                difference=result['smeared_moment_sums'][2]-F(4*n,3)
                require(result['radial_energy_lower']==F(9,2)*difference*difference/S**4,'matched moment Schur coefficient')
            prev=result['radial_energy_lower'];nested+=1
        examples.append({'n':n,'degree':6,'radial_energy_lower':prev})
    # Eight-point center+regular-seven radii make the fourth radial moment defect zero.
    target=radial_data(8,[F(0)]+[F(1)]*7,2,F(151,50))
    require(target['radial_energy_lower']==0,'eight-point target fourth moment cancellation')
    invalid=0
    for call in (lambda:radial_gram(0,3),lambda:radial_gram(13,3),
                 lambda:radial_data(8,[F(1)]*8,2,4),lambda:radial_data(8,[F(0)]+[F(1)]*7,2,2)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid radial hierarchy input accepted')
    return {'passed':True,'independent_expanded_Gram_checks':gram,
        'independent_Cartesian_disk_moment_checks':disks,'nested_projection_checks':nested,
        'radial_examples':examples,'eight_original_fourth_moment_coefficient':eight_radial_coefficient(),
        'invalid_inputs_rejected':invalid,'scope':'Radial moment projection, not an extra copy of full radial energy. General formula and orthogonality are proved in RADIAL_HIERARCHY.md.'}


if __name__=='__main__':
    out=run()
    (HERE/'radial_hierarchy_verification.json').write_text(json.dumps(out,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,default=str,ensure_ascii=False,indent=2))
