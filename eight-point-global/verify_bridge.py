"""Exact algebra checks for the existing-subset/global-pair-kernel bridge."""
from fractions import Fraction as Q
from itertools import combinations
import json

def need(ok,label):
    if not ok:raise ArithmeticError(label)

def z(a=0,b=0):return Q(a),Q(b)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,s):return a[0]*s,a[1]*s
def conj(a):return a[0],-a[1]
def norm(a):return a[0]**2+a[1]**2
def total(values):
    out=z()
    for value in values:out=add(out,value)
    return out
def power(a,k):
    result=z(1)
    for _ in range(k):result=mul(result,a)
    return result
def bilinear(a,M,b):return sum(a[i]*M[i][j]*b[j] for i in range(4) for j in range(4))
def features(u,k):return [u-Q(7,8),u*(u-1),u*u*(u-1),u**3*(u-1)] if k==0 else [1,u-1,(u-1)**2,(u-1)**3]

def check():
    need(Q(49,32)<Q(965267,524288) and Q(1,2)>Q(2,7),'existing strict-subset bounds imply the entire required pair domain')
    need((49-Q(147,4))/8==Q(49,32),'seven-point subset gives the existing radius bound')
    need(Q(2**2*(2-1),8)==Q(1,2),'two-point subset gives the existing separation')
    points=[z(),z(1),z(0,1),z(-1),z(0,-1),z(2,1),z(-2,1),z(0,-2)]
    need(total(points)==z(),'the independent rational example is centred')
    subsets=0;c=z(Q(2,7),Q(-3,11))
    for m in range(2,8):
        for members in combinations(range(8),m):
            W=sum(norm(add(points[i],scale(points[j],-1))) for i,j in combinations(members,2))
            need(W==m*sum(norm(points[i]) for i in members)-norm(total(points[i] for i in members)),
                'exact common-centroid subset identity')
            need(W==m*sum(norm(add(points[i],scale(c,-1))) for i in members)
                -norm(add(total(points[i] for i in members),scale(c,-m))),'exact translated subset identity')
            subsets+=1
    u=[norm(p) for p in points];w=[1,-1,2,0]
    matrices=[[[Q((k+2+i) if i==j else 0)+Q(w[i]*w[j],7) for j in range(4)] for i in range(4)] for k in range(7)]
    def kernel(i,j):
        t=mul(points[i],conj(points[j]))[0];H=[Q(1),t]
        cross=mul(points[i],conj(points[j]))[1]
        need(u[i]*u[j]-t*t==cross*cross,'the actual two-point domain obeys exact Cauchy identity')
        for k in range(2,7):H.append(2*t*H[-1]-u[i]*u[j]*H[-2])
        for k in range(7):need(H[k]==power(mul(points[i],conj(points[j])),k)[0],'exact angular recurrence')
        return sum(H[k]*bilinear(features(u[i],k),matrices[k],features(u[j],k)) for k in range(7))
    K=[[kernel(i,j) for j in range(8)] for i in range(8)]
    radial=[sum(features(t,0)[i] for t in u) for i in range(4)]
    moment_energy=bilinear(radial,matrices[0],radial)
    for k in range(1,7):
        L=[total(scale(power(points[j],k),features(u[j],k)[i]) for j in range(8)) for i in range(4)]
        value=total(scale(mul(conj(L[i]),L[j]),matrices[k][i][j]) for i in range(4) for j in range(4))
        need(value[1]==0 and value[0]>=0,'Hermitian moment form is real and nonnegative')
        moment_energy+=value[0]
    need(sum(map(sum,K))==moment_energy,'full double sum equals the moment quadratic forms')
    b=Q(5,11);p=[(4*u[i]+b-K[i][i]/2)/7 for i in range(8)]
    pair_sum=sum(p[i]+p[j]-K[i][j] for i,j in combinations(range(8),2))
    need(pair_sum==4*sum(u)+8*b-moment_energy/2,'diagonal correction gives the exact pair-sum identity')
    # u*(u-1)^2 = (u^3-u^2) - (u^2-u).
    for t in [Q(0),Q(1),Q(2,3),Q(7,5)]:
        need(t*(t-1)**2==(t**3-t*t)-(t*t-t) and t*(t-1)**2>=0,'radial equality has exactly the required nonnegative form')
    e=[Q(1)];powers=[Q(0)]*7+[Q(7)]
    for k in range(1,8):e.append(sum((-1)**(j-1)*e[k-j]*powers[j] for j in range(1,k+1))/k)
    need(e==[Q(1)]+[Q(0)]*6+[Q(1)],'Newton identities recover the unit heptagon polynomial')
    need(Q(8,7)**28*7**7==Q(8**28,7**21),'normalizations match exactly')
    return {'passed':True,'subset_identities_checked':subsets,'kernel_sum_and_complex_moment_identity_checked':True,
        'domain_inclusion_checked':True,'equality_radial_and_Newton_identities_checked':True,
        'seven_point_classification_required_by_bridge':False,'finite_examples_replace_analytic_proof':False}

if __name__=='__main__':print(json.dumps(check(),indent=2))
