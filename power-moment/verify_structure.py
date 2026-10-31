"""Lacunary Ward reduction and finite power-sum symmetry checks."""
from fractions import Fraction as F
from acceptance import require,write_report
from complex_interval import C
from power_coordinates import coefficients


def add(a,b):
    out=dict(a)
    for key,value in b.items():out[key]=out.get(key,F(0))+value
    return {k:v for k,v in out.items() if v}
def mul(a,b):
    out={}
    for u,x in a.items():
        for v,y in b.items():
            key=tuple(i+j for i,j in zip(u,v));out[key]=out.get(key,F(0))+x*y
    return {k:v for k,v in out.items() if v}
def scalar(a,x):return {k:v*x for k,v in a.items() if v*x}


def run():
    # Independent exact expansion of z*F-n(1-t)*P=G in Q[z,t,a,b].
    one={(0,0,0,0):F(1)};z={(1,0,0,0):F(1)};t={(0,1,0,0):F(1)}
    a={(0,0,1,0):F(1)};b={(0,0,0,1):F(1)}
    identity_checks=0
    for n in range(4,65):
        P=add({(n,0,0,0):F(1)},add(mul(a,z),b))
        residual=add(scalar(mul(add(one,scalar(t,-1)),{(n-1,0,0,0):F(1)}),n),scalar(mul(a,t),-1))
        G=add(mul(mul(a,z),add(scalar(t,n-1),scalar(one,-n))),scalar(mul(b,add(t,scalar(one,-1))),n))
        require(add(mul(z,residual),scalar(mul(add(one,scalar(t,-1)),P),-n))==G,'lacunary physical Ward polynomial identity')
        identity_checks+=1
    regular=[];center=[]
    for n in range(8,101):
        ratio=F(2*(n-4),n-1)
        require(ratio-1==F(n-7,n-1)>0,'regular n-gon fourth-power instability')
        if n in (8,10,11,100):regular.append({'n':n,'norm_ratio_lower':ratio})
    for n in range(11,101):
        ratio=3-F(20,n)
        require(ratio-1==F(2*(n-10),n)>0,'centered polygon fifth-power instability')
        if n in (11,12,100):center.append({'n':n,'norm_ratio_lower':ratio})
    symmetry=0
    # Formal power vectors with the required zeros; Newton must preserve them.
    for n in range(2,13):
        for m in range(2,n+1):
            p=[C(n)]+[C(k+1,2-k) if k%m==0 else C(0) for k in range(1,n+1)]
            c=coefficients(p)
            require(all(not c[k] for k in range(1,n+1) if k%m),'finite symmetry zeros propagate to every coefficient')
            symmetry+=1
    return {'passed':True,'lacunary_identity_checks':identity_checks,'cyclic_symmetry_Newton_checks':symmetry,
        'regular_polygon_instability_examples':regular,'centered_polygon_instability_examples':center,
        'general_n_structure':'For n>=4, any normalized distinct stationary configuration with p2=...=p_(n-2)=0 is a regular n-gon or origin plus regular (n-1)-gon. For n>=11 neither can be a local maximum.',
        'scope':'General theorems are proved analytically. These exact checks validate the polynomial reduction and constants; no general p2 or p3 vanishing assertion is made.'}


if __name__=='__main__':write_report('structure_verification.json',run())
