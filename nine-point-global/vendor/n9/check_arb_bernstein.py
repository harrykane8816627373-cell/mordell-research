"""Independent direct-evaluation checks of tensor indexing and basis conversion."""
from pathlib import Path
import sys,math,json,time
ROOT=Path(__file__).resolve().parent
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
from flint import arb,ctx,fmpq
ctx.prec=192
from arb_bernstein import coefficients,residual_coefficients
from arb_bernstein import cubic_log_upper_power
from cubic_log_bernstein import CubicLogUpper
start=time.monotonic()
def direct(poly,point):return sum((c*point[0]**e[0]*point[1]**e[1]*point[2]**e[2]for e,c in poly.items()),arb(0))
def bernstein(arr,degrees,point):
    answer=arb(0)
    for i in range(degrees[0]+1):
        for j in range(degrees[1]+1):
            for k in range(degrees[2]+1):
                v=arr[i,j,k]
                for index,n,t in zip((i,j,k),degrees,point):v*=math.comb(n,index)*t**index*(1-t)**(n-index)
                answer+=v
    return answer
polys=[{(0,0,0):arb(2),(1,0,0):arb(-3),(0,1,0):arb(5),(0,0,1):arb(7),
        (2,1,3):arb(11),(3,2,1):arb(-13)}]
raw=json.loads((ROOT/'global_pair_bound_exact_contacts.json').read_text());a=arb(2).sqrt();logs=[(2-a).log(),arb(2).log()]
def fld(co):
    out=arb(0)
    for q in co:out=out*a+arb(fmpq(q))
    return out
polys.append({tuple(e['powers']):fld(e['value'][0])+sum((fld(q)*l for q,l in zip(e['value'][1:],logs)),arb(0))for e in raw['polynomial_terms']})
checks=0
for poly in polys:
    degrees=tuple(max(3,max(e[k]for e in poly))for k in range(3))
    for box in [(0.,1.,0.,2.,-1.,1.),(.999999,1.000001,.999998,1.000002,.707105,.707108)]:
        arr=coefficients(poly,degrees,box);res=residual_coefficients(poly,degrees,box,order=1)
        cubic=residual_coefficients(poly,degrees,box,order=3)
        upper=coefficients(cubic_log_upper_power(box),degrees,(0.,1.,0.,1.,0.,1.))
        floated=CubicLogUpper(degrees).upper(box)
        assert all(arb(float(v))>=b for v,b in zip(floated.flat,upper.flat))
        center=[(box[2*k]+box[2*k+1])/2 for k in range(3)];d0=arb(max(center[0]+center[1]-2*center[2],.25))
        for fractions in [(0,0,0),(1,1,1),(fmpq(1,3),fmpq(2,5),fmpq(3,7))]:
            pp=list(map(arb,fractions));xx=[arb(box[2*k])+(arb(box[2*k+1])-arb(box[2*k]))*pp[k]for k in range(3)]
            value=direct(poly,xx)
            assert (bernstein(arr,degrees,pp)-value).contains(0)
            tangent=d0.log()-1+(xx[0]+xx[1]-2*xx[2])/d0
            assert (bernstein(res,degrees,pp)-(value-tangent)).contains(0)
            delta=(xx[0]+xx[1]-2*xx[2]-d0)/d0
            taylor=d0.log()+delta-delta*delta/2+delta*delta*delta/3
            assert (bernstein(cubic,degrees,pp)-(value-taylor)).contains(0)
            checks+=3
print('Tensor basis and tangent checks passed:',checks,'seconds',time.monotonic()-start)
