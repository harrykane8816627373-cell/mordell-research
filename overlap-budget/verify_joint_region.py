"""Whole low-moment/short-edge region and a nonempty four-parameter box."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from overlap_cost import pair_cost_lower,short_edge_capacity

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'stationary-screening'))
sys.path.insert(0,str(HERE.parent/'local-defect'))
sys.path.insert(0,str(HERE.parent/'symmetry-breaking'))
from stationary_moments import low_mode_coefficient,eight_point_certificate
from radial_hierarchy import eight_radial_coefficient
from pair_stability import pair_slack


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)
def scale(a,q):return mul(a,(F(q),F(q)))
def square(a):
    v=[x*x for x in a]
    return (F(0) if a[0]<=0<=a[1] else min(v)),max(v)


def run():
    S=F(151,50);gap=F(523,5000);cr=eight_radial_coefficient(S)
    c2=low_mode_coefficient(8,2,S);d4=eight_point_certificate()['p2_fourth_power_coefficient']
    radial=cr*(F(51,5)-F(64,7))**2;angular=c2*F(4,5)**2
    simple_pair=F(3,500)
    require(pair_cost_lower(F(8,5))>simple_pair,'simple per-edge lower cost')
    require(radial+angular<gap<radial+angular+3*simple_pair,'new three-way whole-class exclusion')
    cap=short_edge_capacity(8,gap,dirichlet_lower=radial+angular)
    require(cap['maximum_short_edges']==2,'conditional edge capacity after joint Dirichlet budget')
    i,j,k,l=616,760,953,780
    points=[(i,0),(-i,0),(0,j),(0,-j),(k,l),(k,-l),(-k,l),(-k,-l)]
    Q=sum(x*x+y*y for x,y in points);f=F(8,Q)
    q=[f*(x*x+y*y) for x,y in points];U=[sum((v**p for v in q),F(0)) for p in range(1,5)]
    p2=f*sum(x*x-y*y for x,y in points);p4=f*f*sum(x**4-6*x*x*y*y+y**4 for x,y in points)
    distances={(a,b):f*((points[a][0]-points[b][0])**2+(points[a][1]-points[b][1])**2)
               for a in range(8) for b in range(a+1,8)}
    edges=sum(d<=F(128,175) for d in distances.values())
    previous=cr*(U[1]-F(64,7))**2+c2*p2*p2+d4*p2**4
    require(U[1]>F(51,5) and p2>F(4,5) and edges==4,'actual normalized input in new joint class')
    require(previous<gap<previous+simple_pair*edges,'input passes displayed previous low-order joint region but fails overlap-aware budget')
    pair_min=min(pair_slack(8,q[a],q[b],d) for (a,b),d in distances.items())
    require(min(distances.values())>F(4,7) and max(q)<F(7,4) and pair_min>0,'listed inherited pair/radius tests')
    lam=F(7,2)
    fourth=32-3*U[1]-p2*p2/8
    sixth=-(8*lam**3*U[2]-78*lam*lam*U[1]+169*28)
    eighth=-(15*lam**4*U[3]-144*lam**3*U[2]+324*lam*lam*U[1]-8*lam*lam*p2*p2+lam**4*p4*p4/8)
    require(fourth>0 and sixth>0,'listed lower scalar cuts pass')
    require(eighth<0,'record prior eighth scalar failure')
    # Full independent +/-1/100 intervals in four reflection-preserving raw parameters.
    intervals=[(F(v)-F(1,100),F(v)+F(1,100)) for v in (i,j,k,l)]
    I,J,K,L=intervals;I2,J2,K2,L2=map(square,intervals)
    rawQ=add(scale(add(I2,J2),2),scale(add(K2,L2),4))
    fbox=(8/rawQ[1],8/rawQ[0])
    rawU=add(scale(add(square(I2),square(J2)),2),scale(square(add(K2,L2)),4))
    ubox=mul(square(fbox),rawU)
    pbox=mul(fbox,add(scale(sub(I2,J2),2),scale(sub(K2,L2),4)))
    edgebox=mul(fbox,add(square(sub(K,I)),L2))
    require(ubox[0]>F(51,5) and pbox[0]>F(4,5) and edgebox[1]<F(128,175),'entire four-parameter box lies in new class')
    old_upper=cr*(ubox[1]-F(64,7))**2+c2*pbox[1]**2+d4*pbox[1]**4
    require(old_upper<gap,'entire box passes displayed previous moment budget')
    return {'passed':True,'original_fourth_moment_class_lower':'51/5','original_p2_modulus_class_lower':'4/5',
        'force_short_edge_cutoff':'8/5','minimum_excluded_edge_count':3,
        'Dirichlet_class_lower':radial+angular,'simple_pair_cost':simple_pair,
        'combined_class_lower':radial+angular+3*simple_pair,'eight_global_gap_upper':gap,
        'conditional_maximum_short_edges':2,'integer_input':points,'normalization_squared':f,
        'input_fourth_moment':U[1],'input_p2':p2,'input_short_edge_count':edges,
        'input_previous_displayed_budget':previous,'minimum_previous_pair_slack':pair_min,
        'previous_fourth_scalar_slack':fourth,'previous_sixth_scalar_slack':sixth,
        'previous_eighth_scalar_slack':eighth,
        'four_parameter_box_radius':'1/100','box_fourth_moment_interval':ubox,
        'box_p2_interval':pbox,'box_short_edge_squared_interval':edgebox,
        'box_previous_displayed_budget_upper':old_upper,
        'scope':'An entire reflection-preserving four-parameter box passes the displayed prior low-order budget and is excluded by retained pair loss. Its sample fails the earlier eighth scalar cut, so no claim is made against all prior conditions.'}


if __name__=='__main__':
    result=run()
    (HERE/'joint_region_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
