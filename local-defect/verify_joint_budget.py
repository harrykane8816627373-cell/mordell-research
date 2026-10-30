"""Nonempty eight-point input and a genuinely combined radial/angular budget."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from radial_hierarchy import radial_data,eight_radial_coefficient

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'stationary-screening'))
sys.path.insert(0,str(HERE.parent/'symmetry-breaking'))
from stationary_moments import low_mode_coefficient,eight_point_certificate
from pair_stability import pair_slack


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    S=F(151,50);gap=F(523,5000)
    cr=eight_radial_coefficient(S);c2=low_mode_coefficient(8,2,S)
    radial_class=cr*(F(52,5)-F(64,7))**2
    angular_class=c2*F(9,10)**2
    require(radial_class<gap and angular_class<gap<radial_class+angular_class,
            'two individually insufficient budgets jointly exclude the whole class')
    points=[(588,0),(-588,0),(0,739),(0,-739),(976,776),(976,-776),(-976,776),(-976,-776)]
    n=8;Q=sum(x*x+y*y for x,y in points);f=F(n,Q)
    q=[f*(x*x+y*y) for x,y in points]
    U=[sum((v**k for v in q),F(0)) for k in range(1,5)]
    p2=f*sum(x*x-y*y for x,y in points)
    p3=0
    p4=f*f*sum(x**4-6*x*x*y*y+y**4 for x,y in points)
    distance=[f*((x-a)**2+(y-b)**2) for i,(x,y) in enumerate(points) for a,b in points[i+1:]]
    pairs=[pair_slack(n,q[i],q[j],f*((points[i][0]-points[j][0])**2+(points[i][1]-points[j][1])**2))
           for i in range(n) for j in range(i+1,n)]
    require(U[0]==8 and U[1]>F(52,5) and p2>F(9,10),'actual input inside excluded joint class')
    require(p2<F(8,7),'actual input passes previous separate p2 bound')
    require(min(distance)>F(4,7) and max(q)<F(7,4) and min(pairs)>0,'listed previous distance/radius/pair conditions pass')
    l=F(7,2)
    fourth_slack=4*n-3*U[1]-p2*p2/n
    sixth_slack=-(8*l**3*U[2]-6*(2*n-3)*l*l*U[1]+(2*n-3)**2*F(n*(n-1),2))
    eighth_slack=-(15*l**4*U[3]-24*(n-2)*l**3*U[2]+9*(n-2)**2*l*l*U[1]-n*l*l*p2*p2+l**4*p4*p4/n)
    require(fourth_slack>0 and sixth_slack>0,'listed previous fourth/sixth scalar cuts pass')
    require(eighth_slack<0,'record stronger prior eighth cut failure honestly')
    actual_radial=radial_data(8,[F(7,8)*v for v in q],2,S)['radial_energy_lower']
    require(actual_radial==cr*(U[1]-F(64,7))**2,'independent radial data vs original fourth-moment formula')
    require(actual_radial+c2*p2*p2>gap,'actual nonempty geometry in globally excluded joint region')
    return {'passed':True,'class_original_fourth_moment_lower':'52/5','class_original_p2_modulus_lower':'9/10',
        'radial_class_energy_lower':radial_class,'angular_class_energy_lower':angular_class,
        'combined_class_energy_lower':radial_class+angular_class,'eight_global_gap_strict_upper':gap,
        'integer_input':points,'normalization_squared':f,'original_fourth_moment':U[1],
        'original_p2':p2,'original_p3':p3,'minimum_squared_distance':min(distance),
        'maximum_squared_radius':max(q),'minimum_previous_pair_slack':min(pairs),
        'previous_fourth_scalar_slack':fourth_slack,'previous_sixth_scalar_slack':sixth_slack,
        'previous_eighth_scalar_slack':eighth_slack,
        'scope':'The whole joint moment class is excluded from global optimizers. The displayed nonstationary input passes only the listed prior tests; its prior eighth-order scalar cut fails, so this is not a comparison against all earlier conditions.'}


if __name__=='__main__':
    out=run()
    (HERE/'joint_budget_verification.json').write_text(json.dumps(out,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,default=str,ensure_ascii=False,indent=2))
