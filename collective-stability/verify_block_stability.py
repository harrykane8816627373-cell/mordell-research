"""Independent exact checks and two entire-coordinate-box exclusions."""
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import json
import sys

from block_stability import *

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'symmetry-breaking'))
from pair_stability import pair_slack


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def direct_determinant(matrix,d=1):
    n=len(matrix)
    total=ZERO
    for order in permutations(range(n)):
        inversions=sum(order[i]>order[j] for i in range(n) for j in range(i+1,n))
        term=ONE
        for i,j in enumerate(order):
            term=mul(term,matrix[i][j],d)
        total=add(total,scale(term,(-1)**inversions))
    return total


def normalized_pairs(points,d=1):
    data=geometry(points,d)
    n,Q=data['n'],data['Q']
    factor=F(n)/Q
    radii=[factor*norm(p,d) for p in data['points']]
    rows=[]
    for i,j in combinations(range(n),2):
        s=factor*norm(sub(data['points'][i],data['points'][j]),d)
        slack=pair_slack(n,radii[i],radii[j],s)
        cij=residual_cross(data,i,j)
        hii=data['lambda']**2*projection(data,i,i)[0]*projection(data,j,j)[0]-norm(cij,d)
        require(hii==data['lambda']*slack/(n*s*s),'old scalar slack equals Schur diagonal')
        rows.append((i,j,s,slack))
    return {'all_old_pairs_pass':all(s>F(2*n,(n-1)**2) and slack>0 for i,j,s,slack in rows),
            'minimum_distance_squared':min(s for i,j,s,t in rows),
            'maximum_radius_squared':max(radii),
            'minimum_pair_slack':min(t for i,j,s,t in rows),'pair_count':len(rows)}


def full_stability(points,d=1):
    """Build original B directly; verify its identities without using blocks."""
    data=geometry(points,d)
    n,Q,l=data['n'],data['Q'],data['lambda']
    v=data['points']
    B=[[ZERO for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                continue
            w=sub(v[i],v[j])
            B[i][j]=scale(divide(ONE,mul(w,w,d),d),-Q/n)
        B[i][i]=neg(sum_complex(B[i]))
    for i in range(n):
        require(sum_complex(B[i])==ZERO,'direct original B annihilates constants')
        require(sum_complex(mul(B[i][j],v[j],d) for j in range(n))==scale(conj(v[i]),l),
                'direct original B satisfies saturated stationary equation')
    A=[[sub((l*l if i==j else F(0),F(0)),
             sum_complex(mul(conj(B[k][i]),B[k][j],d) for k in range(n)))
        for j in range(n)] for i in range(n)]
    require(psd_certificate(A,d)['passes'],'direct original complete stability norm')
    # Check the saturated decomposition by explicit matrix multiplication.
    for i in range(n):
        for j in range(n):
            projected=sum_complex(mul(conj(projection(data,i,a)),
                mul(B[a][b],projection(data,b,j),d),d)
                for a in range(n) for b in range(n))
            rank=scale(mul(conj(v[i]),conj(v[j]),d),l/Q)
            require(B[i][j]==add(rank,projected),'independent saturated projection identity')
    return data


def boxes(points,epsilon=F(1,100)):
    return [((F(a)-epsilon,F(a)+epsilon),(F(b)-epsilon,F(b)+epsilon)) for a,b in points]


def concise(certificate):
    return {k:v for k,v in certificate.items() if k!='matrix'}


def run():
    determinant_checks=0
    for d in (1,3):
        for n in range(1,5):
            matrix=[[(F((i+1)*(j+2)+int(i==j)),F(i-j)) for j in range(n)] for i in range(n)]
            require(determinant(matrix,d)==direct_determinant(matrix,d),'Gaussian determinant versus independent permutation expansion')
            determinant_checks+=1
    stable_examples=[(1,[(1,0),(-1,0)]),
                     (1,[(1,0),(0,1),(-1,0),(0,-1)]),
                     (1,[(0,0),(1,0),(0,1),(-1,0),(0,-1)]),
                     (3,[(0,0),(1,0),(F(1,2),F(1,2)),(F(-1,2),F(1,2)),
                         (-1,0),(F(-1,2),F(-1,2)),(F(1,2),F(-1,2))])]
    triangle_pass_checks=0
    block_pass_checks=0
    scalar_identity_checks=0
    for d,points in stable_examples:
        data=full_stability(points,d)
        scalar_identity_checks+=normalized_pairs(points,d)['pair_count']
        n=len(points)
        for i in range(n):
            for j,k in combinations((a for a in range(n) if a!=i),2):
                require(triangle_certificate(points,i,j,k,d)['passes'],'known stable point passes every three-point Schur test')
                triangle_pass_checks+=1
        if n>=4:
            rows,columns=(0,1),tuple(range(2,min(n,5)))
            require(psd_certificate(cross_block_matrix(data,rows,columns),d)['passes'],'known stable point passes collective cross block')
            block_pass_checks+=1

    triangle_points=[(7,0),(5,6),(0,8),(-4,4),(-7,1),(-3,-6),(1,-8),(4,-3)]
    pair_info=normalized_pairs(triangle_points)
    scalar_identity_checks+=pair_info['pair_count']
    require(pair_info['all_old_pairs_pass'],'all old local pairs pass three-point example')
    triangle=triangle_certificate(triangle_points,3,2,4)
    require(not triangle['passes'] and triangle['schur_determinant']<0,'new joint determinant excludes three-point example')
    require(all(triangle['schur_matrix'][i][i][0]>0 for i in range(2)),
            'three-point failure is not a repeated scalar failure')
    tboxes=boxes(triangle_points)
    old_triangle_box=all_pair_box_certificate(tboxes)
    triangle_box=triangle_box_certificate(tboxes,3,2,4)
    require(old_triangle_box['all_old_local_pairs_pass'],'entire coordinate box passes every old pair')
    require(triangle_box['determinant_interval'][1]<-5,'entire three-point box has determinant upper below minus five')

    six_points=[(70,0),(49,49),(0,70),(-49,49),(-70,0),(-49,-49),(0,-70),(49,-49)]
    six_pairs=normalized_pairs(six_points)
    scalar_identity_checks+=six_pairs['pair_count']
    require(six_pairs['all_old_pairs_pass'],'all old local pairs pass six-point block example')
    require(six_pairs['minimum_distance_squared']==F(2842,4851)>F(4,7),
            'strong old global separation also passes')
    require(six_pairs['maximum_radius_squared']==F(4900,4851)<F(7,4),
            'strong old global radius also passes')
    rows,columns=(0,4,6),(3,5,7)
    collective=cross_block_certificate(six_points,rows,columns)
    require(not collective['passes'],'new six-point collective block fails')
    witness=collective['negative_quadratic_witness']
    require(witness is not None and witness['quadratic_value']<F(-19,10),
            'explicit rational negative collective quadratic direction')
    sboxes=boxes(six_points)
    old_six_box=all_pair_box_certificate(sboxes)
    six_box=cross_block_box_quadratic(sboxes,rows,columns,witness['vector'])
    require(old_six_box['all_old_local_pairs_pass'] and old_six_box['all_old_global_geometry_pass'],
            'entire six-point box passes old pairs and stronger global geometry')
    require(old_six_box['minimum_distance_lower']>F(29,50) and old_six_box['maximum_radius_upper']<F(51,50),
            'simple rational whole-box global geometric margins')
    require(six_box['quadratic_interval'][1]<F(-3,2),'entire six-point box quadratic upper below minus three halves')
    # Check interval containment independently at all coordinate endpoints
    # along one-coordinate edges; whole-box validity follows interval algebra.
    interval_sample_checks=0
    for base,certificate,is_triangle in ((triangle_points,triangle_box,True),(six_points,six_box,False)):
        for i in range(len(base)):
            for axis in (0,1):
                for sign in (-1,1):
                    altered=[list(p) for p in base]
                    altered[i][axis]+=sign*F(1,100)
                    if is_triangle:
                        value=triangle_certificate(altered,3,2,4)['schur_determinant']
                        bound=certificate['determinant_interval']
                    else:
                        matrix=cross_block_matrix(geometry(altered),rows,columns)
                        value=quadratic(matrix,witness['vector'])
                        bound=certificate['quadratic_interval']
                    require(bound[0]<=value<=bound[1],'independent endpoint substitution is inside whole-box interval')
                    interval_sample_checks+=1
    invalid=[lambda:geometry([(0,0),(0,0)]),lambda:geometry([(0,0),(1,0)],0),
             lambda:triangle_certificate(triangle_points,0,0,1),
             lambda:cross_block_certificate(six_points,(0,),(0,)),
             lambda:psd_certificate([[ONE,ZERO],[ZERO,ONE]],max_minors=1),
             lambda:triangle_box_certificate([((0,0),(0,0)),((0,0),(0,0))],0,1,1)]
    for call in invalid:
        try:call()
        except ValueError:pass
        else:raise ArithmeticError('invalid block input was accepted')
    return {'passed':True,'arithmetic':'integers/Fraction; no numerical eigensolver',
            'independent_determinant_checks':determinant_checks,
            'stable_original_operator_examples':len(stable_examples),
            'stable_three_point_checks':triangle_pass_checks,'stable_cross_block_checks':block_pass_checks,
            'old_pair_vs_schur_diagonal_identities':scalar_identity_checks,
            'independent_interval_endpoint_checks':interval_sample_checks,'rejected_invalid_cases':len(invalid),
            'three_point_example':{'points':triangle_points,'old_pairs':pair_info,'new_test':triangle,
                'box_coordinate_half_width':'1/100','old_whole_box':old_triangle_box,'new_whole_box':triangle_box,
                'scope':'All old local two-point tests pass; the stronger old global separation need not pass.'},
            'six_point_example':{'points':six_points,'old_pairs_and_geometry':six_pairs,'new_test':concise(collective),
                'box_coordinate_half_width':'1/100','old_whole_box':old_six_box,'new_whole_box':six_box,
                'scope':'All old local pair tests AND stronger old global distance/radius tests pass over this entire box. No comparison with all other existing global conditions is claimed.'},
            'scope':'Universal block necessity follows from the original saturated operator proof. Finite checks validate arithmetic and the stated nonempty thick classes; passing a block is not an optimality certificate.'}


if __name__=='__main__':
    report=run()
    ROOT.joinpath('block_stability_verification.json').write_text(json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('three_point_example','six_point_example')},ensure_ascii=False,indent=2))
