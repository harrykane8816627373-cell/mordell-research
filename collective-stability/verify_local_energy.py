"""Fraction-only replay for localized annular angular energy certificates."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

from local_energy import (hat_norm,center_weight,upper_support_le,
    lower_support_ge,harmonic_region,validate_orthogonal_hats,
    band_local_certificate,shape_local_certificate,cmul,cpow)
from local_energy import sector_population_upper
from local_energy import orthogonal_holomorphic_bonus


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def weighted_global_cross(k,j,a,b,c):
    """Direct polar integral against global r^(k+2j) cos(k theta)."""
    a,b,c=map(F,(a,b,c))
    t=(b**k-a**k)/(c**k-b**k)
    inside=F(k,2)*(b**(k+j)-a**(k+j)+a**k*(b**j-a**j))
    outside=-F(k,2)*t*(c**k*(c**j-b**j)+c**(k+j)-b**(k+j))
    return inside+outside


def run(write_report=False):
    norm_checks=0
    for k,a,b,c in product(range(1,6),(F(0),F(1,10)),
                          (F(1,3),F(1)),(F(3,2),F(2))):
        if not a<b<c:
            continue
        A,B,C=a**k,b**k,c**k
        t=(B-A)/(C-B)
        # Independent polar integration of the two harmonic annuli.
        inside=F(k,2)*(B-A*A/B)
        outside=F(k,2)*t*t*(C*C/B-B)
        require(inside+outside==hat_norm(k,a,b,c),'annular norm by direct polar integration')
        require(B-A==t*(C-B),'continuity at the peak interface')
        require(center_weight(k,a,b,c,a)==0 and center_weight(k,a,b,c,c)==0,
                'zero values at support boundaries')
        require(center_weight(k,a,b,c,b)==1-(a/b)**k,'peak value from inner piece')
        if a==0:
            require(hat_norm(k,a,b,c)==k*b**k/(1-(b/c)**k),
                    'cap norm simplification')
        norm_checks+=1
        # A global harmonic r^k mode has cancelling cross-energy across the
        # two annuli, independent of the global radius enclosing the hat.
        cross_inner=F(k,2)*(B-A)
        cross_outer=-F(k,2)*t*(C-B)
        require(cross_inner+cross_outer==0,'same-mode global harmonic orthogonality')

    support_checks=0
    for R,r,t in product((F(0),F(1,5),F(1),F(2)),
                         (F(1,20),F(1,3)),(F(0),F(1,4),F(1),F(3))):
        require(upper_support_le(R*R,r*r,t*t)==(R+r<=t),
                'exact outer support comparison')
        require(lower_support_ge(R*R,r*r,t*t)==(max(R-r,F(0))>=t),
                'exact inner support comparison')
        support_checks+=2

    # Exact rational coordinates, deliberately using a small declared smear.
    points=[(1,0),(-1,0),(2,0),(-2,0),(0,5),(0,-5),(0,6),(0,-6)]
    shape=shape_local_certificate(points,[(2,0,F(1,3),1)],h=F(1,10))
    require(shape['normalization_factor_squared']==F(7,132),'shape Coulomb scaling')
    require(shape['angular_delta_lower']==F(1225,1089),'exact local moment on rational shape')
    # Symmetry kills the odd local moment on the same inner pairs.
    odd=shape_local_certificate(points,[(1,0,F(1,3),1)],h=F(1,10))
    require(odd['angular_delta_lower']==0,'local antipodal cancellation')
    factor=shape['normalization_factor_squared']
    shape_bands=[(factor*(F(a)**2+F(b)**2),factor*(F(a)**2+F(b)**2),1)
                 for a,b in points]
    global_p2_squared=F(38416,1089)
    bonus=orthogonal_holomorphic_bonus(8,shape_bands,[(2,0,F(1,3),1)],
                                      {2:global_p2_squared},3,h=F(1,10))
    require(bonus['holomorphic_angular_bonus']==F(19208,9801),
            'same-mode global holomorphic bonus exact value')
    require(shape['angular_delta_lower']+bonus['holomorphic_angular_bonus']
            ==F(30233,9801),'local plus same-mode global j0 add without repeated energy')
    # Contrast: global radial-weighted j=1 need not be orthogonal.
    require(weighted_global_cross(1,1,0,1,2)==-2,
            'nonzero j1 cross-energy example k1,a0,b1,c2')

    # Entire nonempty two-radius class, all at the default optimized smear.
    qin,qout=F(7,20),F(347,260)
    bands=[(qin,qin,35),(qout,qout,65)]
    hat=(2,F(0),F(27,50),F(51,50))
    example=band_local_certificate(100,bands,[hat],[121])
    require(35*qin+65*qout==99,'two-radius example total moment')
    require(example['hats'][0]['band_regions']==['inner_harmonic_piece','outside_support'],
            'filter rejects every outer-circle contribution')
    require(example['hats'][0]['norm']==F(210681,260000),'localized cap norm')
    require(example['delta_lower']==F(31460000,210681),'localized class lower bound')
    gap_upper=F(1657,12)
    require(example['delta_lower']>gap_upper,'entire class excluded from global maxima')

    # Explicit inner 35 points: 8 rational unit quadruples and one triangle.
    # Each quadruple is {u,-u,conjugate(u),-conjugate(u)}.  The triangle has
    # first and second power sums zero by its exact roots-of-unity identity.
    raw_p2=F(0)
    unit_directions=[]
    for j in range(1,9):
        t=F(j,1000)
        u=((1-t*t)/(1+t*t),2*t/(1+t*t))
        require(u[0]*u[0]+u[1]*u[1]==1,'rational unit-circle parametrization')
        p=cpow(u,2)
        raw_p2+=4*p[0]
        unit_directions.extend([u,(-u[0],-u[1]),(u[0],-u[1]),(-u[0],u[1])])
    require(len(set(unit_directions))==32,'all rational quadruple directions distinct')
    p2_inner=qin*raw_p2
    require(11<p2_inner<=F(56,5),'nonempty example satisfies local moment class')
    # Their near-real directions cannot coincide with the non-real triangle
    # directions +/-2pi/3, and none equals its remaining direction 1.
    require(all(z[1]!=0 and abs(z[0])>F(9,10) for z in unit_directions),
            'quadruples distinct from the triangle directions')

    # Compare the same mode k=2 with the prior first global mixed projection.
    S=F(7,4)
    require(upper_support_le(qout,F(1,50),S),'valid global support for old projection')
    prior_first_upper=F(56,5)**2*(1/(2*S*S)+4*(qin+F(3,100)-S)**2/S**4)
    require(prior_first_upper==F(92557824,765625)<gap_upper,
            'same-mode prior first global angular projection alone stays below budget')
    require(example['delta_lower']>prior_first_upper,
            'localized angular certificate strictly stronger in comparison example')
    population=sector_population_upper(100,gap_upper,2,0,F(27,50),F(49,50),
                                      F(7,20),F(99,100))
    require(population['maximum_population']==31,'local aligned population cap')

    # Exact orthogonality criterion: same mode can add on disjoint supports.
    validate_orthogonal_hats([(2,0,F(1,4),F(1,2)),(2,F(1,2),1,2)])
    validate_orthogonal_hats([(1,0,1,2),(2,0,1,2)])
    combined=band_local_certificate(100,
        [(F(1,10),F(1,10),50),(F(47,25),F(47,25),50)],
        [(2,0,F(1,3),F(1,2)),(2,1,F(3,2),2)],
        [F(1),F(4)],h=F(1,10))
    require(combined['angular_delta_lower']==sum(row['angular_delta_lower']
                                                for row in combined['hats']),
            'same-mode spatially disjoint angular penalties add')

    invalid_calls=[
        lambda:hat_norm(0,0,1,2), lambda:hat_norm(1,1,1,2),
        lambda:harmonic_region(F(1),F(1),F(1,4),0,F(1),F(2)),
        lambda:validate_orthogonal_hats([(2,0,1,2),(2,1,2,3)]),
        lambda:band_local_certificate(100,bands,[hat],[-1]),
        lambda:band_local_certificate(100,[(qin,qin,34),(qout,qout,65)],[hat],[121]),
        lambda:shape_local_certificate([(0,0),(0,0)],[(1,0,1,2)]),
        lambda:orthogonal_holomorphic_bonus(8,shape_bands,[(2,0,F(1,3),1)],
                                           {2:global_p2_squared},F(1,2),h=F(1,10)),
        lambda:orthogonal_holomorphic_bonus(8,shape_bands,[(2,0,3,4)],
                                           {2:global_p2_squared},3,h=F(1,10)),
    ]
    for callback in invalid_calls:
        try:
            callback()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid local certificate input accepted')

    report={
        'status':'PASS','arithmetic':'Fraction and outward rational logarithms',
        'annular_norm_identity_checks':norm_checks,
        'radical_support_checks':support_checks,
        'rejected_invalid_cases':len(invalid_calls),
        'rational_shape_angular_delta':str(shape['angular_delta_lower']),
        'rational_shape_global_j0_bonus':str(bonus['holomorphic_angular_bonus']),
        'rational_shape_local_plus_j0':str(F(30233,9801)),
        'example_n':100,'example_bands':[[str(a),str(b),m] for a,b,m in bands],
        'example_hat':[hat[0]]+[str(x) for x in hat[1:]],
        'example_local_moment_squared_lower':'121',
        'example_cap_norm':str(example['hats'][0]['norm']),
        'example_local_angular_delta_lower':str(example['delta_lower']),
        'existing_global_optimum_gap_upper':str(gap_upper),
        'explicit_inner_moment':str(p2_inner),
        'prior_same_mode_degree_one_upper':str(prior_first_upper),
        'aligned_population_example':population,
        'example_nonempty_construction':'35 inner points = 8 rational-unit quadruples plus regular triangle at q=7/20; outer regular 65-gon at q=347/260.',
        'scope':'Class exclusion and same-mode degree-one comparison. No failure assertion for all prior radial-plus-angular certificates or higher global mixed degrees; no full-domain cover.',
    }
    if write_report:
        Path(__file__).with_name('local_energy_verification.json').write_text(
            json.dumps(report,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,default=str,ensure_ascii=False,indent=2))
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    run(parser.parse_args().write_report)
