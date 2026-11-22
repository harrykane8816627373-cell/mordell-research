"""Check the bounds and inputs of the six-real-variable half-turn chart."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def read(path):
    return json.loads(path.read_bytes())


def need(ok,message):
    if not ok:
        raise ArithmeticError(message)


def load(path):
    spec=importlib.util.spec_from_file_location('half_turn_saddle_source',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify():
    manifest=read(HERE/'source_manifest.json')
    for base,key in ((HERE,'files'),(ROOT,'dependencies')):
        for name,h in manifest[key].items():
            need(sha256((base/name).read_bytes()).hexdigest()==h, 'fixed input '+name)
    inner=read(ROOT/'ten-point-conditioned-inner/verification_report.json')
    bands=read(ROOT/'ten-point-angular-radii/verification_report.json')
    gaps=read(ROOT/'ten-point-outer-gap/verification_report.json')
    need(all(x['passed'] and x['unresolved']==0 and x['applies_to_every_global_maximizer']
             and x['normalization_mass']==9 for x in (inner,bands,gaps)), 'accepted global domain inputs')
    a0,a1=(F(x)/4 for x in inner['inner_squared_distance_open_interval'])
    r0,r1=map(F,bands['outer_squared_radius_open_interval'])
    c0,c1=map(F,gaps['outer_successive_angle_cosine_open_interval'])
    need((a0,a1)==(F(41,400),F(11,80)) and (c0,c1)==(F(59,100),F(39,50))
         and (bands['exact_inner_count'],bands['exact_outer_count'])==(2,8), 'chart ranges and populations')

    def weight(r):
        return (17-10*r)*r

    need(0<a0<a1<F(17,20) and 0<weight(a1)<F(9,4), 'inner weight below 9/4')
    need(0<r0<r1<F(17,10) and min(weight(r0),weight(r1))>5, 'outer weights above five')
    need(F(7,10)**2<F(1,2)<c1**2 and c0>F(1,2), 'shortest of four doubled angular gaps')
    need((F(24,5)+F(9,4))**2<50, 'nonzero Ward sum: strict radius lower bound')
    low,high=2*c0*c0-1,2*c1*c1-1
    need(1-max(low*low,high*high)>F(19,20)**2, 'strict lower bound for doubled-gap sine')
    need(F(25)*F(19,20)/(2*F(289,40))**2>F(1,9),
         'nondegenerate second intersection from strictly positive orientation')

    saddle_path=ROOT/'ten-point-saddle-exclusion'
    old=read(saddle_path/'verification_report.json')
    data=load(saddle_path/'verify_all.py')
    v=data.DIRECTION
    need(old['passed'] and old['unique_nonsingular_five_variable_root']
         and F(old['positive_second_variation_lower'])>0, 'accepted saddle curvature')
    need(len(v)==10 and all(tuple(-s for s in v[2*i])==v[2*i+1] for i in range(5)),
         'positive second variation preserves all five antipodal pairs')

    records=read(HERE/'search-results.json')
    counts=[row['occurrences'] for row in records['clusters']]
    need(sum(counts)==records['accepted']<=records['settings']['starts'], 'finite search record counts')
    return {'passed':True,'normalization_mass':9,'assumption':'half-turn symmetry of a global maximizer',
            'free_real_coordinates':6,'phase_charts':2,
            'initial_box':[list(map(str,x)) for x in ((a0,a1),(r0,r1),(r0,r1),(r0,r1),
                                                    (F(7,10),c1),(F(-1),F(1)))],
            'doubled_gap_sine_strict_lower':'19/20','ward_sum_modulus_strict_lower':'24/5',
            'intersection_radical_strict_lower':'1/9','known_saddle_positive_direction_preserves_half_turn':True,
            'numerical_search':{'starts':records['settings']['starts'],'accepted':records['accepted'],
                                'clusters':len(counts),
                                'clusters_within_old_global_bands':sum(x['within_global_bands'] for x in records['clusters'])}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    proof=verify()
    path=HERE/'verification_report.json'
    output=json.dumps(proof,indent=2)+'\n'
    if args.record_report:
        path.write_text(output,encoding='ascii',newline='\n')
    else:
        need(proof==read(path),'saved chart report')
    print(output,end='')


if __name__=='__main__':
    main()
