"""Exact full replay of a refined n12 nonsharp global upper bound.

Only Python's standard library is required. The default command performs the
whole continuous cover, all rational Gram checks and exact witness comparison.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,importlib.util,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def replay(data,prior):
    """The tree is only a proposed partition: every leaf is checked anew."""
    prior.selftest()
    crosschecks=prior.polynomial_crosscheck(data)
    B,shape,box,initial=prior.initial_coefficients(data)
    need(shape==(20,20,12),'fixed complete Bernstein tensor')
    packed=prior.Packed(shape)
    raw=(HERE/'cover.txt').read_bytes()
    actions=raw.decode('ascii').strip()
    stack=[(box,B,0)]
    visited=split=accepted=pruned=maximum_depth=0
    for action in actions:
        need(bool(stack),'tree action has a remaining box')
        box,B,depth=stack.pop();visited+=1;maximum_depth=max(maximum_depth,depth)
        if action=='P':
            need(prior.prune(box),'exact conservative necessary-domain cut')
            pruned+=1
        elif action=='A':
            need(packed.dominates(B,box,F(1,20000)),'exact Bernstein/log-tangent domination on entire leaf')
            accepted+=1
        else:
            need(action in '012','valid fixed split axis')
            axis=int(action);left,right=packed.split(B,axis);mid=sum(box[axis])/2
            bl=list(box);br=list(box)
            bl[axis]=(box[axis][0],mid);br[axis]=(mid,box[axis][1])
            stack.append((tuple(br),right,depth+1))
            stack.append((tuple(bl),left,depth+1))
            split+=1
    need(not stack and visited==2*split+1 and accepted+pruned==split+1,
         'complete binary partition conservation and no pending region')
    return {'passed':True,'n':12,'mass':11,'radial_closed_upper':'121/72',
            'pair_squared_distance_closed_lower':'1/3','c':'1/20000',
            'visited':visited,'split':split,'accepted':accepted,'domain_discarded':pruned,
            'pending':0,'maximum_depth':maximum_depth,'initial_coefficients':initial,
            'independent_exact_polynomial_evaluations':crosschecks,
            'tree_sha256':hashlib.sha256(raw).hexdigest(),
            'arithmetic':'Fraction and packed unbounded Python integers; no floating-point proof bounds'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed refinement input '+name)
    for name,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed prior theorem/helper '+name)
    source=(HERE/'kernel.json').read_bytes();data=json.loads(source)
    need((data['n'],data['mass'],data['feature_order'],data['harmonics_through'])==(12,11,5,11),
         'fixed refined twelve-point kernel')
    gram=load(ROOT/'candidate-kernel/verify_psd.py','refinement_gram')
    Qs=gram.decode(data['Q_blocks']);need(len(Qs)==12,'all twelve Gram blocks')
    pivots=[gram.positive_ldl(Q) for Q in Qs]
    need(all(all(d>0 for d in row) for row in pivots),'positive exact LDL with full reconstruction')
    old_cover=load(ROOT/'twelve-point-kernel-upper/verify_cover.py','refinement_old_cover')
    cover=replay(data,old_cover)
    old_bound=load(ROOT/'twelve-point-kernel-upper/verify_all.py','refinement_old_bound')
    V11,V12=old_bound.witness()
    actual_lower=old_bound.actual_root_lower()
    need(actual_lower>V12,'actual certified candidate strictly exceeds feasible witness')
    B=F.from_float(float.fromhex(data['B']));c=F(1,20000);ratio=F(301,300)
    gap=old_bound.log_lower(ratio*V11)-B-66*c
    need(gap>F(1,40000),'strict rational logarithm comparison to feasible value')
    report={'passed':True,'date':'2026-11-08','n':12,
            'arithmetic':'Python standard-library rational and unbounded integer arithmetic',
            'full_continuous_cover_replayed':True,'exact_positive_definite_Gram_blocks':len(pivots),
            'kernel_sha256':hashlib.sha256(source).hexdigest(),'cover':cover,
            'mass_eleven_logarithmic_global_upper':str(B+66*c),
            'global_maximum_over_old_feasible_value_strict_upper':'301/300',
            'old_feasible_value_over_global_maximum_strict_lower':'300/301',
            'true_candidate_over_global_maximum_strict_lower':'300/301',
            'log_comparison_margin_strict_lower':'1/40000',
            'actual_candidate_root_replayed_and_exceeds_witness':True,
            'old_mass_twelve_feasible_value':str(V12),
            'nonsharp_global_upper_proved':True,'sharp_global_optimality_proved':False,
            'all_maximizers_classified':False,
            'Jacobian_or_numerical_SDP_required_for_acceptance':False,
            'input_sha256':manifest['files'],'dependency_sha256':manifest['dependencies']}
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
