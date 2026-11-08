"""Full exact n12 non-sharp global bound, PSD, cover and witness comparison."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import hashlib,importlib.util,json,sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def log_lower(x):
    """Monotone dyadic rounding and positive atanh sums; no float logarithm."""
    need(x>=1,'range reduction above one')
    k=x.numerator.bit_length()-x.denominator.bit_length()
    while F(2)**k>x:k-=1
    while F(2)**(k+1)<=x:k+=1
    need(k>=0,'nonnegative logarithmic exponent')
    w=x/(F(2)**k);s=w*(1<<64)
    wlo=F(s.numerator//s.denominator,1<<64)
    need(1<=wlo<=w<2,'downward rational mantissa')
    def series(t):
        a=(t-1)/(t+1);need(0<=a<=F(1,3),'positive logarithm series')
        return 2*sum((a**(2*j+1)/(2*j+1) for j in range(80)),F(0))
    return k*series(F(2))+series(wlo)


def witness():
    row=next(r for r in json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())['rows'] if r['n']==12)
    points=row['integer_points'];need(len(points)==12,'twelve-point witness')
    sx=sum(x for x,y in points);sy=sum(y for x,y in points)
    D=12*sum(x*x+y*y for x,y in points)-sx*sx-sy*sy
    P=prod((x-a)**2+(y-b)**2 for i,(x,y) in enumerate(points) for a,b in points[:i])
    need(D>0 and P>0,'nondegenerate collision-free integer witness')
    return F(132**66*P,D**66),F(144**66*P,D**66)


def actual_root_lower():
    dyadic=load(ROOT/'validated-maxima/dyadic.py','upper_dyadic')
    old=sys.modules.get('dyadic');sys.modules['dyadic']=dyadic
    try:prior=load(ROOT/'validated-maxima/verify_all.py','upper_prior_root')
    finally:
        if old is None:del sys.modules['dyadic']
        else:sys.modules['dyadic']=old
    common=json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    row=next(r for r in common['rows'] if r['n']==12)
    accepted=prior.verify_row(row,common)
    need(accepted['classification']=='strict_local_maximum_mod_rotation','actual candidate root accepted again')
    return F(accepted['delta_interval'][0])


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed own input '+name)
    for name,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed theorem/kernel dependency '+name)
    data=json.loads((ROOT/'candidate-kernel/data/n12_exact_interpretation.json').read_bytes())
    gram=load(ROOT/'candidate-kernel/verify_psd.py','upper_gram')
    Qs=gram.decode(data['Q_blocks']);need(len(Qs)==12,'all twelve blocks')
    pivots=[gram.positive_ldl(Q) for Q in Qs]
    need(all(all(d>0 for d in row) for row in pivots),'exact PSD plus complete LDL reconstruction')
    cover=load(HERE/'verify_cover.py','upper_continuous').verify()
    need(cover['passed'] and cover['pending']==0,'full continuous majorant proof')
    B=F.from_float(float.fromhex(data['B']));c=F(21,10000);ratio=F(23,20)
    V11,V12=witness();root_lower=actual_root_lower()
    need(root_lower>V12,'true candidate value strictly exceeds old feasible witness')
    log_gap=log_lower(ratio*V11)-B-66*c
    need(log_gap>F(1,1000),'strict global/witness logarithm comparison')
    report={'passed':True,'date':'2026-11-08','n':12,
            'arithmetic':'Python standard-library exact rational and unbounded integer arithmetic',
            'full_continuous_cover_replayed':True,'exact_positive_definite_blocks':len(pivots),
            'cover':cover,'mass_eleven_logarithmic_global_upper':str(B+66*c),
            'global_maximum_over_old_feasible_value_strict_upper':'23/20',
            'old_feasible_value_over_global_maximum_strict_lower':'20/23',
            'true_candidate_over_global_maximum_strict_lower':'20/23',
            'log_comparison_margin_strict_lower':'1/1000',
            'old_mass_twelve_feasible_value':str(V12),
            'actual_candidate_root_replayed_and_exceeds_witness':True,
            'nonsharp_global_upper_proved':True,'sharp_global_optimality_proved':False,
            'all_maximizers_classified':False,
            'input_sha256':manifest['files'],'dependency_sha256':manifest['dependencies']}
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
