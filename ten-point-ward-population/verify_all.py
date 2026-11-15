"""Check exact radial Gram bounds for the population of the inner band."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UPPER = ROOT/'ten-point-ward-upper'
RADII = ROOT/'ten-point-ward-radii'
BANDS = [(F(0), F(3,10)), (F(41,60), F(71,50)), (F(149,100), F(1539,1000))]
BUDGET = F(27001,1000000)
LOWERS = {0: F(14,10), 1: F(3,100), 4: F(43,1000), 5: F(7,10)}


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bernstein_positive(coefficients, interval, transform):
    """Finite rational midpoint subdivision, every positive hull rechecked."""
    pending = [(interval,0)]
    visited = 0
    while pending:
        (lo,hi),depth = pending.pop()
        visited += 1
        M,denominator = transform(len(coefficients)-1,lo,hi)
        need(denominator>0,'positive Bernstein conversion denominator')
        if all(sum(a*b for a,b in zip(row,coefficients))>0 for row in M):
            continue
        need(depth<20,'one dimensional minorant closes at finite depth')
        m = (lo+hi)/2
        pending.extend([((m,hi),depth+1),((lo,m),depth+1)])
    return visited


def verify():
    manifest = read(HERE/'source_manifest.json')
    for section,base in (('files',HERE),('dependencies',ROOT)):
        for name,value in manifest[section].items():
            need(sha(base/name)==value,'fixed '+section+' '+name)
    spec = importlib.util.spec_from_file_location('n10_population_upper',UPPER/'verify_all.py')
    upper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upper)
    data = read(UPPER/'kernel.json')
    gr = read(UPPER/'verification_report.json')
    rr = read(RADII/'verification_report.json')
    need(gr['passed'] is True and gr['unresolved']==0
         and gr['kernel_sha256']==sha(UPPER/'kernel.json') and gr['tree_sha256']==sha(UPPER/'cover.txt')
         and F(gr['c'])==upper.C==F(3,10000) and F(gr['radius_strict_upper'])==BANDS[-1][1],
         'accepted complete global pair bound and strict radius cap')
    need(rr['passed'] is True and rr['normalization_mass']==9 and rr['applies_to_every_global_maximizer'] is True
         and F(rr['global_pair_slack'])==upper.C,'accepted radius exclusions at all global maximizers')
    by_interval = {tuple(map(F,row['excluded_closed_interval'])): row for row in rr['covers']}
    for left,right in zip(BANDS,BANDS[1:]):
        need(by_interval[(left[1],right[0])]['unresolved']==0,'matching complete radial exclusion')
    need(2*(upper.hx(data['B'])+45*upper.C-upper.log.log_lower(upper.witness()))<BUDGET,
         'strict exact upper budget on every global maximizer radial Gram energy')
    need(6*BANDS[0][1]+4*BANDS[-1][1]<9,'mass excludes inner counts six through ten')

    Q = upper.gram.decode(data['Q_blocks'][0])
    upper.gram.positive_ldl(Q)
    features = upper.old.features(data,0)
    results = []
    entries = read(HERE/'dual.json')
    need([row['inner_population'] for row in entries]==sorted(LOWERS),'all four fixed inner-count alternatives')
    for entry in entries:
        m = entry['inner_population']
        a,b_in,b_out = [upper.hx(entry[k]) for k in ('a','b_in','b_out')]
        s0 = list(map(upper.hx,entry['s0']))
        need(len(s0)==len(Q),'dimension of Gram comparison vector')
        vector = [2*sum(q*s for q,s in zip(row,s0)) for row in Q]
        coeff = [sum(q*f for q,f in zip(vector,row)) for row in features]
        coeff[1] -= a
        checks = []
        for i,interval in enumerate(BANDS):
            if m==0 and i==0:
                continue
            poly = coeff.copy()
            poly[0] -= b_in if i==0 else b_out
            checks.append(bernstein_positive(poly,interval,upper.old.transform))
        lower = 9*a+m*b_in+(10-m)*b_out-upper.base.bilinear(s0,Q,s0)
        need(lower>LOWERS[m]>BUDGET,'strict Gram lower bound excludes the selected count')
        results.append({'inner_count_excluded':m,'strict_Gram_energy_lower':str(LOWERS[m]),
                        'interval_Bernstein_node_counts':checks})
    return {'passed':True,'normalization_mass':9,'applies_to_every_global_maximizer':True,
            'allowed_closed_band_supersets':[list(map(str,row)) for row in BANDS],
            'strict_Gram_energy_upper':str(BUDGET),'minorant_exclusions':results,
            'mass_excludes_inner_population':[6,7,8,9,10],'possible_inner_population':[2,3],
            'sharp_global_optimality_proved':False}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--record-report',action='store_true')
    args = p.parse_args()
    out = verify()
    text = json.dumps(out,indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_text(text,encoding='utf-8',newline='\n')
    else:
        need(out==read(HERE/'verification_report.json'),'saved exact population report')
    print(text,end='')


if __name__=='__main__':
    main()
