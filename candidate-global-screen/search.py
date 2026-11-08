"""Bounded deterministic multistart screening; every result is exploratory."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import time
import numpy as np
from geometry import normalized,diagnostics,force_jacobian

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
spec=importlib.util.spec_from_file_location('screen_prior_optimizer',ROOT/'symmetry-breaking/explore.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)


def log_value(points):
    z=normalized(points);n=len(z);delta=z[:,None,:]-z[None,:,:]
    d2=np.sum(delta*delta,axis=2)
    return float(np.log(d2[np.triu_indices(n,1)]).sum())


def input_data():
    old={r['n']:r for r in json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_text())['rows']}
    best={n:normalized(r['integer_points']) for n,r in old.items()}
    for r in json.loads((ROOT/'finite-search-nov5/integer_witnesses.json').read_text())['rows']:
        best[r['n']]=normalized(r['integer_points'])
    certified={}
    data=json.loads((ROOT/'validated-maxima/certificates.json').read_text());e=data['center_exponent']
    for r in data['rows']:
        v=[q/2**e for q in r['center_numerators']]
        certified[r['n']]=normalized(np.array(v[:1]+[0.0]+v[1:]).reshape(r['n'],2))
    rings={r['n']:r for r in json.loads((ROOT/'radial-angular/population_search_verification.json').read_text())['rows']}
    return old,best,certified,rings


def seeds_for(n,old,best,certified,rings,rng):
    base=best[n];index=np.arange(n)
    spiral=lambda a:np.column_stack((np.sqrt(2*(index+.5)/n)*np.cos(a*index),np.sqrt(2*(index+.5)/n)*np.sin(a*index)))
    seeds=[('anchors','old_integer',normalized(old[n]['integer_points'])),
           ('anchors','best_integer',base),
           ('anchors','certified_center' if n in certified else 'no_certificate_spiral',certified.get(n,spiral(2.399963229728653))),
           ('anchors','population_rings',prior.regular_seed(rings[n])),
           ('anchors','golden_spiral',spiral(2.399963229728653))]
    for amplitude in (.02,.08,.2,.45,.8):
        for k in range(6):seeds.append(('perturbations',f'noise_{amplitude}_{k}',base+amplitude*rng.normal(size=base.shape)))
    for k in range(30):
        r=np.sqrt(2*rng.uniform(size=n));a=rng.uniform(0,2*np.pi,n)
        seeds.append(('uniform_disk',str(k),np.column_stack((r*np.cos(a),r*np.sin(a)))))
    for k in range(20):seeds.append(('gaussian',str(k),rng.normal(size=(n,2))))
    for k in range(4):
        a=2.399963229728653*n+k*np.pi/2
        seeds.append(('continuation',f'previous_add_{k}',np.vstack((best[n-1],(1.4*np.cos(a),1.4*np.sin(a))))))
    larger=best[n+1];order=np.argsort(np.sum(larger*larger,axis=1))
    choices=[int(order[0]),int(order[-1])]+[int(q) for q in rng.choice(n+1,size=4,replace=False)]
    for k,index_to_remove in enumerate(choices):
        seeds.append(('continuation',f'next_remove_{k}_{index_to_remove}',np.delete(larger,index_to_remove,axis=0)))
    for k,a in enumerate((1.3,1.7,2.1,2.6,3.1)):
        seeds.append(('spirals',str(k),spiral(a)+.03*rng.normal(size=(n,2))))
    if len(seeds)!=100:raise ArithmeticError('exact100 seed budget')
    return seeds


def summarize(rows,pinned_log,tolerance=1e-7):
    valid=[r for r in rows if r['gradient_max']<1e-5]
    values=sorted(valid,key=lambda r:-r['normalized_log_value'])
    clusters=[]
    for r in values:
        if not clusters or abs(r['normalized_log_value']-clusters[-1]['highest_log_value'])>tolerance:
            clusters.append({'highest_log_value':r['normalized_log_value'],'lowest_log_value':r['normalized_log_value'],'count':0,'strata':{},'seed_indices':[]})
        c=clusters[-1];c['count']+=1;c['lowest_log_value']=min(c['lowest_log_value'],r['normalized_log_value'])
        c['strata'][r['stratum']]=c['strata'].get(r['stratum'],0)+1;c['seed_indices'].append(r['seed_index'])
    strata={}
    for name in sorted(set(r['stratum'] for r in rows)):
        sample=[r for r in rows if r['stratum']==name]
        strata[name]={'runs':len(sample),'gradient_below_1e_5':sum(r['gradient_max']<1e-5 for r in sample),
                      'at_pinned_value_within_1e_7':sum(abs(r['normalized_log_value']-pinned_log)<=tolerance and r['gradient_max']<1e-5 for r in sample)}
    return {'cluster_log_tolerance':tolerance,'gradient_screen_tolerance':1e-5,
            'pinned_log_value':pinned_log,'runs':len(rows),'gradient_screened_runs':len(valid),
            'best_value_repeat_count':clusters[0]['count'] if clusters else 0,
            'clusters':clusters,'stratified_frequencies':strata,
            'runnerup_log_gap':clusters[0]['highest_log_value']-clusters[1]['highest_log_value'] if len(clusters)>1 else None,
            'runnerup_gap_scope':'distinct observed gradient-screened logvalueclusters only; absence is not an exclusion theorem'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,required=True,choices=range(10,20));p.add_argument('--steps',type=int,default=500);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    old,best,certified,rings=input_data()
    # n=9 is known and supplies an unconstrained continuation seed for n=10.
    if 9 not in best:
        angles=2*np.pi*np.arange(8)/8;best[9]=np.vstack(((0,0),np.column_stack((np.sqrt(9/8)*np.cos(angles),np.sqrt(9/8)*np.sin(angles)))))
    rng=np.random.default_rng(2026110800+args.n);seeds=seeds_for(args.n,old,best,certified,rings,rng)
    pinned=certified.get(args.n,best[args.n]);pinned_log=log_value(pinned)
    runs=[];start=time.monotonic();selected=pinned;selected_log=pinned_log
    for index,(stratum,label,seed) in enumerate(seeds):
        z,info=prior.optimize(normalized(seed),args.steps);value=log_value(z)
        info.update(seed_index=index,stratum=stratum,seed_label=label,normalized_log_value=value)
        radial=z.ravel();rotation=np.column_stack((-z[:,1],z[:,0])).ravel()
        _,_,vh=np.linalg.svd(np.vstack((radial,rotation)),full_matrices=True)
        Q=vh[2:].T
        info['sphere_tangent_min_curvature_float']=float(np.linalg.eigvalsh(-Q.T@force_jacobian(z)@Q)[0])
        info['terminal_integer_points_grid_1e10']=np.rint(z*10**10).astype(np.int64).tolist()
        runs.append(info)
        if value>selected_log:selected,selected_log=z,value
    files=['symmetry-breaking/explore.py','symmetry-breaking/integer_witnesses.json','finite-search-nov5/integer_witnesses.json','validated-maxima/certificates.json','radial-angular/population_search_verification.json']
    result={'date':'2026-11-08','n':args.n,'numpy_version':np.__version__,'random_seed':2026110800+args.n,'max_steps_per_run':args.steps,
            'scope':'finite floating exploration; neither counts nor Hessian eigenvalues prove globaloptimality',
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
            'summary':summarize(runs,pinned_log),'runs':runs,'elapsed_seconds':round(time.monotonic()-start,3),
            'best_integer_points':np.rint(selected*10**12).astype(np.int64).tolist(),'best_log_value':selected_log,
            'best_log_gain_over_pinned':selected_log-pinned_log,'best_geometry':diagnostics(selected)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(args.n,'bestrepeat',result['summary']['best_value_repeat_count'],'clusters',len(result['summary']['clusters']),
          'runnerupgap',result['summary']['runnerup_log_gap'],'gainoverpinned',result['best_log_gain_over_pinned'],
          'seconds',result['elapsed_seconds'],flush=True)


if __name__=='__main__':main()
