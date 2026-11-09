"""Exact parallel Bernstein replay for a near-sharp bound on the proven radial bands."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import os
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
LEVEL=F(1,10**6)
RATIO=F(15000,14999)
INNER=(F(3,25),F(6,25))
OUTER=(F(47,50),F(71,50))
BOXES=((INNER,INNER,(F(-6,25),F(6,25))),
       (INNER,OUTER,(F(-3,5),F(3,5))),
       (OUTER,OUTER,(F(-71,50),F(71,50))))
p=ROOT/'twelve-point-ward-radii/verify_all.py'
spec=importlib.util.spec_from_file_location('radii',p)
radii=importlib.util.module_from_spec(spec);spec.loader.exec_module(radii)
ward=radii.load('twelve-point-ward-upper/verify_all.py','ward')
cover=ward.load('twelve-point-kernel-upper/verify_cover.py','cover')
SHAPE=(20,20,12)
packed=cover.Packed(SHAPE)


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def split(box,control,axis):
    left,right=packed.split(control,axis)
    lo,hi=box[axis]
    middle=(lo+hi)/2
    lbox,rbox=list(box),list(box)
    lbox[axis],rbox[axis]=(lo,middle),(middle,hi)
    return (tuple(lbox),left),(tuple(rbox),right)


def replay(args):
    actions,box,control,depth=args
    stack=[(box,control,depth)]
    accepted=pruned=splits=0
    deepest=depth
    for index,action in enumerate(actions):
        need(bool(stack),'action has remaining box')
        box,control,depth=stack.pop()
        deepest=max(depth,deepest)
        if action==80:
            need(cover.prune(box),f'exact prune {index} {box}')
            pruned+=1
        elif action==65:
            if not packed.dominates(control,box,LEVEL):
                raise ArithmeticError(f'exact leaf fails {index} {box}')
            accepted+=1
        else:
            need(action in b'012','valid split axis')
            (lbox,left),(rbox,right)=split(box,control,action-48)
            stack.extend(((rbox,right,depth+1),(lbox,left,depth+1)))
            splits+=1
    need(not stack and len(actions)==2*splits+1 and accepted+pruned==splits+1,'complete subtree')
    return accepted,pruned,splits,deepest


def partition(actions,box,control,cutoff=4000):
    ends=[0]*len(actions)
    for index in reversed(range(len(actions))):
        action=actions[index]
        if action in b'PA':
            ends[index]=index+1
        else:
            need(action in b'012','valid split action')
            left_end=ends[index+1]
            need(index+1<left_end<len(actions),'complete two-child structure')
            ends[index]=ends[left_end]
    need(ends[0]==len(actions),'one complete prefix tree')
    stack=[(0,box,control,0)]
    jobs=[]
    internal=0
    while stack:
        index,box,control,depth=stack.pop()
        end=ends[index]
        if end-index<=cutoff:
            jobs.append((actions[index:end],box,control,depth))
            continue
        (lbox,left),(rbox,right)=split(box,control,actions[index]-48)
        stack.extend(((ends[index+1],rbox,right,depth+1),(index+1,lbox,left,depth+1)))
        internal+=1
    need(sum(len(a) for a,*_ in jobs)+internal==len(actions),'every prefix node assigned exactly once')
    return internal,jobs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',type=int,default=min(12,os.cpu_count() or 1))
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    need(1 <= args.jobs <= 32,'worker count between one and thirty-two')
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed input '+name)
    radial_report=json.loads((ROOT/'twelve-point-ward-radii/verification_report.json').read_bytes())
    need(radial_report['passed'] and radial_report['row_cut_cover']['pending']==0 and
         radial_report['mass']==11 and
         radial_report['squared_radius_open_bands']==[[str(x) for x in b] for b in (INNER,OUTER)],
         'fixed radial domain theorem at all global maximizers')
    for u,v,t in BOXES:
        need(t[0]==-t[1] and t[1]**2>=u[1]*v[1],'root angular bounds cover the Gram domain')
    data=json.loads((HERE/'kernel.json').read_bytes())
    need((data['n'],data['mass'],data['feature_order'],data['harmonics_through'])==(12,11,5,11),
         'fixed kernel parameterization')
    gram=ward.load('candidate-kernel/verify_psd.py','gram')
    Qs=gram.decode(data['Q_blocks'])
    need(len(Qs)==12 and all(all(d>0 for d in gram.positive_ldl(Q)) for Q in Qs),'all exact Gram blocks')
    A=ward.ward_matrices(data)
    direct=ward.load('candidate-kernel/verify_all.py','direct')
    restrictions=ward.verify_restrictions(A,direct)
    h,shape=ward.combined_polynomial(data,A,cover)
    need(shape==SHAPE,'fixed complete tensor shape')
    evaluations=ward.check_evaluations(data,A,h,direct)
    bound=ward.load('twelve-point-kernel-upper/verify_all.py','witness')
    V,V12=bound.witness()
    need(bound.actual_root_lower()>V12,'true candidate strictly exceeds old witness')
    log_upper=ward.hex_rational(data['B'])+66*LEVEL
    need(bound.log_lower(RATIO*V)-log_upper>F(1,10**7),'strict ratio logarithm comparison')
    cover.selftest()
    report=[]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures={}
        for i,box in enumerate(BOXES):
            actions=(HERE/f'{i}.txt').read_bytes().strip()
            control=radii.initial(h,shape,box,cover)
            internal,jobs=partition(actions,box,control)
            for job in jobs:
                futures[pool.submit(replay,job)]=i
            report.append(dict(nodes=len(actions),accepted=0,pruned=0,splits=internal,maximum_depth=0))
            if args.record_report:
                print('root',i,len(actions),'chunks',len(jobs),file=sys.stderr,flush=True)
        for j,future in enumerate(as_completed(futures),1):
            i=futures[future]
            a,p,s,d=future.result()
            row=report[i]
            row['accepted']+=a
            row['pruned']+=p
            row['splits']+=s
            row['maximum_depth']=max(row['maximum_depth'],d)
            if args.record_report and j%50==0:
                print('exact chunks',j,'/',len(futures),'seconds',round(time.monotonic()-start,2),file=sys.stderr,flush=True)
    for row in report:
        need(row['nodes']==2*row['splits']+1 and row['accepted']+row['pruned']==row['splits']+1,
             'complete exact partition on each radial type')
    result=dict(passed=True,n=12,mass=11,point_pair_slack=str(LEVEL),
                global_maximum_over_old_feasible_value_strict_upper=str(RATIO),
                true_candidate_over_global_maximum_strict_lower=str(1/RATIO),
                log_comparison_margin_strict_lower='1/10000000',
                mass_eleven_logarithmic_global_upper=str(log_upper),
                root_covers=report,nodes=sum(r['nodes'] for r in report),pending=0,
                exact_positive_definite_Gram_blocks=len(Qs),
                exact_zero_bilinear_restrictions=restrictions,
                independent_exact_polynomial_evaluations=evaluations,
                radial_domain_used_as_fixed_theorem=True,
                unknown_maximum_symmetry_assumed=False,sharp_global_value_proved=False,
                input_sha256=manifest)
    target=HERE/'verification_report.json'
    rendered=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        target.write_bytes(rendered.encode('utf-8'))
    else:
        need(json.loads(target.read_bytes())==result,'fresh replay agrees with the saved report')
    print(rendered,end='')


if __name__=='__main__':
    main()
