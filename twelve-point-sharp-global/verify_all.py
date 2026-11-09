"""Exact zero-slack replay of the complete twelve-point sharp certificate."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import lcm, prod
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
spec=importlib.util.spec_from_file_location('sharp_integer_cover',ROOT/'twelve-point-kernel-upper/verify_cover.py')
cover=importlib.util.module_from_spec(spec)
spec.loader.exec_module(cover)
cover.S=1<<62
cover.BIAS=1<<90
cover.W=96
SHAPE=(20,20,12)
PACK=cover.Packed(SHAPE)
D=F(1,3)
INNER=(F(3,25),F(6,25))
OUTER=(F(47,50),F(71,50))
ROOTS=[(u,v,t) for u,v,r in ((INNER,INNER,F(6,25)),(INNER,OUTER,F(3,5)),(OUTER,OUTER,F(71,50)))
       for t in ((-r,F(0)),(F(0),r))]
LOCAL_ACTIONS=b'abcdefghijklmn'


def need(condition,label):
    if not condition:
        raise ArithmeticError(label)


def initial(coefficients,box):
    denominator=lcm(*(a.denominator for a in coefficients.values()))
    size=prod(SHAPE)
    tensor=[0]*size
    for (i,j,k),a in coefficients.items():
        tensor[(i*20+j)*12+k]=a.numerator*(denominator//a.denominator)
    for axis,dimension in enumerate(SHAPE):
        matrix,divisor=cover.transform(dimension-1,*box[axis])
        stride=prod(SHAPE[axis+1:])
        group=dimension*stride
        transformed=[0]*size
        for base in range(0,size,group):
            for offset in range(stride):
                values=tensor[base+offset:base+offset+group:stride]
                for i,row in enumerate(matrix):
                    transformed[base+i*stride+offset]=sum(x*y for x,y in zip(row,values))
        tensor=transformed
        denominator*=divisor
    packed=0
    for value in reversed(tensor):
        lane=cover.S*value//denominator+cover.BIAS
        need(0<=lane<1<<(cover.W-2),'initial coefficient lane stays in range')
        packed=(packed<<cover.W)|lane
    return packed


LOG_GRID=1<<96


def ceil_div(p,q):
    return -(-p//q)


def log_series_lanes(y):
    n,d=y.numerator,y.denominator
    pl=LOG_GRID*(n-d)//(n+d)
    pu=ceil_div(LOG_GRID*(n-d),n+d)
    zl=pl*pl//LOG_GRID
    zu=ceil_div(pu*pu,LOG_GRID)
    lo=hi=0
    for j in range(24):
        lo+=2*pl//(2*j+1)
        hi+=ceil_div(2*pu,2*j+1)
        pl=pl*zl//LOG_GRID
        pu=ceil_div(pu*zu,LOG_GRID)
    return lo,hi+ceil_div(3*pu,49)


LOG2_LANES=log_series_lanes(F(2))


def log_upper_lanes(y):
    k=0
    while y<1:
        y*=2
        k-=1
    while y>=2:
        y/=2
        k+=1
    _,hi=log_series_lanes(y)
    hi+=k*LOG2_LANES[1 if k>=0 else 0]
    return ceil_div(hi*cover.S,LOG_GRID)


def dominates(control,box):
    (ul,uh),(vl,vh),(tl,th)=box
    a=(max(D,ul+vl-2*th)+uh+vh-2*tl)/2
    base=cover.BIAS+cover.ceilq(cover.S*((ul+vl-2*tl)/a-1))+log_upper_lanes(a)
    T=base*PACK.ones
    lo=hi=base
    for axis,slope in enumerate((uh-ul,vh-vl,-2*(th-tl))):
        r=cover.S*slope/(a*(SHAPE[axis]-1))
        vals=[-((-r.numerator*i)//r.denominator) for i in range(SHAPE[axis])]
        lo+=min(vals)
        hi+=max(vals)
        T+=sum(v*o for v,o in zip(vals,PACK.axisones[axis]))
    need(0<=lo<=hi<1<<(cover.W-2),'tangent lane stays in range')
    return ((control+PACK.high-T)&PACK.high)==PACK.high


def quadratic_dominates(control,box):
    """Use log(d) <= log(dhi) - s - s**2/2, s=(dhi-d)/dhi."""
    (ul,uh),(vl,vh),(tl,th)=box
    a=uh+vh-2*tl
    need(a>=D,'positive distance upper endpoint for quadratic log bound')
    axes=((uh-ul)/a,(vh-vl)/a,2*(th-tl)/a)
    linear=[]
    half_squares=[]
    pair_scales=[]
    for axis,(alpha,dim) in enumerate(zip(axes,SHAPE)):
        # The distance deficit is nonnegative on the entire box. Its u and
        # v terms decrease across Bernstein indices, while its t term grows.
        indices=list(range(dim)) if axis==2 else list(reversed(range(dim)))
        degree=dim-1
        step=alpha/degree
        linear.append([cover.S*step.numerator*i//step.denominator for i in indices])
        pair_scales.append([(1<<31)*step.numerator*i//step.denominator for i in indices])
        square=alpha*alpha/(2*degree*(degree-1))
        half_squares.append([cover.S*square.numerator*i*(i-1)//square.denominator for i in indices])
    base=cover.BIAS+log_upper_lanes(a)
    T=base*PACK.ones
    total_upper=0
    for axis in range(3):
        lanes=[x+y for x,y in zip(linear[axis],half_squares[axis])]
        T-=sum(x*mask for x,mask in zip(lanes,PACK.axisones[axis]))
        total_upper+=max(lanes)
    for i,j in ((0,1),(0,2),(1,2)):
        largest=max(pair_scales[i])*max(pair_scales[j])
        need(largest<1<<cover.W,'packed multiplication cannot carry')
        packed=sum(value*mask for value,mask in zip(pair_scales[i],PACK.axisones[i]))
        for scale,mask in zip(pair_scales[j],PACK.full[j]):
            T-=(packed*scale)&mask
        total_upper+=largest
    need(0<=base-total_upper<=base<1<<(cover.W-2),'quadratic threshold lane bounds')
    return ((control+PACK.high-T)&PACK.high)==PACK.high



def split(box,control,axis):
    left,right=PACK.split(control,axis)
    lo,hi=box[axis]
    middle=(lo+hi)/2
    lbox,rbox=list(box),list(box)
    lbox[axis],rbox[axis]=(lo,middle),(middle,hi)
    return (tuple(lbox),left),(tuple(rbox),right)


def replay(job):
    actions,box,control,depth,local_boxes=job
    stack=[(box,control,depth)]
    counts=dict(accepted=0,quadratic=0,pruned=0,local=0,splits=0,maximum_depth=depth)
    local_counts=[0]*len(local_boxes)
    for action in actions:
        need(bool(stack),'action has a remaining box')
        box,control,depth=stack.pop()
        counts['maximum_depth']=max(counts['maximum_depth'],depth)
        if action==ord('P'):
            need(cover.prune(box),'valid Gram, separation, or symmetry pruning')
            counts['pruned']+=1
        elif action==ord('A'):
            need(dominates(control,box),'zero-slack Bernstein/log domination')
            counts['accepted']+=1
        elif action==ord('Q'):
            need(quadratic_dominates(control,box),'zero-slack quadratic log domination')
            counts['accepted']+=1
            counts['quadratic']+=1
        elif action in LOCAL_ACTIONS:
            index=action-ord('a')
            need(all(a<=c<=d<=b for (a,b),(c,d) in zip(local_boxes[index],box)),
                 'entire leaf inside its certified contact box')
            counts['local']+=1
            local_counts[index]+=1
        else:
            need(action in b'012','valid split action')
            (lbox,left),(rbox,right)=split(box,control,action-ord('0'))
            stack.extend(((rbox,right,depth+1),(lbox,left,depth+1)))
            counts['splits']+=1
    need(not stack and len(actions)==2*counts['splits']+1,'complete sharp subtree')
    need(counts['accepted']+counts['pruned']+counts['local']==counts['splits']+1,'sharp leaf conservation')
    return counts,local_counts


def partition(actions,box,control,local_boxes):
    ends=[0]*len(actions)
    for index in reversed(range(len(actions))):
        action=actions[index]
        if action in b'APQ'+LOCAL_ACTIONS:
            ends[index]=index+1
        else:
            need(action in b'012' and index+1<len(actions),'valid binary split')
            right=ends[index+1]
            need(index+1<right<len(actions),'valid pair of children')
            ends[index]=ends[right]
    need(ends[0]==len(actions),'one complete sharp prefix tree')
    todo=[(0,box,control,0)]
    jobs=[]
    internal=0
    while todo:
        index,box,control,depth=todo.pop()
        end=ends[index]
        if end-index<=4000:
            jobs.append((actions[index:end],box,control,depth,local_boxes))
            continue
        (lbox,left),(rbox,right)=split(box,control,actions[index]-ord('0'))
        todo.extend(((ends[index+1],rbox,right,depth+1),(index+1,lbox,left,depth+1)))
        internal+=1
    need(sum(len(job[0]) for job in jobs)+internal==len(actions),'every sharp node assigned once')
    return internal,jobs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs',type=int,default=min(12,os.cpu_count() or 1))
    parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args()
    need(1<=args.jobs<=32,'worker count between one and thirty-two')
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest.items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed input '+name)
    radial=json.loads((ROOT/'twelve-point-ward-radii/verification_report.json').read_bytes())
    need(radial['passed'] and radial['mass']==11 and radial['row_cut_cover']['pending']==0
         and radial['squared_radius_open_bands']==[[str(x) for x in b] for b in (INNER,OUTER)],
         'fixed necessary radius domain of every full maximizer')
    import kernel_bridge
    coefficients,contacts,contact_report=kernel_bridge.reconstruct()
    need(all(all(0<=power<dim for power,dim in zip(p,SHAPE)) for p in coefficients),'complete tensor shape')
    local_boxes,local_report=kernel_bridge.certify_local_boxes(coefficients,contacts)
    cover.selftest()
    report=[]
    contact_leaves=[0]*len(local_boxes)
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures={}
        for index,(u,v,t) in enumerate(ROOTS):
            need(max(abs(t[0]),abs(t[1]))**2>=u[1]*v[1] and t.count(F(0))==1,
                 'signed angular roots cover the Gram domain')
            # For t <= 0, odd powers have negative sign. Choose the upper
            # coefficient endpoint for those monomials to get a lower h.
            low={p:F(a.hi if t[0]<0 and p[2]%2 else a.lo,kernel_bridge.SCALE)
                 for p,a in coefficients.items()}
            box=(u,v,t)
            control=initial(low,box)
            actions=(HERE/f'{index}.txt').read_bytes().strip()
            internal,jobs=partition(actions,box,control,local_boxes)
            report.append(dict(nodes=len(actions),accepted=0,quadratic=0,pruned=0,local=0,
                               splits=internal,maximum_depth=0))
            for job in jobs:
                futures[pool.submit(replay,job)]=index
            if args.record_report:
                print('sharp root',index,len(actions),'chunks',len(jobs),file=sys.stderr,flush=True)
        for j,future in enumerate(as_completed(futures),1):
            counts,leaves=future.result()
            row=report[futures[future]]
            for name in ('accepted','quadratic','pruned','local','splits'):
                row[name]+=counts[name]
            row['maximum_depth']=max(row['maximum_depth'],counts['maximum_depth'])
            contact_leaves=[a+b for a,b in zip(contact_leaves,leaves)]
            if args.record_report and j%100==0:
                print('sharp exact chunks',j,'/',len(futures),'seconds',round(time.monotonic()-start,2),
                      file=sys.stderr,flush=True)
    for row in report:
        need(row['nodes']==2*row['splits']+1 and row['accepted']+row['pruned']+row['local']==row['splits']+1,
             'complete sharp signed radial partition')
    result=dict(passed=True,n=12,mass=11,point_pair_slack='0',
                node_count=sum(row['nodes'] for row in report),pending=0,
                signed_root_covers=report,local_contact_leaves=contact_leaves,
                exact_contact_reconstruction=contact_report,local_boxes=local_report,
                radial_theorem_used_as_fixed_input=True,unknown_maximum_symmetry_assumed=False,
                sharp_global_value_proved=True,all_maximizers_classified=True,input_sha256=manifest)
    target=HERE/'verification_report.json'
    rendered=json.dumps(result,indent=2)+'\n'
    if args.record_report:
        target.write_bytes(rendered.encode('utf-8'))
    else:
        need(json.loads(target.read_bytes())==result,'fresh sharp replay agrees with the saved report')
    print(rendered,end='')


if __name__=='__main__':
    main()
