"""Replay exact stationary row exclusions for all ten-point global maxima."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
UPPER=ROOT/'ten-point-ward-upper'
EXCLUSIONS=((F(3,10),F(41,60)),(F(71,50),F(149,100)))
FLOOR=F(1,800)


def need(test,label):
    if not test:
        raise ArithmeticError(label)


def prunable(box,U,minimum,L):
    # u is the marked point. No ordering between u and v can be assumed.
    (ul,uh),(vl,vh),(tl,th)=box
    amin=max(tl,-th,0);dhi=uh+vh-2*tl
    if amin**2>uh*vh or dhi<minimum:
        return True
    q=9-F(9,8)*(ul+vl)-tl/4
    if q<F(142,25):
        return True
    delta=5*q-28
    if delta>=2:
        return False
    b=max(F(0),81*uh+vh+18*th)*max(F(0),uh+81*vh+18*th)/4096
    return (delta+2)**2*dhi*b**8*q**28<(2-delta)**2*L*7**21


def rowpoly(C, multipliers):
    a,b,c=multipliers
    H=C.copy()
    D={(1,0,0):F(1),(0,1,0):F(1),(0,0,1):-F(2)}
    V={(0,1,0):a,(0,0,0):FLOOR-a,(1,0,0):(a+b-5*c)/9,(0,0,1):b}
    for (i,j,k),v in D.items():
        for (r,s,t),w in V.items():
            x=i+r,j+s,k+t
            H[x]=H.get(x,F(0))-v*w
    H[1,0,0]=H.get((1,0,0),F(0))-c
    H[0,0,1]=H.get((0,0,1),F(0))+c
    return H


def split(box,axis):
    a,b=box[axis];m=(a+b)/2
    return box[:axis]+((a,m),)+box[axis+1:],box[:axis]+((m,b),)+box[axis+1:]


def verify_tree(w,C,shape,box,proof,packed):
    stack=[(box,w.initial(C,shape,box),0)]
    counts={k:0 for k in 'AP012'}
    depth=0
    for node in proof:
        need(bool(stack) and node in counts,'extra node or unknown label')
        counts[node]+=1
        box,B,k=stack.pop();depth=max(depth,k)
        if node=='P':
            need(prunable(box,w.U,w.DMIN,w.L),'valid unsorted domain pruning')
        elif node=='A':
            need(packed.dominates(B,box,F(0)),'exact row log upper tangent')
        else:
            axis=int(node)
            A,Z=split(box,axis)
            BL,BR=packed.split(B,axis)
            stack.extend([(Z,BR,k+1),(A,BL,k+1)])
    need(not stack and counts['A']+counts['P']==1+sum(counts[x] for x in '012'),'whole row tree covered')
    return counts,depth


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for p,h in manifest['files'].items():
        need(hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h,'fixed row input '+p)
    for p,h in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,'fixed prerequisite '+p)
    spec=importlib.util.spec_from_file_location('n10_ward',UPPER/'verify_all.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    data=json.loads((UPPER/'kernel.json').read_bytes())
    A,forces=w.corrections(data)
    C,shape=w.polynomial(data,A,forces)
    need(bool(forces),'rational force majorant')
    packed=w.Packed(shape,True)
    w.old.selftest()
    cutoff=9*FLOOR-36*w.C
    need(w.C==F(3,10000) and cutoff==F(9,20000),'matching global slack and marked row floor')
    # The witness is a feasible configuration in exactly the same normalization.
    need(w.log.log_lower(w.witness())>w.hx(data['B'])-cutoff,'strict contradiction with witness')
    groups=[]
    for key,band in zip(('mid','high'),EXCLUSIONS):
        saved=json.loads((HERE/f'{key}.json').read_bytes())
        need(tuple(F(saved[k]) for k in ('lo','hi'))==band and F(saved['floor'])==FLOOR,'fixed interval and lower floor')
        rows=saved['rows']
        need(bool(rows),'nonempty slab sequence')
        intervals=[tuple(map(F,row['interval'])) for row in rows]
        need(intervals[0][0]==band[0] and intervals[-1][1]==band[1], 'both endpoints of forbidden interval')
        need(all(a<b for a,b in intervals) and all(b==c for (_,b),(c,_) in zip(intervals,intervals[1:])),
             'slabs partition the forbidden interval without gaps')
        totals={k:0 for k in 'AP012'}
        depth=0
        for row,interval in zip(rows,intervals):
            box=(interval,(F(0),w.U),(-w.U,w.U))
            H=rowpoly(C,tuple(map(F,row['multipliers'])))
            counts,dep=verify_tree(w,H,shape,box,row['tree'],packed)
            for k,v in counts.items():
                totals[k]+=v
            depth=max(depth,dep)
        groups.append({'excluded_closed_interval':list(map(str,band)),'slabs':len(rows),
                       'nodes':sum(totals.values()),'counts':totals,'maximum_depth':depth,'unresolved':0})
    print(json.dumps({'passed':True,'normalization_mass':9,'applies_to_every_global_maximizer':True,
                      'no_symmetry_assumed':True,'global_pair_slack':str(w.C),
                      'marked_row_lower':str(9*FLOOR),'total_deficit_lower_if_excluded':str(cutoff),
                      'covers':groups},indent=2))


if __name__=='__main__':main()
