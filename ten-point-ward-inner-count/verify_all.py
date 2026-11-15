"""Exclude exactly three small radii by a sharpened Ward row and Gram term."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def need(ok,label):
    if not ok:
        raise ArithmeticError(label)


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for prefix,key in ((HERE,'files'),(ROOT,'dependencies')):
        for path,expected in manifest[key].items():
            need(hashlib.sha256((prefix/path).read_bytes()).hexdigest()==expected,'fixed input '+path)
    w=load('n10_ward',ROOT/'ten-point-ward-upper/verify_all.py')
    r=load('n10_radii',ROOT/'ten-point-ward-radii/verify_all.py')
    g=load('n10_radial_gram',HERE/'verify_radial.py')
    up=ROOT/'ten-point-ward-upper'
    data=json.loads((up/'kernel.json').read_bytes())
    Q,forces=w.corrections(data)
    C,shape=w.polynomial(data,Q,forces)
    need(bool(forces) and w.C==F(3,10000),'same accepted global kernel')
    saved=json.loads((HERE/'row.json').read_bytes())
    lo,hi,intercept,slope=g.row_spec('relaxed',w.C)
    need(tuple(F(saved[k]) for k in ('lo','hi','floor','slope'))==(lo,hi,intercept/9,slope/9),
         'row certificate matches required affine lower bound')
    rows=saved['rows']
    intervals=[tuple(map(F,row['interval'])) for row in rows]
    need(bool(intervals) and intervals[0][0]==lo and intervals[-1][1]==hi,'fixed interval endpoints')
    need(all(a<b for a,b in intervals) and all(b==c for (_,b),(c,_) in zip(intervals,intervals[1:])),
         'no gaps between slabs')
    packed=w.Packed(shape,True)
    w.old.selftest()

    def prunable(box):
        (ul,uh),(vl,vh),(tl,th)=box
        if ul>=w.U_NEW or vl>=w.U_NEW or uh+vh<=F(1007,5000):
            return True
        if any(a<=x<=y<=b for a,b in r.EXCLUSIONS for x,y in ((ul,uh),(vl,vh))):
            return True
        if uh+vh-2*tl<=F(203,500) and uh+vh+2*th<=F(1,2500):
            return True
        return r.prunable(box,w.U,w.DMIN,w.L)

    counts={k:0 for k in 'AP012'}
    depth=0
    for row,interval in zip(rows,intervals):
        # r.rowpoly subtracts d times the fixed floor and the exact dual
        # corrections. Add d times (fixed floor - desired affine floor).
        P=r.rowpoly(C,tuple(map(F,row['multipliers'])))
        for idx,c in (((1,0,0),F(1)),((0,1,0),F(1)),((0,0,1),F(-2))):
            P[idx]=P.get(idx,F(0))+c*(r.FLOOR-intercept/9)
            j=(idx[0]+1,idx[1],idx[2])
            P[j]=P.get(j,F(0))-c*slope/9
        box=(interval,(F(0),w.U_NEW),(-w.U_NEW,w.U_NEW))
        stack=[(box,w.initial(P,shape,box),0)]
        local={k:0 for k in counts}
        for node in row['tree']:
            need(bool(stack) and node in local,'valid action in row tree')
            local[node]+=1
            box,B,k=stack.pop();depth=max(depth,k)
            if node=='P':
                need(prunable(box),'exact physical, band or Schwarz prune')
            elif node=='A':
                need(packed.dominates(B,box,F(0)),'exact affine row Bernstein acceptance')
            else:
                axis=int(node);l,h=r.split(box,axis)
                BL,BH=packed.split(B,axis)
                stack.extend([(h,BH,k+1),(l,BL,k+1)])
        need(not stack and local['A']+local['P']==1+sum(local[k] for k in '012'),'complete exact row coverage')
        for k,v in local.items():
            counts[k]+=v
    gap=F(1,10**6)
    need(w.log.log_lower(w.witness())>w.hx(data['B'])-gap,'strict witness gap upper')
    radial=g.verify(data,gap,w.C,profile='relaxed')
    print(json.dumps({'passed':True,'normalization_mass':9,'all_global_maximizers':True,
                      'exact_inner_count':2,'inner_radius_strict_upper':'3/10',
                      'outer_count':8,'row_slabs':len(rows),'row_nodes':sum(counts.values()),
                      'counts':counts,'maximum_depth':depth,'unresolved':0,'radial_gram':radial},indent=2))


if __name__=='__main__':main()
