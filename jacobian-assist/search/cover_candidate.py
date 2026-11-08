"""Floating continuous cover discovery for a new candidate; NOT acceptance."""
from pathlib import Path
import argparse,importlib.util,json,math,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];SCRIPT=Path(__file__).resolve().parent;HERE=ROOT/'work/jacobian-search'
spec=importlib.util.spec_from_file_location('bern',SCRIPT/'bernstein_discovery.py')
bern=importlib.util.module_from_spec(spec);spec.loader.exec_module(bern)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',type=Path,required=True)
    ap.add_argument('--c',type=float,required=True);ap.add_argument('--seconds',type=float,default=35)
    ap.add_argument('--maxnodes',type=int,default=350000);a=ap.parse_args()
    HERE.mkdir(parents=True,exist_ok=True);data=json.loads(a.candidate.read_bytes());start=time.monotonic();U=121/72
    C=bern.coefficients(data);box=((0.,U),(0.,U),(-U,U));B=bern.bernstein(C,box)
    stack=[(box,B,0)];actions=[];visited=accepted=pruned=0;maxdepth=0;worstleaf=1e30;unresolved=[]
    while stack and visited<a.maxnodes and time.monotonic()-start<a.seconds:
        box,B,depth=stack.pop();visited+=1;maxdepth=max(maxdepth,depth)
        (ul,uh),(vl,vh),(tl,th)=box;amin=tl if tl>0 else -th if th<0 else 0
        if ul>vh+1e-12 or amin*amin>uh*vh+1e-12 or uh+vh-2*tl<1/3-1e-12:
            pruned+=1;actions.append('P');continue
        dlo=max(1/3,ul+vl-2*th);dhi=uh+vh-2*tl;point=(dlo+dhi)/2
        R=B-bern.dcoeff(box,B.shape)/point;lower=float(R.min())+1-math.log(point)
        if lower>=-a.c:
            accepted+=1;actions.append('A');worstleaf=min(worstleaf,lower);continue
        if depth>=50:
            unresolved.append({'box':box,'lower':lower,'depth':depth});continue
        axis=int(np.argmax([float(np.ptp(R,axis=j).max()) for j in range(3)]));mid=sum(box[axis])/2
        actions.append(str(axis));L,RB=bern.split(B,axis);bl=list(box);br=list(box)
        bl[axis]=(box[axis][0],mid);br[axis]=(mid,box[axis][1])
        stack.append((tuple(br),RB,depth+1));stack.append((tuple(bl),L,depth+1))
    tag=a.candidate.stem+'-c'+str(a.c);tree=HERE/(tag+'.txt')
    tree.write_text(''.join(actions)+'\n',encoding='utf-8',newline='\n')
    out={'scope':'floating discovery; no accepted global bound','candidate':str(a.candidate),'c':a.c,
         'float_complete':not stack and not unresolved,'visited':visited,'accepted':accepted,'pruned':pruned,
         'pending':len(stack),'unresolved_count':len(unresolved),'unresolved_first':unresolved[:6],
         'maximum_depth':maxdepth,'worst_accepted_lower':worstleaf,'seconds':time.monotonic()-start,
         'tree':str(tree)}
    (HERE/(tag+'.json')).write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
