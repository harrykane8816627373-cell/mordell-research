"""Exact arithmetic replay for piecewise Ward kernels on the known ten-point bands."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import contextlib
import io
import sys
import hashlib
import importlib.util
import json
import math

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
BANDS=((F(43,1000),F(51,200)),(F(37,50),F(139,100)))
BOXES=[(i,j,s) for s in (-1,1) for i in range(len(BANDS)) for j in range(i,len(BANDS))]

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

upper=load('ten-point-ward-upper/verify_all.py','next_n10_old_upper')
collapse=load('ten-point-collapse-domain/verify_all.py','next_n10_collapse')
need=upper.need

def parsed(data):
    out=[]
    for kind,band,power,raw in data['band_parameters']:
        need(kind in ('count','mass','center','force') and band in ('all','inner')
             and isinstance(power,int) and 0<=power<=8,'known stationary multiplication direction')
        out.append((kind,band,power,upper.hx(raw)))
    need(len(out)==len(set(row[:3] for row in out)),'distinct zero-sum directions')
    return out

def root_box(i,j,sign):
    product=BANDS[i][1]*BANDS[j][1]*10000
    w=math.isqrt(product.numerator//product.denominator)
    w=F(w+(w*w<product),100)
    return BANDS[i],BANDS[j],(-w,F(0)) if sign==-1 else (F(0),w)

def weight(band,index):
    return band=='all' or band=='inner' and index==0

def root_polynomial(base,shape,params,i,j):
    result=base.copy()
    d={(1,0,0):F(1),(0,1,0):F(1),(0,0,1):F(-2)}
    for kind,band,power,coefficient in params:
        for axis,b1,b2 in ((0,i,j),(1,j,i)):
            if not weight(band,b1):continue
            if kind=='force':
                numerator={(1,0,0):F(1),(0,0,1):F(-1),(2,0,0):F(-5,9),
                           (1,1,0):F(-5,9),(1,0,1):F(10,9)}
            else:
                if kind=='count':row={(0,0,0):F(b2==0)-F(2-int(b1==0),9)}
                elif kind=='mass':row={(0,0,0):F(-1),(1,0,0):F(1,9),(0,1,0):F(1)}
                else:row={(0,0,1):F(1),(1,0,0):F(1,9)}
                numerator={}
                for a,x in row.items():
                    for b,y in d.items():
                        key=tuple(p+q for p,q in zip(a,b));numerator[key]=numerator.get(key,F(0))+x*y
            for deg in range(power+1):
                scale=coefficient*math.comb(power,deg)*(-1)**(power-deg)
                for (a,b,c),p in numerator.items():
                    key=(a+deg,b,c) if axis==0 else (b,a+deg,c)
                    need(all(k<n for k,n in zip(key,shape)),'expanded polynomial fits tensor')
                    result[key]=result.get(key,F(0))+scale*p
    return {key:value for key,value in result.items() if value}

def check_root(base,result,params,i,j):
    def at(poly,u,v,t):return sum(co*u**a*v**b*t**c for (a,b,c),co in poly.items())
    n=0
    for x,y in ((F(1,2),F(1,3)),(F(0),F(1)),(F(3,4),F(2,3))):
        u=x*BANDS[i][0]+(1-x)*BANDS[i][1]
        v=y*BANDS[j][0]+(1-y)*BANDS[j][1]
        for t in (-(u+v)/3,(u+v)/5):
            d=u+v-2*t;need(d>0,'nonzero crosscheck distance')
            change=F(0)
            for kind,band,power,q in params:
                for a,b,b1,b2 in ((u,v,i,j),(v,u,j,i)):
                    if not weight(band,b1):continue
                    g={'count':F(b2==0)-F(2-int(b1==0),9), 'mass':b-(9-a)/9,
                       'center':t+a/9,'force':(a-t)/d-5*a/9}
                    change+=q*(a-1)**power*g[kind]
            need(at(result,u,v,t)==at(base,u,v,t)+d*change,'exact independent piecewise value crosscheck');n+=1
    return n

def prune(box):
    return upper.prune(box) or collapse.pair_box_excluded(tuple(a for a,b in box),tuple(b for a,b in box))

def split(packer,B,box,depth,j):
    l,r=packer.split(B,j);m=sum(box[j])/2
    left=list(box);right=list(box);left[j]=(box[j][0],m);right[j]=(m,box[j][1])
    return (tuple(right),r,depth+1),(tuple(left),l,depth+1)

def verify(data,tree,c,c_inner):
    need((data['n'],data['mass'])==(10,9),'fixed problem normalization')
    n_q=data['harmonics_through']+1
    need(len(data['Q_blocks'])==n_q,'all harmonic PSD blocks')
    pivots=[upper.gram.positive_ldl(Q) for Q in upper.gram.decode(data['Q_blocks'])]
    A,force=upper.corrections(data)
    need(bool(force),'fixed denominator cleared by multiplying by d')
    C,shape=upper.polynomial(data,A,force);upper.crosscheck(data,A,force,C)
    params=parsed(data);packer=upper.Packed(shape,True);upper.old.selftest()
    need([tuple(row['id']) for row in tree['roots']]==BOXES,'all band pairs in both t signs')
    counts=dict(A=0,P=0,S=0);crosschecks=0;statistics=[];max_depth=0
    for (i,j,sgn),row in zip(BOXES,tree['roots']):
        c_region=c_inner if i==j==0 and c_inner is not None else c
        box=root_box(i,j,sgn);g=root_polynomial(C,shape,params,i,j)
        crosschecks+=check_root(C,g,params,i,j)
        B=upper.initial(g,shape,box);stack=[(box,B,0)];seq=[]
        for node in row['tree']:
            need(bool(stack),'tree node has a remaining box')
            box,B,k=stack.pop();max_depth=max(max_depth,k)
            if node in '012':
                stack.extend(split(packer,B,box,k,int(node)));seq.append(node);continue
            if node=='P' and prune(box) or node=='A' and packer.dominates(B,box,c_region):
                seq.append(node);continue
            raise ArithmeticError(f'invalid {node} leaf at region {i,j,sgn}, box {box}')
        need(not stack,'fully expanded root and exact coverage')
        seq=''.join(seq);row['tree']=seq
        stats={k:seq.count(k) for k in 'AP012'}
        splits=sum(stats[k] for k in '012')
        need(stats['A']+stats['P']==1+splits,'binary tree conservation')
        counts['A']+=stats['A'];counts['P']+=stats['P'];counts['S']+=splits
        statistics.append({'id':[i,j,sgn],'counts':stats,'visited':len(seq)})
        print('exact region',i,j,sgn,len(seq),'nodes',flush=True)
    B=upper.hx(data['B'])
    c_inner=c if c_inner is None else c_inner
    return {'passed':True,'scope':'exact bandwise majorant given prior two-inner theorem','c':str(c),'c_inner':str(c_inner),'B':str(B),
            'visited':sum(counts.values()),'counts':counts,'unresolved':0,'maximum_depth_before_repairs':max_depth,
            'regions':statistics,'tensor_shape':shape,'PSD_blocks':n_q,'min_gram_pivot_display':float(min(map(min,pivots))),
            'polynomial_crosschecks':crosschecks,'piecewise_zero_sum_directions':len(params),
            'logarithmic_gap_upper_display':float(B+44*c+c_inner-upper.log.log_lower(upper.witness()))}

C=F(1,10000)
C_INNER=C
RATIO=F(201,200)


def verify_dependencies():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for group,folder in (('files',HERE),('dependencies',ROOT)):
        for name,digest in manifest[group].items():
            need(hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest,'fixed proof source '+name)
    collapse_report=collapse.verify()
    need(collapse_report['passed'] is True,'exact collapse domain replay')
    import contextlib
    import io
    pair_domain=load('maximizer-domain/verify_all.py','banded_pair_domain')
    capture=io.StringIO()
    with contextlib.redirect_stdout(capture):
        pair_domain.main()
    pair_report=json.loads(capture.getvalue())
    need(pair_report['passed'] is True and pair_report['symmetry_assumed'] is False
         and F(pair_report['pair_squared_distance_strict_lower'])==upper.DMIN
         and F(pair_report['eight_point_complement_centered_mass_strict_lower'])==F(142,25),
         'exact Schwarz pair domain checked for all maximizers')
    def saved(directory):
        report=json.loads((ROOT/directory/'verification_report.json').read_bytes())
        need(report['passed'] is True,'accepted evidence '+directory)
        return report
    inner=saved('ten-point-ward-inner-count')
    tail=saved('ten-point-tail-exclusion')
    radii=saved('ten-point-counted-radii')
    need(inner['exact_inner_count']==tail['exact_inner_count']==2 and inner['outer_count']==tail['outer_count']==8,
         'accepted exact split between the two bands')
    need(inner['all_global_maximizers'] is True and tail['all_global_maximizers'] is True
         and radii['applies_to_every_global_maximizer'] is True,'no symmetry assumed in inherited domain')
    intervals={tuple(map(F,row['excluded_closed_interval'])) for row in radii['row_exclusions'] if row['unresolved']==0}
    need(tuple(map(F,radii['small_radius_exclusion']))==(F(0),F(7,1000))
         and {(F(7,1000),BANDS[0][0]),(BANDS[0][1],F(3,10)),(F(41,60),BANDS[1][0]),(BANDS[1][1],F(71,50))}<=intervals,
         'accepted radius cuts leave exactly the two required bands')
    need(tuple(map(F,tail['outer_squared_radius_open_interval']))==(F(41,60),F(71,50))
         and F(tail['inner_squared_radius_strict_upper'])==F(3,10),'accepted outer tail exclusion')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-report',action='store_true')
    parser.add_argument('--progress',action='store_true')
    args=parser.parse_args()
    verify_dependencies()
    kernel=(HERE/'kernel.json').read_bytes()
    certificate=(HERE/'cover.json').read_bytes()
    data=json.loads(kernel)
    tree=json.loads(certificate)
    previous=json.loads((ROOT/'ten-point-ward-upper/kernel.json').read_bytes())
    for key in ('radial_nodes','radial_moment_centers','transforms','Q_blocks','B'):
        need(data[key]==previous[key],'exact original PSD kernel and budget retained')
    need(len(parsed(data))==19,'all piecewise stationary corrections')
    inner_count,outer_count=2,8
    n_inner=math.comb(inner_count,2)
    n_other=inner_count*outer_count+math.comb(outer_count,2)
    need((n_inner,n_other)==(1,44),'correct edge counts by class')
    V=upper.witness()
    budget=upper.hx(data['B'])+n_inner*C_INNER+n_other*C
    log_margin=upper.log.log_lower(V*RATIO)-budget
    need(log_margin>F(1,10000),'strict global/witness logarithm comparison')
    inner_pair_budget=upper.hx(data['B'])-upper.log.log_lower(V)+n_other*C
    need(inner_pair_budget<F(4401,1000000),'strict inner-pair deficit upper')
    need(upper.actual_root_lower()>V,'independently accepted candidate exceeds integer witness')
    with contextlib.redirect_stdout(sys.stderr if args.progress else io.StringIO()):
        report=verify(data,tree,C,C_INNER)
    report.update(kernel_sha256=hashlib.sha256(kernel).hexdigest(),proof_sha256=hashlib.sha256(certificate).hexdigest(),
                  candidate_over_global_strict_lower=str(1/RATIO),global_over_witness_strict_upper=str(RATIO),
                  log_margin_strict_lower='1/10000',log_upper_mass_nine=str(budget),
                  inner_pair_deficit_strict_upper='4401/1000000',
                  exact_inner_population=inner_count,exact_outer_population=outer_count,sharp_global_optimum_proved=False)
    out=json.dumps(report,indent=2)+'\n'
    if args.record_report:
        (HERE/'verification_report.json').write_text(out,encoding='utf-8',newline='\n')
    print(out,end='')


if __name__=='__main__':
    main()
