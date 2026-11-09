"""Exact new nine-point-subset and nonlinear three-point necessary cuts at n12.

Reuses accepted fixed n9 theorem evidence; it does not rerun that full proof.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,prod
import argparse,hashlib,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
QMIN=F(603,100)


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def accepted_nine():
    base=ROOT/'nine-point-global';raw=(base/'source_manifest.json').read_bytes()
    accepted=json.loads((base/'accepted_replay.json').read_bytes())
    need(hashlib.sha256(raw).hexdigest()==accepted['source_manifest_sha256'],'fixed accepted n9 inputs')
    need(accepted['full_nine_point_exclusion_proved'] is True,'accepted n9 theorem flag')
    for name,sha in accepted['files_sha256'].items():
        need(hashlib.sha256((base/name).read_bytes()).hexdigest()==sha,'accepted n9 evidence '+name)
    manifest=json.loads(raw)
    for row in manifest['files']:
        data=(base/row['path']).read_bytes()
        need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],
             'fixed n9 theorem input '+row['path'])
    replay=json.loads((base/'reports/replay.json').read_bytes())
    need(replay['passed'] is True and replay['proof_replay_performed'] is True
         and replay['full_nine_point_exclusion_proved'] is True and replay['unresolved_boxes']==0,
         'complete accepted n9 proof evidence')


def uniform_upper(q):
    delta=6*q-36
    need(0<=delta<2 and q<11,'nine-point Schwarz parameter range')
    return q**36/F(8**28)*(F(2,5)*(11-q))**30*((delta+2)/(2-delta))**2


def witness():
    rows=json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())['rows']
    points=next(r['integer_points'] for r in rows if r['n']==12)
    sx=sum(x for x,y in points);sy=sum(y for x,y in points)
    D=12*sum(x*x+y*y for x,y in points)-sx*sx-sy*sy
    P=prod((x-a)**2+(y-b)**2 for i,(x,y) in enumerate(points) for a,b in points[:i])
    need(D>0 and P>0,'collision-free fixed feasible witness')
    centered=[(12*x-sx,12*y-sy) for x,y in points]
    return centered,F(11,12*D),F(132**66*P,D**66)


def triple_values(points,scale=F(1)):
    sx=sum(x for x,y in points);sy=sum(y for x,y in points)
    q=11-scale*sum(x*x+y*y for x,y in points)-scale*(sx*sx+sy*sy)/9
    delta=6*q-36
    g=scale**3*prod((x-a)**2+(y-b)**2 for i,(x,y) in enumerate(points) for a,b in points[:i])
    b=scale**3*prod(((9*x+sx)**2+(9*y+sy)**2)/F(81) for x,y in points)
    return q,delta,g,b


def nonlinear_retained(q,delta,g,b,L):
    if not(0<delta<2):return True
    return (delta+2)**2*g*b**9*q**36 >= (2-delta)**2*L*8**28


def polynomial_map(poly):
    need(poly['domain']=='QQ' and poly['variables']==['R','y'],'exact derivative polynomial carrier')
    return {tuple(t['exponents']):F(int(t['coefficient']['num']),int(t['coefficient']['den']))
            for t in poly['polynomial']['terms']}


def polynomial_multiply(a,b):
    result={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():result[(i+k,j+l)]=result.get((i+k,j+l),F(0))+x*y
    return {e:q for e,q in result.items() if q}


def derivative_identity():
    # This identity proves the uniform cross-product maximizer, for all R,y.
    source={(27-j,3+j):F(comb(27,j)*4**(27-j)*(-1)**j) for j in range(28)}
    gradients=[{(i-1,j):i*q for (i,j),q in source.items() if i},
               {(i,j-1):j*q for (i,j),q in source.items() if j}]
    base={(26-j,j):F(comb(26,j)*4**(26-j)*(-1)**j) for j in range(27)}
    factored=polynomial_multiply(polynomial_multiply({(0,2):F(1)},base),
                                {(1,0):F(12),(0,1):F(-30)})
    need(gradients[1]==factored,'dy y^3(4R-y)^27 = y^2(4R-y)^26(12R-30y)')
    record=json.loads((HERE/'jacobian_gradient.json').read_bytes())
    need(record['operation_id']=='polynomial_field.scalar.gradient.compute'
         and record['jacobian_version']=='0.23.0','actual optional Jacobian record')
    need(polynomial_map(record['input']['polynomial'])==source,'Jacobian exact input matched')
    need(polynomial_map(record['output']['source_polynomial'])==source,'Jacobian exact source matched')
    need([polynomial_map(q) for q in record['output']['components']]==gradients,
         'actual Jacobian gradient independently checked coefficient by coefficient')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record-report',action='store_true')
    args=parser.parse_args();manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for n,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/n).read_bytes()).hexdigest()==sha,'fixed own input '+n)
    for n,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==sha,'fixed dependency '+n)
    accepted_nine();need(derivative_identity(),'global polynomial identity')
    centered,scale,L=witness();endpoint=uniform_upper(QMIN)
    need(endpoint<L,'strict endpoint excluded by a true feasible lower value')
    derivative_lower=F(36)/QMIN-F(30)/(11-QMIN)+12
    need(derivative_lower>0,'uniform upper strictly increasing on 6<q<=603/100')
    values=[];active=0
    for outside in combinations(range(12),3):
        q,delta,g,b=triple_values([centered[i] for i in outside],scale)
        need(q>QMIN,'old candidate nine-point subset retained')
        need(g>0 and b>0,'old candidate cross-product nonzero')
        need(nonlinear_retained(q,delta,g,b,L),'old candidate conditional nonlinear test retained')
        values.append(q);active+=0<delta<2
    need(len(values)==220,'all 220 candidate nine-point subsets checked')
    example=[(F(6,5),F(0)),(F(1,3),F(6,5)),(F(1,3),-F(6,5))]
    q,delta,g,b=triple_values(example)
    need(q>QMIN and 0<delta<2,'example passes new scalar complement-mass cut')
    need(all(x*x+y*y<F(121,72) for x,y in example),'example passes old radial cap')
    need(all((x-a)**2+(y-c)**2>F(1,3) for i,(x,y) in enumerate(example) for a,c in example[:i]),
         'example passes all its pair separation cuts')
    example_upper=(delta+2)**2*g*b**9*q**36/((2-delta)**2*8**28)
    need(example_upper<F(3,5)*L and not nonlinear_retained(q,delta,g,b,L),
         'accurate triple cut excludes beyond scalar/radius/separation cuts')
    report={'passed':True,'date':'2026-11-08','n':12,'normalization':'sum |x_i|^2=11',
            'scope':'necessary conditions at every global maximizer; no symmetry imposed',
            'all_nine_point_subsets_centered_mass_strict_lower':'603/100',
            'outside_three_point_linear_cut':'10*sum(u_i)+2*sum(t_ij)<4473/100',
            'Schwarz_parameter':'delta=6*q-36',
            'conditional_nonlinear_cut':'(delta+2)^2*g*b^9*q^36 >= (2-delta)^2*L*8^28 when 0<delta<2',
            'uniform_endpoint_upper_over_feasible_lower':str(endpoint/L),
            'uniform_log_derivative_strict_lower':str(derivative_lower),
            'old_candidate_subsets_retained':len(values),
            'old_candidate_minimum_nine_point_mass':str(min(values)),
            'old_candidate_nonlinear_branch_active_subsets':active,
            'nonredundant_example':{'outside_points':[[str(x),str(y)] for x,y in example],
                'q':str(q),'delta':str(delta),'accurate_upper_over_L':str(example_upper/L),
                'accurate_upper_over_L_strict_upper':'3/5',
                'passes_scalar_radial_and_its_pair_cuts':True,
                'excluded_by_nonlinear_necessary_condition':True,
                'complete_twelve_point_global_configuration_claimed':False},
            'Jacobian_exact_gradient_record_independently_checked':True,
            'prior_nine_point_full_replay_reused_not_rerun':True,
            'new_global_sharp_upper_proved':False,'new_pair_only_domain_strengthening_claimed':False,
            'input_sha256':manifest['files'],'dependency_sha256':manifest['dependencies']}
    if args.record_report:
        (HERE/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
