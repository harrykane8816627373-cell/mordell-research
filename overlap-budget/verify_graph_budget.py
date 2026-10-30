"""Actual exact vertex plans and spectral/energy short-edge comparisons."""
from fractions import Fraction as F
from pathlib import Path
import json
from edge_geometry import vertex_plan
from overlap_cost import short_edge_capacity,frobenius_edge_cap

HERE=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    annulus=[(a,b) for a in range(-8,9) for b in range(-8,9) if 9<=a*a+b*b<=64]+[(0,0)]
    plans=[vertex_plan(annulus)]
    require(plans[0]['quiet_population']>=1 and plans[0]['selected_quiet_centers'],'actual conditional input exercises quiet hole branch')
    translated=[(3*x+7,3*y-9) for x,y in annulus]
    again=vertex_plan(translated)
    for key in ('short_edge_count','cover_size','quiet_population','conditional_energy_lower'):
        require(again[key]==plans[0][key],'translation and scale covariance of exact plan')
    stored=json.loads((HERE.parent/'symmetry-breaking/integer_witnesses.json').read_text(encoding='utf-8'))
    for n in (10,30,100):
        plan=vertex_plan(next(row['integer_points'] for row in stored['rows'] if row['n']==n))
        require(plan['short_edge_count']>=plan['cover_size']>0,'actual prior input uses one-endpoint edge cover')
        plans.append(plan)
    report=json.loads((HERE.parent/'symmetry-breaking/witness_verification.json').read_text(encoding='utf-8'))
    gap100=F(next(r for r in report['rows'] if r['n']==100)['global_optimizer_gap_upper'])
    caps=[]
    for n,gap,want in ((8,F(523,5000),17),(100,gap100,136)):
        cap=short_edge_capacity(n,gap)
        require(cap['maximum_short_edges']==want,'finite energy edge capacity')
        old=frobenius_edge_cap(n)
        require(old>want,'comparison with inherited local Frobenius edge count')
        cap['local_Frobenius_edge_cap']=old;caps.append(cap)
    excluded=short_edge_capacity(8,F(1,10),dirichlet_lower=F(1,5))
    require(excluded['entire_class_excluded'] and excluded['maximum_short_edges'] is None,'negative remaining budget is whole-class exclusion')
    invalid=0
    for call in (lambda:vertex_plan([(0,0)]*4),lambda:vertex_plan(annulus,max_points=10),
                 lambda:short_edge_capacity(100,-1),lambda:short_edge_capacity(True,1)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid graph input accepted')
    return {'passed':True,'conditional_vertex_plans':plans,'short_edge_capacity_examples':caps,
        'invalid_inputs_rejected':invalid,
        'scope':'Finite graph plans test geometric hypotheses only. Global edge caps use independently certified gap bounds; passing the cap is not an optimality certificate.'}


if __name__=='__main__':
    result=run()
    (HERE/'graph_budget_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='conditional_vertex_plans'},default=str,ensure_ascii=False,indent=2))
