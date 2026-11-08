"""Exact ten-point necessary-domain constants; no symmetry premise."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import hashlib,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
L=F(97,80)*9**9


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def accepted(suite):
    base=ROOT/suite;a=json.loads((base/'accepted_replay.json').read_bytes())
    raw=(base/'source_manifest.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==a['source_manifest_sha256'],'accepted source manifest '+suite)
    for n,sha in a['files_sha256'].items():need(hashlib.sha256((base/n).read_bytes()).hexdigest()==sha,'accepted evidence '+n)
    manifest=json.loads(raw)
    for row in manifest['files']:
        b=(base/row['path']).read_bytes()
        need(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'fixed theorem input '+row['path'])
    report=json.loads((base/'reports/replay.json').read_bytes())
    need(report['passed'] is True and report['proof_replay_performed'] is True and report['unresolved_boxes']==0,
         'complete prior theorem replay evidence '+suite)


def single_upper(u):
    q=9-F(10,9)*u;delta=5*q-36
    need(0<=delta<2,'single deletion Schwarz range')
    return q**36/F(8**28)*(F(100,81)*u)**9*((delta+2)/(2-delta))**2


def pair_values(u,v,t):
    q=9-F(9,8)*(u+v)-t/4;delta=5*q-28;d=u+v-2*t
    b=(81*u+v+18*t)*(u+81*v+18*t)/4096
    return q,delta,d,b


def pair_feasible(u,v,t):
    q,delta,d,b=pair_values(u,v,t)
    if not(0<=u<F(15403,10000) and 0<=v<F(15403,10000) and t*t<=u*v and d>F(2,5) and q>F(142,25)):
        return False
    if 0<delta<2:return (delta+2)**2*d*b**8*q**28 >= (2-delta)**2*L*7**21
    return True


def uniform_pair_upper(q):
    r=72-8*q;delta=5*q-28
    return F(5,68)*r*(F(80,17)*r)**16/F(2**96)*q**28/F(7**21)*((delta+2)/(2-delta))**2


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for n,sha in manifest['files'].items():need(hashlib.sha256((HERE/n).read_bytes()).hexdigest()==sha,'fixed own input '+n)
    for n,sha in manifest['dependencies'].items():need(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==sha,'fixed dependency '+n)
    for suite in ('eight-point-global','nine-point-global'):accepted(suite)
    rows=json.loads((ROOT/'symmetry-breaking/integer_witnesses.json').read_bytes())['rows']
    points=next(r['integer_points'] for r in rows if r['n']==10)
    sx,sy=sum(x for x,y in points),sum(y for x,y in points)
    D=10*sum(x*x+y*y for x,y in points)-sx*sx-sy*sy
    P=prod((x-u)**2+(y-v)**2 for i,(x,y) in enumerate(points) for u,v in points[i+1:])
    V=F(90**45*P,D**45);need(V>L,'exact feasible mass-nine lower threshold')
    point=single_upper(F(15403,10000));pair=uniform_pair_upper(F(142,25))
    need(point<L and pair<L,'strict domain endpoints are excluded by the feasible lower threshold')
    need(F(9,1)/F(15403,10000)<F(100,9),'point upper function strictly decreases')
    need(F(55,6)-F(136,1)/F(664,25)>0,'pair upper function strictly increases on the excluded interval')
    # Independent exact derivative numerator for max_y y(5r-4y)^16:
    # derivative=(5r-4y)^15*(5r-68y), so its unique maximum is y=5r/68.
    for r in map(F,(1,2,7,25)):
        y=F(5,68)*r
        need(0<y<r/4 and 5*r-68*y==0,'pair cross-product analytic maximizer')
    centered=[(10*x-sx,10*y-sy) for x,y in points];scale=F(9,10*D)
    for i,(x,y) in enumerate(centered):
        for a,b in centered[i+1:]:
            need(pair_feasible(scale*(x*x+y*y),scale*(a*a+b*b),scale*(x*a+y*b)),
                 'saved candidate pair retained by all necessary cuts')
    u=v=F(3,2);t=F(-1,4);q,delta,d,b=pair_values(u,v,t)
    need(q>F(142,25) and u<F(15403,10000) and d>F(2,5) and t*t<u*v,
         'strict counterpair passes simpler necessary cuts')
    need(not pair_feasible(u,v,t),'nonlinear cut adds information beyond the new linear cuts')
    result={'passed':True,'scope':'necessary conditions at every ten-point global maximizer only',
            'symmetry_assumed':False,'normalization':'sum |x_i|^2=9',
            'feasible_lower_threshold':'(97/80)*9^9','point_squared_radius_strict_upper':'15403/10000',
            'pair_squared_distance_strict_lower':'2/5','eight_point_complement_centered_mass_strict_lower':'142/25',
            'linear_pair_cut':'9*(u+v)+2*t < 664/25',
            'point_endpoint_upper_over_threshold':str(point/L),
            'pair_endpoint_upper_over_threshold':str(pair/L),
            'saved_candidate_pairs_retained':45,'counterpair_to_simpler_domain':['3/2','3/2','-1/4'],
            'all_computations_exact':True,'global_sharp_upper_bound_proved':False}
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
