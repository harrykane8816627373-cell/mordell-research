"""Explicit candidate-value caps on whole neighborhoods, using integer intervals."""
from pathlib import Path
from fractions import Fraction
import hashlib,importlib.util,json,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_text(encoding='utf-8'))
    for path,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha,'fixed own input '+path)
    for path,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,'fixed prior dependency '+path)
    dyadic=load(ROOT/'validated-maxima/dyadic.py','dyadic')
    old=sys.modules.get('dyadic');sys.modules['dyadic']=dyadic
    try:prior=load(ROOT/'validated-maxima/verify_all.py','prior_basin_verifier')
    finally:
        if old is None:del sys.modules['dyadic']
        else:sys.modules['dyadic']=old
    data=json.loads((ROOT/'validated-maxima/certificates.json').read_text(encoding='utf-8'))
    claims=json.loads((HERE/'basins.json').read_text(encoding='utf-8'))
    need([r['n'] for r in claims['rows']]==[10,12,19],'complete listed candidate counts')
    B=dyadic.Ball;P=dyadic.PRECISION;S=dyadic.SCALE;results=[]
    for claim in claims['rows']:
        n,k=claim['n'],claim['radius_exponent'];need(0<k<=70,'certified radius range')
        row=next(r for r in data['rows'] if r['n']==n)
        exact_root=prior.verify_row(row,data)
        need(exact_root['classification']=='strict_local_maximum_mod_rotation','prior true local maximum')
        d=2*n-1;small=1<<(P-data['radius_exponent']);radius=1<<(P-k)
        mid=[B.dyadic(q,data['center_exponent']) for q in row['center_numerators']]
        tiny=[B(q.lo-small,q.hi+small) for q in mid]
        need(tiny[0].lo>=S and 32*small<=radius,'root first point and nested radii')
        all_coords=tiny[:1]+[B.integer(0)]+tiny[1:]
        need(all((all_coords[2*i].square()+all_coords[2*i+1].square()).hi<=4*S for i in range(n)),
             'all root point moduli at most two')
        box=[B(q.lo-radius,q.hi+radius) for q in mid]
        need(box[0].lo>0,'large rotation slice')
        _,J,_=prior.equations(box,n)
        W=prior.point_matrix(row['whitener_numerators'],data['matrix_exponent'],d)
        G=dyadic.matmul(dyadic.transpose(W),dyadic.matmul([[-q for q in r] for r in J],W))
        gap=min(G[i][i].lo-sum(G[i][j].magnitude() for j in range(d) if j!=i) for i in range(d))
        need(gap>0,'complete large-box strict negative Hessian')
        results.append({'n':n,'entire_gauge_cube_radius':'2^-'+str(k),
            'whole_orbit_euclidean_radius':'2^-'+str(k+4),
            'rational_whitened_hessian_lower':str(Fraction(gap//(1<<(P-40)),1<<40)),
            'candidate_value_cap_proved_on_entire_cube':True,
            'sphere_configuration_in_orbit_ball_cannot_exceed_candidate':True,
            'outside_this_basin_covered':False})
    report={'passed':True,'arithmetic':'integer dyadic interval at 384 bits',
            'global_optimality_proved':False,'scope':'complete stated local neighborhoods; no complement coverage',
            'rows':results,'own_input_sha256':manifest['files'],'dependency_sha256':manifest['dependencies']}
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
