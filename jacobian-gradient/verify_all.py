"""Standard-library replay of the exact composed Jacobian gradient records."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
AXES=['u','v','t']


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def rational(raw):
    need(type(raw) is dict and set(raw)=={'num','den'},'canonical rational fields')
    num,den=raw['num'],raw['den']
    need(type(num) is str and type(den) is str,'canonical decimal scalar strings')
    q=F(int(num),int(den))
    need(str(q.numerator)==num and str(q.denominator)==den,'reduced scalar with positive denominator')
    return q


def terms(raw):
    need(type(raw) is list,'canonical term list')
    out={};order=[]
    for term in raw:
        need(type(term) is dict and set(term)=={'exponents','coefficient'},'canonical term fields')
        exps=term['exponents']
        need(type(exps) is list and len(exps)==3
             and all(type(e) is int and 0<=e<=16384 for e in exps),'three bounded polynomial exponents')
        exps=tuple(exps);coefficient=rational(term['coefficient'])
        need(coefficient!=0 and exps not in out,'nonzero distinct monomials')
        out[exps]=coefficient;order.append(exps)
    need(order==sorted(order,reverse=True),'descending lexicographic support')
    return out


def polynomial(raw):
    need(type(raw) is dict and set(raw)=={'domain','variables','polynomial'},'canonical polynomial fields')
    need(raw['domain']=='QQ' and raw['variables']==AXES,'complete ordered rational coefficient ring')
    need(type(raw['polynomial']) is dict and set(raw['polynomial'])=={'terms'},'sparse source fields')
    return terms(raw['polynomial']['terms'])


def boundary_replay():
    record=json.loads((HERE/'boundary_gradient.json').read_bytes())
    need(record['passed'] is True and record['jacobian_version']=='0.23.0'
         and record['transport']=='injected MCP','actual recorded boundary operation')
    result=record['result']
    need(result['operation_id']=='polynomial_field.scalar.gradient.compute','boundary gradient operation')
    output=result['output'];source=polynomial(record['input']['polynomial'])
    need(source=={(1,1,0):F(1),(0,0,2):F(-1)},'physical boundary polynomial uv-t^2')
    need(polynomial(output['source_polynomial'])==source,'boundary gradient retains source')
    expected=[{(0,1,0):F(1)},{(1,0,0):F(1)},{(0,0,1):F(-2)}]
    need([polynomial(p) for p in output['components']]==expected,'exact boundary normal (v,u,-2t)')
    # Independent arithmetic for both signs of the active physical boundary.
    # grad G = lambda grad(uv-t^2), lambda>=0 gives t*G_t<=0.
    for t in (F(-1),F(1)):
        u=v=F(1);lam=F(3);gu=lam*v;gv=lam*u;gt=-2*lam*t
        need(gu+t*gt/(2*u)==0 and gv+t*gt/(2*v)==0 and t*gt<0,
             'both physical boundary signs obey the stated KKT orientation')
    return {'passed':True,'physical_normal_exact':['v','u','-2*t'],
            'normal_contact_condition':'t*G_t <= 0',
            'both_sign_examples_checked':True}


def main():
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for path,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha,'fixed gradient evidence '+path)
    for path,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,'fixed prior kernel/helper '+path)
    analytic=load(HERE/'analytic.py','accepted_composed_gradient')
    helper=load(ROOT/'twelve-point-kernel-upper/verify_cover.py','accepted_gradient_helper')
    raw=(ROOT/'twelve-point-kernel-refinement/kernel.json').read_bytes();data=json.loads(raw)
    need((data['n'],data['mass'],data['feature_order'],data['harmonics_through'])==(12,11,5,11),
         'fixed previously accepted three-percent kernel')
    base,blocks=analytic.block_polynomials(data,helper)
    expected_blocks={'affine_base':base,**{str(k):p for k,p in enumerate(blocks)}}
    reconstructed={label:{} for label in expected_blocks}
    generated_grad=[{}, {}, {}]
    native=json.loads((HERE/'native_blocks.json').read_bytes())
    need(len(native)==29,'all 29 actual native gradient calls')
    maximum_source_terms=maximum_gradient_terms=maximum_input_digits=0
    for row in native:
        need(row['operation']=='polynomial_field.scalar.gradient.compute','native gradient operation')
        label=row['label'];need(label in expected_blocks,'fixed harmonic block label')
        need(type(row['integer_limb']) is int and row['integer_limb']>=0,'nonnegative coefficient limb')
        need(type(row['source_term_offset']) is int and row['source_term_offset']%64==0,'fixed source chunk')
        scale=rational(row['exact_restoring_scale']);need(scale>0,'positive exact restoring scale')
        output=row['output'];source=polynomial(output['source_polynomial'])
        need(output['source_components'] is None and len(output['components'])==3,'scalar source, full gradient')
        need(len(source)<=64,'source admission bound')
        need(all(c.denominator==1 and abs(c.numerator).bit_length()<=250 for c in source.values()),
             'bounded integer coefficient payload, no coefficient-guard bypass')
        need(max(max(e) for e in source)<=19,'fixed bounded source exponents')
        components=[polynomial(p) for p in output['components']]
        need(sum(map(len,components))<=192,'bounded total gradient result')
        maximum_source_terms=max(maximum_source_terms,len(source))
        maximum_gradient_terms=max(maximum_gradient_terms,sum(map(len,components)))
        maximum_input_digits=max(maximum_input_digits,max(len(str(abs(c.numerator))) for c in source.values()))
        for j,g in enumerate(components):
            need(g==analytic.derivative(source,j),'every native output checked by complete exact coefficient identity')
            for e,c in g.items():analytic.add(generated_grad[j],e,scale*c)
        for e,c in source.items():analytic.add(reconstructed[label],e,scale*c)
    for label,p in expected_blocks.items():
        need(analytic.clean(reconstructed[label])==p,'native scaled sources recover full harmonic block '+label)
    whole,shape=helper.polynomial(data)
    rebuilt=dict(base)
    for p in blocks:analytic.add_poly(rebuilt,p)
    need(analytic.clean(rebuilt)==whole and shape==(20,20,12),'all blocks reconstruct complete h')
    generated_grad=list(map(analytic.clean,generated_grad))
    artifact=json.loads((HERE/'gradient.json').read_bytes())
    need(artifact['variables']==AXES and artifact['jacobian_version']=='0.23.0'
         and artifact['operation']=='polynomial_field.scalar.gradient.compute','typed derivative artifact metadata')
    need(artifact['kernel_sha256']==hashlib.sha256(raw).hexdigest(),'accepted derivative source binding')
    need(terms(artifact['h'])==whole,'accepted h matches full independently rebuilt polynomial')
    need([terms(p) for p in artifact['gradient']]==generated_grad,'accepted gradient is exactly composed native output')
    for j,g in enumerate(generated_grad):
        need(g==analytic.derivative(whole,j),'full derivative coefficient identity, not sampled evaluation')
    report=json.loads((HERE/'report.json').read_bytes())
    need(report['passed'] and report['actual_Jacobian_native_gradient_calls']==29
         and report['kernel_sha256']==hashlib.sha256(raw).hexdigest(),'producer report binding')
    need(report['gradient_artifact_sha256']==hashlib.sha256((HERE/'gradient.json').read_bytes()).hexdigest()
         and report['native_results_sha256']==hashlib.sha256((HERE/'native_blocks.json').read_bytes()).hexdigest(),
         'accepted derivative and native result digests')
    need(report['complete_h_support']==len(whole) and report['complete_gradient_supports']==list(map(len,generated_grad)),
         'producer support counts')
    points=[(F(1,4),F(2,3),-F(1,5)),(F(7,6),F(3,2),-F(3,4)),
            (F(329,625),F(121,72),F(9033,10000)),tuple(map(F,data['test_pair']))]
    for point in points:
        d=point[0]+point[1]-2*point[2];need(d>0,'positive test distance')
        for j,a in enumerate((1,1,-2)):
            got=analytic.evaluate(generated_grad[j],point)-F(a)/d
            expected=analytic.evaluate(analytic.derivative(whole,j),point)-F(a)/d
            need(got==expected,'exact residual gradient at rational test point')
    numeric=json.loads((HERE/'numeric_crosscheck.json').read_bytes())
    need(numeric['passed'] is True and numeric['points_compared']==83
         and numeric['source_sha256']==manifest['dependencies']['candidate-kernel/discover_sdp.py'],
         'separate bounded numeric diagnostic source binding')
    boundary=boundary_replay()
    result={'passed':True,'arithmetic':'standard-library Fraction and integer arithmetic',
            'source_kernel_sha256':hashlib.sha256(raw).hexdigest(),'actual_native_calls_reconstructed':len(native),
            'source_harmonic_blocks_reconstructed':len(blocks),'full_h_support':len(whole),
            'full_gradient_supports':list(map(len,generated_grad)),
            'complete_source_and_derivative_coefficient_identities_checked':True,
            'rational_points_checked':len(points),
            'native_admission_bounds':{'maximum_source_terms':maximum_source_terms,
                'maximum_total_gradient_terms':maximum_gradient_terms,
                'maximum_integer_coefficient_decimal_digits':maximum_input_digits},
            'physical_boundary_gradient_and_KKT':boundary,
            'Jacobian_runtime_needed_for_this_replay':False,
            'numeric_diagnostics_rerun':False,'new_global_upper_bound_proved':False,
            'scope':'exact derivative composition and corrected boundary orientation; no new extremum or full-domain inequality'}
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
