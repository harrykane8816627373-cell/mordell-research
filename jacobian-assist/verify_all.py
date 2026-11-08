"""Check fixed actual Jacobian records; optionally re-execute its real backends.

Default needs only the standard library and does not start Jacobian. The global
Mordell certificate remains the separate full integer cover acceptor.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse,hashlib,importlib.util,json,os,subprocess

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def rat(q):return F(int(q['num']),int(q['den']))
def hx(q):return F.from_float(float.fromhex(q))


def check_native(name,expected_path,prior):
    report=json.loads((HERE/name).read_bytes());data=json.loads((ROOT/expected_path).read_bytes())
    need(report['passed'] and report['jacobian_version']=='0.23.0','recorded native execution')
    need(report['source_path']==expected_path and report['source_sha256']==hashlib.sha256((ROOT/expected_path).read_bytes()).hexdigest(),
         'native record binds the exact published kernel')
    rows=report['native_Gram_inertia_crosschecks'];need(len(rows)==12,'all twelve matrix operations')
    for k,row in enumerate(rows):
        Q=[[hx(x) for x in r] for r in data['Q_blocks'][k]]
        need(row['harmonic']==k and row['operation']=='matrix.inertia.compute','actual inertia operation')
        need([[rat(x) for x in r] for r in row['input']['matrix']['entries']]==Q,'actual Gram input')
        out=row['output']
        need((out['n_positive'],out['n_negative'],out['n_zero'])==(len(Q),0,0),'recorded exact positive inertia')
        need([[rat(x) for x in r] for r in out['matrix']['entries']]==Q,'returned matrix unchanged')
        need(all(d>0 for d in prior.positive_ldl(Q)),'independent rational LDL reconstructed again')
    bern=report['native_radial_Bernstein_crosscheck'];source=bern['input']['polynomial']['polynomial']['terms']
    coefficients=[F(0)]*7
    for term in source:coefficients[term['exponents'][0]]=rat(term['coefficient'])
    U=F(121,72)
    expected=[sum((coefficients[i]*U**i*F(comb(j,i),comb(6,i)) for i in range(j+1)),F(0)) for j in range(7)]
    need([rat(x) for x in bern['output']['coefficients']]==expected,'recorded FLINT Bernstein coefficients independently reconstructed')
    jet=report['native_validated_log_second_jet']
    need(jet['output']['status']=='ENCLOSED' and jet['rational_center_gradient_and_Hessian_crosscheck'],
         'actual validated log-distance second jet record')
    need(report['new_global_Mordell_bound_proved'] is False,'native cross-check scope retained')
    return {'kernel':expected_path,'exact_positive_inertia_matrices':12,'exact_Bernstein_coefficients':7,
            'recorded_Arb_log_second_jet':'ENCLOSED'}


def check_mcp():
    report=json.loads((HERE/'mcp_report.json').read_bytes())
    need(report['passed'] and report['server']['version']=='0.23.0','actual local MCP handshake')
    need(set(report['tool_names'])=={'math.find','math.run'} and report['MCP_inertia_matches_independently_checked_native_result'],
         'actual discovery and run with mathematical output checked')
    text=next(b['text'] for b in report['operation_output']['content'] if b['type']=='text')
    output=json.loads(text)['output']
    need((output['n_positive'],output['n_negative'],output['n_zero'])==(5,0,0),'MCP inertia result')
    injected=json.loads((HERE/'injected_mcp_report.json').read_bytes())
    need(injected['passed'] and injected['source_path']=='twelve-point-kernel-refinement/kernel.json','actual injected MCP refined-kernel checks')
    out=injected['matrix_inertia']['result']['output']
    need((out['n_positive'],out['n_negative'],out['n_zero'])==(5,0,0),'injected exact inertia')
    sos=injected['actual_radial_feature_sos']['result']['output']
    need(all(sos[k] is True for k in ('is_valid','is_symmetric','reconstructs_polynomial','is_psd')),
         'injected radial SOS Gram certificate accepted')
    return {'actual_stdio_handshake_record_checked':True,'injected_MCP_inertia_and_radial_SOS_records_checked':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jacobian-python',type=Path,help='also run real native APIs and a fresh stdio MCP call')
    args=parser.parse_args()
    manifest=json.loads((HERE/'source_manifest.json').read_bytes())
    for name,sha in manifest['files'].items():
        need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==sha,'fixed optional Jacobian artifact '+name)
    for name,sha in manifest['dependencies'].items():
        need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'fixed kernel dependency '+name)
    spec=importlib.util.spec_from_file_location('jacobian_record_independent_ldl',ROOT/'candidate-kernel/verify_psd.py')
    prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
    rows=[check_native('native_report.json','candidate-kernel/data/n12_exact_interpretation.json',prior),
          check_native('refined_native_report.json','twelve-point-kernel-refinement/kernel.json',prior)]
    mcp=check_mcp();fresh=False
    if args.jacobian_python:
        env={**os.environ,'PYTHONUTF8':'1'}
        for data in ('candidate-kernel/data/n12_exact_interpretation.json','twelve-point-kernel-refinement/kernel.json'):
            result=subprocess.run([str(args.jacobian_python),'-B',str(HERE/'bridge.py'),'--data',str(ROOT/data)],
                cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',timeout=90)
            need(result.returncode==0,'fresh native backend execution: '+result.stderr[-2000:])
            need(json.loads(result.stdout)['passed'] is True,'fresh native computation accepted')
        result=subprocess.run([str(args.jacobian_python),'-B',str(HERE/'mcp_probe.py')],cwd=ROOT,env=env,
            capture_output=True,text=True,encoding='utf-8',timeout=120)
        need(result.returncode==0,'fresh MCP operation: '+result.stderr[-2000:])
        need(json.loads(result.stdout)['passed'] is True,'fresh MCP computation accepted');fresh=True
    print(json.dumps({'passed':True,'rows':rows,'MCP':mcp,'Jacobian_executed_now':fresh,
        'default_scope':'fixed actual execution evidence plus independent rational checks',
        'new_global_bound_proved_by_this_script':False},indent=2))


if __name__=='__main__':main()
