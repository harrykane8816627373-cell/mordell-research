"""Rebuild all exact and interval certificates from the fixed coefficient file.

Run without -O. No SDP optimizer, random search, or pre-existing cover is used.
Only NumPy, SymPy, and python-flint are required.
"""
from pathlib import Path
import sys,json,subprocess,hashlib,platform,time
if not __debug__:
    raise RuntimeError('Assertions must be enabled: do not use Python -O.')
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';start=time.monotonic()
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
import numpy,sympy,flint
source=DATA/'global_pair_bound_exact_contacts.json';sha=hashlib.sha256(source.read_bytes()).hexdigest()
def run(name,*args):
    print('VERIFY',name,*args,flush=True)
    subprocess.run([sys.executable,str(ROOT/name),*map(str,args)],cwd=ROOT,check=True)
run('verify_radius_bounds.py')
run('verify_exact_definition.py')
run('verify_positive_definiteness.py')
run('verify_contact_neighborhoods.py',0)
run('verify_global_cover.py',180)
statefile=DATA/'global_pair_bound_bernstein_frontier.json'
state=json.loads(statefile.read_text())
while state['frontier']:
    assert not state['global_majorant_proved']
    run('verify_global_cover.py',180,'--resume')
    state=json.loads(statefile.read_text())
assert state['global_majorant_proved'] and state['frontier']==[]
c=state['counts']
assert c['visited']==2*c['split']+1
assert c['positive']+c['local']+c['outside']==c['split']+1
for name in ['global_pair_bound_exact_definition_certificate.json','global_pair_bound_psd_certificate.json','global_pair_bound_local_certificate.json','global_pair_bound_bernstein_frontier.json']:
    assert json.loads((DATA/name).read_text())['coefficient_sha256']==sha
assert json.loads((DATA/'global_pair_bound_exact_definition_certificate.json').read_text())['kernel_definition_exact']
assert json.loads((DATA/'global_pair_bound_psd_certificate.json').read_text())['exact_auxiliary_matrices_positive_definite']
assert json.loads((DATA/'global_pair_bound_local_certificate.json').read_text())['local_positivity_proved']
assert hashlib.sha256(source.read_bytes()).hexdigest()==sha
manifest={'global_majorant_proved':True,'all_required_computational_checks_passed':True,
 'mathematical_conclusion_with_supplied_proof':'M8=8^28/7^21; equality exactly center plus regular heptagon, up to rotation and relabeling',
 'coefficient_sha256':sha,'cover_counts':c,'unresolved_boxes':0,
 'domain':{'u_and_v':['0','965267/524288'],'t_squared_at_most':'u*v','u_plus_v_minus_2t_at_least':'2/7'},
 'seconds':time.monotonic()-start,'versions':{'python':platform.python_version(),'numpy':numpy.__version__,'sympy':sympy.__version__,'python_flint':flint.__version__},
 'files_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(ROOT.glob('*.py'))}}
(DATA/'FINAL_CERTIFICATE.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('ALL REQUIRED CERTIFICATES PASSED; UNRESOLVED BOXES: 0',flush=True)
print(json.dumps(manifest,indent=2),flush=True)
