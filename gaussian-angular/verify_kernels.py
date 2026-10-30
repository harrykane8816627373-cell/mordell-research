"""Kernel intervals, actual disk averages and original review replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,json,zipfile
from gaussian_energy import (exp_negative,laplace_moment,add,sub,scale,
    standard_probe,probe_gram,measurement,gaussian_gram,quarter_probe)
from acceptance import require,write_report,nine_points,HERE


def run():
    imported=HERE/'imported-review';mf=json.loads((imported/'source_manifest.json').read_text(encoding='utf-8'))
    for row in mf['files']:
        b=(imported/row['file']).read_bytes()
        require(len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],'original review bytes')
    archive=imported/'local-smoothing-transfer-2026-10-31.zip'
    require(hashlib.sha256(archive.read_bytes()).hexdigest()=='1af797f37b90fec84f6930268e6d80a7c68da54163a020759e0a09870c23e410','original review ZIP')
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():require(z.read(name)==(imported/name).read_bytes(),'review archive member')
    spec=importlib.util.spec_from_file_location('review_original',imported/'verify_gaussian_bridge.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    accepted=json.loads((imported/'gaussian_bridge_verification.json').read_text(encoding='utf-8'))
    require(old.encode(old.run())==accepted,'original exact review replay')
    for x in (F(0),F(1,2),F(1),F(7),F(16),F(64)):
        lo,hi=exp_negative(x)
        require(lo>0 and old.log_interval(lo)[0]<=-x<=old.log_interval(hi)[1],'exponential versus independent atanh logs')
    recurrences=0
    for y in (F(0),F(1,3),F(4,9),F(15)):
        for m in (0,1,4,8):
            v=laplace_moment(m,y)
            if y==0:require(v[0]<=F(1,m+1)<=v[1],'zero-scale moment')
            else:
                lhs=scale(v,y)
                rhs=sub((F(1),F(1)) if m==0 else scale(laplace_moment(m-1,y),m),exp_negative(y))
                require(max(lhs[0],rhs[0])<=min(lhs[1],rhs[1]),'integration-by-parts recurrence')
            recurrences+=1
    probe=standard_probe();M=measurement(nine_points(),probe);H=probe_gram(probe,probe)
    oldM=list(map(F,accepted['measurement_interval']));oldH=list(map(F,accepted['Dirichlet_norm_interval']))
    require(oldM[0]<=M[0]<=M[1]<=oldM[1] and oldH[0]<=H[0]<=H[1]<=oldH[1],'sharper independent Gaussian implementation')
    require(0<H[0]<=H[1]<F(487,2000),'fixed norm upper for global distribution classes')
    for gamma in (F(1,2),F(2),F(3)):
        scaled=standard_probe(F(1,16)*gamma*gamma,F(1,2)*gamma)
        require(probe_gram(scaled,scaled)==H,'exact two-dimensional scale invariance')
    require(gaussian_gram((0,0),(0,0),F(1,8),F(1,8))==(F(1,2),F(1,2)),'Gaussian self norm')
    require(gaussian_gram((0,0),(2,0),F(1,8),F(1,8))[1]<0,'negative cross Gram entries retained')
    invalid=0
    for fn in [lambda:exp_negative(129),lambda:laplace_moment(0,17),lambda:quarter_probe((1,0),(0,2),1),lambda:quarter_probe((1,0),(0,1),1),
               lambda:gaussian_gram((0,0,0),(1,1,1),1,1),lambda:measurement([(0,0,0),(1,0,0)],probe)]:
        try:fn()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid kernel input accepted')
    return {'passed':True,'imported_review_replayed':True,'laplace_recurrences':recurrences,
            'standard_measurement_interval':M,'standard_norm_interval':H,'global_norm_upper':'487/2000',
            'invalid_inputs_rejected':invalid,'scope':'Exact Gaussian and disk-mean arithmetic. No external wave estimate is an acceptance input.'}


if __name__=='__main__':write_report('kernels_verification.json',run())
