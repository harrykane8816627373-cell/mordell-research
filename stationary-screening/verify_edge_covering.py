"""Exact Chernoff rate, all-n tail margin and finite whole-disk covers."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from edge_covering import boundary_certificate

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'joint-bounds'))
sys.path.insert(0,str(HERE.parent/'distribution-geometry'))
from subset_moment import log_interval,outward_decimal
from microscopic_covering import exp_interval


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def run():
    rates=0
    for j in range(1,20):
        e=F(j,20)
        lo,hi=log_interval(1-e)
        require(-e-hi>=e*e/2,'Chernoff quadratic rate enclosure')
        # Independent positive expansion has terms k>=3 left over.
        lower=sum((e**k/k for k in range(2,35)),F(0))
        require(lower>=e*e/2,'rate from positive power series')
        rates+=1
    require(log_interval(512)[1]<F(64,9),'epsilon(64)<2/3 and later decreasing')
    require(exp_interval(F(23,8),12)[0]>16,'Gaussian separation tail <1/2')
    tails=0
    for n in range(64,201):
        require(F(1,8)+F(1,64*n*n)<F(1,4),'all-n interpolation contradiction margin')
        tails+=1
    rows=[]
    for n in (64,100,512,10000,10**6,10**8,10**12,10**18):
        row=boundary_certificate(n)
        eps=row['epsilon_upper']
        require(n*eps*eps>=4*log_interval(8*n)[1],'rounded epsilon is conservative')
        require(0<row['certified_bulk_squared']<2,'positive shrinking bulk')
        require(row['whole_background_covering_strict_upper']<=row['simple_covering_upper'],'exact radial shrink improves simple bound')
        whole,simple=row['whole_background_covering_strict_upper'],row['simple_covering_upper']
        rows.append({'n':n,'whole_background_covering_upper_outward':outward_decimal((whole,whole),12)[1],
                     'simple_covering_upper_outward':outward_decimal((simple,simple),12)[1]})
    invalid=0
    for call in (lambda:boundary_certificate(63),lambda:boundary_certificate(True),lambda:boundary_certificate(100,scale=10)):
        try:call()
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid boundary cover accepted')
    return {'passed':True,'independent_rate_checks':rates,'finite_tail_checks':tails,
        'whole_background_examples':rows,'invalid_inputs_rejected':invalid,
        'theorem':'For every global optimizer and n>=64, all |u|<=sqrt(2) have min_i|z_i-u| < 6/sqrt(n-1)+sqrt(2)-sqrt(2*(1-2*sqrt(log(8*n)/n))).',
        'rate':'O(sqrt(log(n)/n)), original z scale',
        'scope':'Effective whole closed background disk, via shrinking-bulk interpolation. Does not prove boundary O(n^-1/2), optimal constants or crystallization.'}


if __name__=='__main__':
    result=run()
    (HERE/'edge_covering_verification.json').write_text(json.dumps(result,default=str,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,default=str,ensure_ascii=False,indent=2))
