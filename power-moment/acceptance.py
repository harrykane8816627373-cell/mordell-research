from fractions import Fraction as F
from pathlib import Path
import json
from complex_interval import C,I

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def require(ok,label):
    if not ok:raise ArithmeticError(label)

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,C):return [str(x.re),str(x.im)]
    if isinstance(x,I):return {'re':encode(x.re),'im':encode(x.im)}
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x

def write_report(name,report):
    require(report.get('passed') is True,'passed report')
    (HERE/name).write_text(json.dumps(encode(report),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS '+name,flush=True)
