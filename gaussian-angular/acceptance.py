"""Deterministic exact reports, active checks under python -O."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def write_report(name,report):
    require(report.get('passed') is True,'passed report')
    (HERE/name).write_text(json.dumps(encode(report),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS '+name)


def nine_points():
    return [(F(0),F(0)),(F(7,9),F(0)),(F(-7,9),F(0)),(F(0),F(7,9)),(F(0),F(-7,9))]+[(F(s*4,9),F(r*4,9)) for s in (-1,1) for r in (-1,1)]
