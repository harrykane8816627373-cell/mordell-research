"""Exact reports and checks that remain active under -O."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value


def write_report(name,report):
    require(report.get('passed') is True,'passed report')
    (HERE/name).write_text(json.dumps(encode(report),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS '+name,flush=True)
