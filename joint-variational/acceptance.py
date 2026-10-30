"""Deterministic exact reports; checks remain active under -O."""
from fractions import Fraction as F
from pathlib import Path
import json


def require(ok,label):
    if not ok:raise ArithmeticError(label)


def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    return value


def write_report(name,result):
    require(result.get('passed') is True,'passed report')
    (Path(__file__).resolve().parent/name).write_text(json.dumps(encode(result),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS '+name)
