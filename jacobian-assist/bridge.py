"""Actual Jacobian native cross-checks for the existing Mordell kernel data.

Optional discovery dependency only. The old standard-library proof acceptors
remain independent. Every reported backend result comes from Jacobian 0.23.0.
"""
from pathlib import Path
from fractions import Fraction as F
from importlib.metadata import version
from math import comb
import argparse,hashlib,importlib.util,json

from jacobian.math.analysis._models import IntervalExpressionNode,RationalIntervalBox
from jacobian.math.analysis.intervals import RationalBox
from jacobian.math.analysis.operations import second_jet_enclosure
from jacobian.math.matrices.analysis import compute_inertia
from jacobian.math.matrices.values import RationalMatrix
from jacobian.math.polynomials.bernstein import bernstein_coefficients
from jacobian.math.polynomials.values import RationalPolynomial

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'candidate-kernel/data/n12_exact_interpretation.json'


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def rat(q):
    q=F(q)
    # Native Python validation uses integer fields; JSON serialization projects
    # them to the transport's canonical decimal strings.
    return {'num':q.numerator,'den':q.denominator}


def hx(q):return F.from_float(float.fromhex(q))


def polynomial(coefficients):
    return RationalPolynomial.model_validate({'domain':'QQ','variables':['u'],
        'polynomial':{'terms':[{'coefficient':rat(coefficients[i]),'exponents':[i]}
         for i in range(len(coefficients)-1,-1,-1) if coefficients[i]]}})


def matrix(entries):
    return RationalMatrix.model_validate({'domain':'QQ',
        'entries':[[rat(q) for q in row] for row in entries]})


def multiply(a,b):
    result=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]+=x*y
    return result


def gram_crosscheck(data):
    spec=importlib.util.spec_from_file_location('jacobian_bridge_prior_psd',ROOT/'candidate-kernel/verify_psd.py')
    prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
    rows=[]
    for k,block in enumerate(data['Q_blocks']):
        entries=[[hx(q) for q in row] for row in block]
        pivots=prior.positive_ldl(entries)
        need(all(q>0 for q in pivots),'independent standard-library LDL positive')
        source=matrix(entries)
        result=compute_inertia(source)
        need((result.n_positive,result.n_negative,result.n_zero)==(len(entries),0,0),
             'Jacobian exact inertia agrees with full LDL')
        need(result.definiteness=='positive_definite','Jacobian positive definiteness')
        rows.append({'harmonic':k,'operation':'matrix.inertia.compute',
                     'input':{'matrix':source.model_dump(mode='json')},
                     'output':result.model_dump(mode='json'),
                     'independent_LDL_reconstruction_and_positive_pivots':True})
    need(len(rows)==12,'all twelve stored Gram blocks checked')
    return rows


def radial_bernstein_crosscheck(data):
    # This is precisely the final radial feature 10*q(u)^2 of the n12 fit.
    coefficients=[F(1)]
    for node in data['radial_nodes']:coefficients=multiply(coefficients,[-hx(node),F(1)])
    coefficients=[10*q for q in multiply(coefficients,coefficients)]
    source=polynomial(coefficients);U=F(121,72);degree=len(coefficients)-1
    box=RationalBox.model_validate({'domain':'QQ','variables':['u'],
        'intervals':[{'lower':rat(0),'upper':rat(U)}]})
    result=bernstein_coefficients(source,box,(degree,))
    expected=[sum((coefficients[i]*U**i*F(comb(k,i),comb(degree,i))
                   for i in range(k+1)),F(0)) for k in range(degree+1)]
    actual=[q.as_fraction() for q in result.coefficients]
    need(actual==expected,'Jacobian FLINT tensor coefficients agree with independent exact formula')
    return {'operation':'polynomial.bernstein.coefficients.compute',
            'feature':'actual n12 radial feature 10*prod(u-a_r)^2',
            'input':{'polynomial':source.model_dump(mode='json'),
                     'box':box.model_dump(mode='json'),'multidegree':[degree]},
            'output':result.model_dump(mode='json'),
            'independent_exact_coefficients':[str(q) for q in expected],
            'coefficient_identity_checked':True,
            'nonpositive_coefficients_do_not_disprove_positivity':True}


def logarithm_second_jet():
    # This box surrounds a low-margin region missed by the original SLSQP
    # violation search. It lies in d>0; no claim that it is globally minimizing.
    expression={'op':'log','children':[{'op':'sub','children':[
        {'op':'add','children':[{'op':'var','variable':'u'},{'op':'var','variable':'v'}]},
        {'op':'mul','children':[{'op':'const','value':rat(2)},{'op':'var','variable':'t'}]}]}]}
    limits=[(F(5263,10000),F(5265,10000)),(F(16805,10000),F(121,72)),
            (F(9032,10000),F(9034,10000))]
    box_data={'variables':['u','v','t'],
        'intervals':[{'lower':rat(a),'upper':rat(b)} for a,b in limits]}
    source=IntervalExpressionNode.model_validate(expression)
    box=RationalIntervalBox.model_validate(box_data)
    result=second_jet_enclosure(source,box,128)
    need(result.status=='ENCLOSED','Arb validates the complete log-distance second jet')
    center=[(a+b)/2 for a,b in limits];d=center[0]+center[1]-2*center[2]
    need(d>0,'positive rational center distance')
    weights={'u':F(1),'v':F(1),'t':F(-2)}
    def contains(interval,q):return interval.lower.as_fraction()<=q<=interval.upper.as_fraction()
    for item in result.gradient:
        need(contains(item.enclosure,weights[item.variable]/d),'analytic gradient at rational center contained')
    for item in result.hessian:
        q=-weights[item.first_variable]*weights[item.second_variable]/d**2
        need(contains(item.enclosure,q),'analytic Hessian at rational center contained')
    return {'operation':'interval.expression.second_jet_enclosure.compute',
            'input':{'expression':source.model_dump(mode='json'),
                     'box':box.model_dump(mode='json'),'precision_bits':128},
            'output':result.model_dump(mode='json'),
            'rational_center_gradient_and_Hessian_crosscheck':True,
            'scope':'validated complete stated box for log-distance only; no new bound on full kernel residual'}


def run(data_path=DATA):
    installed=version('jacobian');need(installed=='0.23.0','pinned Jacobian version')
    data_path=Path(data_path).resolve()
    need(data_path.is_relative_to(ROOT),'recorded kernel input remains repository-relative')
    raw=data_path.read_bytes();data=json.loads(raw)
    need((data['n'],data['mass'],data['feature_order'],data['harmonics_through'])==(12,11,5,11),
         'fixed existing twelve-point fit')
    return {'passed':True,'jacobian_version':installed,
            'source_path':data_path.relative_to(ROOT).as_posix(),
            'source_sha256':hashlib.sha256(raw).hexdigest(),
            'native_Gram_inertia_crosschecks':gram_crosscheck(data),
            'native_radial_Bernstein_crosscheck':radial_bernstein_crosscheck(data),
            'native_validated_log_second_jet':logarithm_second_jet(),
            'new_global_Mordell_bound_proved':False,
            'scope':'actual optional Jacobian native calculations; existing independent proof acceptors retained'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record',action='store_true')
    parser.add_argument('--data',type=Path,default=DATA)
    parser.add_argument('--output',type=Path,default=HERE/'native_report.json')
    args=parser.parse_args();report=run(args.data)
    if args.record:
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
