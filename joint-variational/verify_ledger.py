"""Fractional edge credits reject double use of the same unordered pair."""
from fractions import Fraction as F
from shared_budget import credit_ledger
from acceptance import require,write_report


def run():
    claims=[{'edges':[(0,2),(2,3)],'weight':F(1,2),'per_pair_lower':1},
            {'edges':[(0,2),(3,4)],'weight':F(1,2),'per_pair_lower':1}]
    row=credit_ledger(5,claims,[(0,1)])
    require(row['disjoint_pair_lower']==2,'shared edge credited fractionally')
    naive=4;true_distinct=3
    require(naive>true_distinct>=row['disjoint_pair_lower'],'explicit naive double-count counterexample')
    invalid=0
    cases=[(5,[{'edges':[(0,1)],'weight':1,'per_pair_lower':1}],[(1,0)]),
           (5,[{'edges':[(0,2)],'weight':1,'per_pair_lower':1}]*2,[]),
           (5,[{'edges':[(0,2),(2,0)],'weight':F(1,2),'per_pair_lower':1}],[]),
           (5,[{'edges':[(0,5)],'weight':1,'per_pair_lower':1}],[]),
           (5,[{'edges':[(0,2)],'weight':-1,'per_pair_lower':1}],[])]
    for args in cases:
        try:credit_ledger(*args)
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid pair credit accepted')
    return {'passed':True,'fractional_credit':row,'naive_sum':naive,'true_distinct_pair_budget':true_distinct,
            'invalid_credits_rejected':invalid,
            'scope':'Validation of pair allocation only. Claimed distance/cost bounds and membership outside the reserved cycle edges remain geometric hypotheses.'}


if __name__=='__main__':write_report('ledger_verification.json',run())
