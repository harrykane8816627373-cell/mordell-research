"""Exact radial Gram implications for a ten-point Ward inner-count proof.

The caller must establish the two radial exclusion bands, the global Ward
pair slack, the upper bound on B-log(V), and the selected row inequality.
This module rechecks the Gram arithmetic and 1D polynomial inequalities.
"""
from fractions import Fraction as F
from math import comb

INNER = (F(0), F(3, 10))
OUTER = ((F(41, 60), F(71, 50)), (F(149, 100), F(1539, 1000)))
GLOBAL_SLACK = F(3, 10000)
U_MASS = F(9, 10)
# (intercept at u=.2044, required lower bound on the combined inner summand).
ROW_LEVELS = {"standard": (F(".004"), F(".0299")),
              "relaxed": (F(".0031"), F(".0290"))}
ROW_SLOPE = F(1, 10)
ROW_CENTER = F(".2044")

B01 = tuple(map(F, ("0", ".92045", "3.69876", "9.08955", "17.91077", "2.26600")))
A01 = F("-.62053")
B3 = tuple(map(F, ("0", ".57025", "1.92520", "3.81069", "6.18625", ".52160")))
A3 = F(".08514")


def need(ok, why):
    if not ok:
        raise ArithmeticError(why)


def mul(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):
            out[i+j] += x*y
    return out


def bernstein(p, lo, hi):
    n=len(p)-1
    q=[sum(p[j]*comb(j,i)*lo**(j-i)*(hi-lo)**i for j in range(i,n+1))
       for i in range(n+1)]
    return [sum(q[j]*F(comb(i,j),comb(n,j)) for j in range(i+1))
            for i in range(n+1)]


def strictly_above(p, interval, floor):
    todo=[(*interval,0)]
    nodes=splits=leaves=maximum_depth=0
    while todo:
        lo,hi,depth=todo.pop()
        maximum_depth=max(maximum_depth,depth)
        nodes+=1
        if min(bernstein(p,lo,hi))>floor:
            leaves+=1
        else:
            need(depth<20, "unresolved exact univariate Gram inequality")
            mid=(lo+hi)/2
            todo.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
            splits+=1
    need(nodes==2*splits+1 and leaves==splits+1,"complete Gram cover")
    return {"interval":list(map(str,interval)),"strict_floor":str(floor),
            "visited":nodes,"maximum_depth":maximum_depth}


def radial_data(data):
    hx=lambda t:F(float.fromhex(t))
    Q=[list(map(hx,row)) for row in data["Q_blocks"][0]]
    means=list(map(hx,data["radial_moment_centers"]))
    nodes=list(map(hx,data["radial_nodes"]))
    need(len(Q)==6 and all(len(row)==6 for row in Q) and len(means)==5
         and len(nodes)==3,"fixed radial feature basis")
    need(Q==[list(row) for row in zip(*Q)],"symmetric radial Gram block")
    L=[row[:] for row in Q]
    for k in range(6):
        need(L[k][k]>0,"strict exact Gram LDL pivot")
        for i in range(k+1,6):
            for j in range(k+1,i+1):
                L[i][j]-=L[i][k]*L[j][k]/L[k][k]
    return Q,means,nodes


def polynomial(radial, b, a):
    Q,means,nodes=radial
    qb=[sum(x*y for x,y in zip(row,b)) for row in Q]
    q=[F(1)]
    for x in nodes:
        q=mul(q,[-x,F(1)])
    p=[10*qb[-1]*c for c in mul(q,q)]
    for j in range(1,6):
        p[0]-=qb[j-1]*means[j-1]
        p[j]+=qb[j-1]
    p[0]+=a*U_MASS
    p[1]-=a
    return p,sum(x*y for x,y in zip(qb,b))/4


def shift(p, intercept, slope=F(0)):
    q=p[:]
    q[0]+=intercept
    q[1]+=slope
    return q


def row_spec(profile="standard", slack=GLOBAL_SLACK):
    row_level,_=ROW_LEVELS[profile]
    baseline=-9*slack
    row_intercept=row_level-ROW_SLOPE*ROW_CENTER
    lo=(baseline-row_intercept)/ROW_SLOPE
    need(INNER[0]<lo<INNER[1],"row cut transition")
    return lo,INNER[1],row_intercept,ROW_SLOPE


def verify(data, gap_upper, slack=GLOBAL_SLACK, profile="standard"):
    """Return exact inequalities, conditional on the stated global/row inputs.

    gap_upper must be a proven upper endpoint of B-log(V) in mass nine.
    The extra row inequality must hold on row_spec(profile,slack)[:2],
    with full row floor intercept+slope*u from the last two fields.
    """
    gap_upper,slack=F(gap_upper),F(slack)
    need(0<slack<=GLOBAL_SLACK,"verified Ward slack range")
    radial=radial_data(data)
    p,penalty=polynomial(radial,B01,A01)
    floors01=(F("-.2667"),F(".03649"))
    cov01=[strictly_above(p,INNER,floors01[0])]+[
        strictly_above(p,I,floors01[1]) for I in OUTER]
    gram_budget=2*(gap_upper+45*slack)
    lower01={n:n*floors01[0]+(10-n)*floors01[1]-penalty for n in (0,1)}
    need(min(lower01.values())>gram_budget,"radial Gram excludes zero and one inner")

    p,penalty3=polynomial(radial,B3,A3)
    inner_floor=ROW_LEVELS[profile][1]
    outer_floor=F("-.0122")
    lo,hi,intercept,slope=row_spec(profile,slack)
    baseline=shift(p,-9*slack)
    cov3=[strictly_above(baseline,(INNER[0],lo),inner_floor),
          strictly_above(shift(p,intercept,slope),(lo,hi),inner_floor)]+[
          strictly_above(baseline,I,outer_floor) for I in OUTER]
    lower3=3*inner_floor+7*outer_floor-penalty3
    need(lower3>2*gap_upper,"radial Gram and row bounds exclude three inner")

    return {"passed":True,"scope":"radial count arithmetic conditional on supplied Ward and row bounds",
            "inner_count_excluded":[0,1,3],"profile":profile,
            "inner_band":list(map(str,INNER)),"outer_bands":[list(map(str,I)) for I in OUTER],
            "gap_upper":str(gap_upper),"gram_budget_upper":str(gram_budget),
            "row_required":[str(x) for x in (lo,hi,intercept,slope)],
            "zero_or_one":{"b":list(map(str,B01)),"a":str(A01),
                           "strict_gram_lower_decimal_display":{str(n):float(v) for n,v in lower01.items()},
                           "one_variable_covers":cov01},
            "three":{"b":list(map(str,B3)),"a":str(A3),
                     "strict_combined_lower_decimal_display":float(lower3),
                     "one_variable_covers":cov3}}

