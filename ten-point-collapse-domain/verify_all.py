"""Exact necessary domain certificates for ten points at mass nine."""
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial, prod
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
L = F(97, 80)*9**9
U_OLD = F(15403, 10000)
U = F(1539, 1000)
D0 = F(2, 5)
EPS = F(1, 500)
POP_CUTOFFS = dict(zip(range(2, 10), map(lambda r: F(r, 10000),
                  [1007, 2044, 3133, 4248, 5296, 6283, 7253, 8195])))
# The sharp k=m(m-1)/2 homogeneous bound is Delta_m <= C_m q^k.
SHARP = {2:F(2), 3:F(1), 4:F(1,16), 5:F(1,5**5), 6:F(1,5**10),
         7:F(1,6**15), 8:F(1,7**21), 9:F(1,8**28)}


def need(ok, label):
    if not ok:
        raise ArithmeticError(label)


def mul(p, q):
    b = [F(0)]*(len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            b[i+j] += x*y
    return b


def bernstein(p, lo, hi):
    n=len(p)-1
    q=[sum((p[j]*comb(j,i)*lo**(j-i)*(hi-lo)**i for j in range(i,n+1)), F(0))
       for i in range(n+1)]
    return [sum((q[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),F(0))
            for i in range(n+1)]


def envelope_two():
    # The odd Taylor polynomial exp_23(x) <= exp(x) for every real x.
    # Verify exp_23(11/2*(s-1)+eps) - (s^2+...+s^9)/8 > 0 on [0,2].
    n=23
    a,b=F(11,2), EPS-F(11,2)
    p=[sum((F(comb(j,k),factorial(j))*a**k*b**(j-k) for j in range(k,n+1)),F(0))
       - (F(1,8) if 2<=k<=9 else F(0)) for k in range(n+1)]
    stack=[(F(0),F(2),0)]
    nodes=leaves=depthmax=0
    least=None
    while stack:
        lo,hi,depth=stack.pop()
        nodes+=1
        depthmax=max(depthmax,depth)
        coeff=bernstein(p,lo,hi)
        if min(coeff)>0:
            leaves+=1
            least=min(coeff) if least is None else min(least,min(coeff))
        else:
            need(depth<20, "weighted pair envelope unresolved")
            mid=(lo+hi)/2
            stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
    need(9/F(2)<a, "decreasing log deficit beyond 2 via max degree")
    need(0<8*EPS<1, "geometric exponential upper bound")
    return {"degree":n, "nodes":nodes, "leaves":leaves, "depth":depthmax,
            "min_Bernstein_coefficient_decimal_display":float(least)}


def envelope_general():
    # (2P'-(m+9)P) = -s^(m-1)(s-1)R_m(s), R_m coefficients nonnegative.
    out={}
    for m in range(3,10):
        p=[F(int(k>=m)) for k in range(10)]
        lhs=[2*(j+1)*(p[j+1] if j+1<=9 else F(0))-(m+9)*p[j]
             for j in range(10)]
        q=[F(2*m+j*(m-9)+j*(j+1)) for j in range(10-m)]
        need(q[0]>0 and min(q)>=0, "log derivative sign")
        need(lhs==mul([F(0)]*(m-1)+[F(1),F(-1)],q), "exact factor identity")
        out[str(m)]=[int(c) for c in q]
    return out


def exp_upper(x, k=25):
    # Tail starting at k+1 is geometric dominated.
    need(0<=x<k+2, "exp comparison range")
    return sum((x**j/factorial(j) for j in range(k+1)),F(0))+(
           x**(k+1)/factorial(k+1)/(1-x/(k+2)))


def radius():
    need(U_OLD-D0/10 < U < U_OLD, "radius activates zero-free case")
    q=9-F(10,9)*U
    delta=5*q-36
    need(0<delta<2, "endpoint range")
    # d/du log(F9) = -40/q+9/u-100/9; first term is negative.
    need(9/U<F(100,9), "decreasing log upper in final interval")
    endpoint=q**36*F(1,8**28)*(F(100,81)*U)**9*exp_upper(2*delta)
    need(endpoint<L, "zero-free Schwarz radius endpoint")
    return float(endpoint/L)


def upper_weight(m):
    return 1/(1-8*EPS) if m==2 else F(1)


def schwarz_upper(m, s, S):
    k=m*(m-1)//2
    h=10-m
    delta=5*s-k
    need(0<=delta<2 and S>=0, "collapse range")
    return SHARP[m]*s**k * upper_weight(m)*h**h*(S/h)**(45-k)*(
           (delta+2)/(2-delta))**2


def exclude_population():
    out={}
    for m,r in POP_CUTOFFS.items():
        k=m*(m-1)//2
        lo,hi=F(k,5),m*r
        need(lo<hi<lo+F(2,5), "population interval")
        worst=F(0)
        stack=[(lo,hi,0)]
        nodes=0
        # s^k and (delta+2)/(2-delta) increase, S=9+10r-(1+10/m)s decreases.
        while stack:
            s0,s1,depth=stack.pop()
            nodes+=1
            S0=9+10*r-(1+F(10,m))*s0
            upper=schwarz_upper(m,s1,S0)
            if upper<L:
                worst=max(worst,upper)
            else:
                need(depth<20, f"population cell upper m={m}, ratio={float(upper/L)}")
                mid=(s0+s1)/2
                stack.extend([(s0,mid,depth+1),(mid,s1,depth+1)])
        out[str(m)]={"at_most":m-1, "u_cutoff":str(r),"strict_uncentered_mass_lower":str(m*r),"verified_nodes":nodes,
                "largest_upper_over_L_decimal_display":float(worst/L)}
    return out


def pair_excluded(d_lo,d_hi,a2_hi):
    if not (D0<=d_lo<=d_hi<F(6,5) and a2_hi>=0):
        return False
    return schwarz_upper(2, d_hi/2, 9-d_lo/2+10*a2_hi)<L


def pair_examples():
    a2max=F(1,10000)
    dmax=F(203,500)
    need(1/dmax-F(22,1)/(9-dmax/2)+5>0, "fixed midpoint monotonicity in d")
    endpoint=schwarz_upper(2,dmax/2,9-dmax/2+10*a2max)
    need(endpoint<L,"short-distance midpoint rectangle")
    u=v=F(1013,10000)
    t=F(-253,2500)
    r=F(1,100000)
    d_lo=u+v-2*t-4*r
    d_hi=u+v-2*t+4*r
    a2_hi=(u+v+2*t+4*r)/4
    q_lo=9-F(9,8)*(u+v+2*r)-(t+r)/4
    need(u-r>0 and u+r<U_OLD and abs(t)+r<u-r and d_lo>D0
         and q_lo>=F(6),"whole witness cube passes previous domain")
    need(pair_excluded(d_lo,d_hi,a2_hi) and pair_box_excluded((u-r,v-r,t-r),(u+r,v+r,t+r)), "whole witness cube excluded")
    return {"midpoint_squared_max":str(a2max),"distance_squared_max":str(dmax),
            "endpoint_over_L_decimal_display":float(endpoint/L),
            "strict_example_center":[str(x) for x in (u,v,t)],"example_cube_radius":str(r)}


def old_input_witness():
    path=ROOT/"symmetry-breaking/integer_witnesses.json"
    raw=path.read_bytes()
    expected="73e5f42ef78db59dd0ff3e05227caf75eb3834a02486c509767c8534c3c954d5"
    need(hashlib.sha256(raw).hexdigest()==expected,"fixed integer witness")
    p=next(row["integer_points"] for row in json.loads(raw)["rows"] if row["n"]==10)
    sx,sy=sum(x for x,y in p),sum(y for x,y in p)
    q=10*sum(x*x+y*y for x,y in p)-sx*sx-sy*sy
    P=prod((x-u)**2+(y-v)**2 for i,(x,y) in enumerate(p) for u,v in p[i+1:])
    V=F(90**45*P,q**45)
    need(V>L,"lower threshold witness")
    centered=[(10*x-sx,10*y-sy) for x,y in p]
    scale=F(9,10*q)
    radii=[scale*(x*x+y*y) for x,y in centered]
    need(sum(radii)==9 and max(radii)<U,"saved witness radii retained")
    need(all(sum(sorted(radii)[:m])>m*r for m,r in POP_CUTOFFS.items()),
         "all saved witness subsets retained")
    count=0
    for i,(x,y) in enumerate(centered):
        for j,(w,z) in enumerate(centered[i+1:],i+1):
            pair=(radii[i],radii[j],scale*(x*w+y*z))
            need(not pair_box_excluded(pair,pair),"saved witness pair retained")
            count+=1
    return {"V_over_L_decimal_display":float(V/L),"source_sha256":expected,
            "point_pairs_retained":count,"subsets_retained":sum(comb(10,m) for m in POP_CUTOFFS)}


def pinned_inputs():
    manifest=json.loads((HERE/"source_manifest.json").read_bytes())
    for base,key in [(HERE,"files"),(ROOT,"dependencies")]:
        for name,sha in manifest[key].items():
            need(hashlib.sha256((base/name).read_bytes()).hexdigest()==sha,
                 f"fixed {key} input: {name}")
    for suite in ("eight-point-global","nine-point-global"):
        base=ROOT/suite
        report=json.loads((base/"reports/replay.json").read_bytes())
        manifest_sha=hashlib.sha256((base/"source_manifest.json").read_bytes()).hexdigest()
        accept=json.loads((base/"accepted_replay.json").read_bytes())
        need(accept["source_manifest_sha256"]==manifest_sha and report["passed"] is True
             and report["proof_replay_performed"] is True and report["unresolved_boxes"]==0,
             f"accepted prior full replay: {suite}")
        for name,sha in accept["files_sha256"].items():
            need(hashlib.sha256((base/name).read_bytes()).hexdigest()==sha,
                 f"fixed theorem acceptance evidence: {suite}/{name}")
    return len(manifest["dependencies"])


def pair_box_excluded(lo,hi):
    """Sufficient exact pruning on (u,v,t) closed boxes, only by the new cuts."""
    lo,hi=tuple(map(F,lo)),tuple(map(F,hi))
    need(len(lo)==len(hi)==3 and all(a<=b for a,b in zip(lo,hi)), "closed pair box")
    u,v,t=lo
    U1,V1,T1=hi
    if u>=U or v>=U or U1+V1<=2*POP_CUTOFFS[2]:
        return True
    return pair_excluded(max(D0,u+v-2*T1),U1+V1-2*t,(U1+V1+2*T1)/4)


def verify():
    dep_count=pinned_inputs()
    return {"passed":True,
         "scope":"exact necessary conditions for each global ten-point maximizer; no symmetry",
         "normalization":"sum |x_i|^2=9", "fixed_dependency_count":dep_count,
         "input_witness":old_input_witness(),"radius_strict_upper":str(U),
         "radius_endpoint_over_L_decimal_display":radius(),
         "pair_weighted_envelope":envelope_two(),"higher_weighted_envelope_coefficients":envelope_general(),
         "subset_mass_and_population":exclude_population(),"pair_prune":pair_examples(),
         "full_global_sharp_bound_proved":False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--record-report",action="store_true")
    args=parser.parse_args()
    out=verify()
    path=HERE/"verification_report.json"
    if args.record_report:
        path.write_bytes((json.dumps(out,indent=2)+"\n").encode("utf-8"))
    else:
        need(out==json.loads(path.read_bytes()), "saved verification report")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
