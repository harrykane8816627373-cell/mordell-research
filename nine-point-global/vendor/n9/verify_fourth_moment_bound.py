"""Exact and rigorous verification of the improved all-local-maxima radius bound."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for ancestor in ROOT.parents:
    dependency=ancestor/'friedman_full_audit_20261023/deps'
    if dependency.is_dir():sys.path.insert(0,str(dependency));break
import sympy as S
from flint import arb,ctx,fmpq
ctx.prec=192
x,y,s,W,Y,T=S.symbols('x y s W Y T')
assert S.cancel((x**4-y**4)/(x-y)**2-(4*x**3/(x-y)-3*x*x-2*x*y-y*y))==0
assert 9-2*(16*S.Rational(15,8)-21)==-9
assert 256-16==240
assert 441+9*27==684 and 9*324==2916
assert [240*S.Rational(9,8)**4*S.Rational(256,243),-672*S.Rational(9,8)**3*S.Rational(256,243),684*S.Rational(9,8)**2*S.Rational(256,243),-2916*S.Rational(256,243)]==[405,-1008,912,-3072]
V=S.Rational(256,27);s0=S.Rational(17,10)
gap=s**3-S.Rational(56,45)*s*s+S.Rational(100,9)-S.Rational(113,90)*V
assert gap.subs(s,s0)>0
assert S.diff(gap,s).subs(s,s0)>0 and S.diff(gap,s,2).subs(s,s0)>0
aa=S.Rational(5,2)*s*s-s**3-S.Rational(100,9)
assert aa.subs(s,s0)<0 and S.diff(aa,s).subs(s,s0)<0 and S.diff(aa,s,2).subs(s,s0)<0
H=405*(s**4+(Y-s**3)**2/(W-s*s))-1008*Y+912*W-3072
G=S.cancel(H.subs(Y,S.Rational(5,2)*W-S.Rational(100,9)))
expanded=405*s**4+2025*aa-1608*s*s+8128+S.Rational(3693,4)*T+405*aa**2/T
assert S.cancel(G.subs(W,T+s*s)-expanded)==0
c=2025-9*S.sqrt(18465)
rho=c/810;kappa=504-405*rho;L=S.simplify((3072-S.Rational(200,9)*kappa)/405)
assert S.simplify(-1008+2*kappa+810*rho)==0
assert S.simplify(912-5*kappa-405*rho**2)==0
assert S.simplify(-3072+S.Rational(200,9)*kappa+405*L)==0
assert (9*S.sqrt(18465))**2==4*405*S.Rational(3693,4)
poly=405*s**4+c*aa-1608*s*s+8128
assert S.expand(S.diff(poly,s)-s*(1620*s*s-3*c*s+5*c-3216))==0
cb=arb(2025)-9*arb(18465).sqrt();assert cb>800 and cb<803
rhob=cb/810;kb=504-405*rhob;Lb=(3072-arb(200)/9*kb)/405
assert rhob>0 and kb>0 and Lb>0
beta_ball=(rhob+(rhob*rhob+4*Lb.sqrt()).sqrt())/2
assert 1620*s0*s0-3*803*s0+5*800-3216>0
assert 3240*s0-3*803>0
def ball(q):return arb(fmpq(str(q)))
def evaluate(q):
    b=ball(q);a=arb(5)/2*b*b-b*b*b-arb(100)/9
    return 405*b*b*b*b+cb*a-1608*b*b+8128
lo=S.Rational('1.77410189294');hi=S.Rational('1.77410189295')
assert evaluate(s0)<0 and evaluate(lo)<0 and evaluate(hi)>0
assert beta_ball>ball(lo) and beta_ball<ball(hi)
U=S.Rational(S.ceiling(hi*2**20),2**20)
assert U>=hi and evaluate(U)>0 and U<S.Rational(9,5)
rec={'global_M9_determined':False,'all_positive_local_maxima_bound_proved':True,
     'scaling':'x=sqrt(8/9) z; u=|x|²','fourth_moment_necessary_inequality':'405 T4 - 1008 T3 + 912 T2 - 3072 <= 0',
     'c':'2025-9*sqrt(18465)','c_ball':str(cb),'increasing_polynomial':str(S.expand(poly)),
     'radial_rho':str(rho),'radial_kappa':str(kappa),'radial_L':str(L),
     'radial_nonnegative_sum_inequality':'sum r_i^2*(r_i-rho)^2 <= L',
     'beta_closed_form':'(rho+sqrt(rho^2+4*sqrt(L)))/2','beta_ball':str(beta_ball),
     'beta_unique_root_above_17_over_10_interval':[str(lo),str(hi)],
     'beta_lower_polynomial_ball':str(evaluate(lo)),'beta_upper_polynomial_ball':str(evaluate(hi)),
     'upper_squared_radius_x':str(U),'upper_squared_radius_x_float':float(U),
     'upper_squared_radius_z':str(S.Rational(9,8)*U),'upper_squared_radius_z_float':float(S.Rational(9,8)*U),
     'mathematical_proof':'FOURTH_MOMENT_RADIUS_PROOF.txt','working_precision_bits':ctx.prec}
(ROOT/'fourth_moment_radius_certificate.json').write_text(json.dumps(rec,indent=2))
print(json.dumps(rec,indent=2))
