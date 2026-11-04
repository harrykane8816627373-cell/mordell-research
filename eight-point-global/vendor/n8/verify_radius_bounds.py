"""Exact checks for two new necessary/exclusion bounds; not a global certificate."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data'
sys.path.insert(0,str(ROOT.parent.parent/'friedman_full_audit_20261023/deps'))
import sympy as s
x,y,t=s.symbols('x y t')
assert s.cancel((x**3-y**3)/(x-y)**2-3*x*x/(x-y)+2*x+y)==0
assert s.cancel((x**4-y**4)/(x-y)**2-4*x**3/(x-y)+3*x*x+2*x*y+y*y)==0
# Candidate minimum squared distance > 6/7.
trig=8*x**3+4*x**2-4*x-1
assert trig.subs(x,s.Rational(5,8))==s.Rational(1,64)
assert s.diff(trig,x).subs(x,s.Rational(1,2))>0
assert s.diff(trig,x,2).subs(x,s.Rational(1,2))>0
# Spectral gap gives at least one unit of quadratic loss.
assert F(10,7)**2>2
radius=F(1,16);tau=F(1,10)
assert 2*radius**2/F(6,7)<tau**2
remainder=F(2,3)*tau/(1-tau)*8/F(6,7)
assert remainder==F(56,81)
gap=1-remainder
assert gap==F(25,81)>0
# Cubic and quartic moment inequalities from ||A v|| <= (7/2)||v||.
assert s.expand(t*(s.Rational(21,2)*t-13)**2-s.Rational(49,4)*t**3)==98*t**3-273*t**2+169*t
assert s.expand(t**2*(14*t-18)**2-s.Rational(49,4)*t**4)==s.Rational(735,4)*t**4-504*t**3+324*t**2
assert s.expand(98*t*(t-s.Rational(39,28))**2-s.Rational(169,8)*t)==98*t**3-273*t**2+169*t
P=s.Poly(t*(t-s.Rational(39,28))**2-s.Rational(169,98),t)
roots=P.intervals(eps=s.Rational(1,10**20))
upper_roots=[(a,b) for (a,b),mult in roots if a>s.Rational(39,28)]
assert len(upper_roots)==1
root_interval=upper_roots[0]
assert root_interval[1]<s.Rational(227,100)
# Combining the fourth and sixth radial moments is substantially stronger.
u=s.symbols('u');U=s.Rational(32,3)
H=98*(t**3+(u-t*t)**2/(8-t))-273*u+1352
Pstrong=s.Poly(882*t**3-2352*t**2+1755*t-1496,t)
assert s.cancel((8-t)*H.subs(u,U)-s.Rational(8,9)*Pstrong.as_expr())==0
threshold=t*t+s.Rational(39,28)*(8-t)-U
assert threshold.subs(t,1)==s.Rational(1,12)
assert s.diff(threshold,t).subs(t,1)>0 and s.diff(threshold,t,2)>0
assert Pstrong.count_roots(1,8)==1
strong_interval=[(a,b)for (a,b),mult in Pstrong.intervals(eps=s.Rational(1,10**20))if a>=1 and b<=8][0]
assert all(Pstrong.eval(endpoint)!=0 for endpoint in strong_interval)
assert strong_interval[1]<s.Rational(421,200)
assert F(1,1)/(F(7,2)*F(7,8))==F(16,49)
# Exact sanity checks on the established candidate and the unstable octagon.
R4=F(64,7);R6=F(512,49);R8=F(4096,343)
assert 98*R6-273*R4+1352<0
assert F(735,4)*R8-504*R6+324*R4<0
assert F(735,4)*8-504*8+324*8==30>0
# Shifted inner point has the advertised distance ratio but a smaller discriminant.
ratio=F(9999999,10000000)**2*F(800,801)**28
assert 0<ratio<1
assert 7*(80**2+1)+7**2==8*5607
out={'global_optimality_proved':False,
     'candidate_orbit_radius':str(radius),'log_discriminant_loss_coefficient':str(gap),
     'candidate_distance_squared_strict_lower_bound':'6/7',
     'all_local_maxima_squared_radius_strict_upper_bound':'421/200',
     'strong_squared_radius_algebraic_bound_polynomial':str(Pstrong.as_expr()),
     'strong_squared_radius_algebraic_bound_interval':list(map(str,strong_interval)),
     'strong_squared_radius_algebraic_bound_approx':float(sum(strong_interval)/2),
     'all_local_maxima_minimum_distance_squared_lower_bound':'16/49',
     'cubic_envelope_root_interval':list(map(str,root_interval)),
     'shifted_inner_point_discriminant_ratio':str(ratio),
     'shifted_inner_point_discriminant_ratio_approx':float(ratio)}
(OUT/'quantitative_bounds_exact_certificate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items()if k not in ['remaining_radius_polynomials','shifted_inner_point_discriminant_ratio']},indent=2))
