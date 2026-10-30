"""Cross integrals by weak Laplacians, and shifted-Gram algebra in Hilbert space."""
from fractions import Fraction as F
from math import comb
from random import Random
from shared_budget import radial_packet,add_centered_holes,checked_gram,dot,solve,exact_certificate,interval_certificate
from affine_pair import fixed_line,pi_bounds
from acceptance import require,write_report


def independent_cross(j,a,T):
    S=T*T;first=second=F(0)
    for l in range(j+1):
        coeff=(-1)**l*comb(j,l)*S**(j-l)
        first+=coeff*a**(2*l)/F(2*l+1)
        integral=(1/a-1/T) if l==0 else (T**(2*l-1)-a**(2*l-1))/F(2*l-1)
        second+=coeff*integral
    return first-second+(1/a-1)*(S-a*a)**j


def run():
    crosses=0;norms=[]
    for a in (F(1,8),F(5,16),F(1,2)):
        last=None
        for degree in range(1,7):
            p=radial_packet(100,a,F(11,10),[0]*degree)
            require(p['c']==[independent_cross(j,a,F(11,10)) for j in range(1,degree+1)],'cross terms via independent distributional Laplacian')
            crosses+=degree
            h=checked_gram(p['H'],p['m'],p['c'],p['Hf'])[-1]
            require(last is None or h<=last,'nested radial tests reduce residual norm');last=h
        norms.append(last)
    rng=Random(20261030);identities=0
    phi=[[F(1),0,1,0],[0,1,1,1]]
    H=[[dot(a,b) for b in phi] for a in phi]
    for j in range(80):
        u=[F(rng.randrange(-8,9),5) for _ in range(4)]
        f=[F(rng.randrange(-8,9),7) for _ in range(4)]
        m=[dot(u,v) for v in phi];c=[dot(f,v) for v in phi];Hf=dot(f,f)
        J,p,h=checked_gram(H,m,c,Hf)[-3:];theta=F(j-40,37)
        shifted=[x+theta*y/2 for x,y in zip(m,c)]
        shifted_value=dot(shifted,solve(H,shifted))-theta*theta*Hf/4
        require(shifted_value==J+theta*p-theta*theta*h/4,'shifted Schur identity')
        require(dot(u,u)+theta*dot(u,f)>=shifted_value,'independent actual Hilbert energy');identities+=1
    line=fixed_line();packet=radial_packet(100,F(5,16),F(11,10),[0,0])
    exact=exact_certificate(packet,4,line)['shared_gap_lower']
    box=interval_certificate(packet,4,line,[(0,0),(F(-1,100),F(1,100))])
    for moment in (F(-1,100),0,F(1,100)):
        point=radial_packet(100,F(5,16),F(11,10),[0,moment])
        require(exact_certificate(point,4,line)['shared_gap_lower']>=box['shared_gap_lower'],'one dual covers entire moment interval')
    require(box['shared_gap_lower']>F(73,1000),'nontrivial whole-box bound')
    base=radial_packet(100,F(5,16),F(11,10),[0])
    hole=add_centered_holes(base,[F(1,10)])
    J=checked_gram(hole['H'],hole['m'],hole['c'],hole['Hf'])[-3]
    b=F(1,10);T=F(11,10)
    require(J==F(100**2,4)*b**4/(1-b**4/T**4),'centered-hole Schur matches independent closed form')
    hole_value=exact_certificate(hole,4,line)['shared_gap_lower']
    nested=add_centered_holes(base,[F(1,20),F(1,10)])
    require(exact_certificate(nested,4,line)['shared_gap_lower']>=hole_value,'overlapping-hole full Gram nesting')
    # Actual normalized square: q-center radius^2=3/8, smear radius^2=1/4.
    # Its centered energy hole b=1/10 is genuine since (b+1/2)^2<3/8.
    import sys
    from affine_pair import ROOT
    sys.path.insert(0,str(ROOT/'joint-bounds'))
    from subset_moment import log_interval
    square_m2=4*(F(3,8)**2+2*F(1,4)*F(3,8)+F(1,4)**2/3)-F(4,3)
    require(square_m2==F(1,16) and F(3,5)**2<F(3,8),'actual square moments and true smoothed hole')
    square_packet=add_centered_holes(radial_packet(4,F(1,2),F(13,10),[0,square_m2]),[F(1,10)])
    square_value=exact_certificate(square_packet,1,line)['shared_gap_lower']
    l2,u2=log_interval(2);l3,u3=log_interval(3)
    actual_gap_lower=14*l2-6*u3-3
    require(square_value<actual_gap_lower,'certificate below independently computed true square gap')
    invalid=0
    for args in [([[1,2],[0,1]],[0,0],[0,0],1),([[1,1],[1,1]],[0,0],[0,0],1),([[1]],[0],[2],1)]:
        try:checked_gram(*args)
        except ValueError:invalid+=1
        else:raise ArithmeticError('invalid augmented Gram accepted')
    return {'passed':True,'independent_radial_cross_checks':crosses,'hilbert_shift_identities':identities,
            'degree_six_residual_norms':norms,'zero_moment_exact_bound':exact,'whole_moment_box_certificate':box,
            'centered_hole_projection':J,'centered_hole_shared_bound':hole_value,
            'actual_square_shared_bound':square_value,'actual_square_gap_lower':actual_gap_lower,
            'invalid_grams_rejected':invalid,
            'scope':'The box is a conditional measurement box, not a claimed family of physical stationary point configurations.'}


if __name__=='__main__':write_report('projection_verification.json',run())
