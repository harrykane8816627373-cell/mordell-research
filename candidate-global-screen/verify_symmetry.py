"""Exact symmetry, radial-layer and contact certification from unique root boxes."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,lcm
import argparse,hashlib,importlib.util,json,sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DEPENDENCIES={
    'validated-maxima/verify_all.py':'4875024b8eacdc48fefdb83b447655270c2fbf8bcf785ddf53337526d5dcc6e0',
    'validated-maxima/dyadic.py':'94ab5e20011cf2742566e7aaf0af71c41e7fc4254e2320dfb4b3af97fcf07ccc',
    'validated-maxima/certificates.json':'d7c5348ea57a65dc49be8437993aa0b23a14a0f0bb5e2e5f0e7c9d84a959df7d',
}


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load_prior():
    for name,expected in DEPENDENCIES.items():need(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,'pinned dependency '+name)
    spec=importlib.util.spec_from_file_location('dyadic',ROOT/'validated-maxima/dyadic.py')
    dyadic=importlib.util.module_from_spec(spec);spec.loader.exec_module(dyadic)
    old=sys.modules.get('dyadic');sys.modules['dyadic']=dyadic
    try:
        spec=importlib.util.spec_from_file_location('screen_prior_exact_verifier',ROOT/'validated-maxima/verify_all.py')
        prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
    finally:
        if old is None:del sys.modules['dyadic']
        else:sys.modules['dyadic']=old
    return prior


def sqrt_ball(q,Ball,SCALE):
    need(q.lo>=0,'nonnegative square-root interval')
    low=isqrt(q.lo*SCALE);high=isqrt(q.hi*SCALE)
    if high*high<q.hi*SCALE:high+=1
    return Ball(low,high)


def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))


def permutation_order(p):
    need(sorted(p)==list(range(len(p))),'bijective permutation')
    seen=set();order=1
    for i in range(len(p)):
        if i in seen:continue
        j=i;length=0
        while j not in seen:seen.add(j);length+=1;j=p[j]
        order=lcm(order,length)
    return order


def pow_bound(text):
    if text=='0':return F(0)
    need(text.startswith('2^'),'certified power upper bound')
    return F(2)**int(text[2:])


def group_orbits(group,n,pairs=False):
    objects={(i,j) for i in range(n) for j in range(i)} if pairs else set(range(n))
    orbits=[]
    while objects:
        obj=min(objects)
        orbit={tuple(sorted((g[obj[0]],g[obj[1]]),reverse=True)) for g in group} if pairs else {g[obj] for g in group}
        need(orbit<=objects,'disjoint complete group orbit')
        objects-=orbit;orbits.append(sorted(orbit))
    return orbits


def axis_moment(points,orbit,axis,power,Ball):
    """Interval for sum (z_i*conj(w))**power, without square roots."""
    a,b=axis;real=Ball.integer(0);imag=Ball.integer(0)
    for i in orbit:
        x,y=points[i];u=x*a+y*b;v=y*a-x*b
        p,q=Ball.integer(1),Ball.integer(0)
        for _ in range(power):p,q=p*u-q*v,p*v+q*u
        real=real+p;imag=imag+q
    return real,imag


def axis_reflection(tau,psi,w):
    """A certified geometric reflection fixing the chosen nonzero point."""
    rotation=tuple(range(len(tau)))
    for _ in range(permutation_order(tau)):
        reflection=compose(rotation,psi)
        if reflection[w]==w:
            need(permutation_order(reflection)==2,'nonidentity reflection fixing harmonic axis')
            return reflection
        rotation=compose(tau,rotation)
    raise ArithmeticError('no certified reflection fixes the harmonic axis')


def verify_row(entry,prior,common):
    n=entry['n'];k=entry['rotation_order'];tau=tuple(entry['rotation_permutation']);psi=tuple(entry['reflection_permutation'])
    need(len(tau)==len(psi)==n and permutation_order(tau)==k,'exact rotation permutation order')
    need(permutation_order(psi)==2,'reflection involution order')
    inverse=tuple(tau.index(i) for i in range(n))
    need(compose(psi,compose(tau,psi))==inverse,'dihedral permutation relation')
    row=next(r for r in common['rows'] if r['n']==n)
    accepted=prior.verify_row(row,common)
    Ball=prior.Ball;SCALE=prior.SCALE;precision=prior.PRECISION
    eta=pow_bound(accepted['preconditioned_residual_upper']);q=pow_bound(accepted['contraction_norm_upper'])
    need(q<1,'strict contraction')
    error=eta/(1-q);err=-((-error.numerator*SCALE)//error.denominator)
    radius=1<<(precision-common['radius_exponent'])
    need(err<radius,'a posteriori root enclosure fits the original uniqueness box')
    mid=[Ball.dyadic(v,common['center_exponent']) for v in row['center_numerators']]
    refined=[Ball(v.lo-err,v.hi+err) for v in mid]
    flat=refined[:1]+[Ball.integer(0)]+refined[1:]
    points=[(flat[2*i],flat[2*i+1]) for i in range(n)]
    images=[]
    for kind,p in (('rotation',tau),('reflection',psi)):
        a,b=points[p[0]];rho=sqrt_ball(a.square()+b.square(),Ball,SCALE)
        need(rho.lo>0,'nonzero point for exact regauging')
        mapped=[]
        for i in range(n):
            x,y=points[p[i]]
            real=(a*x+b*y)/rho
            imag=(a*y-b*x)/rho if kind=='rotation' else (b*x-a*y)/rho
            mapped.extend((real,imag))
        # The first imaginary coordinate is algebraically zero; omitting it is
        # exact, and avoids treating dependent interval products as independent.
        retained=mapped[:1]+mapped[2:]
        need(retained[0].lo>0,'regauged first coordinate positive')
        need(all(v.lo>c.lo-radius and v.hi<c.hi+radius for v,c in zip(retained,mid)),
             kind+' maps the refined actual-root enclosure inside the certified uniqueness slice')
        distance=max((v-c).magnitude() for v,c in zip(retained,mid))
        images.append({'kind':kind,'permutation':list(p),'image_distance_from_old_midpoint_upper':prior.power_upper(distance),
                       'image_in_original_uniqueness_box':True})
    identity=tuple(range(n));group={identity};pending=[identity]
    while pending:
        g=pending.pop()
        for h in (tau,psi):
            c=compose(g,h)
            if c not in group:group.add(c);pending.append(c)
    need(len(group)==2*k,'faithful dihedral permutation group')
    orbits=group_orbits(group,n)
    radii=[x.square()+y.square() for x,y in points]
    radial=[(radii[o[0]],o) for o in orbits];radial.sort(key=lambda item:item[0].lo+item[0].hi)
    for (a,_),(b,_) in zip(radial,radial[1:]):need(a.hi<b.lo,'distinct exact radius layers between group orbits')
    pairorbits=group_orbits(group,n,True);contacts=[]
    for orbit in pairorbits:
        i,j=orbit[0];u,v=radii[i],radii[j]
        if u.lo+u.hi>v.lo+v.hi:u,v=v,u
        x,y=points[i];a,b=points[j];t=x*a+y*b
        contacts.append((u,v,t))
    for i,c in enumerate(contacts):
        for d in contacts[:i]:need(any(a.hi<b.lo or b.hi<a.lo for a,b in zip(c,d)),
                                     'distinct O2 contact triples between dihedral pair orbits')
    zero_fixed=[i for i in range(n) if tau[i]==i]
    result={'n':n,'passed':True,'original_root_classification':accepted['classification'],
            'rotation_order_exact':k,'reflection_exact':True,'dihedral_group_order_exact':len(group),
            'root_refinement_error_upper':str(error),'root_refinement_error_power_upper':prior.power_upper(err),
            'regauged_images':images,'exact_radial_layer_populations':[len(o) for _,o in radial],
            'radial_layers':[{'point_indices':o,'squared_radius_interval':[str(F(r.lo,SCALE)),str(F(r.hi,SCALE))]} for r,o in radial],
            'exact_origin_point_indices':zero_fixed,'O2_pair_contact_types_exact':len(pairorbits),
            'pair_contact_orbits':[{'representative':list(o[0]),'multiplicity':len(o)} for o in pairorbits],
            'vanishing_complex_moment_degrees_1_through_18':[j for j in range(1,19) if j%k],
            'radially_weighted_nonmultiple_harmonics_vanish_for_every_radial_weight':True,
            'symmetry_of_all_global_maximizers_claimed':False}
    if n==10:
        need([len(o) for _,o in radial]==[2,4,4],'D2 rectangle layer populations')
        w=radial[0][1][0];a,b=points[w];a0=radii[w]
        need(a0.lo>0 and tau[w]!=w,'nonzero inner antipodal axis')
        # One of the two certified reflections fixes w. Its geometric axis
        # is therefore the line through w, so all layer second moments in
        # that gauge are real, as an exact consequence of the group action.
        fixing=psi if psi[w]==w else compose(tau,psi)
        need(fixing[w]==w and permutation_order(fixing)==2,'reflection fixes inner axis')
        rectangle=[]
        for radius,orbit in radial[1:]:
            numerator=sum(((x.square()-y.square())*(a.square()-b.square())+4*x*y*a*b
                           for x,y in (points[i] for i in orbit)),Ball.integer(0))
            c=numerator/a0
            need(c.lo>0 or c.hi<0,'outer-layer real second moment nonzero')
            need(c.square().hi<(4*radius).square().lo,'outer rectangle second moment strictly below collinear endpoint')
            rectangle.append({'point_indices':orbit,
                'real_layer_second_moment_interval':[str(F(c.lo,SCALE)),str(F(c.hi,SCALE))],
                'nonzero':True,'absolute_value_strictly_less_than_four_squared_radius':True})
        result['D2_rectangle_nondegeneracy']={'normalization':'sum |z_i|^2=10',
            'inner_axis_point_index':w,'inner_axis_squared_modulus_positive':True,
            'reflection_fixing_inner_axis':list(fixing),
            'inner_layer_second_moment_exact':'2*inner_squared_radius',
            'all_layer_second_moments_real_in_inner_axis_gauge':True,
            'outer_layers':rectangle}
    if n==12:
        need([len(o) for _,o in radial]==[3,3,6],'D3 harmonic layer populations')
        w=radial[0][1][0];a0=radii[w]
        need(a0.lo>0,'nonzero inner triangular axis')
        fixing=axis_reflection(tau,psi,w);moments=[]
        for layer,(radius,orbit) in enumerate(radial):
            need(radius.lo>0,'positive D3 layer radius')
            c,imag=axis_moment(points,orbit,points[w],3,Ball)
            need(imag.lo<=0<=imag.hi,'real third moment implied by certified axis reflection')
            need(c.lo>0 or c.hi<0,'layer third moment numerator nonzero')
            endpoint=radius*a0
            if layer==2:
                need(c.square().hi<(36*endpoint.square()*endpoint).lo,
                     'outer third moment numerator strictly below coincident-triangle endpoint')
            moment={'point_indices':orbit,
                'real_axis_weighted_third_moment_interval':[str(F(c.lo,SCALE)),str(F(c.hi,SCALE))],
                'nonzero':True,
                'regular_triangle_by_certified_rotation':len(orbit)==3}
            if layer==2:moment['square_strictly_below_36_times_radius_product_cubed']=True
            else:moment['regular_triangle_moment_modulus_squared_exact']='9*(layer_squared_radius*axis_squared_radius)^3'
            moments.append(moment)
        result['D3_third_moment_nondegeneracy']={'normalization':'sum |z_i|^2=12',
            'axis_point_index':w,'axis_squared_modulus_positive':True,
            'reflection_fixing_axis':list(fixing),
            'moment_definition':'sum_i_in_layer (z_i*conj(w))^3',
            'all_layer_third_moments_real_in_axis_gauge':True,
            'layers':moments}
    if n==19:
        need([len(o) for _,o in radial]==[1,6,6,6],'D6 harmonic layer populations')
        w=radial[1][1][0];a0=radii[w]
        need(a0.lo>0,'nonzero hexagonal axis')
        fixing=axis_reflection(tau,psi,w);moments=[]
        for radius,orbit in radial[1:]:
            need(radius.lo>0,'positive hexagonal layer radius')
            c,imag=axis_moment(points,orbit,points[w],6,Ball)
            need(imag.lo<=0<=imag.hi,'real sixth moment implied by certified axis reflection')
            need(c.lo>0 or c.hi<0,'hexagonal layer sixth moment numerator nonzero')
            moments.append({'point_indices':orbit,
                'real_axis_weighted_sixth_moment_interval':[str(F(c.lo,SCALE)),str(F(c.hi,SCALE))],
                'nonzero':True,'regular_hexagon_by_certified_rotation':True})
        result['D6_sixth_moment_structure']={'normalization':'sum |z_i|^2=19',
            'axis_point_index':w,'axis_squared_modulus_positive':True,
            'reflection_fixing_axis':list(fixing),
            'moment_definition':'sum_i_in_layer (z_i*conj(w))^6',
            'all_positive_layer_sixth_moments_real_in_axis_gauge':True,
            'regular_hexagon_layer_sixth_moment_modulus_squared_exact':'36*(layer_squared_radius*axis_squared_radius)^6',
            'positive_layers':moments}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record-report',action='store_true');a=p.parse_args()
    prior=load_prior();common=json.loads((ROOT/'validated-maxima/certificates.json').read_bytes())
    entries=json.loads((HERE/'symmetry_certificates.json').read_text())['rows']
    need([(r['n'],r['rotation_order']) for r in entries]==[(10,2),(12,3),(19,6)],'fixed exact symmetry targets')
    results=[verify_row(e,prior,common) for e in entries]
    report={'passed':True,'date':'2026-11-08','arithmetic':'standard-library integer dyadic intervals and integer isqrt',
            'prior_root_acceptance_recomputed':True,'pinned_dependency_sha256':DEPENDENCIES,
            'certificate_sha256':hashlib.sha256((HERE/'symmetry_certificates.json').read_bytes()).hexdigest(),
            'rows':results,'global_optimality_claimed':False}
    target=HERE/'symmetry_report.json'
    if a.record_report:target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    else:need(report==json.loads(target.read_text()),'fresh exact symmetry report reproduced')
    print(json.dumps({'passed':True,'rows':[{k:r[k] for k in ('n','rotation_order_exact','dihedral_group_order_exact','exact_radial_layer_populations','O2_pair_contact_types_exact')} for r in results]},indent=2))


if __name__=='__main__':main()
