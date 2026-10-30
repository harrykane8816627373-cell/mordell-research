"""Translation/rotation-invariant thin-strip population capacity.

The author's complete stability norm is the input, not a new result.
SUBSET_GEOMETRY.md proves the new odd-skip bound. All decisions are rational.
"""
from fractions import Fraction as F


def require_population(n,m,max_terms=100000):
    if type(n) is not int or n<2 or type(m) is not int or not 2<=m<=n:
        raise ValueError('integer n>=2 and 2<=population<=n required')
    planned=m//2
    if type(max_terms) is not int or max_terms<1 or planned>max_terms:
        raise ValueError('complete odd-skip sum exceeds term budget')


def odd_skip_sum(population,length_squared,width_squared=0,max_terms=100000):
    m=population
    require_population(m,m,max_terms)
    L2,W2=F(length_squared),F(width_squared)
    if L2<=0 or W2<0:
        raise ValueError('positive squared length and nonnegative squared width required')
    return sum((F((m-k)**3)/(min(k,m-k)**2*L2+W2*(m-k)**2)
                for k in range(1,m,2)),F(0))


def strip_capacity(n,population,length,width=0,use_global_distance=False,max_terms=100000):
    """Necessary condition for at least m points in ANY oriented rectangle.

    L is length in the ordering direction, W is transverse width. The
    stronger distance option is valid only for GLOBAL maxima. The default
    is valid for every local maximum. The thin-width hypothesis is needed.
    """
    m=population
    require_population(n,m,max_terms)
    if type(use_global_distance) is not bool:
        raise ValueError('global-distance flag must be boolean')
    L,W=F(length),F(width)
    if L<=0 or W<0:
        raise ValueError('positive length and nonnegative width required')
    separation=F(4,n-1) if use_global_distance else F(2*n,(n-1)**2)
    if 2*W*W>=separation:
        raise ValueError('thin-strip hypothesis 2*width^2<separation_squared is required')
    eta=1-2*W*W/separation
    total=odd_skip_sum(m,L*L,W*W,max_terms)
    h=m//2
    l=F(n-1,2)
    lower=eta*total
    upper_squared=l*l*h*(m-h)
    slack=upper_squared-lower*lower
    return {'n':n,'population':m,'length':L,'width':W,
            'separation_squared_lower':separation,'phase_real_factor_lower':eta,
            'odd_skip_reciprocal_lower':total,'cross_operator_lower':lower,
            'cross_operator_upper_squared':upper_squared,
            'squared_stability_slack':slack,'excluded':slack<0,
            'scope':('All global maxima' if use_global_distance else 'All local maxima')+
                    ' with at least the stated population in a translated/rotated rectangle; no regular spacing.'}


def line_length_squared_lower(n,population,max_terms=100000):
    """A rational lower bound for the span of m collinear points.

    Exact for even population in this criterion. For odd m the actual
    strip certificate uses sqrt(floor(m/2)*ceil(m/2)) and is stronger.
    """
    require_population(n,population,max_terms)
    m=population
    A=odd_skip_sum(m,1,0,max_terms)
    return 4*A/((n-1)*m)


def microscopic_strip_capacity(population,scaled_length_squared,scaled_width_squared,max_terms=100000):
    """n-independent GLOBAL capacity after lengths are multiplied by sqrt(n-1)."""
    m=population
    require_population(m,m,max_terms)
    L2,W2=F(scaled_length_squared),F(scaled_width_squared)
    if not 0<=W2<2:
        raise ValueError('global microscopic width must have squared width below two')
    total=odd_skip_sum(m,L2,W2,max_terms)
    eta=1-W2/2
    h=m//2
    slack=F(h*(m-h),4)-(eta*total)**2
    return {'population':m,'scaled_length_squared':L2,'scaled_width_squared':W2,
            'squared_stability_slack':slack,'excluded_from_global_maxima':slack<0,
            'scope':'For every n>=population, physical lengths are these scaled lengths divided by sqrt(n-1).'}


def explicit_two_hundred_points():
    """A deterministic nonempty witness; not an optimum construction.

    Twenty collinear points plus 90 centrally symmetric grid pairs.
    Choosing the grid norms fixes the total moment extremely close to200.
    The subsequent exact normalization preserves the old geometric bounds.
    """
    pool=sorted(((a,b) for b in range(1,10) for a in range(-9,10)
                 if a*a+b*b<=97),key=lambda p:(p[0]*p[0]+p[1]*p[1],p))
    norms=[a*a+b*b for a,b in pool]
    selected=set(range(90))
    target=4571
    current=sum(norms[i] for i in selected)
    swaps=[]
    while current<target:
        choices=[(norms[j]-norms[i],i,j) for i in selected for j in range(len(pool))
                 if j not in selected and 0<norms[j]-norms[i]<=target-current]
        if not choices:
            raise ArithmeticError('deterministic grid norm target was not reached')
        increase,i,j=max(choices)
        selected.remove(i);selected.add(j)
        current+=increase
        swaps.append((i,j))
    points=[(F(27*(2*j-19),380),F(0)) for j in range(20)]
    for i in sorted(selected):
        a,b=pool[i]
        points.extend(((F(a,7),F(b,7)),(F(-a,7),F(-b,7))))
    Q=sum(x*x+y*y for x,y in points)
    if len(points)!=200 or sum(p[0] for p in points) or sum(p[1] for p in points):
        raise ArithmeticError('centered two-hundred-point witness failed')
    if Q!=F(3724007,18620):
        raise ArithmeticError('exact raw moment of deterministic witness failed')
    return {'points':tuple(points),'Q':Q,'normalization_squared':F(200)/Q,
            'selected_grid_indices':tuple(sorted(selected)),'grid_swaps':swaps,
            'grid_pair_norm_sum':current}
