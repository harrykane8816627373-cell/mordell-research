"""Exact checks for full-plane W1 constants and the nonradial class example."""
from fractions import Fraction as F
from pathlib import Path
import json

from distribution_metric import (sqrt_bounds,pi_interval,smear_penalty_upper,
    quantization_lower,transport_bound,effective_global_gap,global_transport_bound,
    asymptotic_transport_certificate,cross_axes_background_interval,
    lipschitz_class_check,box_population,outward_decimal,log_interval)


def require(ok,label):
    if not ok:
        raise ArithmeticError(label)


def run(write_report=False):
    root_checks=0
    for a in range(61):
        for b in range(1,10):
            q=F(a,b)
            lo,hi=sqrt_bounds(q,10000)
            require(lo*lo<=q<=hi*hi and hi-lo<=F(1,10000),'outward square-root interval')
            root_checks+=1
    plo,phi=pi_interval()
    require(3<plo<phi<F(22,7),'independent simple pi upper')

    # The harmonic Kelvin extension preserves the inside harmonic energy.
    # A monomial r^k cos(k theta), normalized to Lip 1 on B_R, has total
    # Dirichlet norm² R²/k; k=1 saturates the uniform R² test bound.
    kelvin_checks=0
    for k in range(1,9):
        for R in (F(1,3),F(1),F(7,4)):
            inside=F(k,2)*R**(2*k)
            outside=F(k,2)*R**(2*k)
            normalized=(inside+outside)/(k*k*R**(2*k-2))
            require(normalized==R*R/k<=R*R,'harmonic inside/outside normalization')
            kelvin_checks+=1

    quantization_checks=0
    for n in range(2,201):
        # The layer-cake rearrangement integral is R/n-R³/6 with R²=2/n.
        coefficient=F(1,n)-F(2,n)/6
        require(coefficient==F(2,3*n),'equal-atom quantization layer-cake identity')
        low=quantization_lower(n)
        require(low>0 and low*low<=F(8,9*n),'quantization bound points downward')
        quantization_checks+=1

    # Independent finite-energy radial measure: unit circle vs uniform sqrt2
    # disk has D=-1/4+log2/2 and exact W1=2(sqrt2-1)/3.
    l2,h2=log_interval(F(2))
    r2lo,r2hi=sqrt_bounds(F(2))
    Dlo,Dhi=-F(1,4)+l2/2,-F(1,4)+h2/2
    Wlo,Whi=F(2,3)*(r2lo-1),F(2,3)*(r2hi-1)
    require(Dlo>0 and Whi*Whi<2*Dlo,'compact-disk Coulomb transport benchmark')

    require(effective_global_gap(100)==F(1657,12),'existing effective global gap at n100')
    finite_examples=[]
    for n in (6,10,30,100,1000,10000,1000000):
        data=global_transport_bound(n)
        require(data['W1_upper']>=data['universal_quantization_W1_lower'],
                'upper certificate compatible with universal quantization lower')
        for row in data['candidates']:
            require(row['full_support_radius_upper']>=r2hi,
                    'full transport support includes background disk')
        finite_examples.append({'n':n,'h':str(data['chosen_h']),
            'W1_upper_outward':str(outward_decimal((data['W1_upper'],data['W1_upper']))[1]),
            'W1_lower_outward':str(outward_decimal((data['universal_quantization_W1_lower'],
                                                   data['universal_quantization_W1_lower']))[0])})
    small_h=transport_bound(4,1,1,h=F(1,100),original_z_scale=False)
    require(small_h['full_support_radius_upper']==r2hi,
            'small smear must still include the entire background disk')
    default=transport_bound(100,F(1657,12),F(3,2))
    require(default['smear_penalty_upper']==0,'h1 has exactly zero penalty')
    larger_gap=transport_bound(100,F(1657,12)+1,F(3,2))
    require(larger_gap['W1_upper']>=default['W1_upper'],'valid gap enlargement weakens bound')

    asym=asymptotic_transport_certificate()
    require(asym['headline_certified'] and asym['limsup_sqrt_n_W1_upper']<F(981,1000),
            'strict asymptotic .981 upper coefficient')
    require(asym['universal_liminf_sqrt_n_W1_lower']<asym['limsup_sqrt_n_W1_upper'],
            'asymptotic lower/upper constants consistent')

    # A fully two-dimensional statistic: distance to the union of the axes.
    # sqrt2>7/5 and pi<22/7 already give the simple strict mean lower 56/165.
    bglo,bghi=cross_axes_background_interval()
    require(bglo>F(56,165),'nonradial statistic simple rational background lower')
    n100=global_transport_bound(100)
    require(n100['W1_upper']<F(281,1000)<F(191,660),
            'n100 full-plane budget excludes a cross-strip class')
    exclusion=lipschitz_class_check((0,F(1,20)),(bglo,bghi),1,n100['W1_upper'])
    require(exclusion['excluded'],'entire axis-cross strip class excluded')
    # Reuse a previously certified attainable witness; do not relabel its
    # verification as a newly replayed 91-configuration search.
    witness_path=Path(__file__).resolve().parents[1]/'symmetry-breaking'/'witness_verification.json'
    stored=json.loads(witness_path.read_text(encoding='utf-8'))
    require(stored['passed'] is True,'prior attainable-witness verification reports PASS')
    stored_row=next(row for row in stored['rows'] if row['n']==100)
    stored_gap=F(stored_row['global_optimizer_gap_upper'])
    require(stored_gap==F(413923829261,500000000000),'exact stored n100 gap input')
    witness_transport=global_transport_bound(100,stored_gap)
    require(witness_transport['W1_upper']<F(14,125)<F(23,165),
            'actual attainable witness strengthens full-plane distribution')
    wide=lipschitz_class_check((0,F(1,5)),(bglo,bghi),1,witness_transport['W1_upper'])
    wide_old=lipschitz_class_check((0,F(1,5)),(bglo,bghi),1,n100['W1_upper'])
    require(wide['excluded'] and not wide_old['excluded'],
            'stronger proved lower bound excludes a wider nonradial class through W1')

    # Nonempty centered example with very accurate radial distribution:
    # each q_j contributes {+sqrt q,-sqrt q,+i sqrt q,-i sqrt q}.
    qs=[F(2*j-1,25) for j in range(1,26)]
    require(4*sum(qs)==100,'fourfold shell example total moment')
    require(max(qs)==F(49,25)<F(99,50),'example meets the existing radius cap')
    cdf_error=F(0)
    for j,q in enumerate(qs,1):
        cdf_error=max(cdf_error,abs(F(j-1,25)-q/2),abs(F(j,25)-q/2))
    require(cdf_error==F(1,50),'actual radial CDF error is only one fiftieth')

    # Rectangle translation from a certified W1 bound, including degenerate
    # erosions and boxes extending outside the background circle.
    box_checks=0
    for rect in ((-1,1,-1,1),(0,1,0,1),(-3,3,-3,3)):
        for eta in (F(1,10),F(1,4),F(3,2)):
            row=box_population(1000,rect,F(1,100),eta)
            require(0<=row['population_lower']<=row['population_upper']<=1000,
                    'rectangle population interval lies in allowed range')
            box_checks+=1
    # For zero transport, the empirical/reference-measure formula should
    # still enclose the exact mass of any rectangle wholly inside B_sqrt2.
    box=box_population(1000,(F(-1,2),F(1,2),F(-1,2),F(1,2)),0,F(1,100))
    mass_lo,mass_hi=1/(2*phi),1/(2*plo)
    require(F(box['population_lower'],1000)<=mass_hi
            and F(box['population_upper'],1000)>=mass_lo,
            'zero-transport rectangle mass normalization')

    invalid_calls=[
        lambda:global_transport_bound(1),lambda:effective_global_gap(5),
        lambda:transport_bound(100,-1,2),lambda:transport_bound(100,1,F(1,10)),
        lambda:transport_bound(100,1,2,h=0),
        lambda:lipschitz_class_check((1,0),(0,1),1,1),
        lambda:box_population(100,(0,0,0,1),1,F(1,10)),
        lambda:box_population(100,(0,1,0,1),1,0),
    ]
    for callback in invalid_calls:
        try:
            callback()
        except ValueError:
            pass
        else:
            raise ArithmeticError('invalid distribution-metric input accepted')

    report={
        'status':'PASS','arithmetic':'Fraction, outward square roots, logarithm and Machin pi intervals',
        'root_interval_checks':root_checks,'harmonic_kelvin_checks':kelvin_checks,
        'quantization_layer_cake_checks':quantization_checks,'rectangle_population_checks':box_checks,
        'invalid_inputs_rejected':len(invalid_calls),'finite_global_examples':finite_examples,
        'asymptotic_h':str(asym['h']),
        'asymptotic_upper_outward':str(outward_decimal((asym['limsup_sqrt_n_W1_upper'],
                                                      asym['limsup_sqrt_n_W1_upper']))[1]),
        'asymptotic_lower_outward':str(outward_decimal((asym['universal_liminf_sqrt_n_W1_lower'],
                                                      asym['universal_liminf_sqrt_n_W1_lower']))[0]),
        'asymptotic_headline':'limsup sqrt(n) W1 < 981/1000; external linear-energy and triangular trial input retained',
        'cross_strip_width':'1/20','cross_background_mean_simple_strict_lower':'56/165',
        'cross_strip_W1_simple_strict_lower':'191/660','cross_strip_excluded':True,
        'prior_witness_gap_source':'symmetry-breaking/witness_verification.json, rows[n=100].global_optimizer_gap_upper',
        'prior_witness_gap_upper':str(stored_gap),
        'prior_witness_W1_upper_outward':str(outward_decimal((witness_transport['W1_upper'],witness_transport['W1_upper']))[1]),
        'wider_cross_strip_width':'1/5','wider_cross_strip_W1_simple_strict_lower':'23/165',
        'wider_cross_strip_excluded_with_stored_witness':True,
        'wider_cross_strip_excluded_with_old_effective_gap_by_this_statistic':False,
        'fourfold_shell_example_n':100,'fourfold_shell_squared_radii':[str(q) for q in qs],
        'fourfold_shell_radial_CDF_error':'1/50',
        'scope':'Implementation and explicit nonradial class checks. Universal transport and quantization arguments are proved in DISTRIBUTION_METRIC.md; no optimal lattice, unknown sharp value or full-domain classification is established.',
    }
    if write_report:
        Path(__file__).with_name('distribution_metric_verification.json').write_text(
            json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    run(parser.parse_args().write_report)
