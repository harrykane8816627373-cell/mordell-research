"""Exact constants and actual-root witness for a radial-to-angular energy bound.

The analytic proof is in README.md. No global eight-point coverage is claimed.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent


def constants(epsilon=F(1, 1200), log_gap=F(1, 32)):
    labels = []

    def need(ok, label):
        if not ok:
            raise ArithmeticError(label)
        labels.append(label)

    a, b, tau = F(31, 32), F(33, 32), F(1, 32)
    separation = F(27, 250)
    need(0 < epsilon and 0 < log_gap, 'positive parameters')
    need(a+tau == 1 and b-tau == 1, 'radial band endpoints')
    need(epsilon < tau*a*a and epsilon < a*tau*tau and epsilon < b*tau*tau,
         'three cubic barriers')
    # f(t)=t(t-1)^2; f'(t)=(t-1)(3t-1). No sample-only identity check.
    radial_polynomial = [F(0), F(1), F(-2), F(1)]
    derivative = {i-1: i*v for i, v in enumerate(radial_polynomial) if i and v}
    factored = laurent_mul({1: F(1), 0: F(-1)}, {1: F(3), 0: F(-1)})
    need(derivative == factored, 'cubic derivative factorization')
    need(6*b+2*tau < 7 and 8*a > 7, 'exactly one small and seven ring points')
    beta = 1/(a*(1+a)**2)
    need(beta < 1/a**2, 'radial matching distance upper bound')
    denominator = separation-2*beta*epsilon
    need(denominator > 0, 'radial projections cannot collide')
    ring = 7*beta/denominator
    central = (7/a-2)/a**2
    cross_squared = 7/a**4
    young_small_slack = 18-central-12
    young_ring_slack = 18-ring-cross_squared/12
    need(young_small_slack > 0 and young_ring_slack > 0,
         'Young inequality gives total energy error at most 18 E')
    # Full coefficient identity for the complete-graph quadratic form.
    lhs = {}
    for i in range(7):
        for j in range(i+1, 7):
            lhs[(i, i)] = lhs.get((i, i), 0)+1
            lhs[(j, j)] = lhs.get((j, j), 0)+1
            lhs[(i, j)] = lhs.get((i, j), 0)-2
    rhs = {(i, j): (6 if i == j else -2) for i in range(7) for j in range(i, 7)}
    need(lhs == rhs, 'complete-graph quadratic identity, all 28 coefficients')
    gamma = b*(1+b)**2
    need(gamma > 1, 'one gamma bounds both central and ring defects')
    need(F(2, 15)/gamma >= log_gap, 'local basin implies the proposed E log gap')
    angular_slack = F(7, 400)/b-(18+F(7, 4)/(b*a*a)+log_gap)*epsilon
    need(angular_slack > 0, 'outside the local basin, angular loss exceeds the proposed E gap')
    decimation_ratio = separation*(56-separation)**27/F(2**15*3**75*7**7)
    need(separation < 2 and decimation_ratio < 1-log_gap*epsilon,
         'short-edge decimation also implies the proposed E gap')
    need(log_gap*epsilon < 1, 'exponential gap implies a polynomial half-gap')
    distance_factor = 1/a**2+b*F(72, 7)
    need(distance_factor*epsilon < F(1, 100), 'challengers enter the q-distance 1/10 basin')
    return {
        'epsilon': str(epsilon), 'log_gap_coefficient': str(log_gap),
        'radial_band': [str(tau), str(a), str(b)],
        'radial_error_beta': str(beta), 'pair_projection_denominator': str(denominator),
        'central_error_coefficient': str(central), 'cross_coefficient_squared': str(cross_squared),
        'ring_error_coefficient': str(ring),
        'Young_small_slack': str(young_small_slack), 'Young_ring_slack': str(young_ring_slack),
        'radial_to_orbit_gamma': str(gamma),
        'local_gap_slack': str(F(2, 15)/gamma-log_gap),
        'outside_basin_angular_slack': str(angular_slack),
        'short_edge_decimation_ratio': str(decimation_ratio),
        'short_edge_gap_slack': str(1-log_gap*epsilon-decimation_ratio),
        'challenger_q_distance_squared_upper': str(distance_factor*epsilon),
        'improvement_over_1_over_35000': str(epsilon/F(1, 35000)),
        'checks': labels,
    }


def laurent_mul(p, q):
    out = {}
    for i, v in p.items():
        for j, w in q.items():
            out[i+j] = out.get(i+j, F(0))+v*w
    return {k: v for k, v in out.items() if v}


def laurent_power(p, n):
    out = {0: F(1)}
    for _ in range(n):
        out = laurent_mul(out, p)
    return out


def trace_seventh_roots(p):
    return 7*sum((v for k, v in p.items() if k % 7 == 0), F(0))


def actual_witness():
    """One small point h/sqrt(m), seven points (zeta^j-h/7)/sqrt(m)."""
    labels = []

    def need(ok, label):
        if not ok:
            raise ArithmeticError(label)
        labels.append(label)

    h = F(1, 50)
    m = 1+8*h*h/49
    ring = {1: F(1), 0: -h/7}
    conjugate_ring = {-k: v for k, v in ring.items()}
    radius = laurent_mul(ring, conjugate_ring)
    pure = {k: h**k+trace_seventh_roots(laurent_power(ring, k)) for k in range(1, 9)}
    radial = {k: (h**(2*k)+trace_seventh_roots(laurent_power(radius, k)))/m**k
              for k in range(1, 4)}
    need(pure[1] == 0 and radial[1] == 7, 'actual root construction is centered with mass seven')
    need(pure[2] == 8*h*h/7 and pure[3] == F(48, 49)*h**3,
         'actual witness has both low moments nonzero')
    epsilon = radial[3]-2*radial[2]+7
    need(F(1, 35000) < epsilon < F(1, 1200), 'actual witness is in the added radial band')
    s2 = pure[2]/m
    s3_mod_squared = pure[3]**2/m**3
    need(s2*s2+F(54, 175)*s3_mod_squared < 1, 'actual witness obeys the old ellipse')
    for k in range(2, 7):
        need(pure[k]**2/m**k < F(1, 64), 'actual witness meets old s'+str(k)+' modulus cap')
    s8 = pure[8]/m**4
    need(abs(s8) > F(1, 8), 'actual witness is outside the old six-moment tube at s8')
    # sin(pi/7)>=2/7 gives squared ring-edge lower bound 16/(49m).
    ring_distance_lower = F(16, 49)/m
    center_distance_lower = (1-8*h/7)**2/m
    need(ring_distance_lower > F(27, 250) and center_distance_lower > F(27, 250),
         'all actual witness edges meet the separation cut')
    raw_q1 = h**3+trace_seventh_roots(laurent_mul(radius, ring))
    need(raw_q1 == -2*h+F(48, 49)*h**3 and raw_q1 != 0,
         'nonzero Ward Q1 residual confirms this is not a stationary counterexample')
    ratio = (1-(8*h/7)**7)**2/m**28
    need(0 < ratio < 1-epsilon/32,
         'exact pair product independently satisfies the new logarithmic gap')
    need(ratio < 1-epsilon/64, 'exact pair product also satisfies the polynomial gap')
    return {
        'description': 'q_center=h/sqrt(m), seven ring points=(zeta^j-h/7)/sqrt(m), j=0..6, h=1/50',
        'normalization_m': str(m), 'E': str(epsilon), 's2': str(s2),
        's3_modulus_squared': str(s3_mod_squared), 's8': str(s8),
        'raw_Q1_residual': str(raw_q1), 'Delta_over_target': str(ratio),
        'ring_pair_distance_squared_lower': str(ring_distance_lower),
        'center_pair_distance_squared_lower': str(center_distance_lower),
        'is_stationary': False, 'checks': labels,
    }


def verify():
    arithmetic = constants()
    rejected = []
    for label, eps, gap in [('larger_1_over_1100_not_certified', F(1, 1100), F(1, 32)),
                            ('stronger_E_over_31_not_certified', F(1, 1200), F(1, 31))]:
        try:
            constants(eps, gap)
        except ArithmeticError as exc:
            rejected.append({'case': label, 'reason': str(exc)})
        else:
            raise ArithmeticError('unsafe strengthened parameter unexpectedly accepted')
    return {
        'passed': True, 'arithmetic': 'Python standard library; exact rational arithmetic',
        'theorem': 'E<=1/1200 implies log(Delta_q/7^7)<=-E/32, and Delta_q/7^7<=1-E/64',
        'hypotheses': 'Eight distinct actual centered points of squared mass 7; no stationarity, stability, moment ellipse or a priori radius cap is assumed.',
        'inherited_dependencies': [
            'Exact seven-point sharp bound, including its disclosed external classification input.',
            'The certified radius-1/5 force-one objective basin from eight-point/ENERGY_BASIN.md.',
        ],
        'constants': arithmetic, 'actual_nonstationary_witness': actual_witness(),
        'rejected_strengthenings': rejected,
        'analytic_bridges_fully_formalized': False,
        'full_eight_point_exclusion_proved': False,
        'remaining_non_candidate_challengers_must_have': 'E>1/1200',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.write_report:
        args.write_report.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    else:
        expected = json.loads((HERE/'verification.json').read_text(encoding='utf-8'))
        if result != expected:
            raise ArithmeticError('deterministic exact report mismatch')
    print(json.dumps(result, indent=2))
