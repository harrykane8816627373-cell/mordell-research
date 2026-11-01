"""Replay the inherited integer identity and a global quantitative moment gap.

Only the Python standard library is used. Supply the pinned 2026-10-16 ZIP,
or extract the dedicated new offline package, which contains inherited-first3/.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import ast
import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.set_int_max_str_digits(0)


def require(ok, label):
    if not ok:
        raise ArithmeticError(label)


def ceiling(x):
    return (x.numerator + x.denominator - 1) // x.denominator


def parse_polynomial(expression, variables):
    """Read the reconstructor's integer monomials without eval or a CAS."""
    result = {}
    for text in expression.split(' + '):
        powers = [0] * len(variables)

        def visit(node):
            if isinstance(node, ast.Constant) and type(node.value) is int:
                return node.value
            if isinstance(node, ast.Name) and node.id in variables:
                powers[variables.index(node.id)] += 1
                return 1
            if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
                return -visit(node.operand)
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
                return visit(node.left) * visit(node.right)
            if (isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow)
                    and isinstance(node.left, ast.Name) and node.left.id in variables
                    and isinstance(node.right, ast.Constant)
                    and type(node.right.value) is int and node.right.value >= 0):
                powers[variables.index(node.left.id)] += node.right.value
                return 1
            raise ValueError('unsupported polynomial syntax')

        value = visit(ast.parse(text, mode='eval').body)
        monomial = tuple(powers)
        require(monomial not in result, 'unique reconstructed monomial')
        if value:
            result[monomial] = value
    return result


def weighted(monomial):
    return sum(a*b for a, b in zip(monomial, [4, 5, 6, 7, 8]))


def read_sources(config, archive=None, directory=None):
    result = {}
    if archive is not None:
        with zipfile.ZipFile(archive) as source:
            for row in config['inherited_sources']:
                require(source.namelist().count(row['path']) == 1, 'unique inherited ZIP member')
                result[row['path']] = source.read(row['path'])
    else:
        directory = directory or ROOT/'inherited-first3'
        if not directory.is_dir():
            raise FileNotFoundError(
                'Supply --source-archive with the original first3 proof ZIP, or use '
                'the dedicated offline package containing inherited-first3/. '
                'See algebraic-gap/README.md. No large identity was silently skipped.')
        for row in config['inherited_sources']:
            result[row['path']] = (directory/row['path']).read_bytes()
    for row in config['inherited_sources']:
        data = result[row['path']]
        require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                'pinned inherited source '+row['path'])
    return result


def constants(config, norm_bounds):
    caps = list(map(F, config['coefficient_caps']))
    require(caps == list(map(F, ['1/2', '61/100', '49/4', '101/10',
                                '343/54', '3', '2401/4096'])), 'fixed coefficient box')
    # q=w/2: sum |q|^2=7; c2=-s2/2, c3=-s3/3.
    require(F(175, 486) < caps[1]**2, 'old ellipse bounds c3')
    cauchy_checks = []
    for j in range(4, 8):
        radius_squared = F(7*(8-j), 8*j)
        bound_squared = (radius_squared+F(7, 8))**8/radius_squared**(8-j)
        require(bound_squared <= caps[j-2]**2, 'Cauchy coefficient bound c'+str(j))
        cauchy_checks.append({'coefficient': j, 'radius_squared': str(radius_squared),
                              'bound_squared': str(bound_squared)})
    require(caps[6] == F(7, 8)**4, 'constant term AM-GM bound')
    # Descending coefficients of P, P', P''. Centering removes X^7.
    coefficients = [F(1), F(0), *caps]
    a2 = sum(v*v for v in coefficients)
    b2 = sum(((8-j)*v)**2 for j, v in enumerate(coefficients) if j < 8)
    c2 = sum(((8-j)*(7-j)*v)**2 for j, v in enumerate(coefficients) if j < 7)
    A, B, C = config['sylvester_row_caps']
    require([A, B, C] == [18, 60, 171], 'fixed row bounds')
    require(a2 <= A*A and b2 <= B*B and c2 <= C*C, 'three Euclidean row bounds')
    D = B+C
    Ldelta = 7*A**6*B**8 + 48*A**7*B**7
    LR0 = 7*A**6*D**8 + 288*A**7*D**7
    LRt = 64*A**7*D**7*B + 2016*A**8*D**6*B + 48*A**8*D**7
    contents = config['generator_contents']
    require(contents == [4, 24, 16, 32, 64, 128, 4], 'primitive generator contents')
    factors = [0, 30, 20, 12, 6, 2]
    lipschitz = [F(LRt+factor*LR0, content)
                 for factor, content in zip(factors, contents[:6])]
    lipschitz.append(F(LR0, 4))
    K = sum(q*l for q, l in zip(norm_bounds, lipschitz))
    threshold = 7**7
    near = 4*K/threshold**2 + F(64*norm_bounds[6], threshold)
    far = 2*caps[1]*Ldelta/threshold
    exponent = config['force_one_gap_exponent']
    require(type(exponent) is int and exponent == 36, 'fixed output exponent')
    require(3*max(near, far) <= 10**exponent, 'force-one p3 <= 10^36 p2')
    # Reject the tempting stronger decimal statement using this same majorant.
    require(3*max(near, far) > 10**(exponent-1), '10^-35 is not certified by this majorant')
    require(config['full_eight_point_exclusion_proved'] is False, 'global scope remains open')
    return {
        'cauchy_coefficient_checks': cauchy_checks,
        'sylvester_row_norm_squares': list(map(str, [a2, b2, c2])),
        'delta_lipschitz': str(Ldelta), 'R0_coefficient_lipschitz': str(LR0),
        'Rt_coefficient_lipschitz': str(LRt), 'primitive_generator_lipschitz': list(map(str, lipschitz)),
        'generator_difference_budget_K': str(K),
        'near_case_ratio_bound': str(near), 'far_case_ratio_bound': str(far),
        'force_one_p2_over_p3_lower_bound': '1/10^36',
        'stronger_1_over_10_to_35_not_certified_by_this_majorant': True,
    }


def verify(config, sources):
    env = os.environ.copy()
    env['PYTHONUTF8'] = '1'
    with tempfile.TemporaryDirectory(prefix='mordell-gap-', dir=ROOT) as temp:
        stage = Path(temp)
        for name, data in sources.items():
            target = stage/name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        for entry in ['generate_first3.py', 'verify_full_discriminant_squared.py']:
            run = subprocess.run([sys.executable, '-B', '-S', '-O', str(stage/entry)],
                                 cwd=stage, env=env, capture_output=True, text=True,
                                 encoding='utf-8', timeout=180)
            require(run.returncode == 0, entry+': '+run.stdout[-2000:]+run.stderr[-2000:])
        old_report = json.loads((stage/'full_discriminant_squared_verification.json').read_text())
        require(old_report['identity_proved'] is True and old_report['coefficient_checks'] == 14179,
                'full inherited integer identity replayed')
        require(old_report['target_content'] == '1', 'primitive target is exactly delta squared')
        fresh = json.loads((stage/'first3_fresh_generators.json').read_text())
        cert = json.loads(sources['full_fast1/full_discriminant_squared_certificate.json'])
        require(fresh['variables'] == cert['variables'] == ['a', 'b', 'c', 'd', 'u'], 'variable order')
        gweights = [55, 54, 53, 52, 51, 50, 54]
        qweights = [60, 61, 62, 63, 64, 65, 61]
        contents = config['generator_contents']
        require(qweights == config['multiplier_homogeneous_weights'], 'multiplier target weights')
        for index, expression in enumerate(fresh['generators']):
            poly = parse_polynomial(expression, fresh['variables'])
            require(math.gcd(*(abs(v) for v in poly.values())) == contents[index],
                    'fresh generator content '+str(index+1))
            for monomial in poly:
                gap = gweights[index]-weighted(monomial)
                require(gap >= 0 and gap % 3 == 0, 'generator homogenization power')
            require(qweights[index]+gweights[index] == 115, 'homogeneous product degree 115')
        delta = parse_polynomial(fresh['delta'], fresh['variables'])
        for monomial in delta:
            gap = 56-weighted(monomial)
            require(gap >= 0 and gap % 3 == 0, 'discriminant homogeneous degree 56')
        require(3+2*56 == 115, 'left side is c3 times delta squared')
        den = abs(int(cert['denominator']))
        require(den > 0, 'nonzero inherited denominator')
        caps = list(map(F, config['coefficient_caps']))
        norms, term_counts, powers = [], [], []
        for i, rows in enumerate(cert['multipliers']):
            bound = F(0)
            seen = set()
            max_power = 0
            for monomial, value in rows:
                require(len(monomial) == 5 and all(type(v) is int and v >= 0 for v in monomial),
                        'valid multiplier monomial')
                require(tuple(monomial) not in seen, 'unique multiplier monomial')
                seen.add(tuple(monomial))
                gap = qweights[i]-weighted(monomial)
                require(gap >= 0 and gap % 3 == 0, 'multiplier homogenization power')
                power = gap//3
                max_power = max(max_power, power)
                product = caps[1]**power
                for m, radius in zip(monomial, caps[2:]):
                    product *= radius**m
                bound += abs(int(value))*product
            norms.append(ceiling(bound/den))
            term_counts.append(len(rows))
            powers.append(max_power)
        require(norms == config['normalized_multiplier_norm_ceilings'], 'seven exact norm ceilings')
        require(sum(term_counts) == 7123, 'all inherited multiplier terms used')
        proof = constants(config, norms)
        return {
            'passed': True, 'arithmetic': 'Python standard library; exact integers and Fractions',
            'inherited_resultants_freshly_reconstructed': True,
            'inherited_identity_replayed': True, 'inherited_integer_coefficient_checks': 14179,
            'homogeneous_identity': 'c3*delta0^2 = sum((Qj_hom/D)*Gj_at_c2_zero)',
            'homogeneous_identity_weight': 115,
            'multiplier_term_counts': term_counts, 'max_c3_powers': powers,
            'normalized_multiplier_norm_ceilings': norms,
            'rational_constant_checks': proof,
            'scope': 'True stationary configurations in the inherited ellipse, with Delta_w >= 4^28*7^7; no orbit-distance assumption.',
            'full_eight_point_exclusion_proved': False,
            'all_nonzero_moment_stationary_points_excluded': False,
            'complete_stationary_classification_replayed': False,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--source-archive', type=Path)
    group.add_argument('--source-directory', type=Path)
    parser.add_argument('--write-report', type=Path)
    args = parser.parse_args()
    config = json.loads((HERE/'certificate.json').read_text(encoding='utf-8'))
    manifest_path = HERE/'source_manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        for row in manifest['files']:
            data = (HERE/row['path']).read_bytes()
            require(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
                    'fixed current input '+row['path'])
    result = verify(config, read_sources(config, args.source_archive, args.source_directory))
    if args.write_report is not None:
        args.write_report.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    else:
        expected = json.loads((HERE/'verification.json').read_text(encoding='utf-8'))
        require(result == expected, 'deterministic complete report')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
