"""Generate the full Ward-corrected gradient with Jacobian 0.23.0."""

from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import tempfile


ROOT = Path(__file__).resolve().parent.parent
KERNEL = ROOT / 'twelve-point-ward-upper/kernel.json'
WORK = ROOT / 'work'


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(output):
    from importlib.metadata import version
    from jacobian.math.polynomials.values import RationalPolynomial
    from jacobian.math.polynomials.vector_calculus.operations import gradient

    if version('jacobian') != '0.23.0':
        raise RuntimeError('Expected Jacobian 0.23.0')
    output = output.resolve()
    if not output.is_relative_to(WORK.resolve()):
        raise ValueError('Choose a generated output under the repository work directory')

    ward = load_module('twelve-point-ward-upper/verify_all.py', 'ward_gradient_kernel')
    cover = load_module('twelve-point-kernel-upper/verify_cover.py', 'ward_gradient_cover')
    composer = load_module('jacobian-gradient/analytic.py', 'ward_gradient_composer')
    source = KERNEL.read_bytes()
    data = json.loads(source)
    polynomial, _ = ward.combined_polynomial(data, ward.ward_matrices(data), cover)

    native = []
    derivatives, _ = composer.scaled_native_gradient(
        polynomial, 'complete_ward_kernel', native, gradient, RationalPolynomial)
    for axis, derivative in enumerate(derivatives):
        if derivative != composer.derivative(polynomial, axis):
            raise ArithmeticError(f'Incorrect composed coefficient on axis {axis}')

    report = {
        'passed': True,
        'kernel_sha256': hashlib.sha256(source).hexdigest(),
        'actual_Jacobian_native_gradient_calls': len(native),
        'complete_h_support': len(polynomial),
        'complete_gradient_supports': list(map(len, derivatives)),
        'global_coefficient_identity_checked': True,
        'includes_ward_corrections': True,
    }
    artifact = {
        **report,
        'variables': list(composer.AXES),
        'jacobian_version': version('jacobian'),
        'operation': 'polynomial_field.scalar.gradient.compute',
        'h': composer.to_plain(polynomial),
        'gradient': [composer.to_plain(part) for part in derivatives],
        'native_blocks': native,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        staged = Path(temporary) / output.name
        staged.write_text(json.dumps(artifact, indent=2) + '\n', encoding='utf-8', newline='\n')
        staged.replace(output)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=WORK / 'jacobian-gradient-ward/gradient.json')
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2))
