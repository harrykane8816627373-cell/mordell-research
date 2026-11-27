"""Write the rational coefficient data for Mordell/RootData.lean to stdout."""

from collections import defaultdict
from fractions import Fraction
from itertools import product
import json
from math import comb, prod
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "ten-point-contact-geometry"
sys.path.insert(0, str(SOURCE))
from check_root_box import plus, rat, scaled, system


def centered_polynomial(poly, center, radius):
    terms = defaultdict(Fraction)
    for exponents, coefficient in poly.items():
        for powers in product(*(range(exponent + 1) for exponent in exponents)):
            terms[powers] += coefficient * prod(
                comb(degree, power) * midpoint ** (degree - power) * radius ** power
                for degree, power, midpoint in zip(exponents, powers, center)
            )
    return {powers: coefficient for powers, coefficient in terms.items() if coefficient}


def inverse(matrix):
    size = len(matrix)
    augmented = [
        row + [Fraction(i == j) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for i in range(size):
        pivot = next(j for j in range(i, size) if augmented[j][i])
        augmented[i], augmented[pivot] = augmented[pivot], augmented[i]
        divisor = augmented[i][i]
        augmented[i] = [entry / divisor for entry in augmented[i]]
        for j in range(size):
            if j != i:
                factor = augmented[j][i]
                augmented[j] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(augmented[j], augmented[i])
                ]
    return [row[size:] for row in augmented]


def literal(value):
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator} / {value.denominator}"


def render_matrix(name, rows):
    body = ";\n    ".join(", ".join(map(literal, row)) for row in rows)
    return f"def {name} : Matrix (Fin 5) (Fin 5) ℝ :=\n  !![{body}]\n"


def render_terms(poly):
    terms = []
    for powers, coefficient in sorted(poly.items(), key=lambda term: (sum(term[0]), term[0])):
        indices = ", ".join(str(i) for i, power in enumerate(powers) for _ in range(power))
        terms.append(f"    ({literal(coefficient)}, [{indices}])")
    return "[\n" + ",\n".join(terms) + "]"


def main():
    data = json.loads((SOURCE / "root_box_jacobian.json").read_bytes())
    evidence = data["output"]["conclusion"]["evidence"]
    preconditioner = [
        [Fraction(round(rat(value) * 256), 256) for value in row]
        for row in evidence["preconditioner"]["entries"]
    ]
    center = [rat(value) for value in evidence["center"]["values"]]
    radius = Fraction(1, 2**20)
    equations = [centered_polynomial(poly, center, radius) for poly in system()]
    correction = [
        plus(
            {tuple(int(i == j) for j in range(5)): Fraction(1)},
            *(scaled(poly, -weight / radius) for poly, weight in zip(equations, row)),
        )
        for i, row in enumerate(preconditioner)
    ]
    header = """import Mordell.PolynomialBounds
import Mathlib.LinearAlgebra.Matrix.Notation

/-! Rational coefficient data for the five-dimensional root box.
The preconditioner is the Jacobian inverse rounded to multiples of 1/256.
The expanded fixed-point map is checked against the original equations in RootBox. -/

noncomputable section

namespace Mordell

"""
    rows = ",\n   ".join(map(render_terms, correction))
    source = (
        header
        + render_matrix("rootPreconditioner", preconditioner)
        + "\n"
        + render_matrix("rootLeftInverse", inverse(preconditioner))
        + "\nset_option maxHeartbeats 3000000 in\nset_option maxRecDepth 12000 in\n"
        + "def rootCoefficients : Fin 5 → Coefficients 5 :=\n"
        + "  ![" + rows + "]\n\nend Mordell\n"
    )
    sys.stdout.buffer.write(source.encode("utf-8"))


if __name__ == "__main__":
    main()
