"""Export the corrected Gram coefficient boxes and rational congruences for Lean."""

from fractions import Fraction
import importlib.util
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ENDPOINT_BITS = 96
BASIS_BITS = 64


def load(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read(relative_path):
    return json.loads((ROOT / relative_path).read_bytes())


def ten_blocks():
    source = load("ten_contact_export", "ten-point-sharp-contacts/verify_all.py")
    nodes, xy, _, _, contacts = source.geometry.candidate()
    source.classify(contacts, source.angle_classes())
    seed = source.read("ten-point-sharp-contacts/seed.json")
    transforms = source.projections(nodes, xy, seed["harmonic_pivots"])
    equations, target, _ = source.build(contacts, nodes, transforms)
    weights, _ = source.correct(equations, target, seed)
    offset = 0
    for k, n in enumerate([5, 6] * 6 + [6, 6]):
        if k == 12:
            offset = len(weights) - 42
        block = [[source.ZERO] * n for _ in range(n)]
        for i, j in combinations_with_replacement(range(n), 2):
            block[i][j] = block[j][i] = weights[offset]
            offset += 1
        yield block, source.geometry.SCALE
    if offset != len(weights):
        raise ArithmeticError("ten-point Gram block layout changed")


def twelve_blocks():
    source = load("twelve_contact_export", "twelve-point-sharp-contacts/verify_all.py")
    names, row_names, matrix, contacts, _, _, _ = source.build()
    minor = read("twelve-point-sharp-contacts/minor.json")
    seed = read("twelve-point-banded-upper/kernel.json")
    weights = source.seed_parameters(seed, names)
    target, _ = source.right_hand_side(contacts)
    residual = [b - source.dot(row, weights) for row, b in zip(matrix, target)]
    rows, columns = minor["selected_rows"], minor["selected_columns"]
    if [row_names[i] for i in rows] != minor["selected_row_names"]:
        raise ArithmeticError("twelve-point contact row layout changed")
    if [names[i] for i in columns] != minor["selected_column_names"]:
        raise ArithmeticError("twelve-point contact column layout changed")
    _, _, correction = source.certify_minor(matrix, rows, columns, residual)
    for index, value in zip(columns, correction):
        weights[index] += value
    values = dict(zip(names, weights))
    for k in range(12):
        n = 4 if k in (3, 6, 9) else 5
        yield [[values[f"Q{k}_{min(i, j)}{max(i, j)}"] for j in range(n)] for i in range(n)], source.SCALE


def interval_data(block, source_scale):
    scale = 1 << ENDPOINT_BITS
    lower = [[value.lo * scale // source_scale for value in row] for row in block]
    upper = [[-((-value.hi * scale) // source_scale) for value in row] for row in block]
    radius_num = max((hi - lo + 1) // 2 for row, high in zip(lower, upper) for lo, hi in zip(row, high))
    center = [[Fraction(lo + hi, 2 * scale) for lo, hi in zip(row, high)] for row, high in zip(lower, upper)]
    return lower, upper, radius_num, center


def approximate_basis(center):
    n = len(center)
    L = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    D = []
    for i in range(n):
        pivot = center[i][i] - sum(L[i][k] ** 2 * D[k] for k in range(i))
        if pivot <= 0:
            raise ArithmeticError("nonpositive exact midpoint pivot")
        D.append(pivot)
        for j in range(i + 1, n):
            L[j][i] = (center[j][i] - sum(L[j][k] * D[k] * L[i][k] for k in range(i))) / pivot
    # Back substitution in L^T P = I leaves the diagonal of P exactly one.
    P = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        for i in reversed(range(j)):
            P[i][j] = -sum(L[k][i] * P[k][j] for k in range(i + 1, j + 1))
    return [[round((1 << BASIS_BITS) * value) if i < j else 0 for j, value in enumerate(row)]
            for i, row in enumerate(P)]


def lean_matrix(matrix, indent):
    return "![\n" + ",\n".join(indent + "![" + ", ".join(map(str, row)) + "]" for row in matrix) + "]"


def emit(label, blocks, count):
    print(f"def {label}GramData : Fin {count} → Σ n, GramCertificate n := ![")
    entries = []
    for block, source_scale in blocks:
        lower, upper, radius_num, center = interval_data(block, source_scale)
        basis = approximate_basis(center)
        entries.append(f"  ⟨{len(block)}, gramFromDyadic {ENDPOINT_BITS} {BASIS_BITS} {radius_num}\n"
                       + "\n".join("    " + lean_matrix(matrix, "      ") for matrix in (lower, upper, basis))
                       + "⟩")
    if len(entries) != count:
        raise ArithmeticError("unexpected Gram block count")
    print(",\n".join(entries) + "]\n")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    print("import Mordell.GramCertificate\n")
    print("/-! Generated by `tools/generate_gram_data.py`. Entries are outward-rounded dyadics. -/\n")
    print("namespace Mordell\n")
    emit("ten", ten_blocks(), 14)
    emit("twelve", twelve_blocks(), 12)
    print("end Mordell")
