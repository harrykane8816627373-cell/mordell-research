"""Export the positive-angle inner-inner subtree with floor-rounded Bernstein controls."""

from fractions import Fraction as Q
from itertools import product
from math import comb, lcm
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "twelve-point-sharp-global"
DEGREES = (9, 9, 4)
KEYS = tuple(product(*(range(d + 1) for d in DEGREES)))
GRID = 1 << 32
BIAS = 1 << 48
WIDTH = 64
LOG_TERMS = 24
MIN_DISTANCE = Q(1, 3)
SOURCE_GRID = 1 << 44
CAP = (Q(6, 25), Q(6, 25), Q(3, 40))


def pack(values):
    packed = 0
    for index, key in enumerate(KEYS):
        lane = values[key] + BIAS
        if not 0 <= lane < 1 << WIDTH:
            raise ArithmeticError("coefficient outside the export grid")
        packed |= lane << (WIDTH * index)
    return packed


def rounded_source(coefficients, source_scale):
    return {power: Q((value.lo * SOURCE_GRID) // source_scale, SOURCE_GRID)
            for power, value in coefficients.items()}


def lower_power_terms(weights):
    trimmed = {key: weights.get(key, Q(0)) for key in KEYS}
    tail = Q(0)
    for (i, j, k), a in weights.items():
        if a < 0 and (i, j, k) not in trimmed:
            tail += a * CAP[0] ** i * CAP[1] ** j * CAP[2] ** k
    trimmed[(0, 0, 0)] += tail
    return trimmed, tail


def transform_stages(weights, box, transform):
    denominator = lcm(SOURCE_GRID, *(v.denominator for v in weights.values()))
    values = {key: v.numerator * (denominator // v.denominator) for key, v in weights.items()}
    stages = [(values, denominator)]
    axis_matrices = []
    for axis, (d, (lo, hi)) in enumerate(zip(DEGREES, box)):
        matrix, scale = transform(d, lo, hi)
        axis_matrices.append((matrix, scale))
        values = {key: sum(matrix[key[axis]][j] * values[key[:axis] + (j,) + key[axis+1:]]
                           for j in range(d + 1)) for key in KEYS}
        denominator *= scale
        stages.append((values, denominator))
    return stages, axis_matrices


def nested_vec(rows):
    return "![" + ", ".join(map(str, rows)) + "]"


def print_tensor(name, values):
    d0, d1, d2 = (d + 1 for d in DEGREES)
    slabs = [nested_vec(nested_vec(values[a, b, c] for c in range(d2)) for b in range(d1)) for a in range(d0)]
    print(f"def {name} : GridIndex degree → ℤ := tensorOf3 ![\n  " + ",\n  ".join(slabs) + "]\n")


def emit_initial(weights, tail, stages, matrices):
    print("import Mordell.BernsteinInit\nimport Mordell.InnerCoverData\n")
    print("/-! Generated integer witnesses for the three root-coordinate transforms. -/\n")
    print("namespace Mordell.InnerCoverData\n")
    print("set_option maxRecDepth 25000\nset_option maxHeartbeats 0\n")
    print(f"def sourceDen : ℕ := {SOURCE_GRID}")
    inside = {key: int(weights.get(key, Q(0)) * SOURCE_GRID) for key in KEYS}
    print_tensor("sourceInside", inside)
    outside = [(key, a) for key, a in sorted(weights.items()) if key not in inside]
    print("def sourceOutside : PowerTerms 3 ℚ := [\n  " + ",\n  ".join(
        "(" + nested_vec(key) + ", " + lean_rat(a) + ")" for key, a in outside) + "]\n")
    den0 = stages[0][1]
    print(f"def inputMultiplier : ℤ := {den0 // SOURCE_GRID}")
    print(f"def tailNumerator : ℤ := {int(tail * den0)}\n")
    for j, (matrix, scale) in enumerate(matrices):
        print(f"def axisDen{j} : ℕ := {scale}")
        print(f"def axisMatrix{j} : Fin (degree {j} + 1) → Fin (degree {j} + 1) → ℤ :=\n  "
              + nested_vec(nested_vec(row) for row in matrix) + "\n")
    for j, (data, den) in enumerate(stages):
        print(f"def stageDen{j} : ℕ := {den}")
        print_tensor(f"stage{j}", data)
    print("end Mordell.InnerCoverData")


def split(values, axis, right):
    d = DEGREES[axis]
    output = {}
    for key in KEYS:
        i = key[axis]
        r = d - i if right else i
        output[key] = sum(comb(r, d - j if right else j) * values[key[:axis] + (j,) + key[axis+1:]]
                          for j in range(d + 1)) // (1 << r)
    return output


def exact_log_upper(a):
    y = 1 / a if a <= 1 else a
    t = (y - 1) / (y + 1)
    value = 2 * sum(t ** (2 * k + 1) / (2 * k + 1) for k in range(LOG_TERMS))
    return -value if a <= 1 else value + 2 * t ** (2 * LOG_TERMS + 1) / (1 - t*t)


def tangent(box, values):
    (ul, uh), (vl, vh), (tl, th) = box
    anchor = (max(MIN_DISTANCE, ul + vl - 2 * th) + uh + vh - 2 * tl) / 2
    if anchor <= 0:
        raise ArithmeticError("nonpositive tangent anchor")
    upper = Q((exact_log_upper(anchor) * GRID).__ceil__(), GRID)
    for key in KEYS:
        u, v, t = [lo + (hi-lo) * i / d for (lo, hi), i, d in zip(box, key, DEGREES)]
        if Q(values[key], GRID) < upper - 1 + (u + v - 2*t) / anchor:
            raise ArithmeticError("grid controls fail the archived tangent")
    return anchor, upper


def lean_rat(q):
    return f"({q.numerator} / {q.denominator})"


def render(actions, box, root):
    nodes = []

    def walk(window, values, name):
        number = len(nodes)
        action = next(actions)
        node = {"packed": pack(values), "box": name}
        nodes.append(node)
        if action == ord("P"):
            node["tree"] = "GridTree.pruned"
            return number
        if action == ord("A"):
            a, u = tangent(window, values)
            node["tree"] = f"GridTree.tangent ⟨{lean_rat(a)}, {lean_rat(u)}⟩"
            return number
        if action not in b"012":
            raise ArithmeticError("unexpected node in the signed inner tree")
        axis = action - ord("0")
        lo, hi = window[axis]
        mid = (lo + hi) / 2
        w1, w2 = list(window), list(window)
        w1[axis], w2[axis] = (lo, mid), (mid, hi)
        n1 = walk(w1, split(values, axis, False), f"box{number}.left {axis}")
        n2 = walk(w2, split(values, axis, True), f"box{number}.right {axis}")
        node["tree"] = f"GridTree.split {axis} grid{n1} tree{n1} grid{n2} tree{n2}"
        node["children"] = (axis, n1, n2)
        return number

    walk(box, root, "root")
    if next(actions, None) is not None:
        raise ArithmeticError("prefix tree has unconsumed actions")
    return nodes


def emit_checks(nodes):
    print("import Mordell.InnerCoverData\n")
    print("/-! Generated independent kernel decisions for the signed cover. -/\n")
    print("namespace Mordell.InnerCoverData\n")
    print("set_option maxRecDepth 20000\nset_option maxHeartbeats 0\n")
    for i in reversed(range(len(nodes))):
        node = nodes[i]
        print(f"private theorem checked{i} : tree{i}.accept 24 degree (1 / 3) true box{i} grid{i} = true := by")
        if "children" in node:
            axis, left, right = node["children"]
            print(f"  have hl : gridSplitOK degree {axis} false grid{i} grid{left} := by decide +kernel")
            print(f"  have hr : gridSplitOK degree {axis} true grid{i} grid{right} := by decide +kernel")
            print(f"  exact GridTree.accept_split hl checked{left} hr checked{right}\n")
        else:
            print("  decide +kernel\n")
    print("theorem archived_tree_accepted : tree.accept 24 degree (1 / 3) true root grid0 = true := checked0\n")
    print("end Mordell.InnerCoverData")


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.path.insert(0, str(ARCHIVE))
    import kernel_bridge
    import verify_all as archive

    coefficients, _, _ = kernel_bridge.reconstruct()
    box = archive.ROOTS[1]
    weights = rounded_source(coefficients, kernel_bridge.SCALE)
    trimmed, tail = lower_power_terms(weights)
    stages, matrices = transform_stages(trimmed, box, archive.cover.transform)
    if sys.argv[1:] == ["--initial"]:
        emit_initial(weights, tail, stages, matrices)
        return
    controls, denominator = stages[-1]
    coeffs = {key: value * GRID // denominator for key, value in controls.items()}
    tape = (ARCHIVE / "1.txt").read_bytes().strip()
    nodes = render(iter(tape), box, coeffs)
    if sys.argv[1:] == ["--checks"]:
        emit_checks(nodes)
        return
    print("import Mordell.GridBernsteinReplay\n")
    print("/-! Generated by `tools/generate_inner_cover.py` from signed cover `1.txt`. -/\n")
    print("namespace Mordell.InnerCoverData\n")
    print("def root : RationalBox 3 := ⟨![" + ", ".join(lean_rat(b[0]) for b in box) + "], !["
          + ", ".join(lean_rat(b[1]) for b in box) + "]⟩")
    print("def degree : Fin 3 → ℕ := ![9, 9, 4]")
    print(f'def traceCode : String := "{tape.decode("ascii")}"\n')
    for i, node in enumerate(nodes):
        print(f"def grid{i} : ℕ := {node['packed']:#x}")
        print(f"def box{i} : RationalBox 3 := {node['box']}")
    print()
    for i in reversed(range(len(nodes))):
        print(f"def tree{i} : GridTree := {nodes[i]['tree']}")
    print("\ndef tree : GridTree := tree0\n")
    print("end Mordell.InnerCoverData")


if __name__ == "__main__":
    main()
