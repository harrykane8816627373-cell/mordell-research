"""Export the archived n=12 center and preconditioner for Lean."""

import json
from fractions import Fraction
from pathlib import Path
import sys


root = Path(__file__).resolve().parents[2]
data = json.loads((root / "validated-maxima/certificates.json").read_bytes())
row, = (row for row in data["rows"] if row["n"] == 12)
center, matrix = row["center_numerators"], row["inverse_numerators"]
if len(center) != 23 or len(matrix) != 23 or any(len(line) != 23 for line in matrix):
    raise ValueError("expected the archived 23-coordinate certificate")


def vector(values):
    return "![" + ", ".join(map(str, values)) + "]"


sys.stdout.reconfigure(encoding="utf-8", newline="\n")
print("import Mathlib.Data.Fin.VecNotation\n")
print("/-! Generated from `validated-maxima/certificates.json` by `tools/generate_twelve_root_data.py`. -/\n")
print("namespace Mordell")
print("namespace TwelveRootData\n")
print(f"def centerBits : ℕ := {data['center_exponent']}")
print(f"def inverseBits : ℕ := {data['matrix_exponent']}")
print(f"def radiusBits : ℕ := {data['radius_exponent']}\n")
print("def centerNumerators : Fin 23 → ℤ :=\n  " + vector(center) + "\n")
print("def inverseNumerators : Fin 23 → Fin 23 → ℤ := ![\n"
      + ",\n".join("  " + vector(line) for line in matrix) + "]\n")

secant_scale = 2**40
coords = [Fraction(q, 2 ** data["center_exponent"]) for q in [center[0], 0, *center[1:]]]
points = [coords[i:i + 2] for i in range(0, len(coords), 2)]


def pair(i, j):
    if i == j:
        return [[0, 0], [0, 0]]
    d = [points[i][a] - points[j][a] for a in range(2)]
    d2 = sum(x*x for x in d)
    return [[round(secant_scale * (int(a == b) / d2 - 2*d[a]*d[b] / d2**2))
             for b in range(2)] for a in range(2)]


pairs = [[pair(i, j) for j in range(12)] for i in range(12)]
print("def pairSecantNumerators : Fin 12 → Fin 12 → Fin 2 → Fin 2 → ℤ := ![")
print(",\n".join("  ![" + ", ".join("![" + ", ".join(map(vector, entry)) + "]" for entry in row)
                + "]" for row in pairs) + "]\n")
axes = [divmod(k, 2) for k in range(24) if k != 1]
jacobian = [[(sum(pairs[i][k][a][b] for k in range(12)) if i == j else -pairs[i][j][a][b])
             - (11 * secant_scale // 2) * int(p == q)
             for q, (j, b) in enumerate(axes)] for p, (i, a) in enumerate(axes)]
print("def secantNumerators : Fin 23 → Fin 23 → ℤ := ![\n"
      + ",\n".join("  " + vector(line) for line in jacobian) + "]\n")
print("end TwelveRootData")
print("end Mordell")
