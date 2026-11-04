# 9 points: the centered regular octagon

[中文](README.zh-CN.md) · [Home](../../README.md)

## Theorem

For $`\sum|z_i|^2=9`$,

```math
M_9=C_9=9^{36}/8^{28}.
```

Equality holds exactly for

```math
\{0\}\cup\left\{\sqrt{\frac{9}{8}}e^{i\theta+2\pi i k/8}:k=0,\ldots,7\right\},
```

up to rotation and relabeling. There is no additional translation freedom.

## Proof

Every maximizer exists, is centered and collision-free, and satisfies the full first- and second-order necessary conditions. After $`x_i=\sqrt{(n-1)/n}z_i`$, these conditions place every relevant pair in an explicitly proved bounded domain; arbitrary nonsymmetric configurations and degenerate second variations remain included.

An exact positive-definite kernel gives an auxiliary upper bound for the logarithm of each squared distance on that whole domain. Contact neighborhoods are treated by rigorous local positivity, and the complement by outward-rounded Bernstein bounds. Summation yields

```math
\log\Delta(x)\le(n-1)\log(n-1)-\tfrac12\sum_{i,j}K_{ij}\le(n-1)\log(n-1).
```

Rescaling proves the claimed bound, which the stated configuration attains. Equality in the positive-definite terms forces the squared radii to be one zero and $`n-1`$ ones; Newton identities force the nonzero points to form a regular $`(n-1)`$-gon.

The full computational premises have a saved complete certificate with zero unresolved boxes. The normal audit checks the distributed bytes and this record. `check` recomputes the exact algebra, positivity and local contacts; `replay` additionally reruns the full cover. These are ordinary computer-assisted proofs, not proof-assistant formalizations.

[Complete proof in Chinese](../../verification/n9/PROOF_ZH.txt) · [Final certificate](../../verification/n9/FINAL_CERTIFICATE_AUDIT.json)

```sh
python verify.py --case n9
python verify.py --level check --case n9
```
