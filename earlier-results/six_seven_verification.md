# Offline exact verification for n=6 and n=7

This bundle replays the finite rational certificate used in the six- and seven-point global-value proofs. It includes the frozen author coefficient charts, parameter and point-root disks, exact checkers, and the results independently reproduced on 2026-10-06. Python 3.12 was used; verification needs only the standard library and no network or computer algebra system.

For the constraint `sum |z_i|^2 = n`, the resulting values are `M6 = 6^15/5^10` and `M7 = 7^21/6^15`. Equality is the origin plus a regular pentagon or regular hexagon, with radii `sqrt(6/5)` or `sqrt(7/6)`, respectively, up to rotation and permutation.

## What the verification proves

The supplied charts contain 80 six-point parameter roots and 174 seven-point parameter roots. Strict rational Rouché inequalities certify every parameter root and every associated point root, with disjoint disks and the required degree counts. Coefficient uncertainty is propagated from the parameter disks. The physical test uses the original pointwise residual

`Q P''(w_i) - n(n-1) conjugate(w_i) P'(w_i) = 0`, where `Q = sum |w_i|^2`.

Every retained candidate is either excluded by this residual, bounded strictly below the target, or identified as the exact attaining regular polygon with its center. Attainment is also checked directly in the mathematical proof and does not follow merely from surviving a necessary-equation test.

| n | Parameter roots | Point-root disks | Force excluded | Strictly below target | Attaining |
|---|---:|---:|---:|---:|---:|
| 6 | 80 | 480 | 69 | 10 | 1 |
| 7 | 174 | 1218 | 155 | 18 | 1 |

The auxiliary six-point component checker proves that the supplied components satisfy the necessary coefficient equations. It verifies forward containment, not completeness of the original ideal. Its nine former `assert` checks were replaced by explicit `ArithmeticError` checks, and the corrected version passes with `-O`.

## External coverage input and its limits

The proof uses the published exact classification computation of Faugère and Svartz, ISSAC 2012, pp.170–178, [DOI](https://doi.org/10.1145/2442829.2442856), [author manuscript](https://inria.hal.science/hal-00777791/file/FS12.pdf). Their classification concerns all solutions of the original physical equilibrium equations, including collinear equilibria; it is not restricted to stable or previously known equilibria.

The coefficient-chart splitting and collision saturation are in §5.3, the physical filtering in §5.4, both on proceedings p.176. The actual computations are reported in §5.7, p.177. Theorem 5 provides necessary symmetric equations; it is not itself a theorem listing the eleven six-point components.

For n=7, §7 Figure 4, p.178, reports a degree-526 `e2=1` component. Its displayed branches have degrees `6+18+356+82+64`. The degree-356 branch is marked as having no physical solutions. Its parameter polynomial, coordinate maps, and independent arithmetic exclusion certificate are absent from this bundle. The retained degree-18 branch is included and all its 18 candidates are independently excluded here. Consequently, the 174 retained roots are not claimed to enumerate all complex algebraic candidates. The missing 356-branch physical exclusion remains an explicitly cited published input.

The page numbers in the original 9-page PDF are 7 for §5.3–5.4, 8 for §5.7, and 9 for Figure 4. The 10-page HAL PDF adds a cover, making these pages 8, 9, and 10.

The eleven six-point components and the original saturated Gröbner decomposition have not been independently reconstructed here. To claim a fully standalone classification, one would additionally need exact complete decomposition certificates for each original coefficient chart, including the missing seven-point physical exclusion. Replaying the present root counts does not establish that missing reverse-coverage statement.

Frozen author sources:

- [Six-point coefficient data](https://web.archive.org/web/20240423101624id_/https://www-polsys.lip6.fr/~jcf/vortices/vortices_N6.txt). The literal `N=7` header is a metadata error: the file name and variables `e2,...,e6` identify the six-point system; exact substitution confirms this.
- [Seven-point retained coefficient data](https://web.archive.org/web/20240423101625id_/https://www-polsys.lip6.fr/~jcf/vortices/vortices_N7.txt).

## Run the checks offline

From this bundle directory, run:

```text
python -B -S -O verify_manifest.py
python -B -S -O certify_fs6_disks.py --verify
python -B -S -O certify_fs7_disks.py --verify
python -B -S -O verify_fs6_components_exact.py
```

All four commands must exit with code 0. `-B` avoids writing bytecode, `-S` disables site packages, and `-O` tests that certificate decisions remain active in optimized Python. The verifiers rewrite the three result JSON files in this same directory. Their decisive checks use integers and `fractions.Fraction`; no floating-point comparison establishes a root count, force exclusion, or value inequality. The `--generate` paths are not needed for replay and their optional numerical dependencies are not included in this verification bundle.

Expected final counts:

```text
n6: algebraic_solutions=80, force_excluded=69, below_target=10, central_pentagon=1
n7: algebraic_solutions_checked=174, force_excluded=155, strictly_below=18, center_hexagon=1
```

## Files

- `archivepolsys_vortices_N6.txt`, `archivepolsys_vortices_N7.txt`: frozen author chart data, parsed as polynomial expressions and never executed.
- `fs6_disk_certificate.json`, `fs7_disk_certificate.json`: frozen Gaussian-rational disk centers and radii.
- `certify_fs6_disks.py`, `certify_fs7_disks.py`, `parse_fs7_exact.py`: exhaustive retained-chart disk coverage and comparisons.
- `verify_fs6_components_exact.py`, `verify_fs6_dual_standard.py`, `fs6_resultant_data.json`: six-point forward-containment checks.
- `fs6_disk_verification_result.json`, `fs7_disk_verification_result.json`, `fs6_exact_component_checks.json`: independently replayed results.
- `proof_fragment.tex`, `coverage_audit.md`: mathematical proof fragment and detailed coverage audit.
- `SHA256_manifest.json`, `verify_manifest.py`: integrity check for the archived inputs, scripts, proof fragment, audit, this README, and the frozen replay results. Check it before replay; platform-dependent newline handling can change regenerated result-file bytes while preserving their content.

The mathematical explanation is `proof_fragment.tex`, a fragment for the existing theorem environments and bibliography, with the detailed coverage audit in `coverage_audit.md`. The original editable copies are `../n6_n7_improved_proof.tex` and `../n6_n7_improved_proof_notes.md`.
