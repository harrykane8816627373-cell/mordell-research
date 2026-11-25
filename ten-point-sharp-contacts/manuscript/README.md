# The Sharp Ten-Point Bound for Mordell's Distance Product

Fengcan Jiang, 26 November 2026.

[Read the PDF](mordell-n10-proof.pdf) · [LaTeX source](mordell-n10-proof.tex) · [Proof certificate](../README.md) · [Final verification report](../global_verification_report.json)

We determine the maximum product of the 45 squared pairwise distances among ten planar points with total squared radius ten. The exact value is algebraic and is approximately 5.3823025895657 × 10¹⁰. The unique maximizing orbit consists of one antipodal pair and two concentric rectangles. The paper includes the algebraic root definition, the cyclic classes used to construct the majorant, and the proof of the equality conditions.

## Cite the preprint

Jiang, Fengcan. (2026). *The Sharp Ten-Point Bound for Mordell's Distance Product* (Version 1). Zenodo. [doi:10.5281/zenodo.24537890](https://doi.org/10.5281/zenodo.24537890)

## Build the paper

Compile the single LaTeX source file with TeX Live from this directory.

```sh
latexmk -pdf mordell-n10-proof.tex
```

## Replay the proof

Check out the [fixed proof revision](https://github.com/harrykane8816627373-cell/mordell-research/commit/cf631a15c2d62b61b6f40b5e16092a4c5a189978) with newline conversion disabled. From the repository root, run the full dependency replay or run the final certificate after the dependencies have been accepted.

```sh
python -B -S -O scripts/verify_ten_point.py
python -B -S -O ten-point-sharp-contacts/verify_sharp.py
```

The replay uses Python's standard library. The final stage checks all six pair domains with 50,267 subdivision nodes, zero slack, and zero unresolved boxes. [Replay guide](../../docs/REPRODUCING.md) · [Candidate root and exact value](../../ten-point-contact-geometry/README.md)

## File integrity

The seven-page PDF was compiled with pdfTeX 1.40.29 through Latexmk.

| File | SHA-256 |
|---|---|
| PDF | `b74696369c6f563c40efad4d40ab3fb24b5b4e09b4efed754db881c04c5b7f62` |
| LaTeX | `94f7bd798a75376a490070d532f26c608cdf517ae37f8b976816c1ef23e3ece3` |
