# The sharp twelve-point bound in Mordell's distance-product problem

**9 November 2026** · Fengcan Jiang's Mordell research project

[Read the PDF](mordell-n12-proof.pdf) · [Download the LaTeX source](mordell-n12-proof.tex) · [Proof certificate](../README.md) · [Replay guide](../../docs/REPRODUCING.md)

The paper proves the sharp value at twelve points and classifies all maximizers up to rotation and relabeling. The three radial layers contain 3, 3, and 6 points. The outer layer is the union of two regular triangles on a common circle.

The proof combines a Ward deficit budget, radial population bounds, exact interpolation at fourteen pair contacts, and a 728,050-node continuous cover with zero slack. The last section develops constructions for other point counts, including higher-order Ward corrections and moment features that recover relative phases.

## Build the paper

The TeX file contains the complete English manuscript and its TikZ figure. With TeX Live installed, run the following command in this directory.

```sh
latexmk -pdf mordell-n12-proof.tex
```

## Replay the proof

The manuscript cites [the fixed certificate revision](https://github.com/harrykane8816627373-cell/mordell-research/commit/3c922f3446f4f32441d75aef1d7907d38b0544f7). From the repository root, the four commands below replay the relevant checks in dependency order using the Python standard library.

```sh
python -B -S -O twelve-point-removal-domain/verify_all.py
python -B -S -O twelve-point-ward-upper/verify_all.py
python -B -S -O twelve-point-ward-radii/verify_all.py
python -B -S -O twelve-point-sharp-global/verify_all.py
```

[Accepted report](../verification_report.json) · [Fixed input manifest](../source_manifest.json) · [Repository home](../../README.md)

## File integrity

The PDF has 17 pages and was compiled with pdfTeX 1.40.29 through Latexmk.

|File|SHA-256|
|---|---|
|PDF|`7fc35ddbba2c3389d6954b01e70c0b434064c6c9ba3c476f897376ff03c9faa4`|
|LaTeX|`d4fc1f6ec8ad0f856094ea5c89d6649f54bed55cf6cfaacee33c93ed6663bc1a`|
