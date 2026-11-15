# 十点全局最大构型的内点对角度与中点

2026-11-15。蒋丰灿研究项目，OpenAI Codex 协作。

归一化为 \(\sum_i|x_i|^2=9\)。任一全局最大构型的两个内点记为 \(x,y\)，令 \(u=|x|^2,v=|y|^2,t=\operatorname{Re}(x\bar y)\)。精确证书给出

\[
0.045<u,v<0.245,\qquad -0.140<t<-0.094,\qquad
\left|\frac{x+y}{2}\right|^2<\frac{1}{25},\qquad
\frac{t}{\sqrt{uv}}<-\frac{11}{20}.
\]

最后一个不等式给两个内点的夹角大于 \(\arccos(-11/20)\)。由完整构型居中，八个外点的重心是 \(-(x+y)/8\)，其模严格小于 \(1/20\)。

记 \(f=h-\log d\) 为[分带核](../ten-point-banded-kernel/README.md)的点对缺额，\(d=u+v-2t\)。由于恰有两个内点，唯一一条内边之外的 44 条边属于内外边或外外边。每条这类边的缺额至少为 \(-1/10000\)。总缺额至多为 \(B-\log\Delta\)，而核预算与可达见证满足 \(B-\log V_{10}<10^{-6}\)。这将内边缺额限制为

\[
 f(x,y)<44\cdot\frac{1}{10000}+10^{-6}=\frac{4401}{10^6}.
\]

较早的[内点对证书](../ten-point-counted-radii/README.md)给 \(t<0\)。在其允许的内带中覆盖全部负内积点对。每个终止盒或者全部满足几何结论，或者由精确 Bernstein 上界与对数切线证明 \(f\ge4401/10^6\)，或者违反已证必要域。角条件用 \(9t^2>4uv\) 与 \(t<0\) 判定。

按[复现指南](../docs/REPRODUCING.md)重放分带核后，在仓库根目录运行 `python -B -S -O ten-point-inner-pair/verify_all.py`。[复核报告](verification_report.json)
