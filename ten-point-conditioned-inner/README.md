# 联合两条驻点行约束十点内边

2026-11-15。蒋丰灿研究项目，OpenAI Codex 协作。

在总平方半径为 9 的尺度下，任一全局最大构型的两内点 (x,y) 满足

\[
\frac{41}{100}<|x-y|^2<\frac{11}{20},\qquad
\frac{\operatorname{Re}(x\bar y)}{|x||y|}<-\frac{77}{100}.
\]

所用输入是[内外点平方半径与中点的已证范围](../ten-point-refined-radii/README.md)和[分带核](../ten-point-banded-kernel-refined/README.md)。内点无需假定对径，也无需预设外点的排列。

## 共享外点的联合缺额预算

令 (u=|x|^2,v=|y|^2,t=\operatorname{Re}(x\bar y),d=u+v-2t\)。对一个外点 (w)，记 (r=|w|^2,p=\operatorname{Re}(x\bar w),d_x=u+r-2p\)。标记点 (x) 的八条外边满足

\[
\sum r=9-u-v,\qquad
\sum p=-u-t,\qquad
\sum (u-p)/d_x=5u-(u-t)/d.
\]

前两项来自质量与质心，第三项是完整力方程的径向投影去掉内边。对 (y) 则交换 (u,v)。这样同一条内边精确进入两条行等式。

设 (f=h-\log d\) 为分带核的缺额。将已知内带分为相邻小区间，并在每个小区间验证若干固定有理数 (A,B,a,b,c) 的外边下界：

\[
f(x,w)\ge A+Bu+ar+bp+c(u-p)/d_x.
\]

在八个外点上求和就给出显式内点对函数。加上 (y) 的外边下界和唯一内边缺额，可得 17 条边的联合下界。剩下的 28 条外外边各至少为 \(-3/40000\)，全构型的缺额和至多为 (B_0-\log V_{10}\)，其中 (B_0) 是分带核原预算，(V_{10}) 是固定整数见证。精确对数比较给出 17 条边的联合缺额严格小于 (2101/10^6\)。

## 连续覆盖与复现

外边证书将下界乘以 (d_x>0)。内边证书将联合下界乘以 (d>0)，新增修正为低次多项式。逐盒使用精确 Bernstein 系数与对数切线上界，要么验证所述角度与距离，要么排除超过联合预算的盒，要么用此前的半径、Gram、点距与中点等约束排除。两个内点可重排，因此只须覆盖所有无序内带小区间对。

先按[复现指南](../docs/REPRODUCING.md)运行分带核和半径收紧步骤，再执行 `python -B -S -O ten-point-conditioned-inner/verify_all.py`。接受器只使用 Python 标准库。[覆盖报告](verification_report.json)
