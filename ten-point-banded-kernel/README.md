# 十点两径向带的全局上界核

设 $M_{10}$ 为十点最大距离积，$V_{10}$ 为固定整数见证的距离积，$C_{10}$ 为认证的严格局部极大构型的距离积。三者统一归一到同一质量。新证书给出

\[
V_{10}<C_{10}\le M_{10}<\frac{201}{200}V_{10},\qquad
\frac{C_{10}}{M_{10}}>\frac{200}{201}.
\]

采用质量 9 时，每个全局最大构型恰有两个平方半径落在 $I=(43/1000,51/200)$，八个落在 $O=(37/50,139/100)$。计数、两段半径和尾段排除已分别认证。[内点计数](../ten-point-ward-inner-count/README.md) · [计数行半径](../ten-point-counted-radii/README.md) · [尾段排除](../ten-point-tail-exclusion/README.md)

## 逐点恒等式及两点上界

定义 $u=|z_i|^2,v=|z_j|^2,t=\operatorname{Re}(z_i\bar z_j),d=u+v-2t$。当 $u\in I$ 时令 $\chi(u)=1$，否则取 0。对每个标记点 $i$，下列四个表达式分别在九个 $j\ne i$ 上求和为零。

\[
\chi(v)-\frac{2-\chi(u)}9,\qquad
v-\frac{9-u}{9},\qquad
t+\frac u9,\qquad
\frac{u-t}{d}-\frac{5u}{9}.
\]

这四式依次来自已证的内点计数、质量、质心和完整驻点的径向力。将任一式乘以仅依赖标记半径 $u$ 的函数，沿同一行求和仍为零。证书使用 $(u-1)^j$ 或 $\chi(u)(u-1)^j$，并把每项与 $i,j$ 调换所得表达式相加，得到在所有无序点对上总和为零的对称修正。

新函数 $h(u,v,t)$ 沿用旧核的十个正定矩阵和预算 $B$。在原 Ward 与径向力修正中重新选取系数，并加入 19 个逐点恒等式方向。由 Gram 和非负及新旧修正的零和关系，任一全局最大构型满足 $\sum_{i<j}h_{ij}\le B$。[原 Ward 核及正定性](../ten-point-ward-upper/README.md)

精确分区对全部点对给出 $h-\log d\ge-1/10000$。因此在质量 9 下

\[
\log M_{10}\le B+\frac{45}{10000}.
\]

对 $\log((201/200)V_{10})$ 取精确有理下界，与右端比较的严格余量大于 $1/10000$。另外，任一最大构型的内点对缺额满足

\[
h_{II}-\log d_{II}\le B-\log V_{10}+44\cdot\frac{1}{10000}.
\]

这只用到内点对以外的 44 条边的松弛，可用于进一步排除内点对的位置。

## 连续分区和重放

对 $u\le v$，取 $I\times I,I\times O,O\times O$ 的闭包，并分开 $t\le0$ 与 $t\ge0$，得到六个根盒。每个根盒中的带指示函数为常数。将两点不等式乘以正数 $d$，完整展开多项式 $dh$。复核器用精确有理点交叉核对分段展开和原式，然后重算 $21\times21\times11$ Bernstein 系数及整数中点分割。在每个接受叶，使用对数的有理切线上界及向外舍入。排除叶必须通过严格的 Gram、点距、删除 Schwarz 或缩并条件。树覆盖全部六盒。[必要点对域](../maximizer-domain/README.md) · [缩并域](../ten-point-collapse-domain/README.md)

系数由有约束数值搜索发现。冻结证书保留旧正定 Gram 块，用线性规划拟合逐点恒等式和其他零和修正。每个新候选用域内采样、局部优化和 Bernstein 分割寻找违反点。复核器用精确有理数解释所有存储系数，并独立检查每个节点。完整重放统计与输入哈希随目录保存。[重放报告](verification_report.json) · [固定输入](source_manifest.json)

先按十点仓库指引完成上游证书，再从仓库根目录运行标准库重放。

```sh
python -B -S -O ten-point-banded-kernel/verify_all.py --progress
```
