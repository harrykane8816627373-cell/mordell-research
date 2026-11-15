# 十点全局最大构型的内层人口

取质量 $\sum_i|x_i|^2=9$。每个全局最大构型恰有 **2 或 3 点**的平方半径小于 $3/10$，其余 8 或 7 点都位于[已验证的外层允许域](../ten-point-ward-radii/README.md)。该必要条件不假定全局最大点对称。候选的内层人口为 2；这里没有排除全局最大点有 3 个内层点。

## 径向 Gram 能量预算

令 $f=h-\log d$，沿用[全局上界](../ten-point-ward-upper/README.md)的核、常数 $B$ 和两点松弛 $c=3/10000$。在全局最大构型上，Ward 与径向力修正的总和为零，因此

$$
\frac12\sum_{i,j} K_{ij}=B-\log\Delta-\sum_{i<j}f_{ij}
\le B-\log V_{10}^{(9)}+45c.
$$

各谐波的 Gram 能量都非负。令 $S=\sum_i F_0(u_i)$ 为径向特征和，则 $S^TQ_0S\le\sum_{i,j}K_{ij}$。复核器用精确见证和有理对数下界验证 $S^TQ_0S<27001/10^6$。由已验证的半径禁区和严格上端，所有 $u_i$ 都属于以下闭区间的并，端点放宽只增大必要域：

$$
[0,3/10]\quad\cup\quad[41/60,71/50]\quad\cup\quad[149/100,1539/1000].
$$

## 将能量预算变成人口排除

对于拟排除的内层人口 $m$，从[数据](dual.json)取有理向量 $s_0$ 和三常数 $a,b_{\rm in},b_{\rm out}$。由 $Q_0$ 正半定，

$$
S^TQ_0S\ge 2s_0^TQ_0S-s_0^TQ_0s_0.
$$

在内层区间精确验证多项式 $2s_0^TQ_0F_0(u)$ 大于 $au+b_{\rm in}$，在两个外层区间上验证大于 $au+b_{\rm out}$。当 $m=0$ 不需要内层比较。逐点求和并使用 $\sum_i u_i=9$，得到

$$
S^TQ_0S>9a+mb_{\rm in}+(10-m)b_{\rm out}-s_0^TQ_0s_0.
$$

对于 $m=0,1,4,5$，右端分别严格大于 $14/10,3/100,43/1000,7/10$，均大于能量预算。若 $m\ge 6$，即使所有内层半径取 $3/10$、其他点取 $1539/1000$，总质量也小于 $9$。于是仅剩 $m=2,3$。先前[缩并人口条件](../ten-point-collapse-domain/README.md)也独立排除了 $m\ge4$。

## 精确复核

先复核上游[全局上界](../ten-point-ward-upper/README.md)和[半径禁区](../ten-point-ward-radii/README.md)，再运行：

```sh
python -B -S -O ten-point-ward-population/verify_all.py
```

只用 Python 标准库。固定的输入与上游复核记录按 SHA-256 绑定。接受器重算 $Q_0$ 的正定分解、总能量比较和各区间的有理 Bernstein 系数。向量与常数由数值凸优化发现，接受器将其存储字面量解释为精确有理数，并重做完整的连续区间证明。可加 `--record-report` 更新[复核报告](verification_report.json)。
