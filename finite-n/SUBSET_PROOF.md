# 从已有严格子集矩定理到显式有理余量

研究资料与项目发起：蒋丰灿（Fengcan Jiang）；本次证明整理与复核包含 OpenAI Codex 协作。日期：2026-10-27。

已有的[一般 n 结构证明](../earlier-results/general_n_structural.tex)已经建立：在总二阶矩固定的平方距离乘积问题中，每个全局最优构型的真子集，都满足严格的内部二阶矩下界。本节核验这一已有结论，并补充一个不用边界 Hopf 引理的初等证明、完整交互多项式的系数约束，以及对所有允许点数均显式有效的有理余量。

这些结论均为**全局最优构型的必要条件**。本节没有把它们推广为任意局部极大点或驻点的性质，也没有由必要条件直接推出原问题的有限点数下构造。

## 1. 问题、已有结果与条件多项式

所有点位于复平面。记

```math
\Delta_n(z)=\prod_{i\lt j}|z_i-z_j|^2,\qquad
M_n=\max_{\sum_i|z_i|^2=n}\Delta_n(z),\qquad
N=\frac{n(n-1)}2,\quad \lambda=\frac{n-1}2.
```

取一个达到最大值的构型。它没有重合点，并且重心为零：平移到重心不改变距离乘积，却降低非零重心构型的总二阶矩；随后放大恢复规定二阶矩，会严格增大乘积。

由齐次性，同一构型也在整个复配置空间上全局最大化

```math
\Phi(z)=\Delta_n(z)\exp\!\left(-\lambda\sum_i|z_i|^2\right).
```

事实上，对任何固定非零形状，最优缩放后的总二阶矩恰为

```math
\frac N\lambda=n.
```

固定一个大小为 m 的真子集，至少包含两个点。定义

```math
a_S=\frac1m\sum_{i\in S}z_i,\qquad
y_i=z_i-a_S,\qquad
Q_S=\sum_{i\in S}|y_i|^2,\qquad
k=\frac{m(m-1)}2,\quad \alpha=\lambda Q_S.
```

只把子集内各点替换为围绕其重心的同一复相似变换，其余点固定。条件目标函数除去一个与参数无关的正因子后，等于

```math
|t^kq(t)|^2e^{-\alpha|t|^2},\qquad
q(t)=\prod_{i\in S,\ j\notin S}(a_S+t y_i-z_j).
```

原构型对应的参数在复平面上全局最优。因此

```math
q(1)\ne0,\qquad
|q(1)|=\max_{|t|=1}|q(t)|,\qquad
\frac{q'(1)}{q(1)}=\alpha-k\in\mathbb R.
```

最后一个等式同时使用实方向与虚方向的一阶导数。由于子集内的点互异、偏移向量之和为零，至少两个偏移向量非零；没有任何交互因子会成为零多项式。因此 q 非常数，其次数满足

```math
2(n-m)\le d:=\deg q\le m(n-m).
```

已有严格子集矩定理就是

```math
\boxed{\quad
Q_S\gt \frac{m(m-1)}{n-1}\qquad(2\le m\lt n).
\quad}
```

## 2. 初等边界恒等式及严格性的证明

下面的恒等式给出该定理的自含证明，同时保留一个可用于必要约束的非负量。设多项式满足

```math
f(1)=1,\qquad |f(e^{i\theta})|\le1,
\qquad
f(t)=\sum_{j=0}^{d}a_jt^j.
```

定义差商多项式及其系数：

```math
R(t)=\frac{f(t)-1}{t-1}
=\sum_{j=0}^{d-1}r_jt^j,
\qquad r_j=\sum_{k=j+1}^{d}a_k.
```

圆周上的模长在指定边界点达到最大值，角向导数为零，因此

```math
\Im f'(1)=0.
```

有精确恒等式

```math
\boxed{
\Re f'(1)=\sum_{j=0}^{d-1}|r_j|^2
+\frac12\int_{-\pi}^{\pi}
 \frac{1-|f(e^{i\theta})|^2}{1-\cos\theta}\,
 \frac{d\theta}{2\pi}.}
```

积分在原点的商按连续延拓理解：其分子是非负光滑函数，在原点取零值，一阶导数也为零，所以商有有限极限。

证明只需初等三角恒等式和 Parseval 恒等式。先由有限几何级数得到

```math
\int_{-\pi}^{\pi}
\frac{1-\cos(k\theta)}{1-\cos\theta}
\frac{d\theta}{2\pi}=k\qquad(k\ge1).
```

将多项式展开代入下式，余弦项给出导数实部，正弦项的组合是可积奇函数，积分为零：

```math
\int_{-\pi}^{\pi}
\frac{1-\Re f(e^{i\theta})}{1-\cos\theta}
\frac{d\theta}{2\pi}=\Re f'(1).
```

其中正弦组合的线性项因前述角向导数条件消失。再使用

```math
1-\Re f=\frac12|1-f|^2+\frac12(1-|f|^2),
\qquad
|1-e^{i\theta}|^2=2(1-\cos\theta)
```

以及差商多项式的 Parseval 恒等式，即得所述精确分解。

右端两项均非负；非常数多项式的差商不恒为零，故导数实部严格为正。应用于交互多项式的归一化

```math
f(t)=\frac{q(t)}{q(1)}
```

即证明已有的严格子集矩定理。整个严格性证明无需调用 Hopf 引理。

此外，保留差商的完整系数，得到

```math
\boxed{
\lambda Q_S-k\ge
\sum_{j=0}^{d-1}
\left|\sum_{k=j+1}^{d}\frac{[t^k]q(t)}{q(1)}\right|^2.}
```

这里方括号表示提取对应幂次的系数。将等式两边乘以分母模长的平方，可把这个必要条件写成无除法形式：

```math
(\lambda Q_S-k)|q(1)|^2\ge
\sum_{j=0}^{d-1}
\left|\sum_{k=j+1}^{d}[t^k]q(t)\right|^2.
```

## 3. 对所有点数有效的显式有理余量

对允许的整数点数定义

```math
D=m(n-m),\qquad
B=\frac{3n-m-1}{m-1},\qquad
\varepsilon_{n,m}=
\frac{1}{2(1+B)\displaystyle\sum_{j=0}^{D-1}B^j}.
```

则可把已有严格结论进一步量化为

```math
\boxed{
\lambda Q_S-k\gt \varepsilon_{n,m},\qquad
Q_S\gt \frac{m(m-1)+2\varepsilon_{n,m}}{n-1}
\quad(2\le m\lt n).}
```

余量是显式的正有理数。某些子集的余量非常小，因此不声称它对所有分支都有显著的数值改善。

证明的关键是先为交互多项式找到一个位置统一受控的零点。为避免与复相似变换参数混淆，以下用 h 表示补集点数，并记

```math
h=n-m,\qquad s=Q_S,\quad u=Q_{S^c},\qquad
s_0=\frac{m(m-1)}{n-1},\quad
u_0=\frac{h(h-1)}{n-1}.
```

总二阶矩的重心分解为

```math
n=s+u+\frac{mh}{n}|a_S-a_{S^c}|^2.
```

设 A 为补集各点到子集重心的平方距离平均值，则

```math
A=\frac uh+|a_S-a_{S^c}|^2,\qquad
\frac{mA}{s}=\frac{n(n-s)}{hs}-\frac us.
```

选择子集内偏移模长最大的点，以及补集中到子集重心距离最小的点。它们满足

```math
|y_i|^2\ge\frac sm\gt 0,\qquad
|z_j-a_S|^2\le A.
```

于是交互多项式有零点

```math
\tau=\frac{z_j-a_S}{y_i},\qquad
q(\tau)=0,\quad \tau\ne1,
\qquad
|\tau|^2\le\frac{mA}{s}.
```

已证明的严格子集定理，以及单点补集的内部二阶矩为零，统一给出

```math
s\gt s_0,\qquad
u\ge u_0:=\frac{h(h-1)}{n-1}.
```

补集内部二阶矩不超过总二阶矩，而补集点数小于总点数，所以

```math
0\le u\le n,\qquad h\lt n,\qquad n^2/h-u\gt 0.
```

相应表达式随 u 递减，并随 s 严格递减。因此

```math
|\tau|^2
\le \frac{n^2/h-u}{s}-\frac nh
\lt \frac{n(n-s_0)-hu_0}{hs_0}
=\frac{3n-m-1}{m-1}=B.
```

最后的代数恒等式可直接令总点数等于两个子集点数之和后展开核验。将上节差商多项式在该零点取值，再用 Cauchy–Schwarz 不等式：

```math
\frac1{|\tau-1|^2}=|R(\tau)|^2
\le\left(\sum_{j=0}^{d-1}|r_j|^2\right)
\left(\sum_{j=0}^{d-1}|\tau|^{2j}\right).
```

由次数上界及零点模长的严格上界，

```math
|\tau-1|^2\le2(1+|\tau|^2)\lt 2(1+B),\qquad
\sum_{j=0}^{d-1}|\tau|^{2j}
\le\sum_{j=0}^{D-1}B^j.
```

故差商系数平方和严格大于所定义的有理余量。结合上节非负恒等式，完成证明。

这一证明没有循环：先证明已有严格定理，再用它控制一个零点的位置，最后由系数恒等式得到更强的显式余量。

## 4. 子集重心及八点半径的改进

为统一处理单点子集，规定单点情形的余量为零：

```math
\varepsilon_{n,1}=0.
```

对任意非空真子集，同时对子集和补集应用矩下界，得到

```math
\left|\sum_{i\in S}z_i\right|^2
\lt \frac{2m^2(n-m)^2}{n(n-1)}
-\frac{2m(n-m)}{n(n-1)}
\left(\varepsilon_{n,m}+\varepsilon_{n,n-m}\right).
```

特别地，每个点的半径平方满足

```math
|z_i|^2\lt \frac{2(n-1)}n
-\frac{2\varepsilon_{n,n-1}}n
\qquad(n\ge3).
```

八点情形取七点子集时，所有常数都可以精确写成有理数：

```math
B=\frac83,\qquad D=7,\qquad
\varepsilon_{8,7}=\frac{2187}{9217846}.
```

因此七点内部二阶矩与每个点的半径平方分别满足

```math
\boxed{
Q_S\gt 6+\frac{2187}{32262461},\qquad
|z_i|^2\lt \frac74-\frac{2187}{36871384}.}
```

这给已有八点紧化约束一个小但严格的有理改进。它不依赖预设对称性，也没有排除尚未解决的一般非零矩分支。

## 5. 适用范围与有限点数下界的区分

两点情形没有满足条件的真子集。单点内部二阶矩为零；全体点的内部二阶矩等于规定总量，原下界取等号。因此“至少两点”与“真子集”两个条件都不可删除。

用于闭合的有理外松弛时，可以把上述严格不等式改成相同右端的非严格不等式，仍得到所有全局最优构型必须满足的条件。该放松便于代数核验，但不改变其必要条件性质。

这些结构限制本身不是有限点数下界构造。给出有限点数的乘积下界，仍需展示一个具体构型、明确缩放到规定总二阶矩，并从下方计算或估计它的距离乘积。本节没有把必要约束反向当作存在性结论。
