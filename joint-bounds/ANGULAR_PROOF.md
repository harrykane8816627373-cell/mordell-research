# 联合径向、角向矩与全局上下界：圆盘上界的严格有限 n 修正

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**。推导、代码及交叉审阅包含 OpenAI Codex 协作。2026-10-28。

本节补充现有库仑缺额方法的角向部分。主要结果有两项：给出可以与既有径向罚项**直接相加**的复幂和罚项；利用有限 Newton 恒等式，推出每个有限 n 的圆盘上界均有显式正修正。本节没有确定新的未知点数锐值，也没有比较所有历史上界包络的强弱。Dirichlet 分解、Fourier 正交与 Newton 恒等式均为标准工具，不作全球优先权声明。

## 1. 归一化与结论

沿用仓库的定义

```math
M_n=\max_{\sum|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i<j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2,
```

```math
U_n=n^{n/2}\left(\frac{2n}{n-1}\right)^N e^{-N/2},\qquad
\delta(z)=\log\frac{U_n}{\Delta_n(z)}.
```

令库仑尺度点、平滑半径、复幂和为

```math
x_i=\sqrt{\frac{n-1}{n}}z_i,\qquad
r=\sqrt{\frac2n},\qquad
p_k(x)=\sum_{i=1}^n x_i^k.
```

取任何满足 \(R\ge\max_i|x_i|+r\) 的数。设 \(\mu=n^{-1}\sum_i\sigma_r(\cdot-x_i)\)，\(\mu_0\) 为半径 \(\sqrt2\) 均匀圆盘概率，\(\nu=\mu-\mu_0\)。对任意有限正整数 d，均有

```math
\boxed{\displaystyle
\delta(z)\ge n^2D(\overline\nu)
+\sum_{k=1}^d\frac{|p_k(x)|^2}{kR^{2k}},}
```

其中 \(\overline\nu\) 表示旋转平均，\(D\) 是[已有库仑分解](../asymptotics/README.md)中的对数能量。该结论适用于任意归一化无碰撞构型，不要求驻点、稳定或特殊角度。尤其可以将[已有径向单元证书](../configuration-bounds/radial_audit.md)的下界代入第一项。

再定义

```math
B_n=\left(1+\sqrt{\frac2{n-1}}\right)^2,\qquad
c_n=\frac{(\log2)^2}{\displaystyle\sum_{k=2}^n B_n^k/k}.
```

则对于每个整数 n≥2，都有

```math
\boxed{\displaystyle M_n\le U_n e^{-c_n}<U_n.}
```

这里改进的是指定的圆盘上界 \(U_n\)。例如小点数已有的专门锐界仍可远强于这个普适修正。\(c_n>0\) 但很小，并且趋于零；不能把它表述为线性渐近常数 γ 的改善。

## 2. 径向与角向能量的正交性

设 \(g(w)=-\log|w|\)、\(u=g*\nu\)。已有证明给出

```math
D(\nu)=\frac1{2\pi}\int_{\mathbb C}|\nabla u|^2,
\qquad \delta(z)\ge n^2D(\nu).
```

所有平滑电荷均紧支撑、具有有界密度且总质量零，因此这个 Dirichlet 表示以及下述弱分部积分适用。把 u 绕原点作旋转平均得 \(u_0\)，令 \(u_a=u-u_0\)。在每个圆周上，\(u_a\) 的角向平均为零。极坐标积分给

```math
\int\nabla u_0\cdot\nabla u_a=0,\qquad
D(\nu)=D(\overline\nu)+\frac1{2\pi}\int|\nabla u_a|^2.
```

事实上 \(\nu-\overline\nu=\mu-\overline\mu\)，因为背景 \(\mu_0\) 是径向的。因此角向电荷的支撑包含于 \(|w|\le R\)。**不需要背景圆盘也包含在 R 内。** 忽略这一点而把 \(R\ge\sqrt2\) 当作前提，会不必要地削弱证书。

## 3. 角向测试函数、系数与可加性

对整数 k≥1，考虑连续函数

```math
\phi_{k,c}(w)=\begin{cases}
\Re(w^k),&|w|\le R,\\
\Re(R^{2k}/\overline{w}^{k}),&|w|>R,
\end{cases}
```

以及将 Re 换成 Im 的 \(\phi_{k,s}\)。它们在圆内、圆外分别调和，跨圆周连续，在无穷远衰减，梯度平方可积。圆周上法向导数可以跳跃；测试函数的弱分部积分不要求导数连续。

记 Dirichlet 内积为 \(\langle f,h\rangle_D=(2\pi)^{-1}\int\nabla f\cdot\nabla h\)。直接对极坐标积分得到

```math
\|\phi_{k,c}\|_D^2=\|\phi_{k,s}\|_D^2=kR^{2k}.
```

其中圆内和圆外各贡献 \(kR^{2k}/2\)。不同 k 的函数相互正交；同一 k 的正弦、余弦也正交；它们全都与径向函数正交。因此不能只计算圆外一半能量，随后仍使用上述完整系数。

由于 \(-\Delta u_a=2\pi(\nu-\overline\nu)\)，弱分部积分给

```math
\langle u_a,\phi\rangle_D=\int\phi\,d(\nu-\overline\nu).
```

对每个中心位于 \(x_i\) 的均匀圆盘，调和均值性质给

```math
\int w^k\,d\sigma_r(w-x_i)=x_i^k.
```

由于 \(\mu\) 的支撑在 \(|w|\le R\)，测试函数在这部分电荷上等于 Re 或 Im 的幂函数；其在任何径向电荷上的积分为零。故

```math
\langle u_a,\phi_{k,c}\rangle_D=\frac{\Re p_k(x)}n,
\qquad
\langle u_a,\phi_{k,s}\rangle_D=\frac{\Im p_k(x)}n.
```

对这组有限正交函数使用 Bessel 不等式，得到

```math
\frac1{2\pi}\int|\nabla u_a|^2
\ge\frac1{n^2}\sum_{k=1}^d\frac{|p_k(x)|^2}{kR^{2k}}.
```

与第2节相加，即得第1节的联合不等式。径向与角向罚项占用正交的能量，所以能够相加；若使用两个没有正交性证明的能量下界，通常只能取最大值。

## 4. 有限 Newton 恒等式强制至少一个非零幂和

下面的引理只涉及有限组复数，不调用优化条件。设 \(x_1,\ldots,x_n\) 不全为零，\(\rho=\max_i|x_i|>0\)。则

```math
\boxed{\displaystyle
\sum_{k=1}^n\frac{|p_k(x)|}{k\rho^k}\ge\log2.}
```

令 \(e_j\) 为这些复数的基本对称多项式，\(e_0=1\)，写

```math
P(t)=\prod_{i=1}^n(1-x_it)=\sum_{j=0}^n(-1)^je_jt^j.
```

选 \(|x_{i_0}|=\rho\)。在根 \(t=1/x_{i_0}\) 处，由三角不等式有

```math
1\le\sum_{j=1}^n\frac{|e_j|}{\rho^j}.
```

Newton 恒等式等价于形式幂级数关系

```math
P(t)=\exp\left(-\sum_{k\ge1}\frac{p_k(x)}k t^k\right).
```

求到次数 n 的系数，只使用 k≤n 的项。对每个系数逐项取绝对值，得到

```math
\sum_{j=1}^n\frac{|e_j|}{\rho^j}
\le\exp\left(\sum_{k=1}^n\frac{|p_k(x)|}{k\rho^k}\right)-1.
```

右侧是有限多项式的指数在正数处的值：把次数大于 n 的非负系数一起加入，只会增大它。这里没有在 \(|t|=1/\rho\) 处代入可能发散的无限对数级数。于是指数至少为2，引理得证。

## 5. 从必要角向能量得到无条件上界

任一全局最大构型必中心化：如果重心非零，平移至重心不改变距离乘积，却减小总二阶矩，随后放大恢复二阶矩会严格增大乘积。因此在全局最大者处 \(p_1(x)=0\)，且

```math
\sum_i|x_i|^2=n-1,\qquad \rho^2\ge\frac{n-1}n.
```

取 \(R=\rho+r\)，在第4节的引理中去掉为零的 k=1 项，然后用 Cauchy–Schwarz：

```math
(\log2)^2\le
\left(\sum_{k=2}^n\frac{|p_k(x)|}{k\rho^k}\right)^2
\le\left(\sum_{k=2}^n\frac{|p_k(x)|^2}{kR^{2k}}\right)
\left(\sum_{k=2}^n\frac{(R/\rho)^{2k}}k\right).
```

又有

```math
\frac R\rho=1+\frac r\rho
\le1+\sqrt{\frac2{n-1}}.
```

第1节的能量不等式于是给

```math
\log\frac{U_n}{M_n}\ge
\frac{(\log2)^2}{\sum_{k=2}^n(R/\rho)^{2k}/k}
\ge c_n>0.
```

该证明对 n=2 也适用：求和仍含 k=2 一项，没有要求存在至少三个点。所有使用全局最优性的步骤只有其存在、无碰撞、中心化以及固定的二阶矩；没有引入任何未完成的驻点分类。

另一方面，\(B_n^n\) 趋于无穷，且 \(\sum B_n^k/k\ge\sum_{k=2}^n1/k\) 也趋于无穷，所以 \(c_n\to0\)。它不改变 \(\log M_n\) 的二次项、\(n\log n\) 项或线性项常数 γ。

## 6. 同半径数据中角向矩能够排除的八点子类

取 n=8，库仑尺度平方半径为

```math
\{|x_i|^2\}=\{0,1,1,1,1,1,1,1\}.
```

原点加单位正七边形在此径向类内；转换到 z 尺度，恰达到候选值 \(C_8=8^{28}/7^{21}\)。因此只依赖这些径向数据的正确证书，不能排除整个类。

但加入条件 \(|p_2(x)|^2\ge2\) 后，\(r=1/2\)、\(R=3/2\)，角向罚项单独就给

```math
\delta(z)\ge\frac2{2(9/4)^2}=\frac{16}{81}.
```

已有可达构型给全局最优缺额上界

```math
\log(U_8/M_8)\le\log(U_8/C_8)
=40\log2-7\log7-14<\frac19<\frac{16}{81}.
```

所以这个同时规定径向数据与二阶复幂和的整类，不含全局最大者。对数比较可用下节的有理区间验收；上述小数近似不承担证明。

该子类并非空集。取原点，以及外圆上由三个三次单位根、方向 \(\pi/6\) 的一对对径点、方向 \(\pi/2\) 的一对对径点组成的七点。七个外点互异、总和零，且

```math
p_2(x)=2e^{i\pi/3}+2e^{i\pi}=-1+i\sqrt3,
\qquad |p_2(x)|^2=4,\qquad p_3(x)=3.
```

同一径向数据允许候选最优者以及这个可排除的非规则角向构型，后者还有 \(p_2p_3\ne0\)，说明新约束提供了径向数据之外的信息。它没有覆盖一般 \(p_2p_3\ne0\) 的全部八点分支。

## 7. 纯有理实现与范围

[angular_bounds.py](angular_bounds.py) 返回 Fraction 证书。设 \(b\ge B_n\) 为有理外包，\(\ell\le\log2\) 为正有理下端，则

```math
0<\widehat c_n:=\frac{\ell^2}{\sum_{k=2}^n b^k/k}\le c_n,
\qquad M_n\le U_ne^{-\widehat c_n}.
```

根式外包由整数平方根及向上取整构造。取 y=1/3，有

```math
2\sum_{j=0}^{s-1}\frac{y^{2j+1}}{2j+1}
\le\log2\le
2\sum_{j=0}^{s-1}\frac{y^{2j+1}}{2j+1}
+\frac{2y^{2s+1}}{(2s+1)(1-y^2)}.
```

同样的 atanh 公式可用于其他正有理数的对数。八点例可将 \(\log7=2\log2+\log(7/4)\) 先作尺度缩减，再在 y=3/11 处包围后一项。

对联合区域证书，如果平方半径带的最大上端为 q，则程序采用安全的有理数

```math
R^2\le q+\frac2n+2\operatorname{sqrtUpper}\left(\frac{2q}n\right)
```

作为分母中的支撑平方半径上端。这里左边的 R 指实际最小充分支撑半径 \(\sqrt q+\sqrt{2/n}\)，程序传入的是右边的有理上包；增大分母只会减小罚项。输入的各个 \(|p_k(x)|^2\) 下端必须独立证明，程序不会从一个假设自动推出另一个假设。

输入有理坐标时，先中心化为 v，设 \(Q=\sum|v_i|^2\)，再取 \(x=\sqrt{(n-1)/Q}\,v\)。虽然 x 的坐标可能含根式，但

```math
|x_i|^2=\frac{n-1}Q|v_i|^2,\qquad
|p_k(x)|^2=\left(\frac{n-1}Q\right)^k|p_k(v)|^2
```

均为有理数。因此精确构型与新证书可以直接衔接，无需复浮点数。

有限程序核验检查实现、常数和具体证书；任意 n 的量词由本文解析证明承担。小幂和、径向分布良好、满足联合必要条件都不能反推全局最优性。全八点锐界、一般 n 的简洁锐值公式、精确 γ 仍然开放；本节不声称外部同行评审或整篇形式化证明。

## 8. 解析输入与文献来源

本节的外部框架输入只到仓库[渐近专题第3、5节](../asymptotics/README.md)已经独立推导的圆盘平滑分解、非负剩余项与零质量 Dirichlet 表示为止。由这个输入出发，第2至5节完整给出本节有限 n 改善的证明，无需使用外部重整化能量渐近定理。

经典圆盘涂抹上界的现代证明见 A. B. Lauritsen, *Floating Wigner Crystal and Periodic Jellium Configurations*, J. Math. Phys. **62**, 083305 (2021)，[arXiv:2103.07975v2](https://arxiv.org/abs/2103.07975v2)，Appendix B，Proposition B.1。该附录明确把带电点平滑为均匀圆盘并保留总中性电荷能量、重叠代价等非负项；这是经典 Lieb–Narnhofer／Sari–Merlini 界的背景，不归入本项目首次成果。

一般库仑分解与线性渐近项的背景见 E. Sandier and S. Serfaty, *2D Coulomb Gases and the Renormalized Energy*, [arXiv:1201.3503v2](https://arxiv.org/abs/1201.3503v2)。它在本仓库承担已注明的渐近定理输入；本节的有限角向证明、Newton 估计及 c_n>0 不依赖其渐近余项、结晶猜想或三角晶格在全部允许构型中的最优性。
