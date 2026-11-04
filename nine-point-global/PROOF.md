# 九点锐界：必要域、全域证书及全部等号

整理与复核：2026-11-04。本文对上游证明作中文重组与补充，保留原作者署名，改写部分遵循上游CC BY 4.0授权。定理、辅助核及原始证书归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu**，见 [arXiv:2611.03088v1，第4节](https://arxiv.org/pdf/2611.03088v1)及[固定代码版本](https://github.com/Jxxxx-2001/mordell-discriminant/tree/31cc934b7dfa0ef2308566a1bce48cbfb2e01559)。以下是带明示计算步骤的推导，不把上游已保存的成功报告当作本仓库的执行记录。原始文本逐字节保存于 `vendor/`；其中 `BASIC_RESULTS.txt` 的“未证”状态是早期历史记录，最终状态以原始 `PROOF_ZH.txt`、下述完整推导及本地实际重放为准。该历史文本还讨论另一距离比问题，其来源排除不作为本定理前提，本目录不重放那些无关结论。

**定理。** 若 $z_1,\ldots,z_9\in\mathbb C$ 且 $\sum|z_i|^2=9$，则

```math
\Delta(z)=\prod_{i<j}|z_i-z_j|^2\le \frac{9^{36}}{8^{28}}.
```

等号恰为原点与半径 $\sqrt{9/8}$ 的正八边形，允许整体旋转与重排。计算机辅助步骤是第4部分的三变量连续域不等式；其余推导不使用九点驻点完备分类。

## 1. 最大点与完整二阶必要条件

约束球面紧，目标连续，故最大值存在。原点加正八边形给出正值 $8^8(9/8)^{36}$，最大点因而无碰撞。若重心 $c\ne0$，居中后质量降至 $q=9-9|c|^2$，再放大 $\sqrt{9/q}$ 会将目标乘以 $(9/q)^{36}>1$。所以最大点居中。

在任一最大点，实Lagrange方程与对点对求和给出

```math
\sum_{j\ne i}\frac1{z_i-z_j}=4\bar z_i.
```

令 $A_{ij}=-(z_i-z_j)^{-2}$（$i\ne j$），对角元使行和为0。于是 $A1=0$、$Az=4\bar z$，而 $\log\Delta-4\sum|z_i|^2$ 的实二次变分为

```math
Q(v)=-2\operatorname{Re}(v^TAv)-8\|v\|^2.
```

它在球面切空间 $\operatorname{Re}(z^*v)=0$ 上非正。径向 $z$ 与切空间在此二次型下正交，且 $Q(z)=-16\|z\|^2<0$，因此在整个实空间非正。再对 $iv$ 应用同一条件，得到 $|\operatorname{Re}(v^TAv)|\le4\|v\|^2$。写 $A=B+iC$，对应实对称矩阵为 $H=\left(\begin{smallmatrix}B&-C\\-C&-B\end{smallmatrix}\right)$，故 $\|H\|\le4$。乘上实正交矩阵 $\operatorname{diag}(I,-I)$ 即为 $A$ 的实表示，所以 $\|A\|\le4$。此论证允许退化最大点。

令 $P=I-11^T/9$。由 $A=PAP$，$|A_{ij}|\le4\|Pe_i\|\|Pe_j\|=32/9$，从而 $|z_i-z_j|^2\ge9/32$。

## 2. 径向矩与全体最大点的半径界

记 $R_j=\sum|z_i|^j$、$S_j=\sum z_i^j$。直接展开，使用居中及驻点方程，有

```math
(Az^2)_i=8(|z_i|^2-1),\qquad
(Az^3)_i=12z_i|z_i|^2-15z_i,
```

```math
4\sum |z_i|^2z_i^2=\tfrac{15}{2}S_2,\qquad
(Az^4)_i=16z_i^2|z_i|^2-21z_i^2-S_2.
```

对去掉常向量投影的 $z^k-S_k/9$ 使用 $\|A\|\le4$。二次与三次试探向量分别给出

```math
3R_4+|S_2|^2/9\le36,\qquad
128R_6-360R_4+2025+(16/9)|S_3|^2\le0.
```

四次试探向量的平方范数为 $256R_8-672R_6+441R_4-9|S_2|^2$。与 $16(R_8-|S_4|^2/9)$ 比较，再用 $|S_2|^2\le324-27R_4$，丢去非负项，得

```math
240R_8-672R_6+684R_4-2916\le0.
```

现令 $x_i=\sqrt{8/9}z_i$、$u_i=|x_i|^2$、$T_j=\sum u_i^j$，则 $T_1=8$，并有

```math
T_2\le256/27,\quad F_3:=2T_3-5T_2+200/9\le0,
\quad F_4:=405T_4-1008T_3+912T_2-3072\le0.
```

置

```math
\rho=\frac{225-\sqrt{18465}}{90},\quad
\kappa=504-405\rho>0,\quad
L=\frac{14372-100\sqrt{18465}}{405}>0,\quad f(a)=a^2(a-\rho)^2.
```

精确恒等式 $F_4+\kappa F_3=405\sum f(u_i)-405L$ 推出 $\sum f(u_i)\le L$。每项非负，故

```math
0\le u_i\le\beta:=\frac{\rho+\sqrt{\rho^2+4\sqrt L}}2
< U:=\frac{1860281}{1048576}.
```

这包括 $u_i<\rho$ 的情形；当 $u_i\ge\rho$ 时，由 $u_i(u_i-\rho)\le\sqrt L$ 解二次式即得。原始 `verify_fourth_moment_bound.py` 精确核验系数恒等式并以Arb严格隔离 $1.77410189294<\beta<1.77410189295<U$。

## 3. 一对半径与剩余七点的共同约束

对最大点中的一对，记 $u=|x_i|^2$、$v=|x_j|^2$、$t=\operatorname{Re}(x_i\bar x_j)$、$s=u+v$、$a=u^2+v^2$、$V=256/27$。有 $t^2\le uv$、$d=u+v-2t\ge1/4$，且 $8-s>0$（由 $u,v\le U<2$）。

对其余七个非负平方半径用加权Cauchy，$T_3\ge u^3+v^3+(T_2-a)^2/(8-s)$。结合 $F_3\le0$，二次函数

```math
F(W)=2\left[u^3+v^3+\frac{(W-a)^2}{8-s}\right]-5W+200/9
```

在 $W\le V$ 上的最小值为 $F(W_0)$，其中 $W_0=\min\{V,a+5(8-s)/4\}$，所以 $F(W_0)\le0$。把真实 $T_2$ 的下界放宽到 $-\infty$ 只会削弱必要条件。

其次，对任何 $y\ge0$，

```math
f(y)\ge(\rho^2/4)\operatorname{dist}(y,\{0,\rho\})^2.
```

分别在 $y\le\rho/2$、$y\ge\rho/2$ 使用另一因子的平方下界即得。把其余七个平方半径各自舍入到最近的 $0$ 或 $\rho$；若有 $k$ 个舍入到 $\rho$，Cauchy给出

```math
f(u)+f(v)+\frac{\rho^2}{28}\min_{k=0,\ldots,7}(8-u-v-k\rho)^2\le L.
```

记 $D_9$ 为满足上述两条共同约束以及 $0\le u,v\le U$、$t^2\le uv$、$d\ge1/4$ 的域。每个全局最大点的每对点都在 $D_9$；不声称任意可行九点构型的点对都在此域。

## 4. 明示的计算机辅助引理

原始 `global_pair_bound_exact_contacts.json`（SHA256 `1d1cef17749fcb242e6f09dffe3c26183ae672fe4707b2bad70bb0b12b0e9d88`）定义八个 $6\times6$ 实对称矩阵 $Q_0,\ldots,Q_7$；元素属于 $\mathbb Q(\sqrt2)+\mathbb Q(\sqrt2)\log(2-\sqrt2)+\mathbb Q(\sqrt2)\log2$，全部有理系数精确保存。

令 $F(u)=(u-8/9,u(u-1),u^2(u-1),\ldots,u^5(u-1))^T$，$V(u)=(1,u-1,\ldots,(u-1)^5)^T$，$H_0=1$、$H_1=t$、$H_k=2tH_{k-1}-uvH_{k-2}$，并置

```math
K=F(u)^TQ_0F(v)+\sum_{k=1}^7H_kV(u)^TQ_kV(v),
\quad b=(8\log8-36)/9,
```

```math
p(u)=\frac{(9/2)u+b-K(u,u,u)/2}{8},\qquad
h(u,v,t)=p(u)+p(v)-K(u,v,t).
```

**计算引理。** 全部 $Q_k$ 严格正定，且 $h(u,v,t)\ge\log d$ 在 $D_9$ 上成立。重放入口 `verify_global.py` 核验如下完整证明链。

1. `verify_pair_bound_definition.py` 从矩阵独立重建全部多项式系数，并精确核验五类接触点 $(0,1,0)$、$(1,1,\sqrt2/2)$、$(1,1,0)$、$(1,1,-\sqrt2/2)$、$(1,1,-1)$ 的函数值与所需导数；交换 $u,v$ 由精确对称性处理。
2. `verify_pair_bound_local_preconditioned.py` 以Arb区间LDL核验八个矩阵正定及各接触邻域。内点用零梯度与正定Hessian。边界点令 $q=uv-t^2$、$g=h-\log d$，验证 $g-\tau q+\eta q^2$ 的Hessian正定及 $\tau-\eta q>0$。修正函数在接触点的值、梯度为零，故在邻域物理域 $q\ge0$ 内 $g\ge0$，包括正八边形对径点。可逆整数合同变换保持正定性。
3. `verify_pair_bound_cover.py` **从完整初始长方体** $[0,U]^2\times[-U,U]$ 重新二分，不读旧空前沿。仅当整盒在必要域外、落于认证接触邻域，或严格Bernstein系数下界为正时才关闭。$u>v$ 整盒按对称性排除，对角线保留。必要域排除有时故意少用条件，扩大覆盖范围，因而不会漏掉真实点对。
4. 对数使用切线上界或全局三次Taylor上界。对 $d,d_0>0$，$\log d\le\log d_0+\delta-\delta^2/2+\delta^3/3$，$\delta=(d-d_0)/d_0$，因为四阶导数严格为负。跨 $q=0$ 的盒还可检查 $g-\lambda q+\eta q^2\ge0$，条件为 $\lambda,\eta\ge0$ 且 $\eta q\le\lambda$；这在物理域是充分条件。
5. 浮点基本运算逐步向外取整；必要时用Arb重新包围全部Bernstein系数；二分中点必须精确。未决盒保留，时间耗尽或断言失败均不得认证成功。检查两种基转换的独立脚本也实际执行。

保存的上游报告与本地 `reports/replay.json` 分开。只有实际运行、空前沿、树守恒以及全部前置检查通过后，驱动才保存成功报告；`--check-inputs` 仅校验字节。程序依赖SymPy精确代数、python-flint/Arb包含运算、NumPy与IEEE-754运算语义，解析推导仍需阅读本证明。没有声称整篇形式化证明或外部同行评审。

## 5. 求和推出全局上界

任取一个存在的全局最大点，用其 $x$ 坐标。递推式给出 $H_k=\operatorname{Re}((x_i\bar x_j)^k)$。定义 $L_0=\sum_iF(u_i)$、$L_k=\sum_i x_i^kV(u_i)$，则

```math
\sum_{i,j}K_{ij}=L_0^TQ_0L_0+\sum_{k=1}^7L_k^*Q_kL_k\ge0.
```

又由 $\sum_i u_i=8$ 与对角修正的定义，

```math
\sum_{i<j}h(u_i,u_j,t_{ij})=8\log8-\frac12\sum_{i,j}K_{ij}.
```

计算引理适用于每个点对，故 $\log\Delta(x)\le8\log8$，从而 $\Delta(z)\le(9/8)^{36}8^8=9^{36}/8^{28}$。由于已约束任意全局最大点，此界对所有可行构型成立；候选直接达到该界。

## 6. 全部等号

若达到锐界，该构型就是全局最大点，所以上述推导全部适用。等号迫使非负矩二次型之和为零。每个 $Q_k$ 严格正定，因此 $L_0=\cdots=L_7=0$。$L_0$ 的第二、第三分量相减给出

```math
\sum_i u_i(u_i-1)^2=0.
```

每项非负，故 $u_i\in\{0,1\}$；九项之和为8，恰有一项为0，其余八项为1。对 $k=1,\ldots,7$，$L_k$ 的第一分量给出 $\sum_i x_i^k=0$。Newton恒等式迫使八个非零根的首一多项式为 $X^8+c$。根模长均为1，故 $|c|=1$，其根组成一个正八边形。还原尺度得全部等号构型；反向代入直接成立。

## 7. 从小点数到一般结构

本证明显示小点数最优性的一个可用机制：严格正定的矩核把全局能量差转化为径向与角向矩的非负二次型；径向矩差固定支撑在 $\{0,1\}$，角向矩消失再通过Newton恒等式固定正多边形。这个机制的认证域使用完整二阶必要条件，不能预设对称性。对 $n\ge10$，中心正多边形已被严格排除，因而不能把同一等号结构机械推广；后续需针对新的多层候选设计矩核及必要域，或转向一般粒子分布与稳定性结构。
