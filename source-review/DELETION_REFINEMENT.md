# 保留原删点定理，并加强方差损失常数

研究资料与项目发起：蒋丰灿。2026-10-27。本节以此次核对的用户修订源码为基础；原有严格跨点数单调性及删点方差罚项归属于用户原稿。本节补充对数平均的最优统一二次系数。

原稿位置：`sources/revised/anc/all_n_deletion_radius.tex` 第15–31行及其完整证明，原有径向方差罚项位于第87–98行。原稿已经证明 d₂=d₃=2、n≥4时 dₙ<dₙ₋₁；不应把相同结论重新标成本轮新成果。

## 1. 最佳对数二次损失

定义连续函数

```math
h(x)=\frac{x-1-\log x}{(x-1)^2}\quad(x\ne1),
\qquad h(1)=\frac12.
```

h 在正半轴严格递减。事实上，x≠1时

```math
h'(x)=\frac{-x+x^{-1}+2\log x}{(x-1)^3}.
```

分子在1处为0，其导数为 −(x−1)²/x²。因此分子在x<1时为正、x>1时为负，恰使上述导数在两侧均严格为负；在1处按连续延拓处理即可。

设正随机变量 X 满足 E X=μ、X≤L，令 T=L/μ≥1。由 h(X/μ)≥h(T) 得

```math
\log(X/\mu)\le X/\mu-1-h(T)(X/\mu-1)^2.
```

取平均，得到

```math
\boxed{\displaystyle
\mathbb E\log X\le\log\mu-
\frac{h(T)}{\mu^2}\operatorname{Var}(X).}
```

T=1时 X=μ 几乎处处，方差为0，结论仍成立。

当 T>1 时，这个系数是给定上界 T 时最佳的统一二次切线系数：在 x=T 的点态不等式中已经取等，任何更大的系数都会失败。它也是仅固定均值和上界时最佳的统一方差系数。为看后者，取概率 p 在 T，概率1−p在 x_p=(1−pT)/(1−p)>0，令 p↓0。该两点分布均值为1，而其对数 Jensen 损失除以方差趋于

```math
\frac{T-1-\log T}{(T-1)^2}=h(T).
```

这里最佳性只针对这一统一均值、方差及上界不等式，不声称对子集离散分布或原极值问题已经最优。

当 T>1 时，还有严格比较

```math
h(T)\gt\frac1{2T}.
```

证明：令 g(T)=T−1−logT−(T−1)²/(2T)，则 g(1)=0，g′(T)=(T−1)²/(2T²)>0。故新式严格加强先前使用的 Var(X)/(2μL) 损失，只要方差正。

## 2. 一般子集递推的加强

记 N_j=j(j−1)/2。对中心化的全局 n 点最大构型，均匀随机 m 点真子集 S 的内部二阶矩为 Q_S。旧的子集乘积恒等式给

```math
\log M_n\le\frac{N_n}{N_m}\log M_m
 +N_n\bigl(\mathbb E\log Q_S-\log m\bigr),
\qquad\mu=\mathbb E Q_S=\frac{n(m-1)}{n-1}.
```

设 Q_S≤L，T=L/μ。将第1节直接代入得

```math
\boxed{\displaystyle
\log M_n\le\frac{N_n}{N_m}\log M_m
 +N_n\log\frac{n(m-1)}{m(n-1)}
 -\frac{N_n h(L/\mu)}{\mu^2}\operatorname{Var}(Q_S).}
```

此前完整平面方差恒等式可以原样代入；特别地，n≥4、2≤m≤n−2时，其普适下界为

```math
v_{n,m}=\frac{n(m-1)(n-m)(n-m-1)}{m(n-1)^2(n-2)}.
```

因此把上述式中的 Var(Q_S)替成 v_{n,m}仍给全局上界。旧的严格子集有理余量可以减小 L，并且 h 递减，故会进一步增大折扣。这是真正联合两批方法，而不是重复普通 Jensen 递推。

## 3. 删去一个点时的显式径向系数

设 n≥3，t_i=|z_i|²、R₄=Σt_i²。删除 i 后内部二阶矩为 Q_i=n−nt_i/(n−1)，故

```math
\mu=\frac{n(n-2)}{n-1},\qquad
\operatorname{Var}(Q_i)=\frac{n}{(n-1)^2}(R_4-n),
\qquad Q_i\le n.
```

取 L=n，就得到

```math
\boxed{\begin{aligned}
\log M_n\le{}&\frac n{n-2}\log M_{n-1}
 +N_n\log\frac{n(n-2)}{(n-1)^2}\\
&-c_n(R_4-n),\\
c_n={}&\frac{n-1}{2}
\left[\frac1{n-2}-\log\frac{n-1}{n-2}\right].
\end{aligned}}
```

常数严格满足

```math
\boxed{\displaystyle
c_n\gt\frac1{4(n-2)}\gt\frac1{4(n-1)}.}
```

最后一项是用户原稿的系数；中间项来自先前 Var/(2μL) 估计；本节给出更大的完整对数系数。若 R₄=n，这些径向罚项都为0，本节不将其误报为所有构型上的正统一亏损。普适正子集方差仍需 m≤n−2。

## 4. 保留纯有理接受证书

T 是正有理数时，可用有显式尾项的 atanh 级数包围 logT。若 T≥1，令 y=(T−1)/(T+1)，则

```math
2\sum_{j=0}^{r-1}\frac{y^{2j+1}}{2j+1}
\le\log T\le
2\sum_{j=0}^{r-1}\frac{y^{2j+1}}{2j+1}
 +\frac{2y^{2r+1}}{(2r+1)(1-y^2)}.
```

用 logT 的有理上端代入 h(T)即可得到折扣系数的有理下端；或者直接使用严格较弱但纯有理的 1/(2T)。同样，c_n 中用 log((n−1)/(n−2))的有理上端便得 c_n 的安全有理下端。所有不等式方向都保留了所需的上界证明。

本节没有由径向必要条件反推出全局最优性；没有得到一般 n 的闭式锐值或八点无条件锐界。
