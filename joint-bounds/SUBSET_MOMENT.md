# 联合子集构型限制的最佳两矩对数上包

研究资料与项目发起：蒋丰灿（Fengcan Jiang）；本补充推导与程序含 OpenAI Codex 协作。2026-10-28。

本节进一步加强已有的子集递推：保留均值、方差下界及构型所给的支撑上界，直接解出这些信息下的对数平均最大值。它严格改善先前的统一二次 Jensen 损失。**在这里实际比较的 147 个参数对中，新递推仍高于既有圆盘界；这部分结果不刷新无条件全局上包络。** 文中的“最佳”仅针对指定的随机变量信息，不针对真实平面构型或 Mordell 问题。

## 1. 给定均值、方差下界和支撑上端的精确上包

设正随机变量满足

```math
X>0,\qquad \mathbb E X=\mu,\qquad X\le L,\qquad
\operatorname{Var}(X)\ge v,\qquad 0<v<\mu(L-\mu).
```

这些假设保证 \(0<\mu<L\)。定义

```math
a=\mu-\frac{v}{L-\mu}>0,\qquad
p=\frac{(L-\mu)^2}{(L-\mu)^2+v}.
```

则

```math
\boxed{\quad
\mathbb E\log X\le p\log a+(1-p)\log L.
\quad}
```

这是上述信息下的精确最大值：概率 \(p\) 取 \(a\)、概率 \(1-p\) 取 \(L\) 的两点分布满足均值 \(\mu\)、方差 \(v\)，并取等。对于一般正分布，等号要求恰有该两点支撑和方差 \(v\)。两点分布的可达性只是在随机变量类内，不保证它来自平面点构型的子集矩。

证明。取二次多项式

```math
q(x)=\log a+\frac{x-a}{a}+c(x-a)^2,\qquad
c=\frac{\log(L/a)-(L-a)/a}{(L-a)^2}<0.
```

它在 \(a\) 与 \(\log x\) 同值同导数，在 \(L\) 同值。Hermite 插值余项给出，对于 \(0<x\le L\)，

```math
\log x-q(x)
=\frac{(\log)'''(\xi)}{3!}(x-a)^2(x-L)\le0,
\qquad (\log)'''(\xi)=2/\xi^3>0,
```

其中 \(\xi\) 在所涉及的正节点之间。重复节点按连续延拓理解；\(x=a,L\) 时误差为零，其余 \(x\in(0,L)\) 时严格为负。

因为 \(c<0\)，

```math
\mathbb E q(X)
=\log a+\frac{\mu-a}{a}
+c\bigl(\operatorname{Var}(X)+(\mu-a)^2\bigr)
\le\log a+\frac{\mu-a}{a}+c\bigl(v+(\mu-a)^2\bigr).
```

右端只用前两矩，等于上述两点分布下 \(q\) 的期望，因而等于 \(p\log a+(1-p)\log L\)。若 \(\mathbb E\log X=-\infty\)，不等式自动成立；其余情形直接取期望即可。

当 \(v=0\) 或 \(L=\mu\) 时应使用普通上界 \(\log\mu\)，不直接代入此处的严格参数公式。若给定正变量的方差下界达到或超过 \(\mu(L-\mu)\)，这些信息本身就不相容，因为 \(X^2<LX\) 在 \(0<X<L\) 上严格成立，且 \(\mu<L\)。实现拒绝不满足严格参数范围的输入。

## 2. 严格强于原有最佳统一二次损失

此前[最优统一二次系数](../source-review/DELETION_REFINEMENT.md)为

```math
h(T)=\frac{T-1-\log T}{(T-1)^2},\quad T=L/\mu>1,
\qquad
\mathbb E\log X\le\log\mu-\frac{h(T)}{\mu^2}\operatorname{Var}(X).
```

本节两点上包严格满足

```math
p\log a+(1-p)\log L
<\log\mu-\frac{h(L/\mu)}{\mu^2}v.
```

可把旧不等式用于本节两点分布。归一化后的节点 \(a/\mu<1<T\)；旧切线二次式只在 \(1,T\) 处接触 \(\log x\)，而在 \(a/\mu\) 处严格大于它，因为 \(h\) 严格递减。节点 \(a\) 的概率为正，故平均后严格。这没有否定旧系数的最佳性：旧结果要求单一方差系数对所有方差都有效；本节使用指定的正方差下界，允许上包随它非线性变化。

## 3. 代入一般 n 的已有子集方法

仍令

```math
M_n=\max_{\sum|z_i|^2=n}\prod_{i<j}|z_i-z_j|^2,
\qquad N_j=j(j-1)/2.
```

取全局最大构型，均匀选择一个 \(m\) 点真子集 \(S\)，并令 \(Q_S\) 为该子集相对其自身重心的内部二阶矩。对于 \(n\ge4\)、\(2\le m\le n-2\)，[已有子集方差恒等式及递推](../global-n/RECURRENCE_PROOF.md)给出

```math
\mu=\mathbb E Q_S=\frac{n(m-1)}{n-1},\qquad
\operatorname{Var}(Q_S)\ge v_{n,m}
=\frac{n(m-1)(n-m)(n-m-1)}{m(n-1)^2(n-2)}>0,
```

```math
\log M_n\le\frac{N_n}{N_m}\log M_m
+N_n\bigl(\mathbb E\log Q_S-\log m\bigr).
```

设 \(h=n-m\)，并使用[先前严格子集有理余量](../finite-n/SUBSET_PROOF.md)中的 \(\varepsilon_{n,j}\)，则构型必要条件给出

```math
Q_S\le L^+_{n,m}:=
\min\left\{n,\frac{2m(n-1-\varepsilon_{n,n-1})}{n},
n-\frac{h(h-1)+2\varepsilon_{n,h}}{n-1}\right\}.
```

第一项来自总矩，第二项来自每个点的半径上界，第三项来自补集的内部矩下界和重心分解。这里直接使用旧约束，不登记为新证明。闭合上包允许把旧严格不等号放宽为非严格。

代入 \(L=L^+_{n,m}\)、\(v=v_{n,m}\)，并按第 1 节定义 \(a,p\)。只要核对 \(0<v<\mu(L-\mu)\)，就有

```math
\boxed{\quad
\log M_n\le\frac{N_n}{N_m}\log M_m
+N_n\left[p\log a+(1-p)\log L-\log m\right].
\quad}
```

对于合法的全局最大构型，该严格参数条件实际上由 \(0<Q_S\le L\)、\(\mathbb E Q_S=\mu<L\)、\(\operatorname{Var}(Q_S)\ge v>0\) 保证；程序也直接作有理数检查。实际方差若已由更强的四阶、角向矩约束从下方控制，也可用更大的有效 \(v\) 重新应用本节，而不必退回统一方差系数。

此结论是全局最大构型上的子集递推。不把全局构型半径界推广到任意驻点或普通构型；\(m=n-1\) 的普适方差下界为零，也不在正方差比较范围内。若已另证该分支的正径向方差，则一般随机变量定理仍可使用。

## 4. 实际比较范围与限制

[实现](subset_moment.py)使用 Python 标准库 `Fraction`；所有对数由带显式尾项的 atanh 级数包围，并按系数符号向外组合。有限比较使用 \(M_2=4\)、\(M_3=27\)、\(M_4=256\)、\(M_5=3125\) 及仓库给出的六七点锐值。六七点结论保留原有 Faugère–Svartz 完整物理驻点分类的外部输入限制，见[覆盖审计](../earlier-results/six_seven_coverage_audit.md)。

实际检查 \(4\le n\le30\)、\(2\le m\le\min(7,n-2)\)，共 147 组；各组均使用相同的 \(L^+_{n,m}\) 比较新旧公式，并与原圆盘上界 \(U_n\) 对照。147 组全部严格改善旧递推公式，147 组也全部仍高于 \(\log U_n\)。例如 \(n=4,m=2\) 时，

```math
\mu=\frac43,\quad L=\frac{629}{210},\quad v=\frac29,
\quad a=\frac{1256}{1047},\quad p=\frac{121801}{131601}.
```

这组的对数上界较旧递推降低约 0.014232，仍比 \(\log U_4\) 高约 0.002174。所有接受判断使用有理区间，近似数仅作阅读提示。因此这部分贡献是更强的跨点数兼容条件，并为未来更强构型方差约束保留了接口，没有把公式改善误报为已刷新当前最佳整体上界。

[独立重放入口](verify_subset_moment.py)检查对数包围、两点极值律、Hermite 多项式的符号与实际有理子集矩分布；有限检查并不代替上述通用证明。运行

```bash
python -B -S -O verify_subset_moment.py
python -B -S -O verify_subset_moment.py --write-report
```

第二条命令重新生成 [147 组比较数据](subset_moment_comparison.json)。报告中的直接枚举仅涉及明示的 4、5、6、8 点测试族，不是一般 n 的优化域穷尽，也不包含全部真实驻点的枚举。
