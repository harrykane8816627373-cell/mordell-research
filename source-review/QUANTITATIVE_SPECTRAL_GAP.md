# 将原稿的有限谱层级障碍量化

研究资料与项目发起：蒋丰灿。2026-10-27 的来源核对及推导包含 OpenAI Codex 协作。

此处是在用户提供的完整修订稿已有定理上的补充。原稿 `main.tex` 第3122–3535行已经证明乘法算子谱层级收敛，以及对 n≥3 的每个有限层均严格高于 Mₙ。本节不将这些既有成果视为新发现；补充的是把“严格”转换成可计算的有理正缺口。

## 1. 一般齐次多项式的行能量界

沿用原稿的 Bombieri 归一化。设非零齐次 holomorphic 多项式 f 在 Cᵈ 上的次数为 N≥1，记

```math
S=\max_{\|u\|=1}|f(u)|^2,\qquad
\Lambda_0=\|f\|_{B,N}^2,
\qquad C_{N,r}=\binom{N+r}{r}.
```

令 Λᵣ 是原稿第3153–3163行定义的乘法矩阵 Aᵣ 的最大奇异值平方。则对每个整数 r≥0，

```math
\boxed{\displaystyle
\Lambda_r^2\ge S^2+
\frac{S(\Lambda_0-S)}{C_{N,r}}.}
```

证明：取一个达到 S 的单位点，经酉坐标变换送到 e₁。Bombieri 范数、S 和谱层级不因该变换改变。在新坐标写

```math
f(u)=\sum_{|\gamma|=N}a_\gamma u^\gamma,
\qquad \gamma_0=(N,0,\ldots,0),\qquad |a_{\gamma_0}|^2=S.
```

原稿已有 Kᵣ=AᵣAᵣ*≥0。令 e₀ 表示输出单项式 u₁⁽ᴺ⁺ʳ⁾ 的归一化坐标向量。该行只有输入单项式 u₁ʳ 对应的一列非零，故

```math
\|K_re_0\|^2
=S^2+S\sum_{\gamma\ne\gamma_0}
\frac{|a_\gamma|^2}{\binom{N+r}{\gamma+r e_1}}.
```

这里多指标组合数取 (N+r)! 除以各分量阶乘。组合数之比满足

```math
\frac{\binom{N+r}{\gamma+r e_1}}{\binom N\gamma}
=\frac{\binom{N+r}{r}}{\binom{\gamma_1+r}{r}}
\le C_{N,r}.
```

因此该行的非对角能量至少为

```math
\frac{S}{C_{N,r}}
\sum_{\gamma\ne\gamma_0}\frac{|a_\gamma|^2}{\binom N\gamma}
=\frac{S(\Lambda_0-S)}{C_{N,r}}.
```

又 ‖Kᵣ‖=Λᵣ、‖e₀‖=1，所以 Λᵣ²≥‖Kᵣe₀‖²，得证。没有假设最大点非退化，也没有要求形状集合有限。

## 2. 与已有 Mₙ 上下界联合的有理正缺口

令 N=n(n−1)/2，原稿定义给出

```math
S=\frac{M_n}{n^N},\qquad
P_n:=n^N\Lambda_0
=\frac{n^N n!\prod_{j=0}^{n-1}j!}{N!}.
```

Pₙ是谱层级的第零层上界，此处记号独立于旧稿中阶乘连乘 Pₙ，避免跨文档混用。若已有证明给出正有理数 L、U，满足

```math
0\lt L\le M_n\le U\lt P_n,
```

则每个有限 r 的纯 holomorphic 谱上界有显式正缺口

```math
\boxed{\displaystyle
n^N\Lambda_r-M_n
\ge\frac{L(P_n-U)}{2P_n\binom{N+r}{r}}\gt0.}
```

证明：原稿层级单调性给 Λᵣ≤Λ₀，且 S≤Λ₀。因而

```math
\Lambda_r-S
=\frac{\Lambda_r^2-S^2}{\Lambda_r+S}
\ge\frac{S(\Lambda_0-S)}{2\Lambda_0 C_{N,r}}.
```

分别使用 S≥L/nᴺ、Λ₀−S≥(Pₙ−U)/nᴺ并乘 nᴺ，即得所述有理数。上界 U 的严格条件不可省略；程序计算表达式不会自动证明用户给出的 L、U 是有效全局界。

例如利用已知 M₃=27 和 P₃=54，任意有限层都有

```math
27\Lambda_r-M_3\ge\frac{27}{4\binom{3+r}{r}}.
```

原稿实际算出的第1、2层上界为135/4、297/10，其正缺口均大于本节的下界。这个量化界不声称是该层误差的精确渐近式。

若希望某层谱上界达到相对误差 ε，则还必须有

```math
\binom{N+r}{r}\ge\frac{L(P_n-U)}{2P_n\varepsilon U}.
```

这是必要条件而非充分条件。它可用于排除不够的谱层数，不能充当 Mₙ 的改进上界。

## 3. 对现有方法的含义和范围

这说明纯 holomorphic 平方证书在任何有限层都不能直接达到精确 Mₙ；应与原稿的实混合多项式证书、物理迹等价模型或关键值增广方法衔接。原稿第4647–4711行已经给出关键值增广后有限精确实 SOS 存在性，不能将本节理解为所有 SOS 都有同一障碍。

`verify_quantitative_spectral_gap.py`用整数、Fraction检查多指标组合数比较、原稿已执行的三点低层谱值，以及 f=uᴺ+vᴺ 的有理二维主块。有限核验只检查实现和常数，普遍结论由上面的证明承担。该结果量化方法自身的误差，并未求出新的未知 Mₙ，也不声称改善目前全局最强上界。
