# 一般 n 的线性渐近项、圆盘上界与八点缺额

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；本文整理包含 OpenAI Codex 协作。日期：2026-10-27。

阅读[六页英文短稿PDF](proof.pdf)及[TeX源码](proof.tex)。

本页把用户提出的库仑能量方向补成明确的定理接口、完整归一化和可复核常数。线性项存在性使用 Sandier–Serfaty 的既有定理；三角晶格值使用 Bétermin–Sandier；更强上界与经典圆盘平滑界一致，不作普适常数界的原创优先权声明。八点锐界仍未证明。

## 结论

所有对数为自然对数。令 Δₙ 为所有无序点对平方距离的乘积，Mₙ 为在 ∑|zᵢ|² = n 下 Δₙ 的最大值，N = n(n−1)/2。则存在与 n 无关的常数 γ，使

```math
\log M_n=\left(\frac12\log2-\frac14\right)n^2+\frac n2\log n+\gamma n+o(n).
```

在下文明确的单位背景重整化能量规范中，

```math
\gamma=-\frac1\pi\min_{\mathcal A_1}W-\frac12\log(2\pi).
```

更窄的严格界为

```math
\gamma_\triangle\le\gamma\le\gamma_D:=\frac34-\frac12\log2,
\qquad
\gamma_\triangle=3\log\Gamma(1/3)+\frac12\log3-2\log2-\frac32\log\pi.
```

[纯有理区间证书](gamma_bounds_certificate.json)给出严格外包：

```math
0.402178895223365173209539459449
<\gamma<
0.403426409720027345291383939271.
```

这里的十进数为严格包围端点；它们不是 γ 的已知精确值。达到左侧精确三角候选值要求全体允许微观配置中的全局三角极小性，不能用仅在晶格族内最优替代。本项目未独立建立该全局结论。

对每个有限 n ≥ 2，另有不依赖渐近定理的全局上界

```math
M_n\le U_n:=n^{n/2}\left(\frac{2n}{n-1}\right)^N e^{-N/2}.
```

因此全 n 还满足

```math
\log M_n<\left(\frac12\log2-\frac14\right)n^2+\frac n2\log n+
\left(\frac34-\frac12\log2\right)n-\frac14.
```

最后的 −1/4 是这个显式上界的修正，不是已经确定的真实有限点数边界展开项。

## 1. 全局问题与库仑能量的精确等价

使用有序点对的 Hamiltonian，而 Δₙ 仍使用无序点对的平方距离：

```math
H_n(x)=-\sum_{i\ne j}\log|x_i-x_j|+\frac n2\sum_i|x_i|^2,
\qquad E_n=\min_{x\in\mathbb C^n}H_n(x).
```

即 Hₙ = −log Δₙ + (n/2)Q。碰撞处能量趋于无穷，二次势在无穷远控制对数项，因此存在无碰撞极小点。

任意无碰撞形状均可写成 x = sz，s > 0，且 Q(z) = n。由 Δₙ 的 2N 次齐次性，

```math
H_n(sz)=-\log\Delta_n(z)-2N\log s+\frac{n^2s^2}{2}.
```

对 s 的唯一全局最小点是 s² = (n−1)/n。先对尺度、再对形状取极小，得到完整全局恒等式

```math
\log M_n=A_n-E_n,
\qquad A_n=N\left[1-\log\left(1-\frac1n\right)\right].
```

所有库仑极小点恰满足 Q(x) = n−1；它与本问题 Q(z) = n 的差不能省略。不存在先验固定形状或仅比较某些驻点的假设。

桥接项有收敛展开

```math
A_n-\frac{n^2}{2}
=-\frac12\sum_{k=0}^{\infty}\frac{n^{-k}}{(k+1)(k+2)}
=-\frac14-\frac1{12n}+O(n^{-2}).
```

特别地，所有 n ≥ 2 都有

```math
-\frac14-\frac1{12(n-1)}\le A_n-\frac{n^2}{2}
\le-\frac14-\frac1{12n}.
```

所以归一化桥接不引入额外线性项。由于 Eₙ 的余项只知道是 o(n)，不能由这里的常数 −1/4 推断 log Mₙ 的真实常数项。

## 2. 连续平衡测度及外部定理接口

设 g(x) = −log|x|，D(μ) 为 g 的双重积分，

```math
I_V(\mu)=D(\mu)+\int V\,d\mu,\qquad V(x)=\frac{|x|^2}{2}.
```

唯一平衡测度是半径 √2 圆盘上的常密度概率测度

```math
\mu_0=\frac1{2\pi}\mathbf1_{|x|\le\sqrt2}\,dx,
\qquad I_0=I_V(\mu_0)=\frac34-\frac12\log2.
```

直接核验：对半径 R 的均匀圆盘，其对数势在盘内是 log R − 1/2 + |x|²/(2R²)，盘外是 log|x|。取 R = √2，g 势加 V/2 在盘内为常数，盘外不小于该常数，给出平衡条件；圆盘二阶矩为 1，故得到 I₀。第5节的能量分解和零质量 Dirichlet 表示也直接证明 I_V(μₓ) ≥ I₀，所以有限 n 的圆盘上界不需要借用外部渐近定理。

Sandier–Serfaty [arXiv:1201.3503v2](https://arxiv.org/abs/1201.3503v2) Theorem 2(C)、式(1.40)给出最低能量的展开。其假设(1.11)–(1.13)在本例成立：二次势光滑并具有充分增长，平衡支持为光滑圆盘，内部密度为严格正常数。无需交换温度趋零与 n 趋于无穷的两个极限。

```math
E_n=n^2I_0-\frac n2\log n+
\left[\frac1\pi\min_{\mathcal A_1}W-\frac12\int m_0\log m_0\right]n+o(n).
```

这里 A₁ 取原文完整允许类：场满足 div E = 2π(ν−1)、curl E = 0，其中 ν 是离散点集的点质量和，并满足原文的有界平均点密度条件。它不限于周期场或 Bravais 晶格。W 使用原文的半 Dirichlet 能加 π log η 截断，并按大方形面积归一化取上极限的约定。由于 m₀ = 1/(2π)，括号内等于 W*/π + log(2π)/2。与第1节全局恒等式合并，即得所述 γ。

同一接口也明确列在 Bétermin–Sandier [arXiv:1404.4485v5](https://arxiv.org/abs/1404.4485v5) Theorem 1.1、式(1.4)。本项目应用这些既有定理，不将它们改署为本项目原创证明。

## 3. 有限 n 圆盘平滑证明

令 σᵣ 为半径 r 的均匀圆盘概率测度，Tᵣf = f*σᵣ。角向均值恒等式给出，d > 0 时

```math
(T_rg)(d)=\begin{cases}
-\log d,&d\ge r,\\
-\log r+\frac12-\frac{d^2}{2r^2},&0<d<r.
\end{cases}
```

因此 Tᵣg ≤ g。Tᵣ 是正平均算子，再平均一次仍有 Tᵣ²g ≤ Tᵣg ≤ g。这个比较对重叠圆盘同样成立，不需要点对分离条件。

积分上式并直接计算二阶矩，得到

```math
D(\sigma_r)=-\log r+\frac14,
\qquad \int |u|^2\,d\sigma_r(u)=\frac{r^2}{2}.
```

对任意互异点 x，定义概率测度 μₓ = (1/n)∑σᵣ(·−xᵢ)。比较点对能量与平滑点对能量，并准确扣除 n 个自能项：

```math
H_n(x)\ge n^2I_V(\mu_x)-nD(\sigma_r)-\frac{n^2r^2}{4}
\ge n^2I_0+n\log r-\frac n4-\frac{n^2r^2}{4}.
```

对 r 优化得到 r² = 2/n，故

```math
E_n\ge L_n:=n^2I_0-\frac n2\log n+n\left(\frac12\log2-\frac34\right).
```

代入第1节恒等式，整理常数及 log n、log(n−1)、log 2 的系数，正好得到 Mₙ ≤ Uₙ。结合线性渐近项的存在性，得到 γ ≤ 3/4 − log 2/2。

该强上端与经典 Lieb–Narnhofer／Sari–Merlini 圆盘涂抹界一致。准确的现代来源是 A. B. Lauritsen, *Floating Wigner crystal and periodic jellium configurations*, J. Math. Phys. 62, 083305 (2021), [DOI](https://doi.org/10.1063/5.0053494)、[arXiv:2103.07975v2](https://arxiv.org/abs/2103.07975v2)，Appendix B Proposition B.1、Corollary 4.4 与 Remark 4.5。本文提供其在有限矩约束问题中的独立推导，不声称这个普适常数界是新发现。

## 4. 三角晶格下端与严格十进制证书

Bétermin–Sandier §7.3 计算单位密度三角晶格的 W 值：

```math
\frac1\pi W(\Lambda_\triangle)=\log\pi-\frac12\log3
-3\log\Gamma(1/3)+\frac32\log2.
```

W* ≤ W(Λ△)，而 γ 对 W* 有负号，所以上式给 γ 的下界，方向不能颠倒。精确 γ 尚未由本文建立；仅在 Bravais 晶格中最优不够。

[verify_gamma_bounds.py](verify_gamma_bounds.py)不用浮点验收：Machin 公式包围 π，整数平方根包围 √3，正项 Taylor 级数包围指数，atanh 级数包围对数；Dedekind η 乘积取16项并显式控制尾项。它包围的是已有解析公式的数值，不证明重整化能量定理或晶格全局最优。

## 5. 与八点剩余难点连接：精确非负缺额

取 r² = 2/n，Kᵣ = Tᵣ²g。令 νₓ = μₓ−μ₀，定义

```math
\zeta(\rho)=\begin{cases}
0,&\rho\le\sqrt2,\\
\rho^2/4-\log\rho-1/2+(\log2)/2,&\rho>\sqrt2.
\end{cases}
```

展开连续能量即得精确恒等式

```math
H_n(x)-L_n=
 n^2D(\nu_x)+2n^2\int\zeta\,d\mu_x+
 \sum_{i\ne j}\bigl[g(x_i-x_j)-K_r(x_i-x_j)\bigr].
```

三项都非负。第一项由零质量电荷的 Dirichlet 表示 D(ν) = (1/(2π))∫|∇(g*ν)|² 得到；这里 ν 紧支撑、有界密度且总质量为零，分部积分的无穷远边界项消失。第二项是盘外溢出代价；第三项是圆盘重叠导致的点对代价，距离至少 2r 时恰为零。

还能使用一个显式点对下界，令 a₊ = max(a,0)：

```math
g(d)-K_r(d)\ge(1-d/r)_+^2+\frac1{16}(1-d/(2r))_+^4.
```

证明：hᵣ = g−Tᵣg 在 d<r 时为 log(r/d)−1/2+d²/(2r²)，且 hᵣ ≥ (1−d/r)₊²。又 g−Kᵣ = hᵣ+Tᵣhᵣ；在两个半径 r 圆盘的交叠区内，取以线段中点为圆心、半径 r(1−d/(2r))/2 的小圆盘，面积比与 hᵣ 的统一下界相乘即得四次项。无需预设构型对称性。

当 Q(x) = n−1、z = √(n/(n−1))x 时，该非负缺额恰为 log Uₙ − log Δₙ(z)。因此八点锐界等价于要求所有 Q(x) = 7 的配置满足

```math
H_8(x)-L_8\ge\delta_8,
\qquad
\delta_8=\log(U_8/C_8)=40\log2-7\log7-14,
\quad C_8=8^{28}/7^{21}.
```

这里 r = 1/2。这给出一个同时容纳一般非零矩分支的连续非负量，但还没有证明所需的正下界 δ₈，不能将等价改写当作主问题已证。

[finite_bound_certificate.json](finite_bound_certificate.json)严格包围 δ₈ 及 U₈/C₈；后者约 1.1101733552。当前得到的是无条件非锐上界，不宣称它优于本项目每一个既有专门类别的估计。

## 6. 复核与出版范围

```bash
python -B -S verify_all.py
```

统一入口先验证源文件哈希，再在临时副本运行两个标准库接受程序。解析平滑比较、外部渐近定理及其假设仍须阅读本文；程序通过不等于全篇已在证明助理中形式化。

未给出 γ 的本项目精确值；未给出真实有限 n 边界修正或 o(n) 的有效统一速率；八点锐界及一般非零矩系统仍开放。本文不是外部同行评审论文，也不作全球优先权声明。
