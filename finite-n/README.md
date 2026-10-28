# 严格子集方法的有效化与每个有限 n 的下界

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；推导、程序与内部复核含 OpenAI Codex 协作。日期：2026-10-27。

**当前包络入口更新（2026-10-28）：** 本页保留该批三角格、矩及子集余量证书。完整原稿已有的全n有理多环及n6..12保存整数构型现已恢复；相位和中心化改善后的完整比较见[source-review](../source-review/README.md)。进一步的构型和上下界见[configuration-bounds](../configuration-bounds/README.md)。本页旧ZIP是2026-10-27历史快照，当前网页状态说明及source_manifest中的相应字节索引单独更新；旧数据与ZIP内部清单保持原字节。

本专题补全严格性细节，给出全局最优构型的显式有理余量，并把有限点数下界整理为可计算的精确证书。既有严格子集定理、系数矩和矩比收敛性已经在[此前稿件](../earlier-results/general_n_notes.tex)中；本轮增加的是有效余量、三角整数格构造、有理 Rayleigh 改善及可执行复核。不作全球优先权声明。

## 1. 归一化与本轮结论

```math
M_n=\max_{\sum_i|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i\lt j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2.
```

对任意全局最大构型、任意真子集 S，只有在 2 ≤ m = |S| < n 时断言严格不等式。单点子集两边均为零；全体点时为等式。这里没有将全局结论推广到所有驻点或局部极大点。

```math
\bar z_S=\frac1m\sum_{i\in S}z_i,\qquad
Q_S=\sum_{i\in S}|z_i-\bar z_S|^2
\gt \frac{m(m-1)+2\varepsilon_{n,m}}{n-1},
\quad
\varepsilon_{n,m}=\frac1{2(1+B)\sum_{j=0}^{m(n-m)-1}B^j},
\quad B=\frac{3n-m-1}{m-1}.
```

全部余量为正有理数。[严格性与有效余量的完整证明](SUBSET_PROOF.md)使用复缩放多项式、初等 Fourier 恒等式和一个有界零点，没有假定构型对称性。八点七元子集给出：

```math
\varepsilon_{8,7}=\frac{2187}{9217846},\qquad
Q_S\gt 6+\frac{2187}{32262461},\qquad
|z_i|^2\lt \frac74-\frac{2187}{36871384}.
```

这些是排除候选构型的必要约束；不能反向当作 M_n 的下界。下界由下面的构型或精确平均独立建立。

## 2. 每个有限 n 的显式有理下界

对所有整数 n ≥ 2，可以无条件使用以下合并下界。R₂ 定义为 4；其他符号在下文或所链接的完整证明中给出：

```math
\boxed{\displaystyle
M_n\ge L_n:=\max\{n^n,\ C_n,\ G_n,\ S_n,\ T_n,\ R_n\}.}
```

其中单位正多边形、原点加正多边形和中心化高斯平均分别给出：

```math
C_n=\frac{n^N}{(n-1)^{(n-1)(n-2)/2}},\qquad
G_n=\frac{n!(n-2)!\,n^N\prod_{j=0}^{n-1}j!}{(N+n-2)!}.
```

C₂=2 也有效，但被其他项严格改善。G_n 的 Gamma 形状参数是 n−1。两式的[既有完整证明](../earlier-results/general_n_proof.md)继续保留；下界不要求全局最优点具有这些构型。

### 三角格：每个 n 均可直接计算的正有理数

令 ω = exp(iπ/3)，使用整数坐标 v = a+bω：

```math
q(a,b)=|a+b\omega|^2=a^2+ab+b^2.
```

将所有整数对按 (q(a,b),a,b) 从小到大排序，取前 n 个 v₁,…,vₙ。这个规则包含平局处理，不需要浮点数或搜索最优参数。定义：

```math
D_n=n\sum_iq(v_i)-q\left(\sum_i v_i\right)
=\sum_{i\lt j}q(v_i-v_j)\gt 0,
```

```math
\boxed{\displaystyle
T_n=\left(\frac{n^2}{D_n}\right)^N
       \prod_{i\lt j}q(v_i-v_j)\in\mathbb Q_{\gt 0}.}
```

证明：取 zᵢ = n(vᵢ−平均值)/√Dₙ。二次型恒等式给出总二阶矩 n，每个距离平方乘以 n²/Dₙ，因而距离积恰为 Tₙ。取同样规则、改用 q(a,b)=a²+b²，就得到已有的高斯整数格下界 Sₙ。

有限枚举也是严格的：矩形 [−n,n]² 外有 q ≥ (n+1)²/2，而内部已有 n 个横轴点满足 q ≤ floor(n/2)²，后者严格更小。因此代码在这个有限矩形内排序，确实取得整个无限格子的前 n 点。

[精确构型证书](triangular_lattice_values.json)列出 n=2,…,30 的坐标、D、距离平方乘积和完整有理数。恢复 T₂=4、T₃=27、T₇=7²¹/6¹⁵；没有据此宣称其他 n 的全局锐值。对已有下界包络（含高斯平均、两个多边形、方格及已有 n≥10 显式扰动）作精确比较，三角格在 n=13–19、23–27、30 严格改善。[比较表](triangular_comparison.md)中的十进数仅供阅读，决定结果的是整数比较；没有断言三角格对每个 n 都优于方格。

### 有理加权矩：严格改善已有矩比

把 Vₙ 的整数系数幂完全合并同类项，定义球面矩：

```math
V_n(z)^k=\sum_\alpha c_{k,\alpha}z^\alpha,\quad
J_{n,k}=\sum_\alpha c_{k,\alpha}^2\prod_i\alpha_i!,\quad
\mu_k=\frac{n^{kN}(n-2)!J_{n,k}}{(kN+n-2)!}.
```

对 n ≥ 3，令：

```math
A=\mu_1\mu_3-\mu_2^2\gt 0,\quad
B=\mu_3-\mu_1\mu_2\gt 0,\quad a=A/B,
```

```math
\boxed{\displaystyle
R_n=\frac{\mu_3-2a\mu_2+a^2\mu_1}
            {\mu_2-2a\mu_1+a^2}
\gt \frac{\mu_3}{\mu_2}\gt \frac{\mu_2}{\mu_1}\gt \mu_1.}
```

Rₙ 是以 (Δₙ−a)² 为非负权重的精确平均，因此不超过 Mₙ。这一公式给每个有限 n 一个有限整数和所定义的正有理下界。[完整矩方法证明](MOMENT_METHOD.md)还给出 Hankel 矩阵推广，以及只用整数比较即可结束的相对误差证书。高阶系数的成本会迅速增长，不能把有限公式理解成廉价算法。

程序默认只计算快速构型和平均下界；指定 `--moment-order 3` 才加入 Rₙ。本专题实际完整计算并公布 n=2,…,7 的前三阶矩，未宣称已算出所有 n 的高阶矩。

## 3. 复核与使用

仅需 Python 标准库：

```bash
python -B -S -O verify_all.py
python -B -S finite_lower_bounds.py 8 13 19 30
python -B -S finite_lower_bounds.py 3 5 7 --moment-order 3
```

第一项独立展开两两差因子，核对 n≤5 的前三阶系数范数；重新计算已公开 n≤7 的矩、n≤30 的两种格子和全部比较；核对严格子集有理常数及整数停止准则。通用定理由证明文本承担，有限样本不替代全 n 证明。中间输出不是浮点极值搜索。

快速构型包络也保留已有十点反例及 n≥11 扰动的有理下界，来源与参数见[此前完整证明](../earlier-results/general_n_proof.md)。这些额外项可以继续提高 Lₙ，但不影响上面六项公式的普遍有效性。

[代码审阅与范围](REVIEW.md)、[文件哈希](source_manifest.json)、[重放报告](verification_report.json)。`.gitattributes` 固定字节，不转换换行。程序通过和内部 AI 审阅不等于第三方同行评审或整篇形式化证明。

## 4. 与八点主问题及渐近线性项的关系

本轮把真子集严格约束有效化，并提供每个有限 n 的有理下界。八点合并构型下界仍达到中心正七边形候选值；无条件全局上界仍未达到该值，一般非零矩分支尚未排空。以上结论没有确定精确 γ 或有限 n 的真实边界修正。[线性渐近项与圆盘上界专题](../asymptotics/README.md)继续独立保留。

未指定许可证；未将通常的 Gaussian、Parseval、Rayleigh 商理论或外部重整化能量定理宣称为项目原创。

