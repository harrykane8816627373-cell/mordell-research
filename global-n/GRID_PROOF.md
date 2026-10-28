# 有效网格的替代估计：全局二阶条件给出平方收敛证书

研究资料与项目发起：蒋丰灿（Fengcan Jiang）；证明整理与实现包含 OpenAI Codex 协作。2026-10-27。

来源更正（2026-10-28）：此前[一般n补充材料审核第9节](../earlier-results/general_n_proof.md)因未取得完整主稿而保留待补；用户现提供的完整原稿已经明确定义网格并给出证明及完整枚举，修订稿还有更强三阶范数控制。本节是后续另一误差估计及安全筛选的整理，不能再概括为原稿缺证。详见[完整源码核对](../source-review/SOURCE_REVIEW.md)。所用力方程、全局二阶稳定性和严格点距都是项目已有工具。

## 1. 有限整数网格及目标

设 n≥3，N=n(n−1)/2，λ=(n−1)/2。定义

```math
M_n=\max_{\sum_i|z_i|^2=n}\prod_{i\lt j}|z_i-z_j|^2,
\qquad
A_{n,q}=\left\lceil q\sqrt{\frac{2(n-1)}n}\right\rceil,
\qquad
\Lambda_{n,q}=\{-A_{n,q},\ldots,A_{n,q}\}^2.
```

对该有限整数方格中每个含 n 个互异点的无序集合 v，定义

```math
D(v)=\sum_{i\lt j}|v_i-v_j|^2
=n\sum_i|v_i|^2-\left|\sum_i v_i\right|^2,
\qquad
G_n(q)=\max_{\substack{v\subset\Lambda_{n,q}\\|v|=n}}
\left(\frac{n^2}{D(v)}\right)^N
\prod_{i\lt j}|v_i-v_j|^2.
```

每个被比较的数均为正有理数，对应实际构型 n(v−平均值)/√D。故完整枚举给出可达下界 G_n(q)≤M_n。网格单位长度在中心化归一化后抵消，但随 q 增长的有限整数箱保留了近似任意最优形状所需的比例精度。程序枚举全部无序集合，不是从中抽样。

## 2. 二阶条件的准确使用

全局最优构型 z 是中心化、无碰撞的；它也全局最大化 Δ(z)exp(−λΣ|z|²)。记

```math
(Bh)_i=\sum_{j\ne i}\frac{h_i-h_j}{(z_i-z_j)^2}.
```

直接展开一阶、二阶项，得

```math
\sum_{j\ne i}\frac1{z_i-z_j}=\lambda\overline z_i,
\qquad
\left|\sum_{i\lt j}\left(\frac{h_i-h_j}{z_i-z_j}\right)^2\right|
\le\lambda\sum_i|h_i|^2.
```

第二式对任意复方向 h 成立：二阶极大条件给出 Re(hᵀBh)+λ‖h‖²≥0，将 h 乘以任意单位复数，即可把 hᵀBh 的实部取为其模长的负值。无需假设驻点是非退化的。

已有严格真子集定理在 m=2 给出 |z_i−z_j|²>4/(n−1)，在 m=n−1 给出 |z_i|²<2(n−1)/n。本节只把它们用于真实全局最大点及其圆整附近，不对所有格点假设这些必要条件。

## 3. 圆整误差定理

取整数 q，使 2q²>n−1。任选满足下列条件的有理数 τ：

```math
\sqrt{\frac{n-1}{2q^2}}\le\tau\lt1,
\qquad
\eta_{n,q}(\tau)=\frac{n(n-1)}{2q^2}
\left(1+\frac{n\tau}{6(1-\tau)}\right).
```

则有全局证书

```math
\boxed{G_n(q)\le M_n\le G_n(q)\exp\bigl(\eta_{n,q}(\tau)\bigr).}
```

证明：将每个 qz_i 的实部、虚部分别圆整到最近整数 v_i。半径界保证 v_i 位于规定整数箱。令 y_i=v_i/q，再中心化，写 y−平均值=z+h。误差满足

```math
\sum_i h_i=0,\qquad E=\sum_i|h_i|^2\le\frac n{2q^2},
\qquad
r_{ij}=\frac{h_i-h_j}{z_i-z_j},\quad |r_{ij}|\le\tau\lt1.
```

最后一步由每个原始圆整误差模长≤1/(√2q)及点距下界得出；中心化不改变差值。因此没有圆整碰撞，Q=Σ|z+h|²>0。另由零和误差恒等式

```math
\sum_{i\lt j}|h_i-h_j|^2=nE,
\qquad
\sum_{i\lt j}|r_{ij}|^2\le\frac{n(n-1)}4E.
```

归一化圆整构型的乘积记为 F(y)。其对数相对最优值满足精确式

```math
\log\frac{F(y)}{M_n}
=2\Re\sum_{i\lt j}\log(1+r_{ij})-N\log(Q/n).
```

在单位圆内用解析对数；实部等于模长的对数。幂级数和几何尾项给出

```math
\left|2\Re\bigl[\log(1+r)-r+r^2/2\bigr]\right|
\le\frac{2|r|^3}{3(1-|r|)}.
```

又 log(Q/n)≤(Q−n)/n。代入力方程，一次项完全抵消，留下

```math
\log\frac{M_n}{F(y)}
\le\lambda E+\Re\sum_{i\lt j}r_{ij}^2
 +\frac{2\tau}{3(1-\tau)}\sum_{i\lt j}|r_{ij}|^2
\le(n-1)\left(1+\frac{n\tau}{6(1-\tau)}\right)E
\le\eta_{n,q}(\tau).
```

其中二次项正是第2节的稳定性约束，控制所有方向而无须分类形状。F(y)属于完整网格候选，故 F(y)≤G_n(q)，定理得证。

## 4. 显式有理停止与一般 n 的值

令 a_n=ceil√((n−1)/2)。若 q≥(n+1)a_n，取 τ=a_n/q，则

```math
\eta_{n,q}(\tau)\le\frac{7n(n-1)}{12q^2}.
```

给定任意正有理数 ε，预先选整数 q 满足

```math
q\ge(n+1)a_n,\qquad
q^2\ge\frac{7n(n-1)(1+\varepsilon)}{12\varepsilon}.
```

利用 exp(x)≤1/(1−x)（0≤x<1），即可得到全有理区间

```math
\boxed{G_n(q)\le M_n\le(1+\varepsilon)G_n(q).}
```

因此对每个固定有限 n，完整网格由整数比较在预定有限步内给出任意相对精度，并且

```math
M_n=\lim_{q\to\infty}G_n(q),\qquad
0\le\log M_n-\log G_n(q)\le\frac{7n(n-1)}{12q^2}
\quad(q\ge(n+1)a_n).
```

由于整数箱随 q 单调扩大，G_n(q)也单调不减。这是对一般 n 真正全局值的有效包围和可达构型序列，不是新的闭式公式或已计算出任意大 n 的精确代数数。旧稿已经有系数矩的全局近似和有限代数选择；这里补的是另一条直接产出构型的有限网格路线。

例如 n=8、ε=1/100，可预先取 q=58。这个计划只说明数学上的充分网格大小，尚未执行其极大的全枚举，不能将它写成八点1%数值证书。程序默认在搜索数量超过预算时立即拒绝，不把部分搜索报告成全局上界。

两点独立处理：M₂=4，规定网格中任意两个互异点归一化后都取得4。无需调用不存在的两点真子集定理。

## 5. 将已有严格子集成果用于安全删格

可以把旧的显式子集余量和本节误差控制联合，而不把最优点的必要条件错误地直接加到所有格点上。记旧专题的有理下界为

```math
s_{n,m}=\frac{m(m-1)+2\varepsilon_{n,m}}{n-1}\quad(2\le m\lt n),
\qquad s_{n,n}=n.
```

对真实全局最大点的最近整数圆整 v，每个 m 子集的内部圆整误差平方和≤m/(2q²)。由平方范数的初等不等式（q≥1）

```math
\|a+b\|^2\ge(1-1/q)\|a\|^2-(q-1)\|b\|^2,
```

得到全部显式有理条件

```math
\boxed{\displaystyle
D_S(v):=\sum_{i\lt j,\ i,j\in S}|v_i-v_j|^2
\ge m q(q-1)s_{n,m}-\frac{m^2(q-1)}2
\quad(2\le m\le n).}
```

q=1时右端为零；q>1时前述范数不等式来自2 Re⟨a,b⟩≥−‖a‖²/q−q‖b‖²。对真子集使用旧余量，对全体使用精确二阶矩 n。保留非严格号可以安全地用于整数筛选。

若在完整有限网格中删除违反任一上述条件的集合，所得受筛最大值记为 G_n^cut(q)。最近圆整的真正最优点一定保留，所以第3、4节的全部上、下误差证书对 G_n^cut(q)仍成立；受筛值不再声称随 q 单调。它可能小于未经筛选的 G_n(q)，但仍有同一个真实圆整见证。这是旧子集成果与新方向的直接组合。

`passes_rounded_subset_cuts`实现这个可选精确检查；已存全枚举值采用未筛网格。测试用可精确表示的四点最优构型检查它的最近圆整确被保留。没有据此声称八点搜索已经穷尽。

## 6. 实际核验范围

代码给出小规模完整枚举及独立的整数归一化复核。`grid_values.json`明确记录实际执行的 n、q、箱大小、总集合数和见证；`grid_plans.json`只记录未执行的误差计划。停止准则的全 n 有效性由上面证明承担，有限样本核验不能代替该证明。

本节只补证明确义的本节版本；旧资料保持原字节，并链接本节说明待补项的处理。无全球优先权声明，无外部同行评审声明。
