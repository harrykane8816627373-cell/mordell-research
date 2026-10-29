# 截断逆半径运输：以有效 n 阶缺额推出二分之一层数系数

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**。推导、精确程序及交叉审阅包含 OpenAI Codex 协作。2026-10-29。

本节把同圆装点容量与连续径向能量直接连接。新的关键步骤是：截断逆半径函数的 Dirichlet 范数可以精确计算；它在截断圆周上的拉普拉斯跳项为负，故平滑圆盘均值相对中心值的上误差是二阶的。联合本轮独立证明的同圆谱容量，以及用户已有的有效多环下界，即得

~~~math
\boxed{\displaystyle\liminf_{n\to\infty}\frac{K_n}{\sqrt n}\ge\frac12.}
~~~

K_n 是任意全局最优构型相对于重心的不同半径数。此推论只需要一个完全有效的 \(\delta_n=O(n)\) 上包，**不调用 Sandier–Serfaty 线性渐近项、三角晶格常数或其未知有效起始阈值**。本文还给出带有限误差的必要层数下界；它不说明真实层数比值趋于1/2，也没有给出一般 n 的最优构型分类。

## 1. 归一化与运输不等式

沿用

~~~math
\Delta_n(z)=\prod_{i<j}|z_i-z_j|^2,\quad
M_n=\max_{\sum|z_i|^2=n}\Delta_n(z),\quad
U_n=n^{n/2}\left(\frac{2n}{n-1}\right)^{n(n-1)/2}e^{-n(n-1)/4}.
~~~

令

~~~math
\delta(z)=\log\frac{U_n}{\Delta_n(z)},\qquad
x_i=\sqrt{\frac{n-1}n}z_i,\qquad \sum_i|x_i|^2=n-1.
~~~

任取 \(0<a\le\sqrt2\)，定义截断逆半径

~~~math
w_a(R)=\min\{1/a,1/R\},\qquad w_a(0)=1/a.
~~~

取 h>0，令 \(r=h\sqrt{2/n}\)，\(\mu_h=n^{-1}\sum_i\sigma_r(\cdot-x_i)\)，\(\mu_0\) 为半径 \(\sqrt2\) 的均匀圆盘概率，\(\nu_h=\mu_h-\mu_0\)。沿用

~~~math
\psi(h)=\frac{h^2-1}2-\log h,\qquad
\delta(z)+n\psi(h)\ge n^2D(\nu_h).
~~~

如果还独立证明同一个 h 下的角向能量下界

~~~math
n^2D(\nu_h-\overline\nu_h)\ge A\ge0,
~~~

则任意已证 \(\delta(z)\le E\) 给出径向预算

~~~math
n^2D(\overline\nu_h)\le E+n\psi(h)-A=:\mathcal E.
~~~

\(\mathcal E<0\) 时，给定条件类已经被排除。\(\mathcal E\ge0\) 时，本节的通用运输结论为

~~~math
\boxed{\displaystyle
\frac1n\sum_i w_a(|x_i|)
\ge \sqrt2-\frac a2
-\frac{\sqrt{\mathcal E}}{\sqrt2\,n a}
-\frac{h^2}{4na^3}.}
~~~

右侧可与0取最大。没有角向输入时取 A=0。这个运输结论适用于任意归一化、无碰撞且满足所给能量预算的构型；把它转换成层数结论时，再加入全局最优者的同圆容量条件。

## 2. 精确 Dirichlet 范数与平衡测量

把 \(w_a\) 视为平面上的径向函数。其梯度在 R<a 时为0，在 R>a 时大小为 \(R^{-2}\)。因此

~~~math
\|w_a\|_D^2:=\frac1{2\pi}\int_{\mathbb C}|\nabla w_a|^2
=\int_a^\infty\frac{dR}{R^3}=\frac1{2a^2}.
~~~

由于 \(w_a\) 径向，和非径向电势正交。零质量电荷的弱分部积分及 Cauchy–Schwarz 给

~~~math
\left|\int w_a\,d(\mu_h-\mu_0)\right|
\le\sqrt{D(\overline\nu_h)}\,\|w_a\|_D
\le\frac{\sqrt{\mathcal E}}{\sqrt2\,na}.
~~~

\(w_a\) 在原点有界、在无穷远衰减，其梯度平方可积；用光滑截断逼近可以完成上述弱测试。紧支撑零质量电荷的无穷远边界项消失。

平衡圆盘的半径密度为 \(R\,dR\)，所以对 a≤\(\sqrt2\)，

~~~math
\int w_a\,d\mu_0
=\frac1a\int_0^a R\,dR+\int_a^{\sqrt2}dR
=\sqrt2-\frac a2.
~~~

这给平滑测量的下界。剩下必须控制圆盘均值与各中心值的差，不能直接把两者相等，因为 \(1/|w|\) 在平面上不是调和函数。

## 3. 负边界跳项给出二阶圆盘误差

在 R>a 的区域，

~~~math
\Delta(1/R)=R^{-3}.
~~~

在 R=a 处，径向导数从0跳至 \(-a^{-2}\)，因此分布意义下

~~~math
\Delta w_a=\mathbf1_{R>a}R^{-3}\,dA
-a^{-2}\,ds_{\{|w|=a\}}\le a^{-3}\,dA.
~~~

这里 ds 是圆周弧长测度。负跳项不可漏掉；其符号使全平面上的上界成立。于是

~~~math
w_a(w)-\frac{|w|^2}{4a^3}
~~~

是超调和函数，故其圆盘均值不大于中心值。利用均匀半径 r 圆盘的二阶矩 \(\int|v|^2d\sigma_r(v)=r^2/2\)，得到

~~~math
\boxed{\displaystyle
\int w_a(|x+v|)\,d\sigma_r(v)
\le w_a(|x|)+\frac{r^2}{8a^3}.}
~~~

可不借助超调和函数术语直接核对常数。对光滑 f 且 \(\Delta f\le C\)，以 x 为中心的圆周均值 m(t) 满足

~~~math
m'(t)=\frac1{2\pi t}\int_{|v|\le t}\Delta f(x+v)\,dA(v)
\le Ct/2.
~~~

积分后圆周均值至多 \(f(x)+Ct^2/4\)，再对半径作均匀圆盘平均，得到 \(f(x)+Cr^2/8\)。对 \(w_a\) 作光滑卷积，保持分布拉普拉斯上界，再取局部一致极限，即得同一结论。

对所有中心求平均并反向排列，

~~~math
\frac1n\sum_i w_a(|x_i|)
\ge\int w_a\,d\mu_h-\frac{r^2}{8a^3}
=\int w_a\,d\mu_h-\frac{h^2}{4na^3}.
~~~

和第2节合并，即证明第1节的运输不等式。使用单纯 Lipschitz 位移估计会出现一阶 r 误差；这里的二阶误差是有限收敛速度改善的关键。

## 4. 从谱容量得到带误差的层数必要条件

本轮独立证明的[同圆谱容量](ANGULAR_CAPACITY.md)在库仑尺度给出每层人口 m 与半径 R 的约束

~~~math
m-1\le\sqrt{8n}\,R.
~~~

偶人口有更强的 \(m\le\sqrt{8n}\,R\)；本节统一使用上述仿射形式。对非零半径层，因为 \(w_a(R)\le1/R\) 且 \(w_a(R)\le1/a\)，

~~~math
m w_a(R)=(m-1)w_a(R)+w_a(R)
\le\sqrt{8n}+1/a.
~~~

单点层同样成立，包括可能的原点层。故 K 层满足

~~~math
K\ge\frac{\sum_i w_a(|x_i|)}{\sqrt{8n}+1/a}.
~~~

联立运输不等式，得到通用有限结果

~~~math
\boxed{\displaystyle
\frac K{\sqrt n}\ge
\frac{\sqrt2-a/2-\sqrt{\mathcal E}/(\sqrt2\,na)-h^2/(4na^3)}
{2\sqrt2+1/(a\sqrt n)}.}
~~~

分子可与0取最大。本节不把同层容量当作另一项可相加的能量罚项：它先约束每条层的加权人口，再将同一个径向测量转换为 K 的下界。

### 4.1 完全有效的二分之一渐近系数

取 h=1、A=0、\(a=n^{-1/4}\)。用户已有有效多环可达下界联合圆盘上界给

~~~math
\delta_n:=\log(U_n/M_n)
<E_n:=\frac n{100}+\frac{77}6\sqrt n+\frac{35}4
\qquad(n\ge6).
~~~

来源见[此前有效径向层数证明](../radial-angular/RADIAL_LAYERS.md)。于是对所有全局最优者，

~~~math
\boxed{\displaystyle
\frac{K_n}{\sqrt n}\ge
\frac{\sqrt2-n^{-1/4}[\,3/4+\sqrt{E_n/(2n)}\,]}
{2\sqrt2+n^{-1/4}}.}
~~~

E_n/n 有界，右端有效地趋于1/2。更明确地，

~~~math
\frac{K_n}{\sqrt n}\ge\frac12-
n^{-1/4}\frac{5/4+\sqrt{1/200+77/(12\sqrt n)+35/(8n)}}
{2\sqrt2+n^{-1/4}}.
~~~

这给一个完全显式的 \(O(n^{-1/4})\) 误差，并推出第1节的 \(\liminf\ge1/2\)。推导没有使用 \(\delta_n/n\) 收敛或 γ 的任何值。更简洁的有限常数推论见本专题的有效层数文档；通用程序保留完整根式上包。

这里的误差是**已证必要下界的误差**，不是 K_n 的真实第二阶渐近展开，也没有给 K_n 上界。

### 4.2 保留同 h 的角向预算

若同 h 的局部／全局角向 Gram 证书给 A>0，径向运输误差中的 E+nψ(h) 可以安全替成 E+nψ(h)−A。这样会提高截断测量及必要层数下界。只有已证明与径向能量正交的角向部分可以如此扣除；另一个完整径向下界不能重复扣除。

## 5. 不使用谱容量的独立 π 基线

用户已有的全局两点间距，在库仑尺度给 \(|x_i-x_j|^2>4/n\)（n≥3）。结合圆周最短角间隔，一条人口 m≥2 的层满足

~~~math
\frac mR<\pi\sqrt n.
~~~

如果 \(a\ge1/(\pi\sqrt n)\)，单点层也有 \(w_a(R)\le1/a\le\pi\sqrt n\)，因此

~~~math
K\ge\frac{\sum_i w_a(|x_i|)}{\pi\sqrt n}.
~~~

联同第1节运输式，得到不使用新谱容量的独立必要条件。仍取 h=1、\(a=n^{-1/4}\)、E=E_n，得到

~~~math
\liminf\frac{K_n}{\sqrt n}\ge\frac{\sqrt2}{\pi}.
~~~

这个基线使用已有间距、有效缺额和径向运输即可成立，可以独立交叉核对新谱容量接口的尺度。它不与第4节的完整层数下界直接相加。n=2 时可直接核对那一对反向等长点的层容量。

## 6. 原尺度 CDF 运输的独立对照

旧证明在原 z 尺度给平方半径 CDF 误差

~~~math
\epsilon\ge\sup_{s\ge0}\left|
\frac1n\#\{|z_i|^2\le s\}-\min\{s/2,1\}\right|,
~~~

其中可取

~~~math
\epsilon=\min\left\{1,
\left(\frac{6E}{n^2}\right)^{1/3}+\frac2{\sqrt n}+\frac1n\right\}.
~~~

保留最后1/n的尺度修正。只用单调函数总变差，原始估计为

~~~math
\frac1n\sum_iw_a(|z_i|)\ge\sqrt2-a/2-\epsilon/a.
~~~

还可联合全局最优者已有的 \(|z_i|<\sqrt2\) 半径上界，把 CDF 下端与0截断后再积分。若记半径 CDF 为 G_z，有

~~~math
\frac1n\sum_iw_a(|z_i|)
=\frac1{\sqrt2}+\int_a^{\sqrt2}\frac{G_z(R)}{R^2}\,dR,
\qquad G_z(R)\ge\max\{R^2/2-\epsilon,0\}.
~~~

直接积分给更稳的有限下界

~~~math
\boxed{\displaystyle
\frac1n\sum_iw_a(|z_i|)\ge
\begin{cases}
\sqrt2-a/2-\epsilon/a+\epsilon/\sqrt2,&a^2\ge2\epsilon,\\
\sqrt2+\epsilon/\sqrt2-\sqrt{2\epsilon},&a^2\le2\epsilon.
\end{cases}}
~~~

两式在边界吻合，且均可与 \(1/\sqrt2\) 取最大。这里需要原尺度半径上界，不能默认用于任意能量预算下的构型。原尺度仿射容量的每层预算为 \(\sqrt{8(n-1)}+1/a\)，与第4节库仑尺度的 \(\sqrt{8n}+1/a\) 不同。

E_n=O(n) 时，该 CDF 路线通常只给 \(O(n^{-1/6})\) 的误差；第1至4节直接能量路线给 \(O(n^{-1/4})\)。本节将 CDF 方式保留为独立有限比较，没有把两个完整径向估计重复相加。

## 7. 精确实现、复核与有限参数比较

[radial_transport.py](radial_transport.py) 所有接受判断均用 Fraction、整数平方根和已认证的根式／对数外包：

- `energy_mean_lower` 实现任意有理 a、h 及同 h 正交角向下端 A 的中心测量下界。
- `affine_layers` 使用本轮谱容量；`pi_layers` 使用独立两点间距基线，并核对单点层的截断阈值。
- `cdf_mean_lower`、`cdf_affine_layers` 在原 z 尺度实现带0截断的 CDF 对照。
- `finite_family_best` 只在列明的有限有理 a、h 集合中取最强证书，不宣称找到最优参数或最优构型。

仅输入用户已有的 E_n，接受程序给以下整数必要下界：

|n|直接能量，a约为n⁻¹ᐟ⁴、h=1|列明有限参数家族|独立CDF加谱容量|能量加π基线|
|---|---|---|---|---|
|10⁴|45|46|39|42|
|10⁶|485|487|459|442|
|10⁸|4954|4958|4835|4476|
|10¹²|499534|499575|496783|449897|

表中只是这些指定证书得到的必要层数，未与全部历史方法取总包络。报告保留实际有理截断值与平滑参数，未宣称执行相应规模的点构型优化或驻点枚举。

运行

~~~bash
python -B -S -O verify_radial_transport.py
~~~

[验收报告](radial_transport_verification.json)核对范数及分布通量、二阶误差常数、独立径向环电荷对偶、两类同层容量的代数接口、CDF截断积分、角向预算扣除、非法输入及上述有限例。一般的分布超调和比较与全 n 量词由正文解析证明承担，有限采样不替代这些步骤。

本节沿用[已有圆盘分解与可调平滑](../radial-angular/RADIAL_LAYERS.md)、[严格子集间距](../finite-n/SUBSET_PROOF.md)及用户已有有效多环下界；新谱容量的证明和限制见[同圆谱容量](ANGULAR_CAPACITY.md)。原有定理与外部来源不改署为本节成果。本文不作全球优先权、外部同行评审或全篇形式化验证声明；完整八点锐界、一般 n 最优构型分类和精确 γ 仍然开放。
