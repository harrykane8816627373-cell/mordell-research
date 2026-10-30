# 从库仑缺额控制完整二维点分布：显式 W₁ 运输上下界

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；推导、精确程序及交叉审阅含 OpenAI Codex 协作。2026-10-30。承接仓库 `26c8fcd`。

本节把已有能量、可达下界和支撑半径限制转成**全平面**经验测度的运输误差。它同时控制所有 Lipschitz 测试函数，不只观察半径。对一般全局最优构型，得到有效的 \(O(n^{-1/2})\) 二维 W₁ 上界；任意 n 个等权原子的球重排下界也是同一阶。沿用已有线性能量和三角试探输入，常数可夹在约0.942809与0.980240之间。

库仑能控制概率距离是经典框架。相关一手文献为 Chafaï–Hardy–Maïda, *Concentration for Coulomb gases and Coulomb transport inequalities*, J. Functional Analysis 275 (2018), 1447–1483，[arXiv:1610.00980](https://arxiv.org/abs/1610.00980)。本文不把这一框架登记为新发明；下面自含推导圆盘情形的常数、处理本项目归一化，并把现有严格下界输入二维分布证书。不作全球优先权声明。

## 1. 测度、距离与有限 n 主结论

仍令 \(\sum_i|z_i|^2=n\)，并记

```math
\eta_z=\frac1n\sum_i\delta_{z_i},\qquad
x_i=\alpha z_i,\quad \alpha=\sqrt{\frac{n-1}{n}},\qquad
\mu_0=\frac1{2\pi}\mathbf1_{|w|\le\sqrt2}\,dw.
```

\(\mu_0\) 是均匀半径 \(\sqrt2\) 圆盘概率，其二阶矩为1。W₁距离定义为把一个概率分布搬运成另一个时的最小平均欧氏距离：

```math
W_1(\eta,\mu)=\inf_{\Gamma\in\Pi(\eta,\mu)}
\int|u-v|\,d\Gamma(u,v)
=\sup_{\operatorname{Lip}(f)\le1}\left|\int f\,d(\eta-\mu)\right|.
```

最后是标准 Kantorovich–Rubinstein 对偶。紧支撑下所有量都有限。

取任何 h>0，令 \(r=h\sqrt{2/n}\)、\(\mu_h=n^{-1}\sum_i\sigma_r(\cdot-x_i)\)、\(\nu_h=\mu_h-\mu_0\)，并设

```math
\delta(z)=\log\frac{U_n}{\Delta_n(z)},\qquad
\psi(h)=\frac{h^2-1}{2}-\log h.
```

此前[可变平滑分解](../radial-angular/RADIAL_LAYERS.md)给 \(n^2D(\nu_h)\le\delta+n\psi(h)\)。若另知每个库仑中心 \(|x_i|\le\rho\)，取

```math
R=\max\{\sqrt2,\rho+r\}.
```

则对任意归一化无碰撞构型都有

```math
\boxed{\quad
W_1(\eta_z,\mu_0)
\le\frac{2r}{3}
+\frac Rn\sqrt{\delta(z)+n\psi(h)}
+1-\sqrt{\frac{n-1}{n}}.
\quad}
```

这里 R 必须覆盖**背景圆盘和全部平滑点圆盘**。旧角向能量可以忽略径向背景超出支撑的例外，不适用于这里的完整 W₁ 推导。

对全局最大构型，旧严格子集矩结果给

```math
\rho\le\sqrt2\frac{n-1}{n}.
```

将任意已证可达下界 \(L_n\le M_n\) 所给的 \(E\ge\log(U_n/L_n)\) 代替 δ，即得所有全局最优构型的有效二维分布上界。程序不会把任意输入的 E 自动当作已证预算。

## 2. 紧支撑库仑电荷的圆盘运输常数

先证明本节使用的圆盘常数。设零质量电荷 \(\nu=\mu-\mu_0\) 支撑在 \(B_R\)，势为 \(u=g*\nu\)、\(g=-\log|\cdot|\)，并具有有限 Dirichlet 能量

```math
D(\nu)=\frac1{2\pi}\int|\nabla u|^2.
```

本项目的两测度有紧支撑、有界密度，满足这项条件。则

```math
\boxed{\quad W_1(\mu,\mu_0)\le R\sqrt{D(\nu)}.\quad}
```

证明。取任意1-Lipschitz实函数 f。先对光滑近似论证，再按一致近似和弱导数极限回到 Lipschitz 函数。令 C 为 f 在半径R圆周上的平均值。圆内保留 \(f-C\)，圆外使用同一边界值的衰减调和延拓，拼成 F。下面明确控制整个平面的能量。

令 H 是圆内具有同一边界值的调和函数。它在原点的值为边界均值0。调和函数的 Dirichlet 最小性给

```math
\int_{B_R}|\nabla H|^2
\le\int_{B_R}|\nabla f|^2
\le\pi R^2.
```

第一步可由展开 \(f-C=H+v\)、v 的边界值为零，以及 \(\int\nabla H\cdot\nabla v=0\) 直接证明。二维 Kelvin 变换保持 Dirichlet 能量：圆外延拓为 \(H(R^2w/|w|^2)\)，其圆外能量等于 H 的圆内能量；H(0)=0保证它在无穷远衰减。因此

```math
\frac1{2\pi}\int_{\mathbb C}|\nabla F|^2\le R^2.
```

总质量为零且电荷支撑在圆内，故 \(\int Fd\nu=\int fd\nu\)。弱分部积分和 Cauchy–Schwarz 给

```math
\left|\int f\,d\nu\right|
=\left|\frac1{2\pi}\int\nabla F\cdot\nabla u\right|
\le R\sqrt{D(\nu)}.
```

取 Lipschitz 对偶上确界即得结论。延拓在边界只有弱导数要求，不要求法向导数连续。有限 Fourier 和、平滑逼近或齐次 Sobolev 极限均可使上述分部积分严格化；无穷远边界项由衰减消失。

## 3. 平滑和原尺度转换的代价

把每个中心点和其平滑圆盘直接耦合，平均移动距离恰为

```math
\int|v|\,d\sigma_r(v)=\frac{2r}{3}.
```

因此 \(W_1(\eta_x,\mu_h)\le2r/3\)。将第2节与已有 \(D\le(\delta+n\psi)/n^2\) 相接，得到 x 尺度的前两项。

最后，逐点耦合 x 与原尺度 z。由于 \(n^{-1}\sum|z_i|\le1\)，

```math
W_1(\eta_z,\eta_x)
\le(1-\alpha)\frac1n\sum_i|z_i|
\le1-\alpha.
```

三角不等式给第1节主式。尺度修正是必要的，不能把 x 的经验测度直接标成原尺度 z 的分布。

另外，独立耦合及参考圆盘均值为零给简单全局上端 \(W_1(\eta_z,\mu_0)\le\sqrt2\)。x尺度相应为 \(\sqrt{2-1/n}\)。实现始终取这个上端和能量上端中的较小值，小 n 的粗能量预算因而不会给出比这一简单界更差的最终报告。

## 4. 任意 n 个等权原子都有同阶量化下界

不管这些点是否最优、是否中心化或满足矩约束，任何 n 个等权原子与参考圆盘都满足

```math
\boxed{\quad
W_1\left(\frac1n\sum_i\delta_{z_i},\mu_0\right)
\ge\frac{2\sqrt2}{3\sqrt n}.
\quad}
```

证明。取任意运输耦合。向第 i 个有标签原子运输的参考子测度记为 \(\gamma_i\)，质量为1/n且 \(\gamma_i\le\mu_0\)。背景密度不超过 \(1/(2\pi)\)，故任何半径t圆盘中的 \(\gamma_i\) 质量至多 \(t^2/2\)。层蛋糕积分给

```math
\int|w-z_i|\,d\gamma_i(w)
\ge\int_0^{\sqrt{2/n}}\left(\frac1n-\frac{t^2}{2}\right)dt
=\frac{2\sqrt2}{3n^{3/2}}.
```

对 i 求和并对运输取下确界即得。它是球重排的普适下界，不保证这些最优圆形运输区域能同时平铺；因而没有声称该常数是真实最优量化常数，更没有确定三角格点。

## 5. 全局最优者的有效速率及渐近窄常数夹逼

此前初等可达下界给 n≥6 的有效预算

```math
\delta_n<\frac n{100}+\frac{77}{6}\sqrt n+\frac{35}{4}.
```

配合全局半径上端，第1节立即给所有全局最优构型的有效 \(W_1=O(n^{-1/2})\)。这个结论只需上述有效预算，不依赖未知有效起始阈值的外部渐近余项。

若进一步沿用[已有线性能量与三角试探常数](../asymptotics/README.md)，则

```math
\frac{\delta_n}{n}\longrightarrow\gamma_D-\gamma
\le\gamma_D-\gamma_\triangle<\frac1{800}.
```

因此对任意固定 h>0，

```math
\limsup_{n\to\infty}\sqrt n\,W_1(\eta_{z^{(n)}},\mu_0)
\le\frac{2\sqrt2}{3}h
+\sqrt2\sqrt{\frac1{800}+\psi(h)}.
```

取有理参数 h=97/100，以带余量对数和向上平方根包围，右端严格小于981/1000。由第4节，下列完整二维分布结论成立：

```math
\boxed{\quad
\frac{2\sqrt2}{3}
\le\liminf\sqrt n\,W_1(\eta_{z^{(n)}},\mu_0)
\le\limsup\sqrt n\,W_1(\eta_{z^{(n)}},\mu_0)
<\frac{981}{1000}.
\quad}
```

接受程序把右端实际外包为0.980239418395以内，普适下端约为0.942809041582。这里上端保留既有外部重整化能量定理及三角试探输入；不是不调用外部输入的有限 n 常数。本文没有给出从哪个有限 n 开始可以直接使用0.981/√n，也未证明真正的 W₁ 常数存在。

W₁速率与旧径向CDF速率是不同类型的陈述，不能直接比较两个指数后声称所有指标都更强。W₁同时控制全部二维 Lipschitz 统计量；圆盘、方格的人口指示函数不连续，仍需边界层估计。

## 6. 已有实际可达下界能显著加强二维有限分布

[此前已认证的100点整数见证报告](../symmetry-breaking/witness_verification.json)给

```math
\delta_{100}\le\frac{413923829261}{500000000000}.
```

它来自一个严格可达构型，不是100点锐值。将该预算直接输入本节，取 h=97/100，得到对**任意100点全局最优构型**的严格上端

```math
W_1(\eta_z,\mu_0)<0.111199952674.
```

相比只使用初等全 n 预算时本接受器给出的0.277165020682，这正好展示“更强可达下界→更小能量预算→更紧二维分布”的闭合链条。这两个数都是必要上界，不是真实 W₁ 的计算值，也不是 M₁₀₀ 已确定的证明。

本次重用旧报告中的已证有理预算，没有重新执行其91项坐标搜索或把它们登记为新计算。该旧报告字节作为明确依赖保留。

## 7. 一个径向分布良好却角向集中的非空类

函数

```math
f(w)=\min\{|\Re w|,|\Im w|\}
```

是到两坐标轴并集的距离，因而1-Lipschitz。均匀参考圆盘中

```math
\int f\,d\mu_0=\frac{8(\sqrt2-1)}{3\pi}>\frac{56}{165}.
```

极坐标计算即可得到等式：平均半径是 \(2\sqrt2/3\)，角向平均为 \(4(1-\sqrt2/2)/\pi\)。简单严格下端只用了 \(\sqrt2>7/5\)、\(\pi<22/7\)。

若100个点都距两轴并集不超过1/20，则 \(\int f d\eta_z\le1/20\)，所以

```math
W_1(\eta_z,\mu_0)>\frac{191}{660}>0.289.
```

即使只用初等旧预算，本节已给全局最优者 \(W_1<0.281\)，故排除这一整个厚交叉带类。使用第6节实际见证预算后，连距两轴并集不超过1/5的更宽整类也被排除，因为

```math
W_1>\frac{23}{165}>0.139
\quad\hbox{而全局最优者 }W_1<0.112.
```

同一个统计量在粗的旧预算下不能排除1/5宽类；验收报告明确比较了这两条必要上端。这不是宣称其他旧几何或稳定性方法都无法排除该类。

类确实非空，且可以同时具有很好的径向数据。取25个平方半径

```math
q_j=\frac{2j-1}{25},\qquad j=1,\ldots,25,
```

每个半径放四点 \(\{\sqrt{q_j},-\sqrt{q_j},i\sqrt{q_j},-i\sqrt{q_j}\}\)。它们全部互异、中心为零、二阶矩为100，\(p_2=0\)，且最大平方半径 \(49/25<99/50\) 满足旧全局半径限制。其相对于均匀圆盘的实际径向CDF误差只有1/50，但所有点都在坐标轴上，二维统计量 f 全为零。

因此小径向CDF误差不代表二维分布充分均匀。这个例子展示本节指标捕获径向数据之外的信息；它仍可能被已有点距或其他条件排除，不登记为所有旧方法遗漏的新最优候选。

## 8. 任意紧支撑 Lipschitz 统计与矩形人口接口

若 f 为任意已证 Lipschitz 常数 L 的函数，包括连续分片线性紧支撑函数，则

```math
\left|\frac1n\sum_i f(z_i)-\int f\,d\mu_0\right|
\le L\,W_1(\eta_z,\mu_0).
```

这同时控制任意有限族测试；不能把每个误差预算相加当作另一套独立能量。`lipschitz_class_check` 接收另行证明的经验、背景积分区间和 L，检查整类是否与同一个 W₁ 上端矛盾。程序不会从任意输入数值自动证明 f 的 Lipschitz 常数或积分区间。

对任意集合 A，令 \(A^{+\varepsilon}\) 为欧氏扩张、\(A^{-\varepsilon}\) 为距补集至少 ε 的内核。距离斜坡函数给

```math
\mu_0(A^{-\varepsilon})-\frac{W_1}{\varepsilon}
\le\eta_z(A)
\le\mu_0(A^{+\varepsilon})+\frac{W_1}{\varepsilon}.
```

对轴对齐矩形，欧氏扩张包含于四边各扩张 ε 的矩形，面积给背景质量安全上端；内缩矩形若四角都在背景圆盘内，其面积给精确背景质量，使用 \(\pi\) 上端可得有理下端。否则当前实现取较弱的0下端。乘 n 后向外取整即给人口区间。小格或小 n 时可能退化为[0,n]；它不是已计算的完整网格运输方案。

## 9. 程序、复核与范围

[distribution_metric.py](distribution_metric.py)使用标准库 `Fraction` 和向外取整平方根，并复用旧 `joint-bounds/subset_moment.py` 的带余量对数包围。默认 h 只在1与97/100两个指定值中选较小的已证上界，没有声称优化了所有平滑参数。

```bash
python -B -S -O distribution-geometry/verify_distribution_metric.py
python -B -S -O distribution-geometry/verify_distribution_metric.py --write-report
python -B -S distribution-geometry/distribution_metric.py 100
```

[报告](distribution_metric_verification.json)核对549个根包围、24个 Kelvin 单项式能量恒等式、199个量化层蛋糕恒等式、9个矩形人口接口及8项非法输入；另直接核对旧100点预算字面值、非径向交叉带和25层四向构型的精确矩与径向误差。有限程序不替代正文的运输对偶、弱分部积分和球重排证明。

本节说明全局最优点集在二维平均运输意义下接近均匀圆盘，并给可算的有限限制。它没有确定最近邻图、真实微观晶格、点集唯一性、全八点锐界或一般 n 的完整最优构型；小 W₁ 也不能反向认证全局最优。
