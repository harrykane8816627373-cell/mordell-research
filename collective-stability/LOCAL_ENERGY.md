# 环带局部角向能量：调和帽函数与不抵消人口限制

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；推导、精确程序及交叉审阅含 OpenAI Codex 协作。2026-10-29。

此前全局混合矩用共同外支撑半径作多项式投影。外层的同一 Fourier 模式可能抵消内层，较大的外支撑半径也会削弱内层证书。本节直接采用有紧径向支撑的分段调和测试函数，把指定角向模式局部化到圆盘或环带。局部矩测量只在**整个平滑圆盘都位于单个调和片内**时由中心值精确给出。

结论是必要能量约束和条件构型类上界；未确定新的未知点数锐值，也未覆盖全部构型域。不声称分段调和投影或 Dirichlet 方法的全球优先权。

## 1. 一般环带帽函数及精确范数

取整数 k≥1 和三个平方半径 \(0\le a<b<c\)，令

```math
t=\frac{b^k-a^k}{c^k-b^k},\qquad q=|w|^2.
```

定义复值函数

```math
\Phi(w)=\begin{cases}
0,&q\le a,\\
w^k-a^k/\overline w^{\,k},&a<q\le b,\\
t(c^k/\overline w^{\,k}-w^k),&b<q<c,\\
0,&q\ge c.
\end{cases}
```

当 a=0 时，内片就是 \(w^k\)，没有原点奇点。取 \(\Phi(0)=0\)。函数在三个界面连续，各片调和，梯度平方可积。只用其实部和虚部作为实 Dirichlet 测试函数。沿用

```math
\langle f,g\rangle_D=\frac1{2\pi}\int_{\mathbb C}\nabla f\cdot\nabla g.
```

两者正交，并且各自的范数平方都恰为

```math
\boxed{\quad
\mathcal N_k(a,b,c)
=\frac{k(b^k-a^k)(c^k-a^k)}{c^k-b^k}>0.
\quad}
```

证明。对于 \(f(R)\cos(k\theta)\)，Dirichlet 范数为

```math
\frac12\int\left(R f'(R)^2+\frac{k^2f(R)^2}{R}\right)dR.
```

在调和片上等于 \(\tfrac12[Rff']\) 的端点差。记 \(A=a^k,B=b^k,C=c^k\)。内外两个片的贡献分别为

```math
\frac k2\left(B-\frac{A^2}{B}\right),\qquad
\frac k2\left(\frac{B-A}{C-B}\right)^2
\left(\frac{C^2}{B}-B\right).
```

相加化简得所述范数。正弦、余弦角积分相同且相互正交。

不同 k 的帽函数正交。**同一个 k，只要支撑环带 \((a,c)\) 的内部互不相交，也正交**，因为梯度支撑互不重叠。边界相接不影响面积积分。相同 k 的重叠帽函数则必须保留完整 Gram 矩阵，不能独立相加。

## 2. 局部测量的可加能量定理

仍取 \(\sum_i|z_i|^2=n\)、\(x_i=\sqrt{(n-1)/n}\,z_i\)，并使用[既有可调平滑分解](../radial-angular/RADIAL_LAYERS.md)：

```math
r=h\sqrt{2/n},\quad
\mu_h=\frac1n\sum_i\sigma_r(\cdot-x_i),\quad
\psi(h)=\frac{h^2-1}{2}-\log h,\qquad
\delta(z)+n\psi(h)\ge n^2D(\nu_h).
```

其中 \(\nu_h=\mu_h-\mu_0\)，\(\mu_0\) 为半径 \(\sqrt2\) 的均匀平衡圆盘，\(\delta=\log(U_n/\Delta_n)\)。对一组按第1节规则正交的帽函数，定义**实际平滑测量**

```math
\mathcal A_\ell=
n\int\Phi_\ell(w)\,d\mu_h(w)
=\sum_i\int\Phi_\ell(w)\,d\sigma_r(w-x_i).
```

则对所有归一化无碰撞构型均有

```math
\boxed{\quad
\delta(z)\ge n^2D(\overline{\nu_h})
+\sum_\ell\frac{|\mathcal A_\ell|^2}{\mathcal N_\ell}
-n\psi(h).
\quad}
```

这里不要求驻点、稳定或规则角度。测试函数的角向平均为零，因此与径向势正交，背景圆盘的测量也是零。弱分部积分给

```math
\langle u_a,\Re\Phi_\ell\rangle_D
=\Re\mathcal A_\ell/n,
\qquad
\langle u_a,\Im\Phi_\ell\rangle_D
=\Im\mathcal A_\ell/n.
```

分段调和界面上的法向导数可以跳跃；连续性和有限 Dirichlet 范数足以作弱测试。Bessel 不等式及径向正交分解即得定理。因为局部帽函数有紧支撑，不需要共同的全局外半径分母。

旧的同 k 全局**径向加权 j≥1** 混合矩测试通常与帽函数不正交，不能直接把两项完整能量相加；可以取两者较强值，或合并到一个 Gram 矩阵。例如 k=1、a=0、b=1、c=2，帽函数与圆内全局测试 \(w|w|^2\) 的交叉 Dirichlet 内积为−2，确实不为零。

但旧的 **j=0 全纯测试有一个额外正交例外**。设其全局支撑平方半径为 S，它覆盖全部平滑圆盘，并且每个帽的 c≤S。全局测试在整个帽支撑内恰为 \(w^k\)，实虚部均调和。紧支撑帽函数在其边界取零，弱分部积分给交叉 Dirichlet 内积为零。即使角向 k 相同也成立。用直接极坐标积分同样可见两片的交叉项分别为

```math
\frac k2(b^k-a^k),\qquad
-\frac k2\frac{b^k-a^k}{c^k-b^k}(c^k-b^k),
```

恰好抵消。因此，在这些支撑条件下，可以把

```math
\sum_{k\in\mathcal K}\frac{|p_k(x)|^2}{kS^k}
```

加到局部帽罚项之和中，再与径向项相加，并且只扣除一次共同的 \(n\psi(h)\)。这里每个 k 的全局全纯测试只计一次。这个正交性不自动推广到 j≥1，也不允许重复加入不同全局支撑半径下的同模式测试。

## 3. 何时能用中心值精确测量

一般情况下 \(\mathcal A\) 是整个圆盘的平均，不能自动换成 \(\sum_i\Phi(x_i)\)。若每个闭圆盘的内部完整位于以下四个调和区域之一：\(q<a\)、\(a<q<b\)、\(b<q<c\)、\(q>c\)，调和均值性质才给

```math
\mathcal A=\sum_i\Phi(x_i).
```

允许圆盘与片边界相切；其开圆盘仍处于单个调和区域，连续性处理边界。a=0 时，原点不构成内片的调和障碍。

对已知平方半径带 \(q_i\in[\alpha,\beta]\)，程序用

```math
(\sqrt\beta+r)^2\le B,
\qquad
(\sqrt\alpha-r)_+^2\ge A
```

认证整盘位于平方半径区间 \([A,B]\) 中，两式分别通过有理平方比较完成。若圆盘可能穿过界面，程序拒绝中心值公式；这不是把可能的测量误差默认为零。

在内片中，中心测量为

```math
\Phi(x)=x^k\left(1-\frac{a^k}{|x|^{2k}}\right);
```

外片中是 \(t x^k(c^k/|x|^{2k}-1)\)，支撑外则为零。有理输入坐标中心化后，即使归一化含公共根式，每个测量都有公共实因子 \(f^{k/2}\)，平方罚项只含有理数 \(f^k\)。

## 4. 一条有角向条件的局部人口上界

假设 m 个被计数点的平滑圆盘完整处于内片，中心平方半径至少为 \(q_0>a\)，并在某个方向 \(\theta_0\) 上满足

```math
\cos(k(\arg x_i-\theta_0))\ge\gamma>0.
```

其余点的帽测量均为零，或沿该方向的实投影非负。幅值

```math
q^{k/2}-a^k q^{-k/2}
```

随 q 严格递增，所以

```math
|\mathcal A|^2\ge
m^2\gamma^2q_0^k\left(1-(a/q_0)^k\right)^2.
```

若已独立证明 \(\delta\le E\) 及同一 h 下的径向罚项下界 \(n^2D(\bar\nu_h)\ge d\)，就有

```math
\boxed{\quad
m^2\le
\frac{\mathcal N_k(a,b,c)\,[E+n\psi(h)-d]}
{\gamma^2q_0^k(1-(a/q_0)^k)^2}.
\quad}
```

它约束带有角向集中条件的局部人口，不是任意圆盘的无条件人口上界。方向投影的不抵消条件、整盘调和条件均不可删除。实现对右端向上包围后，以整数平方根给安全的整数人口上端。

例如 n=100、h=1、\(E=1657/12\)、k=2、\(a=0,b=27/50,c=49/50,q_0=7/20,\gamma=99/100\)，右端为 \(24357900/25289\)，故上述条件下 m≤31。

## 5. 一个非空构型类的局部角向排除

取 n=100，在库仑尺度下有35点满足 \(|x_i|^2=7/20\)，其余65点满足 \(|x_i|^2=347/260\)。总矩为

```math
35\frac7{20}+65\frac{347}{260}=99.
```

加入局部矩条件

```math
\left|\sum_{|x_i|^2=7/20}x_i^2\right|^2\ge121.
```

取帽参数 \(k=2,a=0,b=27/50,c=51/50\)、h=1。所有内圈平滑圆盘均处于内片；所有外圈平滑圆盘均在帽支撑外，平方比较全部严格。因此帽测量只包含内圈。精确范数和罚项为

```math
\mathcal N=\frac{210681}{260000},\qquad
\delta(z)\ge\frac{121}{\mathcal N}
=\frac{31460000}{210681}
>\frac{1657}{12}.
```

右边的既有全局最优缺额预算来自一般 n 可达下界，见[径向人口证明](../configuration-bounds/radial_audit.md)。故整个带局部矩条件的类不含全局最大者。

该类非空且可中心化。令 \(t_j=j/1000\)、\(u_j=(1-t_j^2+2it_j)/(1+t_j^2)\)，j=1,…,8。在平方半径7/20的内圈取8组 \(\{u_j,-u_j,\bar u_j,-\bar u_j\}\) 的缩放点，共32点，再加一个同半径的正三角形。各组重心为零，正三角形的一、二次幂和均为零。内圈二次幂和是实数，且

```math
11<\frac7{20}\,4\sum_{j=1}^8\Re(u_j^2)\le\frac{56}{5}.
```

程序用有理数直接核对这些比较和互异性。外圈取正六十五边形；它的低阶幂和均为零，整体仍中心化。

在这个具体非空例中，旧的同模式 k=2、径向次数0和1的全局混合矩证书，取安全全局支撑平方半径 \(S=7/4\)，其角向罚项至多

```math
\left(\frac{56}{5}\right)^2
\left[\frac1{2S^2}+\frac4{S^4}
\left(\frac7{20}+\frac3{100}-S\right)^2\right]
=\frac{92557824}{765625}
<\frac{1657}{12}.
```

新局部角向罚项单独已经超过预算；这里比较的是明示的同模式第一级全局角向投影。**没有声称旧方法的更高次数、其他模式、完整径向与角向组合、或已有点距及稳定性约束均无法排除该类。** 本例用于展示新的局部投影能力，不是声称找到此前所有方法都不能处理的真实最优候选。

事实上，本轮[同圆人口预算](RING_POPULATION.md)也排除这组35、65人口：\(35^3+65^3=317500>\pi^2\cdot100\cdot99\)。上述具体正六十五边形见证还受此前正规大分量障碍限制。它们不影响本节一般局部能量定理，但本例不应登记为其他方法遗漏的新增最优候选排空。

## 6. 复核与接口范围

[local_energy.py](local_energy.py)实现范数、几何片归属、非重叠局部帽相加、已证局部矩下端的整类证书、有理构型直接测量以及人口上端。`orthogonal_holomorphic_bonus`检查共同外支撑覆盖全部平滑圆盘和全部帽，接收另行已证的 \(|p_k|^2\) 下端映射，返回可正交相加的 j=0 全纯矩能量；它不再次扣除平滑代价。对数包围复用 `joint-bounds/subset_moment.py`；所有接受判断均用 `Fraction`。

```bash
python -B -S -O collective-stability/verify_local_energy.py
python -B -S -O collective-stability/verify_local_energy.py --write-report
```

[核验报告](local_energy_verification.json)包含40组范数及界面恒等式、64项根式支撑比较、有理八点输入的精确测量、上述非空类及局部人口例。接受程序拒绝可能跨界的平滑圆盘和同模式重叠帽的直接相加。有限核验不替代正文的调和均值、弱分部积分和正交投影证明；全八点锐界及一般 n 最优构型仍未解决。
