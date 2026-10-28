# 任意有理构型的显式上升证书与有限驻点残差停止

研究资料与项目发起：蒋丰灿；本次通用上升推导及实现包含 OpenAI Codex 协作。2026-10-28。本节不声称梯度上升框架的原创优先权；补充的是适用于当前齐次距离乘积的明确有理步长、严格增益及适用边界。

## 1. 归一化目标与准确残差

令 v=(v₁,…,vₙ)为无碰撞、中心化的复平面构型，Q=Σ|vᵢ|²>0，N=n(n−1)/2。定义不依赖尺度的可达值

```math
F(v)=\left(\frac nQ\right)^N\prod_{i\lt j}|v_i-v_j|^2\le M_n,
\qquad
f_i=\sum_{j\ne i}\frac1{v_i-v_j},\quad
r_i=\overline{f_i}-\frac NQv_i.
```

以下全部量可用有理数精确表示，只需 vᵢ∈Q+iQ。更一般地，可取 vᵢ=aᵢ+i√d bᵢ，aᵢ,bᵢ∈Q、d为固定正整数；这同时包括正方整数格和三角格（d=3，a+bω=(a+b/2)+i√3 b/2）。更新后的平方距离、二阶矩和F仍为有理数。

成对求和给Σfᵢ=0和Σvᵢfᵢ=N。因此

```math
\sum_i r_i=0,\qquad \sum_i\overline v_i r_i=0,
\qquad R:=\sum_i|r_i|^2
=\sum_i|f_i|^2-\frac{N^2}Q,
\qquad \sum_i f_ir_i=R.
```

于是沿实参数 t，有 Q(v+tr)=Q+t²R，且 logF的一阶方向导数恰为2R。r=0等价于当前形状是约束驻点，**不等价于局部或全局最优**。

## 2. 不需线搜索的有理步长

记

```math
d_{ij}=\frac{r_i-r_j}{v_i-v_j},\qquad
A=\sum_{i\lt j}|d_{ij}|^2,
\qquad C\ge\max_{i\lt j}|d_{ij}|,
\qquad H=2A+\frac NQ R.
```

C取任一正有理上界。实现中，对 dᵢⱼ=α+i√dβ，取 C=max(|α|+d|β|)，这是合法上界且无需平方根。若 R>0，则 C>0、H>0，取

```math
\boxed{\displaystyle t=\frac{R}{H+2CR}\in\mathbb Q_{>0}.}
```

它满足 tC<1/2，因此更新点互异。归一化更新的完整对数增益满足

```math
\boxed{\displaystyle
\log\frac{F(v+tr)}{F(v)}
\ge tR=\frac{R^2}{H+2CR}>0.}
```

证明：对 |w|<1，幂级数给 |log(1+w)−w|≤|w|²/[2(1−|w|)]。利用 tC<1/2以及 log(1+x)≤x，得到

```math
\begin{aligned}
\log\frac{F(v+tr)}{F(v)}
&=2\Re\sum_{i\lt j}\log(1+td_{ij})
-N\log\left(1+\frac{t^2R}Q\right)\\
&\ge2tR-t^2\left(2A+\frac NQ R\right)
=2tR-Ht^2\ge tR.
\end{aligned}
```

最后一步用 Ht≤R。故不必反复尝试步长或采用浮点优化，每个非驻点的有理构型都给一个严格更大的可达有理下界。此外

```math
\log M_n-\log F(v)\ge\frac{R^2}{H+2CR}.
```

这是证明当前值尚距全局至少多少的**缺口下界**，不能倒读成全局误差上界。

## 3. 迭代在有限步内达到指定驻点残差

归一化残差 S(v)=Q R在整体缩放下不变。给定初始F₀>0，沿上述上升过程F≥F₀。因为每条距离平方≤2Q，而所有距离平方的乘积为F(Q/n)ᴺ，得到一致归一化点距下界

```math
\frac{\min_{i\lt j}|v_i-v_j|^2}Q\ge
a:=\frac{F_0}{n^N2^{N-1}}>0.
```

由柯西不等式及残差恒等式，S≤n(n−1)²/a。零和残差还给Σpairs|rᵢ−rⱼ|²=nR，因此A≤nR/(aQ)。对实现采用的C，有

```math
C^2\le(d+1)A,\qquad
CQ\le\frac{n(n-1)(d+1)}a.
```

后式用前述S上界和√(d+1)≤d+1。定义明确的正有理数

```math
B=N+\frac{2n+2n(n-1)(d+1)}a.
```

则 H+2CR≤B R/Q，每步对数增益至少为S/B。另一方面，Σpairs|vᵢ−vⱼ|²=nQ，AMGM给

```math
F(v)\le U^{\rm AMGM}_n=\left(\frac{2n}{n-1}\right)^N.
```

所以对任意正有理ε，最多

```math
\left\lceil\frac B\varepsilon
\left(\frac{U_n^{\rm AMGM}}{F_0}-1\right)\right\rceil
```

次更新内，必定出现 S≤ε。否则每步log增益严格大于ε/B，而log(U/F₀)≤U/F₀−1，矛盾。这是有限驻点残差保证；常数可能极大，不声称实际效率低廉。事实上ΣS有限，归一化形状的每个聚点都是无碰撞驻点。

## 4. 必须保留的边界

梯度残差为零时，上升步骤停止，但可能仍是严格次优局部极大或鞍点。例如三点共线 v=(−1,0,1) 已满足 r=0，其归一化值27/2仍小于M₃=27；五点原点加正方形也是已知次优严格局部最大。不能把“驻点容差已达到”写成“全局误差已达到”。

本方法适合将既有方格、三角格和保存整数候选转为无随机、可重放的严格下界改进，再与当前完整包络取最大。它不能取代全局上界、稳定性检查、关键值覆盖或精确全局证明，也没有排空尚未解决的一般分支。
