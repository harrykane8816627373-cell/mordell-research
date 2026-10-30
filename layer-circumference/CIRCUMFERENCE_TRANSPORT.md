# 用径向运输控制全部层的周长总和

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；推导与精确程序包含 OpenAI Codex 协作。2026-10-30。

本节把此前的单层谱容量与新的凸点对预算连接。它控制的是不同层半径之和ΣR_j，不是按人口加权的点半径之和；混淆两者会使新层数推论失效。

## 1. 尺度、已知输入与测试函数

原约束Σ|z_i|²=n，N=n(n−1)/2，w_i=√((n−1)/2)z_i，Σ|w_i|²=N。记

~~~math
q_i=|w_i|/\sqrt n,\qquad q_j^{\rm layer}=R_j/\sqrt n.
~~~

沿用作者已证全局半径界q_i≤1，以及[同层谱容量](../population-transport/ANGULAR_CAPACITY.md)

~~~math
m_j-1\le4R_j=4\sqrt n\,q_j^{\rm layer}.
~~~

该容量覆盖单点层和原点。q平面上的平衡背景是单位均匀盘，点盘平滑半径为1/√n。对零质量电荷，库仑Dirichlet能量在该尺度变化下不变。若δ_n≤E，旧精确分解仍给n²D≤E。

对0<a≤1取整个平面上的径向函数

~~~math
f_a(q)=(1-q)\min\{1/a,1/q\},\qquad f_a(0)=1/a.
~~~

q>1时它可以为负；没有在1处截零。它趋于−1，添加常数1即可使其在无穷远趋零，不改变零质量测试。

单位背景平均与Dirichlet范数分别是

~~~math
B_a=1-a+a^2/3,\qquad
H_a=\frac1{2\pi}\int|\nabla f_a|^2=\frac{1+a^{-2}}2.
~~~

前者由2∫_0^1q f_a(q)dq分成[0,a]和[a,1]直接得到。后者是∫_0^a q/a²dq+∫_a^∞q⁻³dq。

## 2. 平滑误差的方向

q<a时f=(1−q)/a，分布拉普拉斯的面积部分为−1/(aq)；q>a时f=1/q−1，其拉普拉斯为q⁻³。在q=a处径向导数跳为−a⁻²+a⁻¹≤0。所以全平面分布意义下Δf≤a⁻³。

沿用超调和圆盘均值论证，半径1/√n平滑满足

~~~math
T_nf_a(q)\le f_a(q)+\frac1{8na^3}.
~~~

径向对偶及Cauchy–Schwarz遂给

~~~math
\frac1n\sum_i f_a(q_i)\ge
T_n(a,E):=1-a+a^2/3
-\frac{\sqrt E}{n}\sqrt{\frac{1+a^{-2}}2}
-\frac1{8na^3}.
~~~

常数项和负跳项均已保留。它与此前截断逆半径测试相关，但背景量、范数和几何用途不同。

## 3. 从人口容量转到层半径总和

因为实际点q≤1，f_a(q)≥0，且q f_a(q)≤1−q、f_a≤1/a。于是每层

~~~math
m_jf_a(q_j)\le4\sqrt n(1-q_j)+1/a.
~~~

对各层求和，再用第2节，得到

~~~math
\boxed{\frac{\sum_jR_j}{n}
\le\frac K{\sqrt n}+\frac K{4an}-\frac{T_n(a,E)}4.}
~~~

这正是[相邻边凸预算](PAIR_CONVEXITY.md)所需的周长上端。若K/√n有界、E=O(n)，取a=n⁻¹ᐟ⁴，则T_n→1，故

~~~math
\frac{\sum_jR_j}{n}\le\frac K{\sqrt n}-\frac14+o(1).
~~~

## 4. 有限n接受接口

给定已证E、a及假设“层数至多K”，以本节上端估计ΣR_j，用L≥n−2K估计边数。若所得平均边长上端≤一个已认证截止c，且(n−2K)P_lower(c)>E，则整个至多K层类被排除。

[circumference.py](circumference.py)以向外平方根和有理π上端完成该过程。它也处理周长上端已为负的直接矛盾。所有能量输入必须独立证明；该接口不会从一组输入点或小驻点残差认定全局最优。

这里分别用了n²D≤E和点对损失≤E两条必要条件，**没有将两份E相加**。更精细地分配同一缺额中的Dirichlet和点对份额仍是可继续优化的方向。
