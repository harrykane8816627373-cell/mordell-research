# 十点半转对称分支的六变量方程

本节假设一个全局最大构型具有半转对称。取总平方半径 9，旋转后内点为 \(\pm\sqrt a\)。每一对的平方映射为一个点，五个像记作 \(w_0=a,w_1,\ldots,w_4\)，并写 \(w_j=r_j u_j\)、\(|u_j|=1\)。[内边范围](../ten-point-conditioned-inner/README.md)给 \(41/400<a<11/80\)，[逐点径向范围](../ten-point-angular-radii/README.md)控制四个外部 \(r_j\)。

把四个外像按正向角序编号。由于原构型具有半转对称，四个相邻外像角差恰为原外角差的两倍，且都小于 \(\pi\)。用[相邻外角约束](../ten-point-outer-gap/README.md)选最短的外像间隙作为 \(w_1\) 到 \(w_2\) 的间隙，记其半角余弦为 \(c\)。最短间隙至多 \(\pi/2\)，故 \(1/\sqrt2\le c<39/50\)。以两张有理参数图写

\[
u_1=\sigma\frac{1-t^2+2it}{1+t^2},\quad |t|\le1,\ \sigma=\pm1,
\qquad
u_2=u_1\left(2c^2-1+2ic\sqrt{1-c^2}\right).
\]

## 精确消去后两个角

质量恒等式给 \(r_4=9/2-a-r_1-r_2-r_3\)。令 \(b_j=(17-10r_j)r_j\)、\(b_0=(17-10a)a\)。二阶 Ward 恒等式给

\[
b_3u_3+b_4u_4=V,
\qquad V=-b_0-b_1u_1-b_2u_2.
\]

由 \(b_j>5\)、\(0<b_0<9/4\) 和最短间隙至多 \(\pi/2\)，有 \(|V|>\sqrt{50}-9/4>24/5\)。因此可用两圆交点公式表示

\[
A=\frac{|V|^2+b_3^2-b_4^2}{2|V|^2},\qquad
B=\sqrt{b_3^2/|V|^2-A^2},\qquad
b_3u_3=(A-iB)V,\quad b_4u_4=V-b_3u_3.
\]

符号由角序固定。该交点满足 \(\operatorname{Im}(\overline{b_3u_3}\,b_4u_4)=B|V|^2>0\)。由加倍角差的余弦在 \((2(59/100)^2-1,2(39/50)^2-1)\) 内，正弦严格大于 \(19/20\)。又有 \(b_j\le289/40\)，故 \(B>1/9\)。在所有可能的半转对称全局最大构型上，消元公式的分母和根式均远离零。

## 方程和数值搜索

余下连续参数为 \((a,r_1,r_2,r_3,c,t)\)，以及两张离散参数图。设 \(S_k=\sum_{j=0}^4 w_j^k\)、\(T_k=\sum_{j=0}^4 |w_j|w_j^k\)。完整力等价于

\[
2\sum_{j=0}^k S_jS_{k-j}-(2k+1)S_k-10T_k=0,\qquad k=0,\dots,4.
\]

将差商按五个平方像展开，上式是逐点力乘以 Vandermonde 矩阵所得。前两条已经用于质量和 Ward 消元，后三条是六个实方程。可用 \(\Delta=4^5\prod_j |w_j|\prod_{i<j}|w_i-w_j|^4\) 作目标值筛选。一张完整的六变量根覆盖仍待完成。

另外，仓库中[近优鞍点的正二阶变分](../ten-point-saddle-exclusion/README.md)使用的有理方向在每一对上都有相反的分量，因此保持半转对称。该鞍点在半转对称的变分类中也无法取局部极大。

[有限数值搜索](search-results.json)从 26,000 个半转对称起点出发，同时允许失去反射。24,211 个终点通过完整力残差门槛，按对数乘积分成八个数值簇。已有全局分带内只观察到候选与近优鞍点两簇。该记录用于选择区间覆盖的起点。

从仓库根目录运行 `python -B -S -O ten-point-half-turn-chart/verify_all.py` 可复核精确常数及已认证鞍点的方向。另可用 `python -B ten-point-half-turn-chart/search.py --starts 26000 --optimize-first 2600 --seed 1011202616 --output PATH` 重跑探索性搜索，需 NumPy 和 SciPy。[复核报告](verification_report.json)
