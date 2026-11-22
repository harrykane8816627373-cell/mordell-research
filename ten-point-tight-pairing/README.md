# 十点外点近对径配对的九阶奇矩界

取总平方半径为 \(9\)。对任意十点全局最大构型，依逆时针方向记八个外点为 \(x_1,\ldots,x_8\)。每组循环配对 \((x_i,x_{i+4})\) 满足
\[
\frac{\operatorname{Re}(x_i\bar x_{i+4})}{|x_i|\,|x_{i+4}|}<-\frac{31}{32},
\qquad |x_i-x_{i+4}|^2>\frac{10647}{3200},
\qquad \left|\frac{x_i+x_{i+4}}2\right|^2<\frac{117}{4000}.
\]
其夹角偏离正对径小于 \(2\arcsin(1/8)\approx14.36^\circ\)。把两个内点配成第五对，五个中点构成的复向量 \(m\) 满足 \(\|m\|_2<363/1000\)。旧证书的角度界约为 \(24.02^\circ\)，每个外点配对中点平方模小于 \(73/1250\)。[旧配对证书](../ten-point-outer-moments/README.md)

## 九阶联合投影

记 \(u=|x|^2\)、\(H=\{1,3,5,7,9\}\)。在曲率上界中，令 \(G_k\ge0\) 为第 \(k\) 个 Gram 罚项。核的无序点对和为 \(\sum_{i<j}h_{ij}=B-\tfrac12\sum_k G_k\)。内内、内外、外外分别有 \(1,16,28\) 对，其缺额下界分别为 \(-c_{II},-c_{IO},-c_{OO}\)。设 \(V\) 是可行整数见证的归一化距离积，因全局最大构型的距离积至少为 \(V\)，得到
\[
\sum_k G_k\le 2(B-\log V+c_{II}+16c_{IO}+28c_{OO})<\frac3{625}.
\]
[曲率核](../ten-point-curved-kernel/README.md)

取 \(F(\theta)=\sum_{k\in H}w_k\cos(k\theta)\)，其中 \(w_1=1\)。为每个 \(k\) 选取五次径向多项式 \(g_k\)，使 \(u^{k/2}g_k(u)\) 逼近 \(w_k I_{\mathrm{out}}(u)\)。写 \(R=\sqrt u\)、\(c=\cos(\theta-\omega)\)。对每个径向带和两个误差符号分别作 Bernstein 覆盖，验证
\[
\left|\sum_{k\in H}\bigl(R^k g_k(R^2)-w_k I_{\mathrm{out}}(R^2)\bigr)T_k(c)\right|
< e_{\mathrm{band}}+s(R^2-9/10).
\]
这里 \(T_k\) 是 Chebyshev 多项式。由于所有 \(T_k\) 都是奇函数，只需验证 \(c\in[0,1]\)；求和时质量守恒消去仿射斜率。[固定有理系数](projection.json)

设 \(y_k\) 是 \(g_k\) 在第 \(k\) 个变换后特征中的系数，并用一阶、三阶 Ward 恒等式吸收允许的零矩方向。应用 Cauchy–Schwarz 后，对任意方向 \(\omega\) 有
\[
\left|\sum_{\mathrm{out}}F(\theta_i-\omega)\right|
<2e_{\mathrm{in}}+8e_{\mathrm{out}}+
\sqrt{\frac3{625}\sum_{k\in H}y_k^\mathsf{T}Q_k^{-1}y_k}
<\frac{21}{40}.
\]
复核器重算五组有理 Gram 方程的解和连续近似域的 52 个 Bernstein 节点。[Ward 恒等式](../stationary-kernel/README.md)

## 五点短弧的排除

设 \(\alpha=\arccos(106/125)\)、\(\beta=\arccos(1/8)\)。外点方向的两两最短夹角大于 \(\alpha\)。选定的投影函数满足 \(F''(\theta)=-\cos\theta\,P(\cos^2\theta)\)，有理四次式 \(P\) 在 \([0,1]\) 上正。Jacobian 在该闭区间上的 Sturm 零点计数为零，复核器还用 13 个 Bernstein 节点独立验证正性。由偶性、凹性及 \(F(\pm\pi/2)=0\)，\(F\) 在向前半圆非负，并在正半区间递减。[Sturm 结果](jacobian_concavity.json) · [外角分离](../ten-point-outer-moments/README.md)

若长度为 \(2\beta\) 的闭弧包含五个方向，把弧旋至 \([-\beta,\beta]\)。从四个内部间隙各扣除 \(\alpha\)，剩余长度在六个间隙间组成单纯形。五个投影和是凹函数，其最小值可在六个顶点上取得。第 \(j\) 个顶点把 \(j\) 个方向以间隔 \(\alpha\) 靠左排列，其余 \(5-j\) 个靠右排列。其余三个方向的投影和至少为 \(-F(0)-2F(\alpha)\)。为证此式，把后三个方向反向，利用 \(F(\theta+\pi)=-F(\theta)\)。若有一个投影非正，反向投影和至多为 \(2F(0)<F(0)+2F(\alpha)\)。若三个投影皆正，最宽的一对方向相差至少 \(2\alpha\)，由偶性、凹性和单调性，该对投影和至多为 \(2F(\alpha)\)。

轮流检查六个顶点。用 \(\cos(\beta-j\alpha)\) 的有理区间计算得出，五个投影和的最小值减去 \(F(0)+2F(\alpha)\) 仍大于 \(21/40\)，与联合投影界矛盾。按 \(\theta_{i+8}=\theta_i+2\pi\) 延拓圆周次序，得到
\[
2\beta<\theta_{i+4}-\theta_i<2\pi-2\beta,
\qquad
\cos(\theta_{i+4}-\theta_i)<2\cos^2\beta-1=-\frac{31}{32}.
\]
外带平方半径在 \((169/200,164/125)\)。平方点距大于 \(2(169/200)(1+31/32)\)。若 \(u,v\) 为配对的平方半径，中点平方模小于 \((u+v-2(31/32)\sqrt{uv})/4\)。该上界对 \(u,v\) 分别凸，四个端点给出严格上界 \(117/4000\)。再加入内点中点平方模上界 \(57121/4000000\)，得到五个中点的联合范数界。[细化半径带和内点中点界](../ten-point-refined-radii/README.md)

## 复核

先按复现指南接受上游曲率核、细化半径和外角分离证书。从仓库根目录运行 `python -B -S -O ten-point-tight-pairing/verify_all.py`。复核只用 Python 标准库。[复现指南](../docs/REPRODUCING.md) · [固定输入](source_manifest.json) · [复核报告](verification_report.json)
