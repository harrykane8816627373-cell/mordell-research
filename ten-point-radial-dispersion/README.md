# 十点最优构型的外点径向离散

2026-11-16。蒋丰灿研究项目，OpenAI Codex 协作。

取总平方半径为 9。设任一全局最大构型的八个外点为 \(w_i\)，令 \(r_i=|w_i|^2\)，\(v=\sum_i r_i/8\)。那么
\[
\sum_i|r_i-v|>\frac1{25}.
\]
因此至少一个平方半径与均值相差超过 \(1/200\)，平方偏差和大于 \(1/5000\)。若正偏差有 \(k\) 个，正、负偏差绝对值的和各为总绝对偏差的一半。由 \(1/k+1/(8-k)\ge1/2\)，八个外点的平方半径跨度严格大于 \(1/100\)。所用上游条件适用于全部全局最大构型。[径向带与中点证书](../ten-point-refined-radii/README.md) · [内边距离证书](../ten-point-conditioned-inner/README.md)

## 驻点方程与外点多项式递推

反设 \(E=\sum_i|e_i|\le\varepsilon=1/25\)，其中 \(e_i=r_i-v\)。因 \(\sum e_i=0\)，每个 \(|e_i|\le d_0=\varepsilon/2\)。记两内点为 \(x,y\)，令 \(S=x+y,T=xy\)，\(a=(x-y)/2\)。旋转使 \(a\) 为正实数，设 \(u=a^2\in(41/400,11/80)\)。质量等式 \(8v+2u+|S|^2/2=9\) 及 \(|S|^2<1/10\) 给出 \(347/320<v<11/10\)。从而 \(|w_i|<R_0=53/50\)。

记 \(u_1=|x|^2,u_2=|y|^2\)。质心和一阶 Ward 矩给
\[
\Big(v-\frac{u_1+u_2}{2}\Big)S-(u_1-u_2)a=\sum_i e_i w_i.
\]
因 \(u_1-u_2=2\operatorname{Re}(S\bar a)\)，左端的模至少为 \((v-3u-|S|^2/4)|S|\)，右端至多为 \(ER_0\)。有理比较给出 \(|S|<s=5\varepsilon/3=1/15\)。再代入质量等式得到 \(v\ge v_-=(9-2(11/80)-s^2/2)/8\)。此外 \(T=-u+S^2/4\)，故 \(|T|\le T_+=11/80+s^2/4\)。[驻点矩恒等式](../stationary-kernel/README.md)

写 \(P(z)=z^2-Sz+T\)，\(Q(z)=\prod_i(z-w_i)=\sum_j q_jz^j\)，则 \(q_8=1,q_7=S\)。对外点，完整力方程为
\[
\frac{Q''(w_i)}{2Q'(w_i)}+\frac{P'(w_i)}{P(w_i)}=5\bar w_i.
\]
乘以 \(2w_i P(w_i)Q'(w_i)\)。将 \(zP Q''+(2zP'-10vP)Q'\) 除以 \(Q\)，由最高两项得到商 \((88-80v)z+90(v-1)S\)。余式 \(\mathcal R\) 次数至多七，其在 \(w_i\) 处的值等于 \(10e_i P(w_i)Q'(w_i)\)。插值得到
\[
\mathcal R(z)=10\sum_i e_iP(w_i)\frac{Q(z)}{z-w_i}.
\]

令 \(A_j=(j-9)(j+10-10v)\)。比较 \(z^j\) 系数给出
\[
[z^j]\mathcal R=A_j(q_{j-1}-S q_j)+(j+1)(j-10v)Tq_{j+1}.
\]
设 \(U_k\) 为 \(|q_k|\) 的上界。由综合除法，\(Q(z)/(z-w_i)\) 的 \(z^j\) 系数为 \(\sum_{k=j+1}^8 q_kw_i^{k-j-1}\)。所以余式系数的上界
\[
E_j=10\varepsilon(R_0^2+sR_0+T_+)\sum_{k=j+1}^8 U_kR_0^{k-j-1}
\]
只用到已递推过的系数。从 \(U_8=1,U_7=s\) 开始，依次对 \(j=7,6,5,4,3\) 取
\[
U_{j-1}=sU_j+\frac{(j+1)(11-j)}{(9-j)(j-1)}T_+U_{j+1}+\frac{E_j}{(9-j)(j-1)}.
\]
同一个 \(j=7\) 的式子及 \(|T|\ge41/400-s^2/4\) 给 \(|q_6|\) 的下界
\[
L_6=\frac{8(10v_--7)}{2(17-10v_-)}(41/400-s^2/4)-s^2-E_7/12.
\]

## 近似同圆的系数互反关系

由基本对称式，\(q_6=\sum_{i<j}w_iw_j\)，且
\[
\frac{q_2}{q_0}=\sum_{i<j}\frac{1}{w_iw_j}=
\sum_{i<j}\frac{\bar w_i\bar w_j}{(v+e_i)(v+e_j)}.
\]
将倒数分母沿 \(s_0\mapsto(v+s_0e_i)(v+s_0e_j)\) 在 \(s_0=0\) 展开至一次。求和时，一次项的分子为
\[
\sum_{i<j}(e_i+e_j)\bar w_i\bar w_j=-\bar S\sum_i e_i\bar w_i-\sum_i e_i\bar w_i^2.
\]
因此一次误差至多为 \(\varepsilon(sR_0+R_0^2)/v^3\)。二阶余项的分子至多为 \(R_0^2\sum_{i<j}(e_i^2+|e_i e_j|+e_j^2)\)。这个和等于 \(\frac{13}{2}\sum_i e_i^2+\frac12 E^2\)，至多为 \(15\varepsilon^2/4\)，因为 \(\sum_i e_i^2\le d_0 E\)。

再用 \(|q_0|\ge(v-d_0)^4\)，得到
\[
|q_2|\ge (v-d_0)^4\left(\frac{L_6}{v^2}
-\frac{\varepsilon(sR_0+R_0^2)}{v^3}
-\frac{15R_0^2\varepsilon^2}{4(v-d_0)^4}\right).
\]
用 \(v_-\le v<11/10\) 做有理比较，括号下端为正。所得下界大于 \(173/1000\)，而上面的驻点递推给出 \(U_2<171/1000\)，矛盾。

按[复现指南](../docs/REPRODUCING.md)接受上游内点约束，再运行 `python -B -S -O ten-point-radial-dispersion/verify_all.py`。[有理比较报告](verification_report.json)
