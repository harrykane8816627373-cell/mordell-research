# 半转分支的两组纯角度邻域

取质量九。假设一个全局最大构型有半转对称。将五个反向点对平方，旋转使内像位于正实轴。比较四个外像的辐角，并按循环顺序匹配。对四个角分别取最短角差的绝对值。

| 外角参照 | 四个平方像各自的角差上界（弧度） | 结论 |
|---|---:|---|
| 两矩形候选 | \(11/500\) | 最大构型属于候选轨道 |
| 近优鞍点 | \(3/125\) | 排除全局最大构型 |

半径可以取整个[全局必要径向域](../ten-point-half-turn-radial-concavity/README.md)。验证将实际半径从该域定位到两个有理宽度的局部盒，再使用先前认证的局部判据。

## 从角度定位半径

设 \(s_i=\log r_i\)，其中 \(r_i\) 是平方像的模。固定角向量 \(\theta\) 后令 \(F(s,\theta)=\log\Delta-10\sum_i r_i\)。在参照构型 \((s_*,\theta_*)\) 处设 \(P=(-F_{ss})^{-1}(s_*,\theta_*)\)，该逆矩阵存在。对每组满足上表的角，考虑
\[
Q_\theta(h)=h+P F_s(s_*+h,\theta).
\]
若盒 \(X=\{|h_i|\le a_i\}\) 上 \(Q_\theta(X)\subset\operatorname{int}X\)，Brouwer 不动点给出一个局部径向驻点。程序验证整个局部盒落在已证严格凹的径向域内。径向唯一性因此把任何真实全局最大点在该组角处的半径定位到 \(s_*+X\)。

按角变动和径向变动分别积分，从 \(F_s(s_*,\theta_*)=0\) 得到
\[
|Q_\theta(h)_i|\le
\sum_j |(P F_{s\theta})_{ij}|\, b_j+
\sum_j |(I+P F_{ss})_{ij}|\,a_j .
\]
第一项在固定 \(s_*\) 的角盒取上界，第二项在角径直积盒取上界。这里 \(b_0=0\) 且四个外像的 \(b_j\) 取各自容差。若 \(q=s_i-s_j,\ \psi=\theta_i-\theta_j,\ D=\cosh q-\cos\psi\)，每条边在两类矩阵中分别有系数 \(-2\sinh q\sin\psi/D^2\) 和 \(2(1-\cos\psi\cosh q)/D^2\)，乘相应的秩一差分矩阵。径向矩阵还加 \(-10\operatorname{diag}(r_i)\)。验证器将这些权重与区间 \(P\) 合成，再逐行检查上式严格小于 \(a_i\)。

将平方像逆映射回原点对时，可以选取与参照角差不超过 \(b_i/2\) 的代表。每个代表的位移平方至多为 \(r_{*,i}[(e^{a_i/2}-1)^2+e^{a_i/2}b_i^2/4]\)。将五项相加，候选处的半差范数小于 \(7/250\)。鞍点处将平方范数加倍，得到整个构型的位移小于 \(21/500\)。[候选局部锐界](../ten-point-quartic-midpoints/README.md) · [鞍点轨道排除](../ten-point-saddle-exclusion/README.md)

先接受五径向严格凹性，再从仓库根运行 <code>python -B -S -O ten-point-half-turn-angle-neighborhoods/verify_all.py</code>。这一检查只需 Python 标准库。[固定输入](source_manifest.json) · [复核报告](verification_report.json)
