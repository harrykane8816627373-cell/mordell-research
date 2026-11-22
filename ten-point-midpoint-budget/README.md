# 十点最大构型的五个配对中点范数小于 \(149/500\)

归一化为 \(\sum_i|z_i|^2=9\)。在任一全局最大构型中，两个内点配成一对，按圆周顺序相隔四位的八个外点配成另外四对。记五个配对中点的向量为 \(m\)。则
\[
\|m\|_2<\frac{149}{500}=0.298.
\]
外点角度证书给出的原联合上界是 \(363/1000\)。[配对定义和上游约束](../ten-point-tight-pairing/README.md)

## 逐对缺额和全局约束

对一点对 \(z,w\)，令 \(u=|z|^2,v=|w|^2,t=\operatorname{Re}(z\bar w)\)、\(d=u+v-2t\)。沿用曲率核 \(h\)，取缺额 \(G=h-\log d\)。配对的中点平方模为 \(M=(u+v+2t)/4\)。在唯一内边及四条配对外边上分别验证
\[
G\ge \lambda M+\mu(u+v-s_j)-c_j,
\qquad
\lambda=\frac{69}{500},\quad \mu=-\frac{14}{25000},
\]
其中内边使用 \(s_{\mathrm I}=6/25,c_{\mathrm I}=93/200000\)，配对外边使用 \(s_{\mathrm O}=219/100,c_{\mathrm O}=491/200000\)。检查的连续域包含每个全局最大构型的相应边。内边使用既有半径、夹角、距离及中点范围，配对外边使用细化半径和近对径约束。[曲率核](../ten-point-curved-kernel/README.md) · [内边几何](../ten-point-conditioned-inner/README.md)

五条配对边各顶点恰出现一次，\(\sum_{\text{配对边}}(u+v)=9=s_{\mathrm I}+4s_{\mathrm O}\)。因此仿射项在五条配对边求和时消去。剩余的 16 条内外边和 24 条外外边使用[角序分组证书](../ten-point-angular-upper/README.md)，合计补偿额为 \(c_{\mathrm{rest}}=0.0019528\)。在全局最大点，\(\sum_{\text{全部边}}G\le B-\log C_{10}^{(9)}\)，其中 \(C_{10}^{(9)}\) 为既有十点候选在质量九尺度下的距离积。遂得
\[
\lambda\|m\|_2^2
\le B-\log C_{10}^{(9)}+c_{\mathrm I}+4c_{\mathrm O}+c_{\mathrm{rest}}.
\]
复核程序以候选距离积的有理下端向下估计其对数，并验证右侧严格小于 \(\lambda(149/500)^2\)。

## 复核

写 \(P=d[h-\lambda M-\mu(u+v-s_j)]\)，这是有理多项式。内边的复核使用 1,145 个 Bernstein 节点及曲率对数上界。外边先在 \(t\in[-\sqrt{uv},-(31/32)\sqrt{uv}]\) 上证明 \(dP_t+2P+2d\le0\)，共 45 个节点。由于 \(P/d-\log d\) 对 \(t\) 的导数是 \((dP_t+2P+2d)/d^2\)，只需在较大的端点 \(t=-(31/32)\sqrt{uv}\) 检查外边下界。代入 \(u=r^2,v=s^2,t=-(31/32)rs\) 后，对 \(r,s\) 使用 585 个节点核验 \(P\ge d(\log d-c_{\mathrm O})\)。复核用向外舍入的 Bernstein 系数和有理对数切线上界。

先接受上游曲率核、细化半径、内边几何、近对径约束和角序分组证书，再从仓库根目录运行

    python -B -S -O ten-point-midpoint-budget/verify_all.py

[覆盖树](cover.json) · [固定输入](source_manifest.json) · [复核报告](verification_report.json)
