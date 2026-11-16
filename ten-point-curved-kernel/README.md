# 十点分带核的曲率上界

设 $M_{10}$ 为十点最大距离积，$V_{10}$ 为整数构型见证的距离积，$C_{10}$ 为认证的严格局部极大构型的距离积。统一使用同一质量时，本证书证明

\[
 V_{10}< C_{10}\le M_{10}<\frac{1253}{1250}V_{10},\qquad
 \frac{C_{10}}{M_{10}}>\frac{1250}{1253}.
\]

证明使用质量 9。任一全局最大构型有两个平方半径在 $I=(27/500,217/1000)$，八个在 $O=(37/50,139/100)$。唯一内点对的点积在 $(-27/200,-1/10)$，平方距离小于 $3/5$，中点平方模小于 $1/40$，夹角余弦小于 $-2/3$。 [分带与内点对定理](../ten-point-banded-kernel-refined/README.md)

写 $u=|z_i|^2,v=|z_j|^2,t=\operatorname{Re}(z_i\bar z_j)$ 和 $d=u+v-2t$。证书重新拟合基础核的 Gram 矩阵，并用有理 $LDL^T$ 分解验证十个块全部严格正定。对标记点 $i$，内点指示函数 $\chi$、质量、质心和径向力分别给出四个零和行条件。

\[
 \sum_{j\ne i}\left(\chi(v)-\frac{2-\chi(u)}9\right)=0,\quad
 \sum_{j\ne i}\left(v-\frac{9-u}9\right)=0,\quad
 \sum_{j\ne i}\left(t+\frac u9\right)=0,\quad
 \sum_{j\ne i}\left(\frac{u-t}{d}-\frac{5u}9\right)=0.
\]

用仅依赖 $u$ 的函数乘每条恒等式并对称化。另对 $k=1$ 取 $\ell(u)=1,u$，对 $k=2,3$ 分别取 $\ell(u)=u-17/10,u-8/5$，有 $L_k=\sum_i z_i^k\ell(u_i)=0$。取 $g(u)=\chi(u)(u-1)^j$，则核

\[
 J(z_i,z_j)=\operatorname{Re}((z_i\bar z_j)^k)
 \big(\ell(u)g(v)+g(u)\ell(v)\big)
\]

的有序点对和为 $2\operatorname{Re}(L_k\overline{\sum_i z_i^k g(u_i)})=0$。在两点函数中引入 $-J_{ij}-(J_{ii}+J_{jj})/18$，保持无序总和不变。这允许用内带径向函数调整前三阶角谐波。基础核及其他修正沿用原来的求和公式。 [基础 Ward 核](../ten-point-ward-upper/README.md) · [完整驻点的 Ward 恒等式](../stationary-kernel/README.md)

总和界为 $\sum_{i<j}h_{ij}\le B$。内内、内外、外外分别有 $1,16,28$ 对。精确覆盖给出 $h-\log d\ge-c_{\rm class}$，三类松弛分别为 $3/50000,53/1000000,53/1000000$。所以质量 9 下 $\log M_{10}^{(9)}\le B+c_{II}+16c_{IO}+28c_{OO}$。复核器用整数见证的有理对数下端比较该上界，并独立验证实际候选的距离积严格大于整数见证。

在一个分割盒中，物理 $d$ 位于 $[d_{\rm lo},d_{\rm hi}]$，其中 $d_{\rm lo}\ge2/5$。取有理中点 $a$。由于 $(\log d)''=-d^{-2}\le-d_{\rm hi}^{-2}$，Taylor 积分余项给出

\[
 \log d\le \log a+\frac{d-a}{a}-\frac{(d-a)^2}{2d_{\rm hi}^2}.
\]

乘正数 $d$ 后用三次阈值检验多项式 $dh$，并保留普通切线作为另一种接受叶。对数值 $\log a$ 向上取有理端点。三次阈值在盒坐标中展开，每个单项式用下降阶乘比转换为 Bernstein 系数。阈值系数向上舍入，$dh$ 的初始系数和中点分割向下舍入。程序另用直接有理张量转换交叉核对打包的曲率阈值。

五个根盒覆盖内内负点积、内外和外外的两种点积符号。逐盒检查物理 Gram 条件、删除 Schwarz 条件、缩并必要条件以及内点对的几何约束。冻结树有 **24,681** 个节点，其中 **7,374** 个接受叶使用曲率上界，未决盒为零。 [缩并域](../ten-point-collapse-domain/README.md) · [重放报告](verification_report.json)

先重放上游分带与内点对定理，再从仓库根目录运行标准库接受器。

    python -B -S -O ten-point-curved-kernel/verify_all.py --progress

加 --record-report 可刷新报告。数值拟合只用于发现核，接受器对存储系数和连续分区作有理数和任意精度整数复核。十点锐值及最大构型分类仍待证明。
