# 半转对称最大构型的外点径向离散

取总平方半径为 9。设一个全局最大构型具有半转对称，记四个外点对的平方半径为 \(r_j\)，令 \(v=\frac14\sum_j r_j\)。则

\[
\sum_{j=1}^4|r_j-v|>47/1000,
\qquad \max_j r_j-\min_j r_j>47/2000.
\]

每个平方半径对应两个外点。因此八个外点对均值的总绝对偏差大于 \(47/500\)。这适用于半转分支，可用于[半转区间覆盖](../ten-point-half-turn-chart/README.md)。对于任意全局最大构型，仍可使用[一般外点径向离散界](../ten-point-radial-dispersion/README.md)。

## 外点的平方像多项式

令内点平方像为 \(a\)，由[内边证书](../ten-point-inner-edge/README.md)有 \(a\in(43/400,53/400)\)。令四个外点平方像为 \(w_j\)，\(|w_j|=r_j\)，且 \(Q(z)=\prod_{j=1}^4(z-w_j)=\sum_{k=0}^4 q_k z^k\)。写 \(e_j=r_j-v\)，有 \(v=(9/2-a)/4\)。在半转约化力方程中，任一外根满足

\[
w_j\frac{Q''(w_j)}{Q'(w_j)}+\frac{2w_j}{w_j-a}+\frac12-5r_j=0.
\]

乘以 \((w_j-a)Q'(w_j)\) 后作多项式除法，得到余式

\[
\mathcal R=z(z-a)Q''+[2z+(1/2-5v)(z-a)]Q'-(22-20v)Q
=5\sum_j e_j(w_j-a)\frac{Q(z)}{z-w_j}.
\]

右式由四个外根上的 Lagrange 插值得到。比较 \(z^k\) 系数给出

\[
[z^k]\mathcal R=(k-4)(k+11/2-5v)q_k-a(k+1)(k+1/2-5v)q_{k+1}.
\]

反设 \(\sum |e_j|\le\varepsilon=47/1000\)，则 \(|e_j|\le\varepsilon/2\)。令 \(R=v_++\varepsilon/2\)，其中 \(v_+\) 为 \(v\) 的上界。二阶 Ward 等式在平方像中给出

\[
(17-10v)q_3+10\sum_j e_j w_j=(17-10a)a.
\]

先由该式控制 \(|q_3|\)，再对 \(k=2,1\) 顺次控制 \(|q_k|\)。插值系数用 \([z^k]Q(z)/(z-w_j)=\sum_{\ell=k+1}^4q_\ell w_j^{\ell-k-1}\) 估计。

## 在互反式中代入 Ward 等式

展开 \(q_1/q_0=-\sum_j\overline{w_j}/(v+e_j)^2\) 的一次项，再代入上述 Ward 等式，得到

\[
\frac{q_1}{q_0}=\frac{15v-17}{5v^3}\overline{q_3}+\frac{(17-10a)a}{5v^3}+\mathcal E,
\qquad |\mathcal E|\le\frac{3R\varepsilon^2}{2(v_- -\varepsilon/2)^4}.
\]

这里用到 \(\sum e_j^2\le\varepsilon^2/2\)。又有 \(|q_0|\ge(v_--\varepsilon/2)^4\)。据此得到 \(|q_1|\) 的有理下界，与余式递推的上界相差严格超过 \(1/1000\)。这排除了四个外半径接近同一半径的情形。若四个偏差中有 \(k\) 个正值，正偏差之和至多为 \(k(4-k)(\max r_j-\min r_j)/4\)。用 \(k(4-k)\le4\) 即得跨度结论。

从仓库根目录运行 `python -B -S -O ten-point-half-turn-dispersion/verify_all.py` 复核有理比较。[固定输入](source_manifest.json) · [复核报告](verification_report.json)
