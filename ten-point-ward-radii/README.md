# 十点全局最大构型的半径禁区

2026-11-10。蒋丰灿研究项目，OpenAI Codex 协作。

在质量 \(\sum_i|x_i|^2=9\) 的归一化下，每个全局最大构型的平方半径都避开两个闭区间：

\[
[3/10,41/60]\quad\text{和}\quad[71/50,149/100].
\]

结合[缩并必要条件](../ten-point-collapse-domain/README.md)的严格上端，允许的区域为 \([0,3/10)\cup(41/60,71/50)\cup(149/100,1539/1000)\)。任一四点的平方半径和大于 \(4\cdot3133/10000\)，因此最多三个点落在最内区域。以上排除不预设全局最大点的对称性。

## 一行的三个精确平衡

记 \(u=|x_i|^2\)、\(v=|x_j|^2\)、\(t=\operatorname{Re}(x_i\bar x_j)\)、\(d=u+v-2t\)，其中 \(i\) 是标记点。对完整驻点，质量、质心、径向力依次给出

\[
\sum_{j\ne i}v=9-u,\qquad
\sum_{j\ne i}t=-u,\qquad
\sum_{j\ne i}\frac{u-t}{d}=5u.
\]

最后一个式子是 \(\sum_{j\ne i}(x_i-x_j)/d=5x_i\) 与 \(x_i\) 的实内积。令 \(h\) 为[十点全局核](../ten-point-ward-upper/README.md)的两点函数，\(f=h-\log d\)，\(B\) 为该核的总和预算。已证所有全局最大点的点对满足 \(f\ge-3/10000\)，且 \(\sum_{i<j}h\le B\)。

对每个禁区，将标记半径 \(u\) 分成有限个闭区间。每个小区间用一组固定有理数 \((a,b,c)\) 证明以下连续不等式：

\[
f(u,v,t)\ge \frac{1}{800}
 +a\big(v-1+u/9\big)+b(t+u/9)+c\big((u-t)/d-5u/9\big).
\]

在标记点的九条边上求和，三个修正项全部消失，给 \(\sum_{j\ne i}f\ge9/800\)。其余 36 对使用原全局证书，所以 \(\sum_{i<j}f\ge9/800-36(3/10000)=9/20000\)。程序验证可达有理见证值 \(V_{10}\) 满足 \(\log V_{10}>B-9/20000\)。如果全局最大构型有一条半径落在禁区，就会有 \(\log\Delta\le B-9/20000<\log V_{10}\)，矛盾。

## 精确复核

先重放[缩并必要条件](../ten-point-collapse-domain/README.md)与[全局核](../ten-point-ward-upper/README.md)，再从仓库根目录运行：

```sh
python -B -S -O ten-point-ward-radii/verify_all.py
```

接受器读取两个固定的树集合，覆盖 \(u\) 的两个完整闭区间，且不假定标记点比另一点更靠里。将不等式乘以正数 \(d\) 后，接受器用精确 Bernstein 系数与对数切线上界验证每个盒，或用[十点 Schwarz 必要域](../maximizer-domain/README.md)排除该盒。所有区间、系数和比较均为有理或整数运算。[复跑报告](verification_report.json)列出覆盖统计。这个结论尚未将最大构型定位到[局部锐界范围](../ten-point-local-basin/README.md)。
