# 十点近优鞍点及其轨道邻域的排除

记 \(\Delta(z)=\prod_{i<j}|z_i-z_j|^2\)。取质量 \(\sum_i|z_i|^2=9\)，并令
\[
L(z)=\log\Delta(z)-5\sum_i|z_i|^2,\qquad
F_i(z)=\sum_{j\ne i}\frac1{\overline{z_i-z_j}}-5z_i.
\]
实梯度为 \(2F\)。与十点质量约束下的局部极大相比较，归一化给出同一个 \(L\) 的无约束局部极大条件。事实上，沿尺度 \(sz\) 有 \(L(sz)-\log\Delta(z)=90\log s-5s^2\sum_i|z_i|^2\)，在质量九时唯一取极大于 \(s=1\)。

## 代数定义与分类

仓库先前数值筛选中的次高簇对应一个对称群为 \(D_2\) 的完整驻点 \(Z_s\)。令 \(a,b,c,r>0\)、\(|x|<r\)，用
\[
Z_s=\bigl(A,-A,iB,-iB,C,-C,p+iq,-p-iq,p-iq,-p+iq\bigr),\\
A=\sqrt a,\quad B=\sqrt b,\quad C=\sqrt c,\quad
p=\sqrt{(r+x)/2},\quad q=\sqrt{(r-x)/2}
\]
表示构型。五个反向点对在平方映射下的像是
\[
a,\quad -b,\quad c,\quad x+i\sqrt{r^2-x^2},\quad x-i\sqrt{r^2-x^2}.
\]
它们的模长依次为 \(a,b,c,r,r\)。前三个像的第 \(k\) 次幂为 \(A_k=a^k\)、\(B_k=(-b)^k\)、\(C_k=c^k\)，后两个的幂和为 \(R_k\)。计算 \(R_k\) 时取 \(R_0=2,R_1=2x,R_k=2xR_{k-1}-r^2R_{k-2}\)。记 \(S_k=A_k+B_k+C_k+R_k\)、\(T_k=aA_k+bB_k+cC_k+rR_k\)。根系统为
\[
2\sum_{j=0}^{k}S_j S_{k-j}-(2k+1)S_k-10T_k=0\quad (k=0,\ldots,4).
\]
证明沿用候选构型的五变量法。设 \(w_i\) 是五个平方像，令 \(\mathcal E_i=\frac12+2w_i\sum_{j\ne i}(w_i-w_j)^{-1}-5|w_i|\)。对于代表 \(z_i\)，有 \(\mathcal E_i=z_i\overline{F_i(z)}\)。每个无序对上的差商展开说明第 \(k\) 条方程左端等于 \(2\sum_i w_i^k \mathcal E_i\)。五根互异，故 Vandermonde 可逆，全部多项式方程等价于完整驻点方程。\(k=0\) 给出 \(2a+2b+2c+4r=9\)。

Jacobian 对有理盒返回唯一非奇异实根。盒内 \(a<b<r<c\)，且 \(x\) 接近 \(-0.0729405452\)。四种平方半径的人口是 \((2,2,4,2)\)，其中内点的两个平方半径约为 \(0.1118366998\)，外带的三个值约为 \(1.0074790085,1.0962698481,1.1881445955\)。这满足已有的两内八外半径条件。内层的两个点迫使任何旋转对称保持这对点，故旋转群阶恰为二。共轭也是对称，由此得到 \(D_2\)。

在质量十下，\(\Delta(Z_s)\) 位于 \([53745013433,53745013449]\)，比原候选 \(C_{10}\) 低约 \(0.145\%\)。固定一个最外点在实轴正半轴，十九维切片上的实力 Jacobian 具有惯性 \((1,18,0)\)，依次计正、负、零特征值。标准库复核器以 384 位向外舍入构造切片矩阵的有理中心，算子范数误差小于 \(10^{-6}\)。在中心两侧分别加减 \(10^{-6}I\)，精确消元均给出同一惯性。Jacobian 另外复核这两个有理矩阵。

## 排除整个轨道邻域

若质量九的完整驻点与 \(Z_s\) 的旋转和重排轨道的欧氏距离至多 \(21/500\)，则它有严格正的二阶变分，因而不能是局部或全局极大。该结论覆盖此范围内破坏两个对称轴或半转对称性的构型。

令 \(v,w\in\mathbf C^{10}\) 是复核器内两组有理向量。记 \(\delta_{ij}=z_i-z_j\)、\(k_{ij}=v_i-v_j\)、\(\ell_{ij}=w_i-w_j\)。在二阶变分中加入一阶力的线性组合，得到
\[
G(z)=D^2L(z)[v,v]+\operatorname{Re}\sum_i\overline{w_i}F_i(z)
=\operatorname{Re}\sum_{i<j}
\left(-\frac{2k_{ij}^2}{\delta_{ij}^2}
+\frac{\ell_{ij}}{\delta_{ij}}\right)
-10\|v\|_2^2-5\operatorname{Re}\sum_i\overline{w_i}z_i.
\]
在任意完整驻点，\(G(z)=D^2L(z)[v,v]\)。数值构造先选取鞍点的正曲率方向作为 \(v\)，再用该处的力 Jacobian 选取 \(w\)，以削小 \(G\) 的线性变化。接受器直接检验有理向量所得的 \(G(Z_s)>1697/10000\) 及 \(\|\nabla G(Z_s)\|_2<13/2000\)。

取 \(\rho=21/500\)、\(\sigma=(99/70)\rho>\sqrt2\,\rho\)，令 \(\delta^s_{ij}\) 表示 \(Z_s\) 的点差。对 \(z=Z_s+e\) 且 \(\|e\|_2\le\rho\)，逐边有 \(|e_i-e_j|\le\sqrt2\,\rho\)，距离分母恒不为零。保留 Taylor 二次项的符号。令
\[
H_{ij}=\frac{\ell_{ij}\delta^s_{ij}-6k_{ij}^2}{(\delta^s_{ij})^4}.
\]
对 \(H_{ij}=A+iB\)，二次项 \(\operatorname{Re}(H_{ij}(e_i-e_j)^2)\) 对应的实边矩阵为 \(\bigl(\begin{smallmatrix}A&-B\\-B&-A\end{smallmatrix}\bigr)\)。所有边组装成二十维对称矩阵，复核器证明该矩阵严格大于 \(-74I\)。单边复解析项沿位移线段的三阶复导数之模的一半不超过
\[
U_{ij}=\frac{|24k_{ij}^2-3\ell_{ij}\delta^s_{ij}|+3|\ell_{ij}|\sigma}
                 {(|\delta^s_{ij}|-\sigma)^5}.
\]
因此余项的绝对值至多为 \(\tfrac13\sum_{i<j}U_{ij}|e_i-e_j|^3\)，再用 \(|e_i-e_j|<\sigma\) 化为二次式。复核器用有理上端替换各权重后构造实 Laplacian \(L_U\)，精确验证 \(970I-L_U\) 正定。Jacobian 独立复核二次矩阵和 Laplacian 的有理移位。球内
\[
G(Z_s+e)>\frac{1697}{10000}
-\frac{13}{2000}\rho-74\rho^2-\frac{970}{3}\sigma\rho^2
=\frac{156613}{31250000}.
\]
用旋转和重排等变性把结论传至整个轨道邻域。

## 有限数值筛选

固定种子 20261116010。8,400 个无对称约束起点使用 L-BFGS 和完整力的旋转切片 Newton 修正，全部终点的未约束梯度无穷范数低于 \(10^{-8}\)。按对数距离积的 \(10^{-7}\) 容差分簇，7,115 次达到 \(C_{10}\) 的数值簇，1,284 次达到中心加九边形的数值簇，1 次停在上述鞍点的数值簇。另从 2,450 个起点直接解驻点方程，2,372 个通过同一残差筛选。两种搜索合计得到 47 个数值簇，在保存的代表点中只有 \(C_{10}\) 和上述鞍点满足已证的两内八外半径带。数值分簇不能给出全体驻点分类。

## 复现

基础仓库来源为提交 `68f0d9ee885f8e4241a834de5887058e57a5bd31`。从仓库根目录运行：

```text
python -B -S -O ten-point-saddle-exclusion/verify_all.py
python -B ten-point-saddle-exclusion/search.py --output search-replay.json --scale 100
```

第一项只用 Python 标准库，对固定根记录独立重放多项式、导数、反矩阵、Krawczyk 像、截面惯性与邻域排除。加 `--jacobian-inputs <path>` 可输出四个供 Jacobian 复核的有理矩阵。第二项需要 NumPy 2.5.5 和 SciPy 1.18.1，属于有限浮点探索。[根证书](saddle-jacobian.json) · [惯性复核](saddle-inertia-jacobian.json) · [搜索记录](search.json) · [重放报告](verification_report.json)
