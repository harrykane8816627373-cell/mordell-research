# 十点候选附近消去全部配对中点

取质量 \(\sum_i|z_i|^2=9\)，把两个内点配成一对，按角序相隔四位配对八个外点。记五对为 \(m_i\pm a_i\)，其中第零对是内点对。设 \(a^*\) 是已认证的两矩形 \(D_2\) 候选的有向半差向量。允许整体旋转、重排外点对及翻转各半差符号。

**若一个全局最大构型满足 \(\|a-a^*\|_2\le 7/250\)，则它属于候选的旋转、重排轨道。** 这一充分条件仅要求控制半差。采用[五对中点预算](../ten-point-midpoint-budget/README.md)、[外点配对](../ten-point-tight-pairing/README.md)和[内点中点界](../ten-point-refined-radii/README.md)已认证的 \(\|m\|_2<149/500\)、\(|m_0|<239/2000\) 和 \(|m_i|<43/250\)（\(1\le i\le4\)），下文在整个上述范围上证明消去中点会严格增大能量，除非 \(m=0\)。前一版[加权局部证书](../ten-point-weighted-pair-basin/README.md)要求较小的中点范数 \(7/40\)。

## 保留四次项的符号

令 \(L=\log\Delta-5\sum_i|z_i|^2\)。质量 9 下的全局最大点也是 \(L\) 的无约束全局最大点。事实上，将任意非退化构型归一化到质量 9 后，再沿正比例缩放 \(s\) 时，\(L\) 的变化为 \(90\log s-45s^2+45\)，其最大值为零。

对每个 \(i<j\) 和每个符号，记 \(d=a_i^*\pm a_j^*\)、\(v=(a_i-a_i^*)\pm(a_j-a_j^*)\)、\(\eta=m_i-m_j\)。令 \(w=\eta^2/d^2\)、\(w'=\eta^2/(d+v)^2\)。对径配对的四边乘积给出

\[
L(a,m)-L(a,0)=2\sum\log|1-w'|-10\|m\|_2^2.
\]

每个边的已认证候选长度下端记作 \(\underline d\)。取 \(R=7/250\)、\(e=(10/7)R\)、\(\theta=1-e/\underline d\)。取 \(s_{ij}\) 为两点对各自中点范数上界的和，令 \(\rho=s_{ij}^2/(\underline d-e)^2\)。区间复核逐边得到 \(\rho<1/4\)，且 \(|w'|\le\rho\)。由对数及倒四次幂的级数，

\[
2\log|1-w'|\le -2\operatorname{Re}w'-\operatorname{Re}w^2+E|w|^2,
\quad
E=\theta^{-4}\left(1+\frac{2\rho}{3(1-\rho)}\right)-1.
\]

这里四次项来自 \(-\operatorname{Re}(w')^2\)，替换 \((w')^2\) 为 \(w^2\) 的误差至多 \((\theta^{-4}-1)|w|^2\)。剩余项由 \(2|w'|^3/[3(1-\rho)]\) 控制。分别使用内点和外点的中点界可减小各边的 \(\rho\)。

## 四次式的精确矩阵界

记

\[
Q_M(m)=10\|m\|_2^2+2\operatorname{Re}\sum w.
\]

既有[加权局部证书](../ten-point-weighted-pair-basin/README.md)重放正定性 \(Q_M(m)>(7/5)\|m\|_2^2\) 及一次、二次扰动上界。以 \(x=e/(73/100)\)，得到

\[
-2\operatorname{Re}\sum w'-10\|m\|_2^2
\le -\left(1-12R-\frac{26R^2(3-2x)}{(1-x)^2}\right)Q_M(m).
\]

将 \(m\) 按实部与虚部排列成十维实向量 \(p\)，并把单项式 \(p_ap_b\)（\(a\le b\)）排列成 55 维向量。设 \(Q\) 是代表 \(Q_M\) 的实矩阵。二次型 \(Q\otimes I\) 在向量 \(p\otimes p\) 上的值为 \(Q_M(m)\|m\|_2^2\)。对于非对角指标 \(a<b\)，以 \(e_a\otimes e_b+e_b\otimes e_a\) 作为对称张量的基向量。每条边的 \(w\) 是单项式向量的复线性函数。若其系数的实部和虚部向量为 \(s,t\)，则边的四次项 \(-\operatorname{Re}w^2+E|w|^2\) 对应矩阵 \((E-1)ss^T+(E+1)tt^T\)。

把这些区间矩阵求和后，复核器检查其严格小于 \((61/10)(Q\otimes I)\) 在对称张量基下的限制。检验使用向外舍入的区间和有理 \(LDL^T\) 分解。因而

\[
\sum\bigl(-\operatorname{Re}w^2+E|w|^2\bigr)
\le \frac{61}{10}Q_M(m)\|m\|_2^2.
\]

两式合并后，\(\left(1-12R-26R^2(3-2x)/(1-x)^2-(61/10)(149/500)^2\right)>1/20\)，所以 \(L(a,m)\le L(a,0)-Q_M(m)/20\)。在全局最大构型上，\(m\ne0\) 会与全局性矛盾。取 \(m=0\) 后，半差落入已经复核的 \(7/250\) 局部锐界，等号仅为候选轨道。仍需证明任一全局最大构型的半差能进入这个范围，或另行排除其补集。

## 用候选的外角作为固定参考

先把内点半差旋到正实轴。对固定的内点、八个外半径和平均外角，取参考外角 \(\widehat\theta\) 为候选外角依角序平移到同一平均值。候选的三类相邻外角均严格满足新证书的范围，所以参考角属于[相邻间隙凸域](../ten-point-outer-gap/README.md)。[收紧后的半径比较](../ten-point-inner-midpoint/README.md)给出这个凸域上的固定平均角负 Hessian 下界 \(161/50\)。令 \(G=\nabla_\theta\log\Delta(\widehat\theta)\) 在平均值零子空间中的投影，则对全局最大点的真实外角 \(\theta\) 有

\[
\|\theta-\widehat\theta\|_2\le(50/161)\|G\|_2.
\]

在原半径上放置参考角，得到参考半差向量 \(a^{\rm ref}\)，内点半差保留原值。四个外点对取半差的映射具有欧氏算子范数 \(1/\sqrt2\)。弧长支配弦长，且外平方半径小于 \(1.271\)，故 \(\|a-a^{\rm ref}\|_2\le\sqrt{1.271/2}\,\|\theta-\widehat\theta\|_2\le(4/5)\|\theta-\widehat\theta\|_2\)。因此只要验证

\[
\|a^{\rm ref}-a^*\|_2+(40/161)\|G\|_2\le7/250,
\]

即可应用本页的半差判据。式中只需两个内点中点坐标、内点半差长度、七个独立外半径及平均外角，合计十一个连续参数。剩余域可以用这个显式判据逐盒处理。

## 复核

接受上游证书后，从仓库根目录运行

```sh
python -B -S -O ten-point-quartic-midpoints/verify_all.py
```

程序仅使用 Python 标准库。另把有理下界乘以 \(2^{10}\)，逐项向下取整，并减去 \(55I\) 以补偿取整误差的算子范数。Jacobian 的 `matrix.inertia.compute` 独立复核该 55 阶整数下界的惯性为 \((55,0,0)\)。可加 `--export-jacobian PATH` 重新导出同一个整数矩阵输入。 [Jacobian 结果](jacobian_inertia.json) · [固定输入](source_manifest.json) · [复核报告](verification_report.json)
