# 十点全局最大构型的内边约束

取总平方半径为 9。记两个内点为 \(x,y\)，内点中点为 \(m=(x+y)/2\)，内边平方距离为 \(d=|x-y|^2\)。每个全局最大构型满足

\[
0.43<d<0.53,\qquad
\frac{\operatorname{Re}(x\bar y)}{|x|\,|y|}<-0.87,\qquad
\bigl|\operatorname{Re}(m\overline{x-y})\bigr|/|x-y|<0.056.
\]

## 用其余 44 条边控制唯一内边

使用[曲率核](../ten-point-curved-kernel/README.md)，记缺额 \(G=h-\log d\)。全局最大构型的总缺额至多为 \(B-\log C_{10}^{(9)}\)。对唯一内边之外的 16 条内外边及 28 条外外边，应用[角序不同类别的缺额下界](../ten-point-angular-upper/README.md)，得到内边严格满足 \(G<1073/500000\)。

令 \(u=|x|^2,v=|y|^2,t=\operatorname{Re}(x\bar y)\)。以[新内半径带](../ten-point-angular-radii/README.md)和原有内积区间为起始盒。对每个子盒按以下规则验收。可用既有几何约束或[内点中点界](../ten-point-inner-midpoint/README.md)排除整个盒。也可直接核对欲证距离、夹角及轴上投影。投影的模为 \(|u-v|/(2\sqrt d)\)。其余验收盒使用已向外舍入的 Bernstein 系数及对数上界推出 \(G\ge1073/500000\)，与必要条件矛盾。

## 内点半差平方变量

令 \(w=d/4\)，并固定外点、内点中点及连线方向。在新范围 \([43/400,53/400]\) 上，沿用[内点中点证书的割线估计](../ten-point-inner-midpoint/README.md)。以 \(W=53/400\) 为平方半差上界，仍有 \(|m|<87/1000\)，故函数 \(g(u)= (\sqrt u-87/1000)^2-W\) 保持正值。按外带两个端点及外点平方半径总和界计算，得到 \(F''(w)<-53/2\)。给定参考值 \(w_0\)，切片上界为 \(F(w_0)+F'(w_0)^2/53\)。

运行 `python -B -S -O ten-point-inner-edge/verify_all.py` 复核。[固定输入](source_manifest.json) · [单盒覆盖](cover.txt) · [结果](verification_report.json)
