# 9点：中心加正八边形

[English](README.md) · [返回首页](../../README.zh-CN.md)

## 定理

在$`\sum|z_i|^2=9`$下，

```math
M_9=C_9=9^{36}/8^{28}.
```

全部等号构型为

```math
\{0\}\cup\left\{\sqrt{\frac{9}{8}}e^{i\theta+2\pi i k/8}:k=0,\ldots,7\right\},
```

允许整体旋转和重标号，没有额外的平移自由度。

## 证明

最大值存在；任意最大点必居中、无碰撞，并满足完整一阶与二阶必要条件。变换$`x_i=\sqrt{(n-1)/n}z_i`$后，这些条件把全部有关点对放入明确证明的有界域。非对称构型和二阶退化情形均保留在该域中。

用精确正定核在整个必要域上构造点对距离平方对数的上界。取等邻域采用严格局部正性，其余区域采用逐步向外舍入的Bernstein覆盖。求和得到

```math
\log\Delta(x)\le(n-1)\log(n-1)-\tfrac12\sum_{i,j}K_{ij}\le(n-1)\log(n-1).
```

恢复尺度得到所述上界，给定构型达到上界。正定项取等迫使平方半径恰为一个0与$`n-1`$个1；Newton恒等式再迫使非零点组成正$`(n-1)`$边形。

计算前提已有未决盒数为0的完整证书。默认审计检查分发文件与该记录；`check`重算精确代数、正定性和局部取等邻域；`replay`另行重跑全域覆盖。这里是普通计算机辅助证明，没有将其称为整篇证明助理形式化。

[完整证明](../../verification/n9/PROOF_ZH.txt) · [最终证书](../../verification/n9/FINAL_CERTIFICATE_AUDIT.json)

```sh
python verify.py --case n9
python verify.py --level check --case n9
```
