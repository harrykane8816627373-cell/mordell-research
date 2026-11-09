# Mordell 距离乘积研究

蒋丰灿（Fengcan Jiang）研究项目，推导、精确证书与代码包含 OpenAI Codex 协作。更新于 **2026-11-09**。

**最新成果 · 十二点全局锐界与唯一性**

[英文论文 PDF](twelve-point-sharp-global/manuscript/mordell-n12-proof.pdf) · [LaTeX 源文件](twelve-point-sharp-global/manuscript/mordell-n12-proof.tex) · [证明与证书](twelve-point-sharp-global/README.md) · [复现指南](docs/REPRODUCING.md)

## 十二点的最大值与构型

在质量归一化下，研究
```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\prod_{1\le i<j\le n}|z_i-z_j|^2.
```

**M₁₂ = C₁₂ ≈ 2.599464093232166 × 10¹⁴。** 最大构型在整体旋转和点重排下唯一。三层人口为 **3、3、6**，内、中层各为正三角形，外层为两条共圆正三角形。证明用 Ward 缺额锁定径向分带，再使核在十四类点对上精确取等。最终 **728,050 个节点**的连续覆盖补偿为零，未决为零。 [论文与复现入口](twelve-point-sharp-global/manuscript/README.md) · [精确复核报告](twelve-point-sharp-global/verification_report.json)

## 当前结果

|范围|结果与入口|
|---|---|
|2 至 5 点|Mₙ = nⁿ，等号为单位正 n 边形。 [五点与既有证明](earlier-results/README.md)|
|6、7 点|M₆ = 6¹⁵/5¹⁰，M₇ = 7²¹/6¹⁵，等号为原点加相应正多边形。相关证明引用 Faugère–Svarts 物理驻点分类。 [证明与依赖审计](earlier-results/six_seven_coverage_audit.md)|
|8 点|M₈ = 8²⁸/7²¹，等号为原点加半径 √(8/7) 的正七边形。 [证书复现](eight-point-global/README.md) · [解析桥接](eight-point-global/PROOF_BRIDGE.md)|
|9 点|M₉ = 9³⁶/8²⁸，等号为原点加半径 √(9/8) 的正八边形。 [证明与证书复现](nine-point-global/README.md)|
|12 点|真实 D₃ 构型的全局锐值及唯一性。 [英文论文](twelve-point-sharp-global/manuscript/README.md)|
|其他点数及一般结构|一般锐值仍开放。现有结果包括严格可达下界、必要几何条件、环人口刚性、相位相容性与局部稳定性证书。 [证明状态总表](configuration-bounds/PROOF_STATUS.md) · [一般结构](general-structure/README.md)|

八、九点的全局核与计算机辅助证明归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu**。本仓库保留固定来源、许可与署名，并重放其证书。 [原预印本](https://arxiv.org/abs/2611.03088v1)

## 查找资料与复现

|我想找|入口|
|---|---|
|读论文、下载 PDF 或编译 LaTeX|[十二点论文目录](twelve-point-sharp-global/manuscript/README.md)|
|核对定理、适用范围与外部依赖|[证明状态总表](configuration-bounds/PROOF_STATUS.md)|
|运行精确证书、查软件依赖|[复现指南](docs/REPRODUCING.md)|
|了解十二点上界到锐界的推导过程|[候选与核研究记录](docs/CANDIDATE_KERNEL_RESEARCH.md)|
|查一般点数、旧证书及历史发布|[研究索引](docs/RESEARCH_INDEX.md) · [仓库状态](docs/REPOSITORY_STATUS.md)|
|查作者原稿及早期证明|[原稿档案](source-review/README.md) · [既有证明](earlier-results/README.md)|

十二点全域证书用 Python 标准库运行。完整依赖链的顺序，以及八、九点固定依赖和入口，见[复现指南](docs/REPRODUCING.md)。

**English.** An exact computer-assisted proof gives the sharp twelve-point value and the unique maximizing orbit. The archive includes the manuscript, source, fixed certificate inputs, and replay programs. [Read the paper](twelve-point-sharp-global/manuscript/README.md)
