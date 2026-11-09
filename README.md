# Mordell 距离乘积研究

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[GitHub](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。更新于 **2026-11-09**。

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i<j}|z_i-z_j|^2.
```

**十二点全局锐界与全部最大构型已证明。** 八、九点全局锐界均已补齐固定来源的证明与完整证书重放。一般点数的锐值问题仍开放。

## 十二点的全局锐界

**M₁₂=C₁₂**。真实 D₃ 候选达到全局最大值，且最大构型在整体旋转和点重排下唯一。三层人口为 **3、3、6**，由内、中两条正三角形及最外两条共圆正三角形组成。728,050节点的精确两点覆盖补偿为零。等号由基础正定核的径向和角向零矩恢复。[全局证明、精确值与复现](twelve-point-sharp-global/README.md)

## 十二点的上一版分带全局上界

已证 **V₁₂<C₁₂≤M₁₂<(15000/14999)V₁₂**，真实候选值超过全局最优值的 **99.99333%**。三个必要径向域的整数覆盖共707,469个节点，未决0。每对点的补偿压到1/10⁶。[完整证明与复现](twelve-point-banded-upper/README.md)

## 十二点最大构型的内三外九与锐接触

在质量11尺度，任意全局最大构型恰有 **3 个点**的平方半径在 **(0.12, 0.24)**，其余 **9 个点**在 **(0.94, 1.42)**。逐点 Ward 不等式的99,603个节点全部通过精确重放，再用径向 Gram 块排除其他人口。删点法另给任意十点子集的去重心质量严格大于7.573。[分带与人口证书](twelve-point-ward-radii/README.md) · [删点必要域](twelve-point-removal-domain/README.md)

早期锐接触核在真实候选的十四类点对上已实现精确取等，全部十二个基础 Gram 块及接触处曲率通过区间验证。该版本的严格负缺额例保留为发现记录。[精确接触](twelve-point-sharp-contacts/README.md)

## 十二点的 Ward 全局上界

十二点的全局值满足 **V₁₂<C₁₂≤M₁₂<(379/378)V₁₂**。候选值与全局最大值之比严格大于 **378/379（约99.7361%）**，可达值与上界的相对差小于1/378。标准库精确覆盖处理268239个节点，未决0。[完整证书](twelve-point-ward-upper/README.md)

新核使用19个完整驻点 Ward 参数和5个质量参数。三个合并谐波块各有一个负特征值。基核的正半定性与修正项的驻定零和恒等式共同保证全局最大点的核和非负。Jacobian 已逐项核对包含修正的完整解析梯度。[惯性记录](twelve-point-ward-upper/jacobian_inertia.json) · [梯度生成器](scripts/build_ward_gradient.py)

首轮接触数值拟合保留为发现记录。[接触实验](twelve-point-ward-upper/contact_trial.json)

## 十二点的上一版301/300界

已证 **V₁₂<C₁₂≤M₁₂<(301/300)V₁₂**。真实候选的目标值与全局最大值之比严格大于 **300/301（约99.6678%）**。点对补偿降至1/20000，整数覆盖处理229113个节点，未决0。全部12个Gram块、真实根及可达值比较均通过复核。全局锐值与构型唯一性仍待证明。[完整核证书](twelve-point-kernel-refinement-2/README.md)

29份Jacobian原生求导结果组合成了核梯度，完整系数恒等式通过验证。新三点必要式给出九点子集去重心质量q>6.03。Ward恒等式提供19个驻定零和修正参数，质量约束另提供5个。数值全局界由正半定核覆盖给出。[组合梯度](jacobian-gradient/README.md) · [三点必要域](twelve-point-subset-domain/README.md) · [零和核构造](stationary-kernel/README.md)

## 上一轮：Jacobian与十二点3%全局包围

[新连续核证书](twelve-point-kernel-refinement/README.md)已将上一轮15%的范围收紧为 **V₁₂<C₁₂≤M₁₂<103V₁₂/100**。因此已认证真实候选与全局最大值之比严格大于 **100/103≈97.0874%**。新核在质量11尺度只需点对补偿1/2500；完整标准库重放87385个节点、未决0，没有假设未知最大点具有D₃对称。

[Jacobian 0.23.0辅助环境](jacobian-assist/README.md)已配置本地MCP，并实际执行新旧核的精确惯性、真实径向Bernstein系数、Arb二阶包围及平方和证书。新核由广域多谷反例反馈与SDP发现，严格上界由独立整数覆盖承担。原15%证书继续保存；候选恰为全局最优仍未证明。

## 十二点上一轮：15%全局包围（历史证据）

[十二点连续核上界](twelve-point-kernel-upper/README.md)给出原质量12尺度的严格结论

```math
V_{12}\le M_{12}<\frac{23}{20}V_{12},
```

其中V₁₂是仓库旧整数坐标经精确中心化、归一化后的可达值。已认证真实12点候选值还严格大于V₁₂，因此它与全局最大值的比严格大于 **20/23≈86.9565%**。这项结论覆盖全局最大值，没有预设最大构型的D₃对称；尚未证明候选达到最大值。

证明将原核增加点对余量21/10000，精确核验正定矩块，并在整个必要两点域上完成Bernstein与对数切线覆盖：37053个节点、18526次分割、16237个正性叶、2290个域排除叶、未决0。正式接受只需Python标准库，约30秒实际重放完整覆盖。

## 候选全局最优：优先12、10、19点

[1000次确定性多起点搜索](candidate-global-screen/README.md)覆盖10—19点，结合完整局部稳定性及接触条件的复杂度，优先推进12点，其次10点和19点。三个已认证真实候选分别具有精确D₃、D₂、D₆对称与14、15、22类点对；对称性只针对这些候选，不预设所有全局最大点对称。搜索频率是数值证据，未给全局最优概率。

[十点必要两点域](maximizer-domain/README.md)联合八、九点锐界、严格子集及二阶Schwarz信息，证明质量9下任何全局最大点都满足平方半径<1.5403、点距平方>0.4、八点补集中心化质量>5.68，以及额外非线性点对约束。[整区域局部上界](candidate-basin/README.md)则证明：10、12点到候选轨道的欧氏距离≤2⁻¹⁸、19点≤2⁻²⁰时，距离积不能超过候选值。

[候选适配正定核](candidate-kernel/README.md)给出矩罚项与两点上界的严格接合公式，也保留10、12点原始数值核的精确失败点对。十二点的最终证明使用分带、精确接触和连续覆盖。 [候选与核研究记录](docs/CANDIDATE_KERNEL_RESEARCH.md)

## 11月8日：四环共振、非线性缺额与高点数局部认证

[四环人口定理](four-ring-rigidity/THEOREMS.md)继续加强逐点力条件：不等人口驻点至少有三条环小于最大人口；四环必须有唯一最大人口，且它严格小于其余人口总和。外部总人口不超过最大值两倍时，还需要人口重复或加法共振；特别地，(s,s,s,m) 的异人口情形只能保留 m=2s，未声称保留分支存在。

[非线性相位证书](phase-deficit-sharp/README.md)把指定四正方形环的认证缺额从约2.36195提高至5.97446865。保留指数正项后，[完整家族比较](family-exclusion-nov8/README.md)证明：该16点家族及每个半径独立变化±1/10000的整个盒，在所有连续相位下的值均低于一个已知16点可达值的5%，所以整个家族不能包含全局最优者。

[25、35、42、46点的区间认证](validated-maxima-nov8/README.md)全部通过完整力与Hessian检验。结合前轮，已有12个不同点数的严格局部极大轨道和15点的一个不稳定轨道的实际证书。详细范围见[11月8日研究记录](docs/NOV8_RESEARCH.md)；一般锐值仍开放。

## 11月5日：完整力平衡、旋转阶与真实局部极值

本轮得到[三环人口刚性及旋转阶上限](ring-rigidity/THEOREMS.md)：恰由三条互异半径的同心正多边形环组成的驻点，人口必须全部相等，允许任意连续相位与可选中心。更一般地，具有阶 m 旋转对称的局部最大点必须满足 `floor(m²/4)≤2(n−1+h)`，h 表示是否含原点；违反时完整 Hessian 有增益方向。

[自由构型搜索](finite-search-nov5/README.md)在16个点数上严格改善旧保存下界，最大为46点的 **4.931680542211…%**。[整数区间认证](validated-maxima/README.md)进一步给出10、11、12、13、14、17、18、19点的8个严格局部极大轨道，以及15点的一个真实不稳定驻点；局部极大性与全局最优性分开。

[循环分数预算](phase-packing/README.md)把指定四正方形环的全连续相位有理上界降低超过33.1%，并扩展到整个半径盒；这类上界只适用于所列规则环家族。完整结果、证明接口和后续缺口见[本轮研究记录](docs/NOV5_RESEARCH.md)。

## 当前结果

|点数或问题|结论|证明与边界|
|---|---|---|
|2≤n≤5|Mₙ=nⁿ，等号为单位正 n 边形|[五点作者方法与既有证明](earlier-results/README.md)|
|n=6、7|M₆=6¹⁵/5¹⁰，M₇=7²¹/6¹⁵；等号为原点加相应正多边形|[主稿与依赖审计](earlier-results/six_seven_coverage_audit.md)；本项目现有全局证明保留 Faugère–Svartz 分类输入|
|n=8|M₈=8²⁸/7²¹；等号为原点加半径 √(8/7) 的正七边形|[全域证书复现](eight-point-global/README.md) · [解析桥接](eight-point-global/PROOF_BRIDGE.md)；完整重放 386915 个盒，未决 0|
|n=9|M₉=9³⁶/8²⁸；等号为原点加半径 √(9/8) 的正八边形|[完整证明与证书重放](nine-point-global/README.md)；完整初始域遍历 31613 个盒，未决 0|
|n=12|M₁₂=C₁₂，D₃候选在旋转及重排下唯一|728050节点零补偿覆盖，三层人口3、3、6。[全局锐界](twelve-point-sharp-global/README.md)|
|n≥10|中心正 (n−1) 边形不再全局最优|[十点整数比较](docs/TEN_POINT_TRANSITION.md)；n≥11 沿用 Gueron–Shafrir 不稳定性定理；不确定实际最优者|
|10≤n≤100|旧91项保存下界保留；本轮在其中16项再取得严格改善|[本轮自由构型](finite-search-nov5/README.md) · [旧基线](symmetry-breaking/README.md)；均为可达下界|
|10、12、19点候选范围|精确D₂/D₃/D₆候选、整局部区域的候选值上界及十点全局必要域|[候选筛选](candidate-global-screen/README.md) · [局部上界](candidate-basin/README.md) · [十点必要域](maximizer-domain/README.md)；补集未覆盖、未证明全局锐值|
|局部驻点与稳定性|12个不同点数的严格局部极大轨道、1个不稳定轨道有实际区间证书|[25、35、42、46点](validated-maxima-nov8/README.md) · [前轮8点数](validated-maxima/README.md)；不穷尽驻点、不认证全局最优|
|一般 n 的结构|三、四环人口刚性与共振、旋转阶上限、非线性循环缺额及家族排除|[完整状态表](configuration-bounds/PROOF_STATUS.md) · [四环定理](four-ring-rigidity/README.md) · [家族排除](family-exclusion-nov8/README.md)；未证明晶格结晶或简洁锐值公式|

八、九点全局核与其计算机辅助证明归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu** 的[arXiv:2611.03088v1](https://arxiv.org/abs/2611.03088v1)。本仓库固定来源、保留许可与署名、复现证书并写明与已有定理的衔接。

前轮[连续相位定理](general-structure/CONTINUOUS_PHASES.md)与十二点固定家族的精确见证继续保留。本轮将相位预算、完整逐点力和完整二阶条件接入同一研究主线，各方法保留自身作用域。

## 阅读与复现

|需要了解什么|入口|
|---|---|
|十二点可达值与全局最大值的严格比值界|[十二点连续核上界](twelve-point-kernel-upper/README.md)|
|接近全局最优的候选、必要点对域与正定核|[候选与核研究记录](docs/CANDIDATE_KERNEL_RESEARCH.md)|
|11月8日的新结构、家族排除及高点数认证|[本轮研究记录](docs/NOV8_RESEARCH.md)|
|11月5日的新定理、下界与局部认证|[本轮研究记录](docs/NOV5_RESEARCH.md)|
|结论、假设、外部依赖和未解决问题|[证明状态总表](configuration-bounds/PROOF_STATUS.md)|
|九点完整证明与计算范围|[九点专题](nine-point-global/README.md)|
|一般点数与最优构型的数学结构|[结构研究](general-structure/README.md) · [一般 n 研究索引](docs/RESEARCH_INDEX.md)|
|证书运行、固定输入及软件依赖|[复现指南](docs/REPRODUCING.md)|
|专题关系、历史 ZIP 和 Release|[研究索引](docs/RESEARCH_INDEX.md)|
|本次仓库整理及历史状态如何理解|[仓库状态](docs/REPOSITORY_STATUS.md) · [整理前首页](docs/HISTORY_2026-11-04.md)|
|作者原稿与五、六、七点证明|[既有证明](earlier-results/README.md) · [完整原稿](source-review/README.md)|

轻量核验使用 Python 标准库：

```bash
python -B -S -O scripts/verify_current.py
```

该命令重放所列小型证书、核对已验收的全域证据及输入绑定；**完整区间覆盖需按[复现指南](docs/REPRODUCING.md)安装固定数值依赖后另行执行**。哈希核对、局部搜索和小残差分别只说明其实际检查范围。

旧八点局部森林仍保存 28 个开放后代盒，表示这条历史坐标搜索尚未逐盒完成；八点全局锐界由新的全域核证书与解析桥接证明。历史证明、原稿、ZIP、来源和原始报告均保留，不把历史未决数改记为零。

仓库自有材料未指定统一许可证；导入的第三方代码、数据与说明沿用其明确许可。计算机辅助证书和协作审阅不等于外部同行评审或整篇形式化证明。

**English.** Fengcan Jiang's Mordell distance-product research. The n=12 maximum and the unique maximizing orbit have an exact computer-assisted proof. The attributed n=8 and n=9 certificates have been reproduced with fixed provenance and full-domain replay.
