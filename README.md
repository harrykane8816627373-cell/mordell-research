# Mordell 距离乘积研究

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[GitHub](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。更新：**2026-11-08**。

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i<j}|z_i-z_j|^2.
```

**八、九点全局锐界均已补齐固定来源的证明与完整证书重放。一般 n≥10 的锐值与全部最优构型仍未确定。** 定理、带条件的排除、可达下界和历史部分覆盖分别注明范围。

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
|n≥10|中心正 (n−1) 边形不再全局最优|[十点整数比较](docs/TEN_POINT_TRANSITION.md)；n≥11 沿用 Gueron–Shafrir 不稳定性定理；不确定实际最优者|
|10≤n≤100|旧91项保存下界保留；本轮在其中16项再取得严格改善|[本轮自由构型](finite-search-nov5/README.md) · [旧基线](symmetry-breaking/README.md)；均为可达下界|
|局部驻点与稳定性|12个不同点数的严格局部极大轨道、1个不稳定轨道有实际区间证书|[25、35、42、46点](validated-maxima-nov8/README.md) · [前轮8点数](validated-maxima/README.md)；不穷尽驻点、不认证全局最优|
|一般 n 的结构|三、四环人口刚性与共振、旋转阶上限、非线性循环缺额及家族排除|[完整状态表](configuration-bounds/PROOF_STATUS.md) · [四环定理](four-ring-rigidity/README.md) · [家族排除](family-exclusion-nov8/README.md)；未证明晶格结晶或简洁锐值公式|

八、九点全局核与其计算机辅助证明归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu** 的[arXiv:2611.03088v1](https://arxiv.org/abs/2611.03088v1)。本仓库固定来源、保留许可与署名、复现证书并写明与已有定理的衔接。

前轮[连续相位定理](general-structure/CONTINUOUS_PHASES.md)与十二点固定家族的精确见证继续保留。本轮将相位预算、完整逐点力和完整二阶条件接入同一研究主线，各方法保留自身作用域。

## 阅读与复现

|需要了解什么|入口|
|---|---|
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

**English:** Fengcan Jiang's Mordell distance-product research. The n=8 and n=9 global certificates have been reproduced with fixed provenance and full-domain replay. Exact general-n maxima and all maximizing configurations for n≥10 remain open here.
