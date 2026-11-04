# Mordell 距离乘积研究

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[GitHub](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。更新：**2026-11-04**。

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i<j}|z_i-z_j|^2.
```

**八、九点全局锐界均已补齐固定来源的证明与完整证书重放。一般 n≥10 的锐值与全部最优构型仍未确定。** 定理、带条件的排除、可达下界和历史部分覆盖分别注明范围。

## 当前结果

|点数或问题|结论|证明与边界|
|---|---|---|
|2≤n≤5|Mₙ=nⁿ，等号为单位正 n 边形|[五点作者方法与既有证明](earlier-results/README.md)|
|n=6、7|M₆=6¹⁵/5¹⁰，M₇=7²¹/6¹⁵；等号为原点加相应正多边形|[主稿与依赖审计](earlier-results/six_seven_coverage_audit.md)；本项目现有全局证明保留 Faugère–Svartz 分类输入|
|n=8|M₈=8²⁸/7²¹；等号为原点加半径 √(8/7) 的正七边形|[全域证书复现](eight-point-global/README.md) · [解析桥接](eight-point-global/PROOF_BRIDGE.md)；完整重放 386915 个盒，未决 0|
|n=9|M₉=9³⁶/8²⁸；等号为原点加半径 √(9/8) 的正八边形|[完整证明与证书重放](nine-point-global/README.md)；完整初始域遍历 31613 个盒，未决 0|
|n≥10|中心正 (n−1) 边形不再全局最优|[十点整数比较](docs/TEN_POINT_TRANSITION.md)；n≥11 沿用 Gueron–Shafrir 不稳定性定理；不确定实际最优者|
|10≤n≤100|91 项保存整数构型改善此前下界包络|[自由构型](symmetry-breaking/README.md)；可达下界，不是全局最优认证|
|一般 n 的结构|有限上下界、分布与层数约束；新增连续相位相容性与环对缺额|[完整状态表](configuration-bounds/PROOF_STATUS.md) · [一般点数结构](general-structure/README.md)；未证明晶格结晶或简洁锐值公式|

八、九点全局核与其计算机辅助证明归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu** 的[arXiv:2611.03088v1](https://arxiv.org/abs/2611.03088v1)。本仓库固定来源、保留许可与署名、复现证书并写明与已有定理的衔接。

新增结构定理把正规多环的相位竞争化为人口的二进赋值条件，并对不相容循环给出严格正缺额。十二点的三个正方形环提供精确见证：连续相位可超过同一固定家族全部正负符号相位的两倍；这项比较不宣称改进无约束十二点最优下界。

## 阅读与复现

|需要了解什么|入口|
|---|---|
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
