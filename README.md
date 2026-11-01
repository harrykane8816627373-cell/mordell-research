# Mordell 距离乘积研究

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[GitHub](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。更新：**2026-11-01**。

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i\lt j}|z_i-z_j|^2.
```

**八点完整无条件锐界、一般非零矩分支和一般 n 的简洁锐值公式仍开放。** 本仓库分别记录已证定理、带条件的排除、可达下界和未完成覆盖。

## 阅读入口

|需要了解什么|入口|
|---|---|
|目前哪些结论已证，使用了哪些外部输入|[当前证明状态总表](configuration-bounds/PROOF_STATUS.md)|
|八点最新进展与准确剩余分支|[八点攻关](eight-point/README.md) · [剩余覆盖清单](eight-point/OPEN_BRANCHES.md)|
|各专题的关系、历史版本和下载|[研究索引](docs/RESEARCH_INDEX.md)|
|如何实际运行证书|[复现指南](docs/REPRODUCING.md)|
|作者原稿与五、六、七点证明|[既有证明](earlier-results/README.md) · [完整原稿](source-review/README.md)|

## 最新八点进展：消去 p₇ 限制的完整管域排除

取力1尺度 w=√(7/2)z、p_k=Σw_i^k。对中心化、Σ|w_i|²=28的八点构型，若

```math
|p_2|\le1/10,\quad |p_3|\le1/5,\quad |p_4|\le2/5,\quad
|p_5|\le4/5,\quad |p_6|\le8/5,\quad |p_8|\le32/5,
```

则到“原点加半径2正七边形”轨道的距离 d<1/20，并且

```math
\log\frac{\Delta_w}{4^{28}7^7}\le-\frac{d^2}{35}.
```

这六个条件旋转不变，**无须另设 p₇ 条件**。整个区域内只有候选轨道能达到候选值；真实驻点也只能是候选。新区域严格包含此前全矩半径1/10的局部盒，另有精确非空见证。[完整证明与常数复核](eight-point/MOMENT_TUBE.md)

管域外的六类矩条件及旧11盒之外，仍未完成真实逐点力方程、完整稳定性和目标阈值的联合穷尽。[剩余清单](eight-point/OPEN_BRANCHES.md)明确保留这个缺口。

## 主要成果与边界

|问题|当前结论|范围与证明|
|---|---|---|
|五点锐值|M₅=5⁵，单位正五边形|[作者二阶条件与矩方法](earlier-results/five_point_method.md)|
|六、七点锐值|M₆=6¹⁵/5¹⁰；M₇=7²¹/6¹⁵|[主稿](earlier-results/six_seven_proof.tex)；全局证明引用 Faugère–Svartz 完整物理驻点分类，[依赖审计](earlier-results/six_seven_coverage_audit.md)|
|八点候选与部分排除|候选 8²⁸/7²¹ 可达；新管域、旧矩盒及若干几何类已排除|[八点当前状态](eight-point/README.md)；**完整锐界未证**|
|有限 n 下界|10≤n≤100的91项保存整数见证改善此前包络|[自由构型](symmetry-breaking/README.md)；不是全局最优认证|
|最优构型半径层数|n≥10¹⁴时 K>0.51√n；初等渐近系数>0.513|[最新层数证明](annular-layer/README.md)；保留外部渐近输入时系数>0.556，后者无有效有限起点|
|一般 n 上界与分布|严格圆盘修正、线性缺额、二维运输、含边界覆盖|[证明总表](configuration-bounds/PROOF_STATUS.md)；未确定晶格、精确 γ 或所有有限 n 锐值|

所有结论以正文的尺度、假设和外部依赖为准。历史较弱常数集中在[研究索引](docs/RESEARCH_INDEX.md)，不与当前最强结论混列。

## 复核与下载

在仓库根目录运行当前幂和与八点检查（Python 标准库）：

```bash
python -B -S -O scripts/verify_current.py
```

[最新八点离线包](eight-point/mordell-eight-point-2026-11-01.zip) · [所有专题和历史附件](docs/RESEARCH_INDEX.md) · [完整复现说明](docs/REPRODUCING.md)。历史 ZIP 保持原字节；GitHub 自动生成的 Source code 压缩包不包含 Release 中的大型研究附件。

程序认证精确算术、有限盒及所列常数；解析和几何桥梁仍需阅读证明。浮点搜索、小残差、有限目录和验证通过不能代替全域覆盖。没有外部同行评审、整篇形式化证明或全球优先权声明；仓库未指定许可证。

**English:** Jiang Fengcan's research on finite-point Mordell discriminant maxima. Exact certificates and conditional exclusions are indexed with their hypotheses. The complete eight-point sharp bound and a simple general-n sharp formula remain open.
