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
|八点最新进展与准确剩余分支|[中心化切片与质量对偶](joint-refinement/CENTERED_MASS.md) · [当前覆盖](joint-refinement/README.md) · [联合归约](docs/EIGHT_POINT_REDUCTION.md)|
|各专题的关系、历史版本和下载|[研究索引](docs/RESEARCH_INDEX.md)|
|如何实际运行证书|[复现指南](docs/REPRODUCING.md)|
|作者原稿与五、六、七点证明|[既有证明](earlier-results/README.md) · [完整原稿](source-review/README.md)|

## 当前续算：中心化切片与质量对偶

在实际质量7、无碰撞的八点上，完整实能量 Hessian 半正定已经强迫全部逐点驻点方程：旋转向量 iq 的二次型恒为0，而 H(iq)=ig，其中 g 为真实梯度。这一结论需要完整 Hessian；低阶主块不能替代。

本轮在目标值二次子问题中精确保留坐标和为0，并用严格的质量对偶参数收紧上界。只有在交换方向凸性通过核对后，才枚举切片顶点；这一新筛选额外排除3789个传播后包围域，无需待检点驻定。

与旧树合并后有310017个节点、154992个排除叶，仍有17个连续开放区域，部分区域仍较宽。**完整排空仍未完成。** [新解析证明](joint-refinement/CENTERED_MASS.md) · [精确森林与复核](joint-refinement/README.md)

按上一阶段快照追踪，原19个开放盒中6个已完整关闭；余下13个原区域细分后留下当前17个开放盒。[旧开放域进度](joint-refinement/reports/frontier_progress.json)给出坐标盒统计，完整问题仍然开放。

## 上一阶段：等价多项式矩阵与冻结覆盖

在实际质量7、中心化且无碰撞的八点上，真实驻点平方和与完整实能量 Hessian 可等价合并为一个8×8整数系数 Hermitian 多项式矩阵半正定条件；也可保留12个实 Ward 等式和6×6完整稳定性矩阵，最高次数分别为8和14。该模型保留全部逐点物理条件。

冻结的上一阶段覆盖树有60021个节点、29990个排除叶和21个开放坐标区域；新续算见上节。另认证一个旧椭圆内的非零矩真实驻点：其旧四模式子块正定，但完整 Hessian 不稳定，且目标值低于候选。[等价证明、真实驻点盒与历史检查点](joint-closure/README.md)

## 已证径向区域：缺额门槛推进到 1/1200

令 q=√(7/8)z、Σq=0、Σ|q|²=7，记 E=Σ|q_i|²(|q_i|²−1)²。继承已披露依赖的七点锐界与现有局部能量定理，得到

```math
E\le1/1200\quad\Longrightarrow\quad
\log(\Delta_q/7^7)\le-E/32,\qquad
\Delta_q/7^7\le1-E/64.
```

在此区域，只有中心正七边形能达到候选值。无需预设驻点、稳定性、旧椭圆或半径上界。相对上传统一包的1/35000，径向门槛扩大175/6，约29.17倍；**E>1/1200 的剩余物理系统仍未排空**。

新证明用中心化消去径向扰动的一阶项，再以圆上七点的全局角向凹性控制候选轨道距离。原统一包按字节保留，14个标准库入口已重放。[完整证明、真实点见证与复现](radial-gap/README.md)

## 非零矩代数耦合（保留）

在旧椭圆、真实驻点方程和 Δ_w≥4²⁸7⁷ 的交集中，已推出力1尺度下的新必要条件

```math
|p_2|\ge10^{-36}|p_3|.
```

它把旧的“p₂=0 ⇒ p₃=0”整数证书量化到非零矩分支，不要求靠近候选。完整重放重新构造结果式、检查14179个整数系数，并计算七个乘子范数。常数很小，只排除靠近 p₂=0 轴的细楔形区域，**不能完成剩余非零矩分支的排空**。[证明、联合目标值筛选式与复现](algebraic-gap/README.md)

## 已证目标值区域：直接对数余项与五倍矩管域

取力1尺度 w=√(7/2)z、p_k=Σw_i^k。d 为到“原点加半径2正七边形”全部旋转、重排轨道的欧氏距离。对 Σ|w_i|²=28 的八点构型，已证明

```math
d\le\frac15\quad\Longrightarrow\quad
\log\frac{\Delta_w}{4^{28}7^7}\le-\frac{d^2}{30}.
```

目标值排除的欧氏半径从1/20扩大到1/5。进一步，中心化构型只要满足

```math
|p_2|\le1/2,\quad |p_3|\le1,\quad |p_4|\le2,\quad
|p_5|\le4,\quad |p_6|\le8,\quad |p_8|\le32,
```

就必有 d²<3/100，进入上述排除区域。六个模长门槛都是上一管域的**五倍**，仍**无须另设 p₇ 条件**。证明联合实际距离 Laplacian 的完整谱、对数展开余项和 Fourier–Parseval 根误差估计。[完整证明与精确报告](eight-point/ENERGY_BASIN.md)

整个区域内只有候选轨道能达到候选值；这里不据此宣称较大区域内所有真实驻点都只能是候选。管域之外的[六类剩余条件](eight-point/OPEN_BRANCHES.md)及完整物理系统的联合穷尽仍未完成。

原有 p₂、p₃ 椭圆与零矩分类继续作为前置条件，后续只在它们限定的非零矩物理系统中推进。[三个尺度与当前完整归约](docs/EIGHT_POINT_REDUCTION.md)明确说明新旧结果怎样取交集。

## 主要成果与边界

|问题|当前结论|范围与证明|
|---|---|---|
|五点锐值|M₅=5⁵，单位正五边形|[作者二阶条件与矩方法](earlier-results/five_point_method.md)|
|六、七点锐值|M₆=6¹⁵/5¹⁰；M₇=7²¹/6¹⁵|[主稿](earlier-results/six_seven_proof.tex)；全局证明引用 Faugère–Svartz 完整物理驻点分类，[依赖审计](earlier-results/six_seven_coverage_audit.md)|
|八点候选与部分排除|候选 8²⁸/7²¹ 可达；E≤1/1200、矩管域及若干几何类已排除|[最新径向界](radial-gap/README.md) · [当前联合归约](docs/EIGHT_POINT_REDUCTION.md)；**完整锐界未证**|
|有限 n 下界|10≤n≤100的91项保存整数见证改善此前包络|[自由构型](symmetry-breaking/README.md)；不是全局最优认证|
|最优构型半径层数|n≥10¹⁴时 K>0.51√n；初等渐近系数>0.513|[最新层数证明](annular-layer/README.md)；保留外部渐近输入时系数>0.556，后者无有效有限起点|
|一般 n 上界与分布|严格圆盘修正、线性缺额、二维运输、含边界覆盖|[证明总表](configuration-bounds/PROOF_STATUS.md)；未确定晶格、精确 γ 或所有有限 n 锐值|

所有结论以正文的尺度、假设和外部依赖为准。历史较弱常数集中在[研究索引](docs/RESEARCH_INDEX.md)，不与当前最强结论混列。

## 复核与下载

在仓库根目录运行当前幂和与八点检查（Python 标准库）：

```bash
python -B -S -O scripts/verify_current.py
```

[联合系统部分覆盖包](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-joint-partial-2026-11-01) · [径向—角向包](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-radial-2026-11-01) · [全局代数包](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-algebraic-2026-11-01) · [能量阶段包](eight-point/mordell-eight-point-energy-2026-11-01.zip) · [所有专题和历史附件](docs/RESEARCH_INDEX.md) · [完整复现说明](docs/REPRODUCING.md)。历史 ZIP 保持原字节；GitHub 自动生成的 Source code 压缩包不包含 Release 中的大型研究附件。

程序认证精确算术、有限盒及所列常数；解析和几何桥梁仍需阅读证明。浮点搜索、小残差、有限目录和验证通过不能代替全域覆盖。没有外部同行评审、整篇形式化证明或全球优先权声明；仓库未指定许可证。

**English:** Jiang Fengcan's research on finite-point Mordell discriminant maxima. Exact certificates and conditional exclusions are indexed with their hypotheses. The complete eight-point sharp bound and a simple general-n sharp formula remain open.
