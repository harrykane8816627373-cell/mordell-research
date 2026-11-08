# 当前证明状态与统一入口

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**；推导与代码包含 OpenAI Codex 协作。更新于2026-11-05。定理的全部假设和证明以所链接正文为准。

**八点当前结果：[全局证书复现](../eight-point-global/README.md) · [与既有严格子集定理的桥接](../eight-point-global/PROOF_BRIDGE.md) · [实际重放记录](../eight-point-global/reports/replay.json)。** Ji、Xu、Yu预印本的固定三变量核证书已完整重放，并与本仓库已有子集界接合，得到 M₈=8²⁸/7²¹及全部等号。质量7下Δ≤7⁷，等号恰为中心加单位正七边形，所以用户的Δ≥7⁷、s₂s₃≠0联合系统为空。这条桥接不需要七点完整驻点分类；旧局部森林及其28个开放盒保留为历史有限搜索记录。

**九点当前结果：** Ji、Xu、Yu 的[预印本 §4.4](https://arxiv.org/pdf/2611.03088v1)给出 M₉=9³⁶/8²⁸，等号恰为原点加半径 √(9/8) 的正八边形。本仓库已补齐[完整解析证明](../nine-point-global/PROOF.md)与[固定来源的独立重放](../nine-point-global/README.md)：从完整初始域遍历31613个盒，15807个终端盒，未决0；全部8个上游检查实际执行，来源及软件依赖保留。

一般 n 的当前最强结论集中在下方表格；各次推进与较弱历史常数见[研究索引](../docs/RESEARCH_INDEX.md)。六、七点的外部分类依赖、有限层数与渐近层数的不同前提继续保留。

```math
M_n=\max_{\sum_i|z_i|^2=n}\prod_{i\lt j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2.
```

## 0. 11月5日一般点数研究

当前新增四条主线，完整接口见[本轮研究记录](../docs/NOV5_RESEARCH.md)：

|结果|严格结论与范围|
|---|---|
|[三环人口刚性与旋转阶](../ring-rigidity/README.md)|三条互异半径正规环驻点必须等人口；阶m旋转对称局部最大者满足floor(m²/4)≤2(n−1+h)，来自完整逐点力与Hessian|
|[连续相位循环预算](../phase-packing/README.md)|原始/对偶精确分数预算、一般奇加权循环、整个半径盒；指定四正方形环有理上界较旧贪心降低超过33.1%，不是无约束上界|
|[自由构型整数下界](../finite-search-nov5/README.md)|16项严格改善旧保存见证；46点增益约4.93168%，没有认证一般最优值|
|[真实驻点与完整稳定性](../validated-maxima/README.md)|10、11、12、13、14、17、18、19点共8个严格局部极大轨道；15点有真实不稳定驻点；均有整数区间证据，不穷尽分支|

## 1. 小点数、构型与未解决问题

|范围|已取得的结果|证明入口与适用边界|
|---|---|---|
|五点|$`M_5=3125`$，等号单位正五边形，允许旋转、重排|[作者方法](../earlier-results/five_point_method.md)、[PDF](../earlier-results/five_point_proof.pdf)、[TeX](../earlier-results/five_point_proof.tex)；主线为二阶条件、径向矩、投影及四次对数上界|
|六、七点|$`M_6=6^{15}/5^{10}`$、$`M_7=7^{21}/6^{15}`$，等号原点加相应正多边形|[主稿](../earlier-results/six_seven_proof.tex)、[复核](../earlier-results/six_seven_verification.md)、[覆盖审计](../earlier-results/six_seven_coverage_audit.md)；全局证明依赖 Faugère–Svartz 完整物理驻点分类，七点356次枝的物理排空为论文输入，未在本地独立重建|
|六、七点含原点类|不调用上述分类的锐界与方差惩罚|[独立证明](../earlier-results/six_seven_origin_proof.tex)；结论仅适用于含原点类|
|八点|M₈=8²⁸/7²¹，等号恰为原点加半径√(8/7)的正七边形，允许旋转与重排；质量7下非零二三阶矩挑战者系统为空|[Ji–Xu–Yu证书复现](../eight-point-global/README.md)、[既有严格子集定理接合](../eight-point-global/PROOF_BRIDGE.md)；固定外部代码、原作者署名及软件依赖明确保留，未声称外部同行评审或形式化证明|
|九点|M₉=9³⁶/8²⁸；唯一等号形状为中心正八边形|[Ji–Xu–Yu证书复现及完整证明](../nine-point-global/README.md)；8个上游检查、完整初始域31613个盒、未决0；[实际执行记录](../nine-point-global/reports/replay.json)|
|n≥10|中心正 (n−1) 边形严格全局次优|[十点整数比较与转换证明](../docs/TEN_POINT_TRANSITION.md)、[Ji–Xu–Yu §5](https://arxiv.org/pdf/2611.03088v1)；n≥11 使用 Gueron–Shafrir 不稳定性定理；未确定实际最优者|
|中心正多边形类|5≤n≤10模旋转严格局部最大；n=10使用四阶项，n≥11有显式增益扰动|[一般 n 审核第6、7节](../earlier-results/general_n_proof.md)；五点中心正方形、十点中心正九边形仍全局次优|
|一般 n|已有有限上下界、全局网格包围、结构与渐近定理|没有得到所有有限 n 的简洁锐值公式；有限代数刻画不代表已执行巨大消元或枚举|

八点已证的 $`p_2=0`$、$`p_3=0`$ 分类保留真实驻点及相应稳定约束；七点共任意圆、六共圆且两例外等距、两直线等类别各有假设。旧路线中的连续核只连接部分远端带，当时尚有过渡区域未覆盖。这些旧分类继续保留各自条件；新的八点全局结果由上述核证书与子集域桥接承担。

## 2. 既有一般 n 方法与后续联合

新增[一般点数结构专题](../general-structure/README.md)研究连续相位相容性与环对缺额；受限环族中的极值与无约束全局最优须区分。

|结论|当前状态与范围|证明、数据或代码入口|
|---|---|---|
|连续多环相位相容性与严格缺额|所有环对同时取最大，当且仅当人口二进赋值的最小值最多重复两次、其余互异；不相容奇循环有显式正缺额|[完整证明与精确见证](../general-structure/README.md)；固定人口和半径的相位家族上界，不是无约束 Mₙ 上界；十二点连续相位超过同一家族符号最优两倍|
|高斯平均、截断指数上界及唯一最优参数|已证；唯一性针对上界参数，不是构型唯一性|[一般 n 审核第2节](../earlier-results/general_n_proof.md)|
|格点、有理多环和保存整数构型|每个 n 有可计算可达下界；格点在二次指数尺度渐近最优|[原始主稿](../source-review/original_main.tex)、[修订主稿](../source-review/revised_main.tex)、[来源核对](../source-review/SOURCE_REVIEW.md)|
|三角整数格与加权矩|已证有限有理下界；前三阶矩实际计算到n=7，格点比较至n=30|[有限 n 专题](../finite-n/README.md)、[矩证明](../finite-n/MOMENT_METHOD.md)、[格点数据](../finite-n/triangular_lattice_values.json)|
|严格真子集内部矩约束|对全局最优构型及2≤m<n严格成立；已加强为显式正有理余量|[完整证明](../finite-n/SUBSET_PROOF.md)；单点两边为零，全体点为等式，不推广到所有驻点|
|标准化最大值相邻严格下降|原稿已有完整证明；后续加入随机子集四阶方差折扣|[来源核对](../source-review/SOURCE_REVIEW.md)、[全局方法专题](../global-n/README.md)|
|完整有效网格|原稿已有定义与误差证明；后续联合稳定矩阵、点距与较小整数箱改善|[原审核第9节](../earlier-results/general_n_proof.md)、[混合网格](../source-review/HYBRID_GRID.md)；上界要求完整覆盖，抽样只给下界|
|不同半径层数|有效K>0.51√n适用于所有n≥10¹⁴；仅初等预算的渐近系数>0.513，保留此前外部渐近输入时>0.556|[最新环带与可变平滑证明](../annular-layer/README.md)；任意角度，区分有效有限式与外部输入的渐近式；仍不确定真实层数|
|环带人口与层数共同对偶|有限式保留Σα_j(n_j−K_j)、统一层价格、同一h的Dirichlet能量和一次平滑修正；2048个固定环带认证新渐近系数|[环带对偶证明](../annular-layer/ANNULAR_DUAL.md)；固定分区，不是微观人口分类；受限支撑线不能用于任意实际点距|
|层半径自身的分布|q=R/√n尺度，每个固定窗口(a,b)有liminf K_n((a,b))/√n≥(b−a)/2；有限弱极限的额外层有正半径加权质量|[层分布证明](../layer-circumference/LAYER_THEOREM.md)；固定窗口，不扩大到缩小的微观窗口或晶格分类|
|共同变分预算|一般n的移位Gram对偶；百点四层二阶幂和模<47/20，指定中心空洞下<17/10|[主变分证明](../joint-variational/FRAMEWORK.md)、[最新百点证书](../annular-layer/FINITE_JOINT.md)；只给条件必要式，未执行全部构型盒穷尽|
|高斯局部角向密度|任意相似尺度/旋转的平滑加权对比≤√(487E/2000)/n；百点四层对比<0.375，指定低阶能量下<0.073、指定空洞下<0.27|[全尺度证明](../gaussian-angular/DISTRIBUTION.md)、[最新百点加强](../annular-layer/FINITE_JOINT.md)；与旧k=4测试需保留交叉项，非独立相加或硬人口误差改善|
|内盘微观覆盖|所有全局最优者：n≥64、任意原尺度圆心平方半径≤1，最近点距<6/√(n−1)；更大固定内盘有明确起点|[全局插值证明](../distribution-geometry/MICROSCOPIC_COVERING.md)；不声称背景边界的同阶覆盖|
|整个闭背景盘覆盖|所有全局最优者、n≥64：原z尺度覆盖半径O(√(log n/n))，有明确公式|[扩张内盘证明](../stationary-screening/EDGE_COVERING.md)；包括背景边界，仍保留√log n损失|
|驻点混合矩与空洞联合|一般k的卷积闭合、带残差接口、高阶矩消元四次罚项；空洞与守恒矩的完整Schur能量|[驻点与矩校正专题](../stationary-screening/README.md)；保留支撑、驻点和真实测量条件，不独立重复相加|
|完整二维分布与局部人口|W₁有效O(n⁻¹ᐟ²)，百点W₁<0.111199952674；任意中心有界半径圆盘人口质量误差O(n⁻¹ᐟ³)|[二维分布专题](../distribution-geometry/README.md)；窄渐近常数另需已注明外部输入；含边界覆盖为较粗阶|
|平移旋转薄矩形人口|局部最大稳定性的奇跳谱容量；200点局部最大者的任意长2.7宽0.01矩形不能容纳20点|[条带证明与完整非空见证](../distribution-geometry/SUBSET_GEOMETRY.md)；必要几何条件，非锐值证明|
|确定的线性渐近项|沿用已注明展开，最新γ上端为3/4−log(2)/2−37·2⁻²⁴|[渐近专题](../asymptotics/README.md)、[最新参数接口](../overlap-budget/PARAMETER_THEOREM.md)；[近期文献与本地认证边界](../local-defect/LITERATURE_SCOPE.md)另列，未在本地认证精确γ|
|全 n 圆盘上界及最新严格修正|旧Newton保持；线性缺额提高为37·2⁻²⁴(n−11√n)₊|[新参数证明](../overlap-budget/PARAMETER_THEOREM.md)、[旧Newton证明](../joint-bounds/ANGULAR_PROOF.md)；两证书取较强者，允许相加的项另见正文，不声称最优常数|
|短边图人口与共同预算|力为1尺度8/5截止：八点全局≤12对、百点≤101对；八点指定矩条件下≤2对，四阶矩门槛降至101/10|[精确积分与新联合门槛](../layer-circumference/LAYER_THEOREM.md)；点对代价可与一个完整Dirichlet下端相加，禁止重复计数|
|径向矩与角向共同预算|匹配二阶矩的径向Gram；八点四阶矩≥52/5且二阶幂和模≥9/10的整类排除|[联合证明](../local-defect/RADIAL_HIERARCHY.md)；同一径向能量不能再次重复相加，非空输入的旧高阶失败明确保留|

## 3. 有限构型与上下界

|结果|证明与实现|核验范围及含义|
|---|---|---|
|全部正负相位精确分解|[证明](phase_graph.md)、[实现](phase_graph.py)、[比较](phase_graph_comparison.json)|2≤n≤100的99个家族，与完整旧包络相比提高13项；只优化固定人口、半径的正负相位家族，组过大时预算会拒绝计算|
|残差驱动的明确有理上升|[证明](ascent_audit.md)、[实现](rational_ascent.py)、[见证](ascent_envelope_comparison.json)、[独立审阅](ascent_independent_audit.md)|38个测试形状；指定种子在n≤30各做一步，相位之外另提高13项；n=26至少提高98.30%，n=10、11、12仅微小严格增益|
|有效径向CDF与中心空洞限制|[证明](radial_audit.md)、[实现](radial_distribution.py)、[核验](radial_distribution_verification.json)|任意归一化无碰撞构型由缺额控制；输入可达下界后，全局最优者CDF误差为O(n⁻¹ᐟ³)，至少一点距重心O(n⁻¹ᐟ⁴)；不推出精确原点|
|有宽度径向带及构型类上界|[径向证明](radial_audit.md)与有理单元证书|带宽定理要求库仑尺度平方半径落在[0,2]；一般单元证书不要求所有半径≤2。n=100厚环例有δ≥1225/6，而全局最优δ≤1657/12，故该整类无全局最优者|

设 $`\delta(z)=\log(U_n/\Delta_n(z))`$，径向CDF结论为

```math
\sup_{s\ge0}\left|
\frac{\#\{i:|z_i|^2\le s\}}n-\min\{s/2,1\}
\right|
\le\left(\frac{6\delta(z)}{n^2}\right)^{1/3}+\frac2{\sqrt n}+\frac1n.
```

有理上升的严格增益证书给出目标差额的下界，不能误读成对最优差额的上界。残差零、相位家族内最优、满足径向人口限制，均不构成全局最优认证。

## 4. 已直接纠正的历史状态

旧页基于不完整补充附件所列的待补状态现已撤回。下列证明属于恢复的用户既有成果，不登记为本轮首次发现：

|旧条目|现有证据与判断|
|---|---|
|十点四阶项−510及严格局部最大|完整参数规范、横向消元、旋转切片及余项在原始/修订主稿；原入口已精确核验，仍为全局次优构型|
|完整网格定义及误差证明|原稿量词与证明齐全；历史实际全枚举1874376集合，一般 n 大型网格未执行|
|十点超过121/100的保存整数构型|坐标、中心化、二阶矩、全部距离乘积已完整复算；是严格可达下界，没有确定M₁₀|
|六点特定Hermite majorant障碍|原节点、插值及最大修正齐全；本轮重放verify_hermite.py通过，修正后缺口仍大于3/100；仅限制该方法|

以上状态已同步到[原一般 n 审核](../earlier-results/general_n_proof.md)及[既有成果首页](../earlier-results/README.md)，不是只增设本页链接。[状态记录](STATUS_UPDATE.md)列出当前变更；原始主稿、历史ZIP和旧SHA保留原字节，当前网页及来源索引单独更新。

## 5. 如何复核

当前幂和与八点证书统一入口，在仓库根目录运行：

```bash
python -B -S -O scripts/verify_current.py
```

仅重放所列当前证书；完整原有层数累计验证另见[复现指南](../docs/REPRODUCING.md)。本目录下述命令只验证构型专题。

从[本专题首页](README.md)取得离线包，在根目录运行：

```bash
python -B -S -O verify_all.py
```

入口核对公开字节并重放本专题代码、完整旧下界基线及原Hermite入口。执行范围以[verification_report.json](verification_report.json)为准；不追溯改写历史报告的执行计数，没有冒充运行更大的网格、更多上升步或全部驻点消元。

完整八点材料另见[研究Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)。解析、几何及外部分类输入仍须阅读证明；标准库检查、哈希通过和AI交叉审阅不等于外部同行评审或整篇形式化证明。一般n锐值、精确γ及真实边界修正继续开放；八点全局锐界与本页指定的非零矩剩余系统由新全域证书接合解决。
