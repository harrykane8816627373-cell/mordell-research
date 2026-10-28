# 当前证明状态与统一入口

研究资料与项目发起：蒋丰灿（Fengcan Jiang）；推导与代码包含 OpenAI Codex 协作。更新于2026-10-28。此页汇总当前判断，定理的全部假设和证明以所链接正文为准。

```math
M_n=\max_{\sum_i|z_i|^2=n}\prod_{i\lt j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2.
```

## 1. 小点数、构型与未解决问题

|范围|已取得的结果|证明入口与适用边界|
|---|---|---|
|五点|$`M_5=3125`$，等号单位正五边形，允许旋转、重排|[作者方法](../earlier-results/five_point_method.md)、[PDF](../earlier-results/five_point_proof.pdf)、[TeX](../earlier-results/five_point_proof.tex)；主线为二阶条件、径向矩、投影及四次对数上界|
|六、七点|$`M_6=6^{15}/5^{10}`$、$`M_7=7^{21}/6^{15}`$，等号原点加相应正多边形|[主稿](../earlier-results/six_seven_proof.tex)、[复核](../earlier-results/six_seven_verification.md)、[覆盖审计](../earlier-results/six_seven_coverage_audit.md)；全局证明依赖 Faugère–Svartz 完整物理驻点分类，七点356次枝的物理排空为论文输入，未在本地独立重建|
|六、七点含原点类|不调用上述分类的锐界与方差惩罚|[独立证明](../earlier-results/six_seven_origin_proof.tex)；结论仅适用于含原点类|
|八点|候选值 $`8^{28}/7^{21}`$ 可达；若干低矩、几何与稳定分支已证|[完整八点归档](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)；完整无条件上界及一般 $`p_2p_3\ne0`$ 分支未解决|
|中心正多边形类|5≤n≤10模旋转严格局部最大；n=10使用四阶项，n≥11有显式增益扰动|[一般 n 审核第6、7节](../earlier-results/general_n_proof.md)；五点中心正方形、十点中心正九边形仍全局次优|
|一般 n|已有有限上下界、全局网格包围、结构与渐近定理|没有得到所有有限 n 的简洁锐值公式；有限代数刻画不代表已执行巨大消元或枚举|

八点已证的 $`p_2=0`$、$`p_3=0`$ 分类保留真实驻点及相应稳定约束；七点共任意圆、六共圆且两例外等距、两直线等类别各有假设。连续核已连接部分远端带，尚有过渡区域未覆盖。不能去掉条件后把这些成果拼接为完整八点证明。

## 2. 既有一般 n 方法与后续联合

|结论|当前状态与范围|证明、数据或代码入口|
|---|---|---|
|高斯平均、截断指数上界及唯一最优参数|已证；唯一性针对上界参数，不是构型唯一性|[一般 n 审核第2节](../earlier-results/general_n_proof.md)|
|格点、有理多环和保存整数构型|每个 n 有可计算可达下界；格点在二次指数尺度渐近最优|[原始主稿](../source-review/original_main.tex)、[修订主稿](../source-review/revised_main.tex)、[来源核对](../source-review/SOURCE_REVIEW.md)|
|三角整数格与加权矩|已证有限有理下界；前三阶矩实际计算到n=7，格点比较至n=30|[有限 n 专题](../finite-n/README.md)、[矩证明](../finite-n/MOMENT_METHOD.md)、[格点数据](../finite-n/triangular_lattice_values.json)|
|严格真子集内部矩约束|对全局最优构型及2≤m<n严格成立；已加强为显式正有理余量|[完整证明](../finite-n/SUBSET_PROOF.md)；单点两边为零，全体点为等式，不推广到所有驻点|
|标准化最大值相邻严格下降|原稿已有完整证明；后续加入随机子集四阶方差折扣|[来源核对](../source-review/SOURCE_REVIEW.md)、[全局方法专题](../global-n/README.md)|
|完整有效网格|原稿已有定义与误差证明；后续联合稳定矩阵、点距与较小整数箱改善|[原审核第9节](../earlier-results/general_n_proof.md)、[混合网格](../source-review/HYBRID_GRID.md)；上界要求完整覆盖，抽样只给下界|
|不同半径层数|任意全局最优者有 $`K\gt\sqrt n/12`$；渐近下端为11/50|[联合专题](../source-review/README.md)；不要求规则角度|
|确定的线性渐近项|$`\log M_n=(\frac12\log2-\frac14)n^2+\frac n2\log n+\gamma n+o(n)`$|[渐近专题](../asymptotics/README.md)说明 Sandier–Serfaty 外部输入与尺度校正；精确γ、真实有限 n 边界修正未解决|
|全 n 圆盘上界|$`M_n\le U_n=n^{n/2}(2n/(n-1))^N e^{-N/2}`$|[证明与三项非负缺额](../asymptotics/README.md)；经典普适界不登记为项目原创，本轮未刷新无条件Uₙ|

## 3. 本轮有限构型与上下界

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

从[本专题首页](README.md)取得离线包，在根目录运行：

```bash
python -B -S -O verify_all.py
```

入口核对公开字节并重放本专题代码、完整旧下界基线及原Hermite入口。执行范围以[verification_report.json](verification_report.json)为准；不追溯改写历史报告的执行计数，没有冒充运行更大的网格、更多上升步或全部驻点消元。

完整八点材料另见[研究Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)。解析、几何及外部分类输入仍须阅读证明；标准库检查、哈希通过和AI交叉审阅不等于外部同行评审或整篇形式化证明。一般n锐值、八点完整锐界、一般非零矩分支、精确γ及真实边界修正继续开放。
