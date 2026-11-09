# 研究索引：当前主线、方法依赖与历史版本

## 十二点301/300上界、解析梯度与零和核

当前最紧结论为V₁₂<C₁₂≤M₁₂<(301/300)V₁₂，候选比值>300/301（约99.6678%）。见[完整229113节点证书](../twelve-point-kernel-refinement-2/README.md)、[精确组合梯度](../jacobian-gradient/README.md)、[九点子集质量与三点必要式](../twelve-point-subset-domain/README.md)、[Ward零和核构造](../stationary-kernel/README.md)。原3%与15%证书作为历史证据保留。

## Jacobian辅助与十二点3%全局界

上一轮结论为V₁₂<C₁₂≤M₁₂<103V₁₂/100，真实候选比例>100/103。见[新完整证书](../twelve-point-kernel-refinement/README.md)、[Jacobian配置及真实计算](../jacobian-assist/README.md)、[研究记录](CANDIDATE_KERNEL_RESEARCH.md)。87385节点完整覆盖只需标准库；原15%界作为历史证据保留。

2026-11-08。返回[首页](../README.md) · [当前证明状态](../configuration-bounds/PROOF_STATUS.md)。

首页与证明总表保留当前最强结论；本页集中历次专题入口。各专题中的历史比较、较弱常数和运行计数按原日期理解，不能把多条有共同预算的能量界直接相加。

## 候选范围锁定与正定核

十二点已取得全局定量范围V₁₂≤M₁₂<23V₁₂/20；锐界研究继续优先12点，其次10点和19点。完整接口与范围见[候选与核研究记录](CANDIDATE_KERNEL_RESEARCH.md)。

|专题|已接受结果与边界|
|---|---|
|[十二点非锐全局上界](../twelve-point-kernel-upper/README.md)|增加21/10000点对余量后的连续核证书，37053节点完整覆盖；现有真实候选值与全局值之比>20/23，未确定锐值或全部最大构型|
|[候选筛选与精确对称](../candidate-global-screen/README.md)|10—19点共1000次有限搜索记录；真实10、12、19点候选精确D₂/D₃/D₆对称及15、14、22类点对，不推出全局最大者必须对称|
|[十点必要两点域](../maximizer-domain/README.md)|八、九点锐界与二阶Schwarz加强全部全局最大点的半径、补集质量与非线性点对约束；不适用于任意驻点|
|[整区域候选值上界](../candidate-basin/README.md)|对10、12、19点候选轨道周围的明确连续区域，完整区间Hessian证明目标不超过候选；区域补集尚未覆盖|
|[候选适配正定核](../candidate-kernel/README.md)|正定矩块、刚性零集与两点求和接口；保存原始无余量SDP拟合的精确失败点对；十二点加余量的连续证书见上项|

## 11月8日继续研究

|专题|新增结果及范围|
|---|---|
|[四环人口共振](../four-ring-rigidity/README.md)|最大人口唯一、严格人口预算、低次数加法共振及三同一异的加倍规则；仅必要条件|
|[非线性角向缺额](../phase-deficit-sharp/README.md)|完整正弦/对数下端、标量对偶与291个Bernstein叶；指定家族全部相位、整个半径盒|
|[16点完整家族排除](../family-exclusion-nov8/README.md)|利用指数正项与实际整数见证，排除所列四正方形环及半径扰动盒；不覆盖全部16点构型|
|[高点数真实局部极值](../validated-maxima-nov8/README.md)|25、35、42、46点严格局部极大轨道及目标值区间；固定引用前轮精确验证器|

本轮推导和后续边界见[11月8日记录](NOV8_RESEARCH.md)。

## 11月5日一般点数研究

|专题|本轮新增结论与边界|
|---|---|
|[三环刚性与旋转阶](../ring-rigidity/README.md)|三条互异半径正规环驻点必须等人口；阶m旋转对称局部最大者满足floor(m²/4)≤2(n−1+h)，允许任意相位|
|[循环分数预算](../phase-packing/README.md)|精确原始/对偶证书、一般奇加权循环与整个半径盒；指定规则环家族上界，不是无约束Mₙ上界|
|[自由构型下界](../finite-search-nov5/README.md)|16项整数下界严格改善旧基线，最大46点约4.93%；搜索不是全局穷尽|
|[驻点存在与完整稳定性](../validated-maxima/README.md)|8个严格局部极大轨道与15点不稳定轨道；整数区间压缩与完整Hessian，不是全部驻点分类|

各方法的逻辑接口与未解决问题见[本轮记录](NOV5_RESEARCH.md)。

## 小点数全局结果与历史八点路线

[九点专题](../nine-point-global/README.md)补齐 Ji–Xu–Yu 的九点定理、固定来源和实际复现范围。[整理前首页](HISTORY_2026-11-04.md)保留完整阶段叙述。下表局部方法的“开放”仅描述各自当时的覆盖，不再表示八点全局锐界开放。

|专题|作用及其边界|
|---|---|
|[八点全局锐界与非零矩排除](../eight-point-global/README.md)|Ji–Xu–Yu固定核证书复现，联合已有严格子集定理接合；质量7下Δ≤7⁷，等号唯一中心正七边形；不调用七点分类|
|[11月4日局部续算记录](../joint-refinement/NOV4_PROGRESS.md)|全局核接入前的局部阶段：全局最优者平方半径<231/160，联合旧椭圆后<23/16；固定中心平移子集预算；原区域仍0–9，旧14盒没有新增整盒关闭，导入两棵部分树后28开放盒；新增10004步骤已独立重放；当前执行记录区分复用与新算术|
|[七点–Schwarz半径与共同子集预算](../joint-refinement/SUBSET_BUDGET.md)|19dbbd3阶段：平方半径<119/80；完整闭合原编号10；上轮30个开放盒关闭16个，余10个原域、14个开放盒；合并424413节点、212193排除叶|
|[点差、复余项与二列稳定性](../joint-refinement/EDGE_REMAINDERS.md)|8be9f91阶段：新增35001个处理步骤、17491个排除叶；原11域仍未闭合，细分后30个开放叶；总账390821节点、195381排除叶|
|[附件合并与力方程续算](../joint-refinement/FORCE_CLOSURES.md)|768598a阶段：复核附件6项闭合，再完整闭合1282、7062、7954、9814；原15个未解域减为11个，合并355801节点、177890排除叶；合并附件既有729/1000内点门槛|
|[完整四阶稳定性与旧椭圆](../joint-refinement/JOINT_STABILITY.md)|父阶段：3/4内点条件、−1/64负曲率间隙与径向预算；原父证书保留|
|[精确质量范围与参数包络](../joint-refinement/MASS_DUAL_ENVELOPE.md)|父阶段：共同参数包络与中心化质量范围；原证书、原算术与原报告继续保留|
|[中心化切片与质量对偶](../joint-refinement/CENTERED_MASS.md)|前一阶段：保留两个精确坐标和等式，交换凸性保证顶点枚举；原版本节点继续重放|
|[全局对数二次界与角点上包](../joint-refinement/GLOBAL_QUADRATIC.md)|前一阶段：对角二次与角点筛选；原版本节点、旋转饱和及旧 Taylor 结果继续保留|
|[完整物理系统与部分覆盖](../joint-closure/README.md)|前一阶段：S₁₆与完整 Hessian 等价合并为无分母 PSD 多项式矩阵；冻结60021节点、21开放区域；含一个精确非零矩驻点的高模式不稳定性证书|
|[径向—角向局部目标排除](../radial-gap/README.md)|历史局部结果：E≤1/1200时 log(Δ_q/7⁷)≤−E/32；不需要驻点、稳定性或预设半径；保留七点锐界与局部能量输入，该局部方法未覆盖E较大区域|
|[上传统一包的原始快照](../radial-gap/inherited/mordell-eight-point-unified-2026-11-01.zip)|按原字节保留，49条哈希及14个标准库入口重放通过；含径向 Ward、Schur–Ward、投影稳定性及根配对资料，导入解析结论未全部独立复证|
|[非零矩全局代数耦合](../algebraic-gap/README.md)|此前：旧椭圆内的真实驻点挑战者满足 \|p₂\|≥10⁻³⁶\|p₃\|；重放完整旧整数证书，排除细楔形区域；常数很弱，该方法未完成全域覆盖|
|[八点历史能量路线](../eight-point/README.md)|2026-11-01第二阶段：直接对数余项给 d≤1/5 的排除；Parseval 根估计将六矩门槛提高五倍，p₇ 无额外限制；该阶段未完成全域覆盖|
|[旧椭圆与当前完整归约](EIGHT_POINT_REDUCTION.md)|三个尺度的准确换算；记录旧路线中椭圆内 p₂p₃≠0 的联合条件；新全域核已排除达到候选值的剩余系统|
|[有限幂和](../power-moment/README.md)|2026-10-31：Newton 闭合、Hankel 目标、稀疏结构、旧局部刚性与11盒；旧局部盒已包含于新管域|
|[驻点闭合与空洞](../stationary-screening/README.md)|驻点混合矩、八点必要椭圆、完整 Gram 与背景边界覆盖|
|[完整八点 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)|历史零矩、共圆、反射等分类和大型精确对象；各分类保留驻点、稳定性或几何前提|

## 一般 n 的上界、能量与层数

[一般点数结构专题](../general-structure/README.md)给出连续相位相容性与环对缺额的独立研究入口；它与下表已有上下界、层数方法互补。

下表按方法衔接排列；当前有限有效层数及两种渐近输入见最后一项。

|专题|主要接口|
|---|---|
|[线性渐近项](../asymptotics/README.md)|外部展开的尺度核对、圆盘上界、γ 的界|
|[联合上下界](../joint-bounds/README.md)|有限严格圆盘修正、复矩罚项、内环收缩与子集递推|
|[径向—角向](../radial-angular/README.md)|混合矩、半径层数与人口优化；百点至少四层|
|[联合稳定性](../collective-stability/README.md)|全 n 层数、局部角向能量与多点块条件|
|[谱人口与运输](../population-transport/README.md)|局部谱容量、有效二分之一层数、薄环带及跨界测试|
|[二维分布与微观几何](../distribution-geometry/README.md)|W₁、内盘覆盖、任意中心人口与条带容量|
|[线性正缺额](../local-defect/README.md)|近邻与空洞产生初等有效线性缺额；径向—角向共同预算|
|[重叠积分与短边预算](../overlap-budget/README.md)|保留近相切积分，线性缺额系数提高37倍|
|[层周长与层分布](../layer-circumference/README.md)|凸点对损失与真实层周长、窗口层密度|
|[共同变分框架](../joint-variational/README.md)|移位 Gram、层周长和点对资源的统一对偶|
|[高斯角向探针](../gaussian-angular/README.md)|完整交叉 Gram、多尺度角向密度及相似变换|
|[可变平滑与环带人口](../annular-layer/README.md)|当前最强层数与百点条件界；初等有效结论和保留外部输入的渐近结论分开|

## 可达构型、全局网格与作者原稿

[十点转换证明](TEN_POINT_TRANSITION.md)以整数算术复核印刷坐标的严格比较；n≥11 的不稳定性保留文献输入，不把十点见证称为最优构型。

|专题|主要接口|
|---|---|
|[自由构型与规则分量障碍](../symmetry-breaking/README.md)|10—100点保存整数下界；正规分量人口与两点稳定性|
|[构型与上下界](../configuration-bounds/README.md)|相位精确分解、有理上升、径向人口与统一状态表|
|[有限 n 与严格子集](../finite-n/README.md)|有限矩下界和严格真子集的显式余量|
|[跨点数递推与有效网格](../global-n/README.md)|有限穷尽误差界；只有实际完整枚举才给相应上界|
|[完整原稿恢复](../source-review/README.md)|原始及修订主稿、来源核对、原有网格、十点局部证明和保存整数构型|
|[五、六、七点](../earlier-results/README.md)|五点作者方法及六七点证明；后两项全局性保留外部物理驻点分类依赖|

## 下载与历史追溯

各专题首页链接自己的离线 ZIP、固定字节清单和实际运行报告。保留的[能量阶段八点包](../eight-point/mordell-eight-point-energy-2026-11-01.zip)包含能量证明、两阶段管域及重建参考谱所需代码；它不是全部研究资料的替代品。

当前 Taylor 续算的代码、精确森林与报告见[续算专题](../joint-refinement/README.md)。下列离线包保留各自发布阶段的快照，较早包中的21个开放盒计数不会追溯改写。

历史[联合系统部分覆盖 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-joint-partial-2026-11-01)提供 `mordell-eight-point-joint-partial-2026-11-01.zip`，528151字节。包含等价多项式模型、真实驻点盒、60021节点的部分覆盖及必需旧输入；解压后只需标准库，可重放17个旧入口和3个新入口。**该历史树保存21个开放区域。** 全局锐界由新的八点全域核承担。 SHA-256：

```text
5df931d78ee48c95ddfcc4ac9acbeff4aa4ceb2f191fdef3356dfffc36c938f4
```

此前[径向—角向 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-radial-2026-11-01)提供 `mordell-eight-point-radial-2026-11-01.zip`，169952字节，包含1/1200证明、上传统一 ZIP 的原始字节及必需的旧能量输入。解压后一条命令可重放17个所列入口，全部只需标准库。SHA-256：

```text
7fc0ac0d299f3c79cc320611b60b0632ead25ff53f0ff97c17c6e91137dc9f49
```

此前[全局代数耦合 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-algebraic-2026-11-01)提供 `mordell-eight-point-algebraic-2026-11-01.zip`，14998359字节，仅携带对应证明与3个必需旧源文件，可独立重放全部14179个整数系数和有理常数。大型整数对象放在 Release，不在 Git 目录重复提交。SHA-256：

```text
41c4aa95cdd4386f0eac7625e4c5cb4993c63d5f2a3c364002d05c2b1447f49b
```

完整八点历史附件 **mordell-research-complete-2026-10-27.zip** 位于[研究 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)，大小103026838字节，SHA-256：

```text
00970d35ab15ababf310eb53f2e534ca31a2094b9b1c6f3462983956bdaa1654
```

它包含19个历史及新增证明包、1516个文件和3个必需大型精确对象。2026-10-27报告的14个套件、24个入口是当时实际执行范围，不表示本次重跑了那些大型对象。GitHub 的 Source code.zip/tar.gz 仅包含 Git 目录。

[整理前的完整首页](https://github.com/harrykane8816627373-cell/mordell-research/blob/c954da2b1e4a1dbdb9ea698cb9902eeffc85fdfc/README.md)可按原提交追溯。本次合并重复的导航和历次摘要；旧证明、历史 ZIP、来源与精确证书不删除。未发现待合并 PR 或冗余远端分支：整理前只有 main。
