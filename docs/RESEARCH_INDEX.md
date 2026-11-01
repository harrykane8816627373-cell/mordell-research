# 研究索引：当前主线、方法依赖与历史版本

2026-11-01。返回[首页](../README.md) · [当前证明状态](../configuration-bounds/PROOF_STATUS.md)。

首页与证明总表保留当前最强结论；本页集中历次专题入口。各专题中的历史比较、较弱常数和运行计数按原日期理解，不能把多条有共同预算的能量界直接相加。

## 八点与有限幂和

|专题|作用及其边界|
|---|---|
|[全局对数二次界与角点上包](../joint-refinement/GLOBAL_QUADRATIC.md)|当前：新增12531次对角二次排除、7321次角点排除；合并290019节点、144991排除叶，19个连续区域仍开放；保留旋转饱和与旧 Taylor 结果|
|[完整物理系统与部分覆盖](../joint-closure/README.md)|前一阶段：S₁₆与完整 Hessian 等价合并为无分母 PSD 多项式矩阵；冻结60021节点、21开放区域；含一个精确非零矩驻点的高模式不稳定性证书|
|[径向—角向全局排除](../radial-gap/README.md)|最新：E≤1/1200时 log(Δ_q/7⁷)≤−E/32；不需要驻点、稳定性或预设半径；保留七点锐界与局部能量输入，E较大时仍开放|
|[上传统一包的原始快照](../radial-gap/inherited/mordell-eight-point-unified-2026-11-01.zip)|按原字节保留，49条哈希及14个标准库入口重放通过；含径向 Ward、Schur–Ward、投影稳定性及根配对资料，导入解析结论未全部独立复证|
|[非零矩全局代数耦合](../algebraic-gap/README.md)|此前：旧椭圆内的真实驻点挑战者满足 \|p₂\|≥10⁻³⁶\|p₃\|；重放完整旧整数证书，排除细楔形区域；常数很弱、全域仍开放|
|[当前八点攻关](../eight-point/README.md)|2026-11-01第二阶段：直接对数余项给 d≤1/5 的排除；Parseval 根估计将六矩门槛提高五倍，p₇ 无额外限制；全域仍开放|
|[旧椭圆与当前完整归约](EIGHT_POINT_REDUCTION.md)|三个尺度的准确换算；椭圆内的 p₂p₃≠0 分支仍是攻关对象，新门槛继续与全部旧约束取交集|
|[有限幂和](../power-moment/README.md)|2026-10-31：Newton 闭合、Hankel 目标、稀疏结构、旧局部刚性与11盒；旧局部盒已包含于新管域|
|[驻点闭合与空洞](../stationary-screening/README.md)|驻点混合矩、八点必要椭圆、完整 Gram 与背景边界覆盖|
|[完整八点 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)|历史零矩、共圆、反射等分类和大型精确对象；各分类保留驻点、稳定性或几何前提|

## 一般 n 的上界、能量与层数

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

最新[联合系统部分覆盖 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/eight-point-joint-partial-2026-11-01)提供 `mordell-eight-point-joint-partial-2026-11-01.zip`，528151字节。包含等价多项式模型、真实驻点盒、60021节点的部分覆盖及必需旧输入；解压后只需标准库，可重放17个旧入口和3个新入口。**21个开放区域仍未排空。** SHA-256：

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
