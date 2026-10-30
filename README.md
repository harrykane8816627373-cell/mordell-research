# Mordell 距离乘积研究：有限点数构型与一般 n

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[harrykane8816627373-cell](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。当前整理：**2026-10-30**。

统一归一化为

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i\lt j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2.
```

**八点完整无条件锐界、一般非零矩分支，以及一般 n 的简洁全局锐值公式仍未解决。** 最新补充联合全局插值、能量预算与稳定性，将径向约束推进到二维运输、内盘微观覆盖和任意位置人口限制。具体假设见各证明。

## 从这里阅读

- **[当前证明状态总表](configuration-bounds/PROOF_STATUS.md)**：统一查询结论、适用范围、证明和代码；历史待补条目的当前判断也在此汇总。
- **[最新驻点闭合、矩校正空洞与含边界覆盖](stationary-screening/README.md)**：八点必要矩区域收紧到|p₂|²+(27/100)|p₃|²<64/49；守恒矩与空洞用完整Gram联合；整个闭背景盘覆盖推进到O(√(log n/n))，并原字节复核用户新角向包。
- **[此前二维分布、微观覆盖与局部人口](distribution-geometry/README.md)**：单位内盘在n≥64有6/√(n−1)覆盖保证；全二维W₁有效达到n⁻¹ᐟ²阶；任意中心圆盘人口、含边界覆盖与薄条带谱容量均有精确证书。
- **[此前谱人口容量、径向运输与有效二分之一层数](population-transport/README.md)**：全局渐近下端提高到1/2；所有n≥10000有带明确四次根误差的有限式；一万点至少46层；八点薄环带与跨界角向证书。
- **[此前通用层数、局部能量与多点稳定性](collective-stability/README.md)**：全n层数加强为K>n/(π√(n−1))，渐近下端提高至13/30；新增环带角向与多点联合证书。
- **[此前一般 n 构型：规则分量障碍、两点约束与自由移动](symmetry-breaking/README.md)**：10≤n≤100的91项下界全部严格改善；39—42点至少三层，80—85点至少四层，并有任意驻点的正规分量人口限制。
- **[此前径向层数、混合矩能量与人口优化](radial-angular/README.md)**：全n层数下界提高至3√n/22，渐近常数提高至101/250；任意百点全局最优者至少四层；新构造提高87项已存下界。
- **[此前联合上下界与构型约束](joint-bounds/README.md)**：严格圆盘修正、径向与角向复矩可加证书、全 n 内环收缩和最佳两矩子集递推；相对其此前包络提高45项。
- **[最新有限 n 构型与上下界](configuration-bounds/README.md)**：全部正负相位的精确分解、带严格增益证书的有理上升、径向人口限制与构型类上界。
- **[五、六、七点与较早一般 n 证明](earlier-results/README.md)**：作者原稿、全局证明依赖、PDF、TeX及精确复核；[一般 n 审核](earlier-results/general_n_proof.md)已直接更新过期状态。
- **[完整原稿恢复与联合方法](source-review/README.md)**：原始/修订主稿、来源核对、完整旧包络、半径层数、删除余量和混合网格。
- **[全部有限 n 的下界与严格子集余量](finite-n/README.md)** · **[跨点数递推与有效全局网格](global-n/README.md)** · **[线性渐近项与有限圆盘上界](asymptotics/README.md)**。

## 目前证明到了哪里

|问题或结论|当前状态|证明与限制|
|---|---|---|
|五点全局锐值|$`M_5=5^5=3125`$，单位正五边形|[作者二阶条件与矩方法](earlier-results/five_point_method.md)、[PDF](earlier-results/five_point_proof.pdf)|
|六、七点全局锐值|$`M_6=6^{15}/5^{10}`$；$`M_7=7^{21}/6^{15}`$|[主稿](earlier-results/six_seven_proof.tex)、[覆盖审计](earlier-results/six_seven_coverage_audit.md)；全局证明引用 Faugère–Svartz 物理驻点完整分类，七点356次枝未在本地独立重建|
|八点完整锐界|目标为 $`M_8=8^{28}/7^{21}`$；尚未证明|原点加正七边形达到候选值；最新必要椭圆和二阶幂和模<8/7见[驻点闭合](stationary-screening/STATIONARY_MOMENTS.md)；一般非零矩分支仍开放，p_k=Σz_i^k|
|一般 n 的有限上下界与最优层数|严格圆盘修正保持；91项已存下界改善；最新liminf Kₙ/√n≥1/2，并有有效有限式|[严格上界修正](joint-bounds/ANGULAR_PROOF.md)、[最新谱容量与运输](population-transport/README.md)；不确定一般锐值或γ|
|一般 n 最优点的二维分布|固定内盘微观覆盖、有效W₁=O(n⁻¹ᐟ²)、任意中心圆盘人口误差O(n⁻¹ᐟ³)；最新闭背景盘覆盖O(√(log n/n))|[二维分布](distribution-geometry/README.md)、[扩张内盘证明](stationary-screening/EDGE_COVERING.md)；未确定局部晶格|
|严格真子集约束及有效网格|已证；有显式正有理余量与有限全局误差界|[子集证明](finite-n/SUBSET_PROOF.md)、[混合网格](source-review/HYBRID_GRID.md)；网格上界要求完整穷尽规定网格|
|十点中心正九边形及旧缺件条目|四阶严格局部最大、完整网格、保存整数下界和特定 Hermite 障碍现均已核验|[已更新原审核](earlier-results/general_n_proof.md)；十点该局部最大仍全局次优，旧证明确认归用户原有成果|
|确定的线性渐近项|已证 $`\log M_n=cn^2+\frac n2\log n+\gamma n+o(n)`$，$`c=\frac12\log2-\frac14`$|[外部定理、尺度与常数](asymptotics/README.md)；精确 $`\gamma`$ 及真实有限 n 边界修正仍未知|

## 最新驻点闭合与矩校正能量

[新专题](stationary-screening/README.md)联合用户新提供的角向包与仓库现有混合矩方法。驻点恒等式使默认平滑的一级混合矩满足A_(k,1)=2p_k+n⁻¹Σp_ℓp_(k−ℓ)。代入完整Gram后，八点全局最优者在原尺度必有

~~~math
|p_2|^2+\frac{27}{100}|p_3|^2<\frac{64}{49},\qquad |p_2|<\frac87.
~~~

这收紧了提供包中的低阶必要矩区域；消去未定高阶矩还能留下正四次罚项，但没有排除整个非零矩分支。另一条完整Gram将二阶矩匹配、驻点偶极及偏心空洞的相关项联合，得到超过独立空洞项的严格补偿能量。

全局插值则通过随n扩张的内盘，将**整个闭背景盘**的覆盖上界推进到O(√(log n/n))，适用于所有n≥64，不需能量渐近余项。用户原包按字节保存并复算；其中与旧成果重合及通用折扣较弱的部分均在[输入审核](stationary-screening/INPUT_AUDIT.md)注明。

## 此前二维分布与微观几何

[新专题](distribution-geometry/README.md)把用户既有全局最优性、严格子集间距、稳定性矩阵与可达下界直接联合。对任意全局最优者，n≥64时所有|u|≤1均满足

~~~math
\min_i|z_i-u|<\frac6{\sqrt{n-1}}.
~~~

更大固定内盘也有明确起点。加权插值和分离点Gaussian积分使空洞控制达到与点距同阶的微观尺度；没有把这条内盘结论推广至背景边界。

另一条运输路线控制完整二维分布：到均匀半径√2圆盘的W₁距离有有效O(n⁻¹ᐟ²)上界及同阶普适下界。已有百点整数见证还给任意百点全局最优者W₁<0.111199952674，并排除一类径向CDF误差很小但角向集中的构型。

任意圆心的局部Gauss能量另给圆盘人口区间、不交多空洞的可加代价及包含背景边界的O(n⁻¹ᐟ⁴)覆盖；薄矩形谱容量排除任意200点局部最大者在长2.7、宽0.01矩形中容纳20点。四套证明和标准库重放见[专题](distribution-geometry/README.md)。这些结论尚未分类一般n的全局最优者。

## 此前谱容量与有效径向运输

[本轮专题](population-transport/README.md)用角序交错点集提取原稳定性算子的同层容量，再用截断逆半径的Dirichlet运输连接总体径向分布。严格得到：任意全局最优者在全部n≥10000满足

~~~math
K>\frac{\sqrt n}{2}-\frac{11}{20}n^{1/4},\qquad
\liminf_{n\to\infty}\frac{K_n}{\sqrt n}\ge\frac12.
~~~

这将此前渐近下端13/30推进到1/2，且只用已有初等有效可达下界；不需要精确γ或外部线性渐近项的未知有效阈值。明示有限参数族进一步给一万点至少46层、百万点至少487层，均为必要条件，不是精确层数。

同圆谱容量另给所有局部最大适用的全n界K>n/(1+2√(2(n−1)))。允许半径变化的薄带判据排除全部八点半径处于[0.99,1.01]的局部最大，并给正径向方差缺口。相交Green测试及跨界圆盘有理积分扩大了角向能量接受范围，可与同h径向预算联合。证明、四套精确重放和离线包见[本轮专题](population-transport/README.md)。

## 上一轮层数与多点约束

[本轮专题](collective-stability/README.md)联合用户既有严格子集间距、同圆角向装点和总二阶矩，证明任意全局最优者满足K>n/(π√(n−1))>√n/π，适用于全部n≥2。它不要求每层等角，且正确处理单点层。保留整数人口的立方预算，可进一步排除指定人口模式。

新的分箱径向能量方法把各格共享的总层数预算保留下来，以127条有理凸切线严格推出liminf Kₙ/√n≥13/30，将此前0.404提高到约0.433333。其外部线性渐近与三角试探输入保持原适用范围；有限公式保留K+127修正，不能把渐近系数直接用于所有有限n。

局部角向帽函数允许不同环带的偏差分别计入能量，并在覆盖条件下与全纯零级矩正交相加。多点块稳定性则给出两个非空厚区域的整类排除；其中一个区域通过全部旧两点及较强点距、半径条件。四条标准库接受入口、解析证明、来源与离线包见[本轮专题](collective-stability/README.md)。这些必要条件尚未完成一般最优构型的全域覆盖。

## 上一轮自由构型进展

[新专题](symmetry-breaking/README.md)将原稿两环障碍推广为：驻点若包含以重心为中心的正m边形，且另有M>0个非中心点，则m≤M；人口恰好平衡时，剩余点的多项式还须满足稀疏形式。完整余式检验说明上一轮99个所选规则候选中91个并非驻点；其中23个人口模式的全部半径、连续相位均被排除。

从二阶稳定性另得到点距与半径差的联合必要不等式：接近允许最小间距的一对点，其半径差也必须小。它适用于所有局部最大点，并给一个非空八点厚区域的严格排除，但不完成八点全覆盖。

放开每个点自由移动后，整数坐标证书在10≤n≤100全部91项严格提高旧完整包络；27点至少提高104.12%，100点至少提高35.12%，10、11、12点仅有极小改善。新下界联合通用层数定理，证明**39—42点的任意全局最优者至少三层，80—85点至少四层**。候选并未被认证为全局最优；证明、逐项数据、构型图及标准库离线重放均在[新专题](symmetry-breaking/README.md)。

## 上一轮径向—角向进展

[径向—角向新专题](radial-angular/README.md)将任意全局最优构型的不同半径数加强为

```math
K>\frac3{22}\sqrt n\quad(n\ge2),\qquad
\liminf_{n\to\infty}\frac{K_n}{\sqrt n}\ge\frac{101}{250}.
```

证明保留实际半径空隙的能量权重，并引入可调平滑半径；还适用于有总宽度的近似环层。混合矩Gram投影补足普通复幂和的跨层抵消，整数人口与精确空带积分给更强区域排除。

新方法用角平均闭式半径、相对熵损失和整数动态规划联合选择环人口与半径，在2≤n≤100的实际核验中相对完整旧包络提高87项；所选99个固定半径人口家族的连续相位余量均≤百万分之一。由新下界与通用层数界进一步证明：**任意百点全局最大构型至少有4个不同半径**，不要求规则角度。代理目标最优不等于实际全局最优，详见[完整证明和复现](radial-angular/README.md)。

## 更早联合上下界进展

[联合上下界专题](joint-bounds/README.md)证明，对全部n≥2，令

```math
B_n=\left(1+\sqrt{\frac2{n-1}}\right)^2,\qquad
c_n=\frac{(\log2)^2}{\sum_{k=2}^n B_n^k/k},
\qquad M_n\le U_ne^{-c_n}<U_n.
```

这是指定圆盘上界的严格有限点数改善，cₙ很小且趋零，不改变线性渐近常数γ。其证明同时得到径向缺额与各非零复幂和罚项可相加的构型类上界，并给出八点一个p₂p₃≠0子类的严格排除。

下界方面，所有n≥13的原多环构造，无论各环怎样连续旋转，都可通过收缩最内六点环并重新归一化，使乘积提高超过10%。在13≤n≤100的新固定半径正负相位计算中，相对此前完整已存包络提高45项，n=91至少提高24.83%。最佳两矩子集递推另严格改善147组旧公式，但这些递推仍高于圆盘界。证明、范围及标准库重放入口见[专题首页](joint-bounds/README.md)。

## 此前构型与上下界的联合进展

[最新专题](configuration-bounds/README.md)把完整旧多环、格点、矩平均和保存整数构型合并为基线，再加入两类精确改善。三十点的人口6、12、12及平方半径1/5、4/5、8/5，在正负相位家族中给出相对旧两方案的增益

```math
\left(\frac{90155}{84693}\right)^{12}\gt2.
```

全部相位家族在2≤n≤100已核对，相对完整旧包络提高13项；有理一步上升在n≤30又提高相位之外的13项，n=26至少提高98.30%。小残差只控制驻点条件，不能认证全局最优。

此前专题使用的无条件有限圆盘上界为（最新严格修正见上节）

```math
M_n\le U_n:=n^{n/2}\left(\frac{2n}{n-1}\right)^N e^{-N/2}.
```

本轮新增的是径向必要条件和对整类构型的更强上界：库仑缺额给出有效圆盘人口误差、中心空洞限制，以及角度任意的厚环排除证书。此前还已证明任意全局最优构型的不同半径数满足

```math
K\gt\frac{\sqrt n}{12}\quad(n\ge2),\qquad
\liminf_{n\to\infty}\frac{K_n}{\sqrt n}\ge\frac{11}{50}.
```

无需预设规则角度。新旧结果的假设、计数与证明入口统一列于[状态总表](configuration-bounds/PROOF_STATUS.md)。

## 八点研究入口

完整八点研究见 [research-2026-10-27 发布版本](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)及其完整研究附件。已证内容包括 $`p_2=0`$、$`p_3=0`$ 真实驻点分类与值界，两直线、抛物线、四轴反射等类别，七点共圆的稳定驻点归约，以及六点共圆、两例外等距时的排除。连续核连接三个远端径向带，但尚未覆盖全部过渡区域。各自的驻点、稳定性和几何假设必须保留，不能拼接为完整八点定理。

## 下载与复核

每个专题首页均提供证明、Python标准库接受入口、字节清单和离线包。最新[驻点闭合与矩校正研究包](stationary-screening/mordell-stationary-screening-2026-10-30.zip)在解压根目录运行：

```bash
python -B -S -O stationary-screening/verify_all.py
```

其他专题离线包：[此前谱容量与运输](population-transport/mordell-population-transport-2026-10-29.zip)、[此前层数与联合稳定性](collective-stability/mordell-collective-stability-2026-10-29.zip)、[此前自由构型](symmetry-breaking/mordell-symmetry-breaking-2026-10-29.zip)、[此前径向—角向](radial-angular/mordell-radial-angular-2026-10-28.zip)、[更早联合上下界](joint-bounds/mordell-joint-bounds-2026-10-28.zip)、[此前构型与上下界](configuration-bounds/mordell-configuration-bounds-2026-10-28.zip)、[完整原稿与联合方法](source-review/mordell-source-review-2026-10-28.zip)、[有限 n 与子集余量](finite-n/mordell-finite-n-2026-10-27.zip)、[全局网格更正包](global-n/mordell-global-n-corrected-2026-10-28.zip)、[五六七点历史专题](earlier-results/mordell-earlier-proofs-2026-10-27.zip)、[渐近专题 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/general-n-linear-2026-10-27)。历史ZIP保持原字节；当前网页状态及清单另行同步。

完整八点附件为 **mordell-research-complete-2026-10-27.zip**，大小103026838字节，SHA-256：

```text
00970d35ab15ababf310eb53f2e534ca31a2094b9b1c6f3462983956bdaa1654
```

该附件含19个历史及新增证明包、1516个文件和3个必需大型精确对象。2026-10-27实际重放14个套件、24个程序入口全部PASS，另5个由更强包覆盖而跳过。GitHub自动生成的Source code.zip/tar.gz只包含Git目录，不能替代该完整附件。接受程序及覆盖范围见附件内 `docs/REPRODUCING.md`；来源见附件内 `docs/PROVENANCE.md`，引用建议见附件内 `CITATION.md`。

程序验证代数、严格区间及有限覆盖，解析与几何桥梁仍须阅读证明。AI协作交叉审阅不等于外部同行评审或整篇形式化证明；浮点探索和失败搜索不承担不存在性结论。未指定许可证，不作全球优先权声明。

**English:** Jiang Fengcan's ongoing research on finite-point Mordell discriminant maxima. Proved partial theorems, constructive lower bounds, conditional upper bounds and reproducible exact certificates are indexed separately. The complete eight-point sharp bound and a simple general-n sharp formula remain open.
