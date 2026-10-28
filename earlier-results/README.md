# 五、六、七点及一般 n：既有证明专题

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**。本专题从已公开的完整研究归档抽取原稿、证明、精确数据和复核代码，使较早成果能够单独阅读与复现。推导及程序包含AI协作；不作外部同行评审或全球优先权声明。

**状态更新：2026-10-28。** 完整原始/修订源码现已取得并核验，原先因补充梗概缺件而保留的十点四阶论证、完整网格、十点>121/100保存整数构型，以及六点Hermite方法障碍均已更新为已核验。[当前一般n审核](general_n_proof.md)、[完整来源核对](../source-review/SOURCE_REVIEW.md)、[原字节完整修订稿](../source-review/revised_main.tex)。本目录2026-10-27 ZIP仍是历史冻结快照，网页当前状态与来源索引单独更新。

统一定义：Mₙ 为满足 ∑|zᵢ|² = n 的 n 点配置中，距离乘积 Δₙ = ∏ᵢ<ⱼ |zᵢ − zⱼ|² 的最大值。

|类别|结论|证明入口|
|---|---|---|
|五点|$M_5=5^5=3125$；等号仅单位正五边形，允许旋转与重排|[作者原稿PDF](five_point_proof.pdf)、[TeX](five_point_proof.tex)、[方法衔接与审阅](five_point_method.md)|
|六点|$M_6=6^{15}/5^{10}$；等号仅原点加半径 √(6/5) 的正五边形|[六七点主稿](six_seven_proof.tex)、[精确证明片段](six_seven_proof_fragment.tex)|
|七点|$M_7=7^{21}/6^{15}$；等号仅原点加半径 √(7/6) 的正六边形|[精确复核说明](six_seven_verification.md)、[覆盖审计](six_seven_coverage_audit.md)|
|六七点含原点类|不调用驻点分类的锐界与方差惩罚证明|[独立证明](six_seven_origin_proof.tex)|
|一般 n|显式有限上下界、参数最优性、线性渐近项、有效网格、精确代数框架与下述结构结果|[当前状态审核](general_n_proof.md)、[完整原修订稿](../source-review/revised_main.tex)、[后期不分类增补](general_n_notes.tex)|

六、七点等号也允许旋转与重排。**全局证明使用Faugère–Svartz（ISSAC2012）的已发表物理驻点完整覆盖**，附本地有理算术对保留80/174参数根的独立核验；七点356次枝的物理排空仍是明确引用的论文输入。未宣称本地已重建全部分类。[DOI](https://doi.org/10.1145/2442829.2442856)，[作者稿](https://inria.hal.science/hal-00777791/file/FS12.pdf)。不再分发第三方论文全文。

五点主线是作者二阶条件、径向矩、二维投影及四次对数上界方法，不依赖驻点分类；Wang方法的分类交叉检查另保留于five_point/classification_crosscheck/，不能混淆方法归属。

## 一般 n 已证明的具体范围

取 n ≥ 2，N = n(n−1)/2，Pₙ = 0!1!⋯(n−1)!，Sₘ(x) = ∑ₖ₌₀ᵐ xᵏ/k!。对任意 α > 0，

$$\frac{n!(n-2)!n^NP_n}{(N+n-2)!}\le M_n\le\frac{P_nS_{n-1}(\alpha)^n}{\alpha^N}.$$

上界函数有唯一最优正参数，且

$$\log M_n=\left(\frac{\log2}{2}-\frac14\right)n^2+\frac n2\log n+O(n).$$

O(n)是这里保留的历史初等精度，当前已加强为确定的γn+o(n)，见[渐近专题](../asymptotics/README.md)。还包括：高斯整数及全n有理多环构造；十点超过中心正九边形113/100的简例和更强的>121/100保存整数构型；5≤n≤10中心正多边形的严格局部极大性（十点使用四阶横向消元），但十点全局次优；n≥11显式改善；Hermite共线子类精确值及六点特定majorant的障碍；全局真子集内部矩约束；有限代数选择与关键值增广后有限精确证书。巨量枚举没有据此宣称已执行。

旧提纲相应条目已经由完整原稿、原数据及精确重放补齐；此前“仍待补”的当前状态撤回。通用全纯谱层有限不锐、物理迹等价和有限增广实SOS的适用条件见[来源核对](../source-review/SOURCE_REVIEW.md)。一般n简洁全局锐值公式及八点完整锐界仍未得到。

## 下载完整专题代码与证书

[下载2026-10-27历史专题证明包](mordell-earlier-proofs-2026-10-27.zip)。其字节及旧SHA保持不变；当前补齐的原稿、核验和后续方法见[source-review完整离线包](../source-review/mordell-source-review-2026-10-28.zip)。网页提供主要论文、证明及五点脚本；历史ZIP内有原先的依赖、六七点精确数据及目录结构。下述历史复核命令及未设网页链接的相对路径以解压该ZIP为准。

## 离线复核

解压本专题ZIP，在根目录运行：

```bash
python -B -S verify_prior.py
```

该入口先检查source_manifest.json，再在临时副本执行10个标准库程序入口：五点3个，六七点4个，一般n3个。复核记录写入verification/replay.json；解析论证和外部分类输入仍须阅读证明文本。可选五点分类交叉审计另需SymPy，不由默认入口重复执行。

本专题修复了五点作者联系脚注的TeX排版，并以完全相同的原字节补回六七点旧清单要求的proof_fragment.tex别名。未改变数学正文、算法或精确系数。ZIP内source_provenance.json记录来源与变换，ZIP内source_manifest.json绑定实际公开字节。

完整八点问题与一般非零矩分支仍开放。[完整研究版本](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)。本专题没有另行指定许可证。


新五点PDF为15页；编译诊断及首尾页阅读检查通过。专题接受程序运行结果见[verification_report.json](verification_report.json)。

专题ZIP SHA-256：`1faf966c8127b4ffa2c5aa0e439e0c2c1d724afbe3928cb4afaffc7279871e2a`；大小1398131字节；112文件；10个接受入口已从专题副本重放全部PASS。见[SHA256SUMS.txt](SHA256SUMS.txt)。

## 一般 n 的后续加强

[线性渐近项、圆盘上界与八点缺额](../asymptotics/README.md)给出确定的γn+o(n)项、更强的严格常数界和全部有限n的显式上界；本页保留较早证明的历史入口。

[完整原稿恢复与方法联合](../source-review/README.md)保留原稿完整网格和严格递推的正确归属、全n有理多环与保存构型、有效半径层数及后续误差改善。[构型与上下界继续推进](../configuration-bounds/README.md)给出当前的新研究及复核入口；有限样本比较或小残差不替代全局最优性证明。

