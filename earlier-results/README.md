# 五、六、七点及一般 n：既有证明专题

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**。本专题从已公开的完整研究归档抽取原稿、证明、精确数据和复核代码，使较早成果能够单独阅读与复现。推导及程序包含AI协作；不作外部同行评审或全球优先权声明。

统一定义 $M_n=\max_{\sum |z_i|^2=n}\prod_{i<j}|z_i-z_j|^2$。

|类别|结论|证明入口|
|---|---|---|
|五点|$M_5=5^5=3125$；等号仅单位正五边形，允许旋转与重排|[作者原稿PDF](five_point_proof.pdf)、[TeX](five_point_proof.tex)、[方法衔接与审阅](five_point_method.md)|
|六点|$M_6=6^{15}/5^{10}$；等号仅原点加半径$\sqrt{6/5}$的正五边形|[六七点主稿](six_seven_proof.tex)、[精确证明片段](six_seven_proof_fragment.tex)|
|七点|$M_7=7^{21}/6^{15}$；等号仅原点加半径$\sqrt{7/6}$的正六边形|[精确复核说明](six_seven_verification.md)、[覆盖审计](six_seven_coverage_audit.md)|
|六七点含原点类|不调用驻点分类的锐界与方差惩罚证明|[独立证明](six_seven_origin_proof.tex)|
|一般$n$|显式有限上下界、唯一最优上界参数、渐近律及下述结构结果|[完整接受证明](general_n_proof.md)、[综合TeX主稿](general_n_notes.tex)|

六、七点等号也允许旋转与重排。**全局证明使用Faugère–Svartz（ISSAC2012）的已发表物理驻点完整覆盖**，附本地有理算术对保留80/174参数根的独立核验；七点356次枝的物理排空仍是明确引用的论文输入。未宣称本地已重建全部分类。[DOI](https://doi.org/10.1145/2442829.2442856)，[作者稿](https://inria.hal.science/hal-00777791/file/FS12.pdf)。不再分发第三方论文全文。

五点主线是作者二阶条件、径向矩、二维投影及四次对数上界方法，不依赖驻点分类；Wang方法的分类交叉检查另保留于five_point/classification_crosscheck/，不能混淆方法归属。

## 一般 n 已证明的具体范围

$n\ge2$，$N=n(n-1)/2$，$P_n=\prod_{k=0}^{n-1}k!$，$S_m(x)=\sum_{k=0}^m x^k/k!$。对任意$\alpha>0$，

$$\frac{n!(n-2)!n^NP_n}{(N+n-2)!}\le M_n\le\frac{P_nS_{n-1}(\alpha)^n}{\alpha^N}.$$

上界函数有唯一最优正参数，且

$$\log M_n=\left(\frac{\log2}{2}-\frac14\right)n^2+\frac n2\log n+O(n).$$

还包括：高斯整数构造在二次指数尺度的渐近最优性；多环距离乘积公式；十点超过中心正九边形113/100的精确反例；$5\le n\le9$中心正多边形的严格局部极大性；$n\ge11$显式改善；Hermite共线子类精确值；全局极大点真子集内部矩的严格必要条件；依赖JPT高度界的有限代数选择定理（巨量枚举未执行）。

这里没有宣称一般$n$的简单闭式全局解。历史提纲general_n/historical_outline/general_n_proposed.tex中十点四次项、未完整定义的网格误差界及无坐标数值改善仍待补。接受与未接受条目见完整audit.md。

## 下载完整专题代码与证书

[下载专题证明包](mordell-earlier-proofs-2026-10-27.zip)。网页提供主要论文、证明及五点脚本；全部依赖、六七点精确数据和目录结构在专题ZIP内。下述命令及未设网页链接的相对路径以解压ZIP为准。

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
