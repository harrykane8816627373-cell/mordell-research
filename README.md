# Mordell 距离乘积研究：五至八点与一般 n

研究资料与项目发起：**蒋丰灿**（[harrykane8816627373-cell](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。

**完整八点无条件锐界，以及一般 $p_2p_3\ne0$ 分支，仍未解决。** 本项目记录进行中的部分定理及可复核证明，不宣称主问题已经完成。

在归一化 $\sum_{i=1}^8|z_i|^2=8$ 下，目标是

$$\Delta_8=\prod_{i<j}|z_i-z_j|^2\stackrel{?}{\le}\frac{8^{28}}{7^{21}}.$$

目标等号配置为原点加半径平方 $8/7$ 的正七边形；$p_k=\sum_i z_i^k$。

## 一般 n 更新：确定的线性项与有限上界

[查看完整证明、PDF、精确常数与代码](asymptotics/README.md)。将既有 Sandier–Serfaty 重整化能量定理精确应用于矩约束，得到

```math
\log M_n=\left(\frac12\log2-\frac14\right)n^2+\frac n2\log n+\gamma n+o(n).
```

结合三角晶格计算与经典圆盘平滑界，本项目已核验的区间收紧为

```math
0.402178895223365173209539459449<\gamma<0.403426409720027345291383939271.
```

另有对每个 n ≥ 2 的显式有限上界

```math
M_n\le n^{n/2}\left(\frac{2n}{n-1}\right)^{n(n-1)/2}e^{-n(n-1)/4}.
```

新稿明确列出外部定理、精确尺度校正及三项非负缺额；八点目标还需要一个未证明的正缺额下界。这里没有确定本项目的精确 γ 或真实有限 n 边界修正，也没有宣称经典普适界为本项目原创。

## 已有五、六、七点及一般 n 证明

[打开既有证明专页](earlier-results/README.md)，含五点作者原稿PDF、六七点证明、一般n完整接受证明、TeX源码、精确证书和独立离线专题包。

|点数或范围|已证结论|入口|
|---|---|---|
|五点|$M_5=3125$，等号单位正五边形；作者二阶条件与矩方法|[五点PDF](earlier-results/five_point_proof.pdf)|
|六点|$M_6=6^{15}/5^{10}$，等号原点加正五边形|[六七点证明](earlier-results/six_seven_proof.tex)|
|七点|$M_7=7^{21}/6^{15}$，等号原点加正六边形|[覆盖依赖与精确复核](earlier-results/six_seven_verification.md)|
|一般$n$|有限上下界、唯一最优上界参数、渐近律、高斯整数构造、多环与稳定性等定理|[一般n证明](earlier-results/general_n_proof.md)|

此处$M_n$采用$\sum|z_i|^2=n$。六七点全局证明明确使用Faugère–Svartz已发表完整物理分类覆盖，附保留80/174参数根的独立有理算术复核；另有不调用分类的六七点含原点类证明。一般n的历史待补项目逐条区分，没有宣称简单闭式全局解。

## 2026-10-27 的八点研究进展

- **六共圆与等距例外：** 八点真实驻点若六点共任意圆，另两点到该圆圆心等距，则只能是正八边形，因而不稳定。无需预设反射；证明含 69 节点有理多项式组合证书和 624 项标准库检查。两个例外到圆心不等距的主支仍开放。
- **连续核的跨域连接：** 以共享的 12 维实际特征，将既有局部核连接到半径平方 $1/2,4/7,9/10$ 附近的三个远端带，连续角度与径向带宽由 4285 项严格区间检查覆盖；不需要点距分离假设。过渡半径及其他远端区域尚未全部覆盖。
- **既有几何与低矩结果：** 保留 $p_2=0$、$p_3=0$ 分支的真实驻点分类及相应值界，以及两直线、抛物线、四轴反射等已证类别。七点共任意圆的八点稳定驻点只可能是目标构型。各结论的驻点、稳定性或全局极大假设分别见完整证明索引。
- **一般分支的结构归约：** 保留真实点、共轭配对及完整稳定条件的矩阵与多项式约束已经建立；一般系统尚未排空。若干辅助核证书族的不可能性是对证明方法的限制，不能当成原极值问题的结论。

## 完整资料与复核

完整材料见 [research-2026-10-27 发布版本](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)，请下载其中的 **mordell-research-complete-2026-10-27.zip** 附件，并用同页 `SHA256SUMS.txt` 核验。

网页仓库含成果首页、earlier-results既有证明专页及专题源码包；完整八点历史文件和大型数据仍以本Release研究附件为准。GitHub自动生成的Source code.zip / tar.gz只包含Git目录，不能替代完整八点研究附件。

完整研究归档包含 19 个历史及新增证明包、1500 个证明源码文件、3 个必需的大型精确数学对象、统一入口、来源清单及复核报告，共 1516 个文件。大型对象是验收必需数据，未删减数学系数。

2026-10-27 从公开整理副本重放：**实际执行 14 个套件、24 个程序入口全部 PASS**；另有 5 个旧套件由更强包覆盖而跳过。此次包含可选 SymPy 样例。不能将目录中的 19 个套件全部称为已执行。

解压完整归档后，在归档根目录运行：

```bash
python -B -S verify_all.py --list
python -B -S verify_all.py --quick --offline
python -B -S verify_all.py --all --offline
```

默认接受程序使用 Python 标准库；可选符号样例需要 SymPy。大型证书还原及 JSON 解析需要数百 MB 内存；建议至少 4 GB 可用内存及约 1 GB 临时磁盘。完整运行通常需要数分钟。归档内的索引是 `docs/RESULTS.md`、`docs/OPEN_PROBLEMS.md`、`docs/REPRODUCING.md` 和 `verification/replay_20261027.json`。

完整归档文件名：`mordell-research-complete-2026-10-27.zip`；大小：**103026838 字节**；SHA-256：

```text
00970d35ab15ababf310eb53f2e534ca31a2094b9b1c6f3462983956bdaa1654
```

## 来源与边界

程序检查代数、严格区间及覆盖树；解析、几何与紧性桥梁仍须阅读证明文本。这不是由 Lean/Coq 检查整篇论证的形式化证明。协作代理审阅或不同实现的交叉检查不等于独立第三方同行评审。浮点探索、有限域试算和失败搜索不承担不存在性结论。

未指定许可证，不作全球优先权声明。来源与 AI 协作说明在完整归档的 `docs/PROVENANCE.md`，引用建议在 `CITATION.md`。

**English:** Jiang Fengcan's ongoing research archive on the eight-point Mordell discriminant problem, developed with assistance from OpenAI Codex. Partial theorems, exact certificates and reproducible verification code are preserved. The unconditional sharp inequality remains open; AI-assisted cross-review is not external peer review.

