# 八点 Mordell 判别式极值问题

研究资料与项目发起：**蒋丰灿**（[harrykane8816627373-cell](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。

**完整八点无条件锐界，以及一般 $p_2p_3\ne0$ 分支，仍未解决。** 本项目记录进行中的部分定理及可复核证明，不宣称主问题已经完成。

在归一化 $\sum_{i=1}^8|z_i|^2=8$ 下，目标是

$$\Delta_8=\prod_{i<j}|z_i-z_j|^2\stackrel{?}{\le}\frac{8^{28}}{7^{21}}.$$

目标等号配置为原点加半径平方 $8/7$ 的正七边形；$p_k=\sum_i z_i^k$。

## 2026-10-27 的研究进展

- **六共圆与等距例外：** 八点真实驻点若六点共任意圆，另两点到该圆圆心等距，则只能是正八边形，因而不稳定。无需预设反射；证明含 69 节点有理多项式组合证书和 624 项标准库检查。两个例外到圆心不等距的主支仍开放。
- **连续核的跨域连接：** 以共享的 12 维实际特征，将既有局部核连接到半径平方 $1/2,4/7,9/10$ 附近的三个远端带，连续角度与径向带宽由 4285 项严格区间检查覆盖；不需要点距分离假设。过渡半径及其他远端区域尚未全部覆盖。
- **既有几何与低矩结果：** 保留 $p_2=0$、$p_3=0$ 分支的真实驻点分类及相应值界，以及两直线、抛物线、四轴反射等已证类别。七点共任意圆的八点稳定驻点只可能是目标构型。各结论的驻点、稳定性或全局极大假设分别见完整证明索引。
- **一般分支的结构归约：** 保留真实点、共轭配对及完整稳定条件的矩阵与多项式约束已经建立；一般系统尚未排空。若干辅助核证书族的不可能性是对证明方法的限制，不能当成原极值问题的结论。

## 完整资料与复核

完整材料见 [research-2026-10-27 发布版本](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/research-2026-10-27)，请下载其中的 **mordell-research-complete-2026-10-27.zip** 附件，并用同页 `SHA256SUMS.txt` 核验。

网页仓库暂用作公开首页与下载索引；GitHub 自动生成的 Source code.zip / tar.gz 只对应首页提交，**不是完整研究资料**。仅克隆首页仓库也不能替代下载完整研究附件。

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
