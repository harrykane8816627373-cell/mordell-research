# Mordell 距离乘积研究：有限点数构型与一般 n

研究资料与项目发起：**蒋丰灿（Fengcan Jiang）**（[harrykane8816627373-cell](https://github.com/harrykane8816627373-cell)）。推导、精确证书与代码包含 OpenAI Codex 协作。当前整理：**2026-10-28**。

统一归一化为

```math
M_n=\max_{\sum_{i=1}^n|z_i|^2=n}\Delta_n(z),\qquad
\Delta_n(z)=\prod_{i\lt j}|z_i-z_j|^2,\qquad N=\frac{n(n-1)}2.
```

**八点完整无条件锐界、一般非零矩分支，以及一般 n 的简洁全局锐值公式仍未解决。** 本仓库保留已证定理、可达下界、条件构型类上界及可重放证书，具体假设见各证明。

## 从这里阅读

- **[当前证明状态总表](configuration-bounds/PROOF_STATUS.md)**：统一查询结论、适用范围、证明和代码；历史待补条目的当前判断也在此汇总。
- **[最新有限 n 构型与上下界](configuration-bounds/README.md)**：全部正负相位的精确分解、带严格增益证书的有理上升、径向人口限制与构型类上界。
- **[五、六、七点与较早一般 n 证明](earlier-results/README.md)**：作者原稿、全局证明依赖、PDF、TeX及精确复核；[一般 n 审核](earlier-results/general_n_proof.md)已直接更新过期状态。
- **[完整原稿恢复与联合方法](source-review/README.md)**：原始/修订主稿、来源核对、完整旧包络、半径层数、删除余量和混合网格。
- **[全部有限 n 的下界与严格子集余量](finite-n/README.md)** · **[跨点数递推与有效全局网格](global-n/README.md)** · **[线性渐近项与有限圆盘上界](asymptotics/README.md)**。

## 目前证明到了哪里

|问题或结论|当前状态|证明与限制|
|---|---|---|
|五点全局锐值|$`M_5=5^5=3125`$，单位正五边形|[作者二阶条件与矩方法](earlier-results/five_point_method.md)、[PDF](earlier-results/five_point_proof.pdf)|
|六、七点全局锐值|$`M_6=6^{15}/5^{10}`$；$`M_7=7^{21}/6^{15}`$|[主稿](earlier-results/six_seven_proof.tex)、[覆盖审计](earlier-results/six_seven_coverage_audit.md)；全局证明引用 Faugère–Svartz 物理驻点完整分类，七点356次枝未在本地独立重建|
|八点完整锐界|目标为 $`M_8=8^{28}/7^{21}`$；尚未证明|原点加正七边形达到候选值；低矩和若干几何分支已证，一般 $`p_2p_3\ne0`$ 分支仍开放，$`p_k=\sum_i z_i^k`$|
|一般 n 的有限上下界|已证，且可用指定构型继续提高下界|[有限下界](finite-n/README.md)、[圆盘上界](asymptotics/README.md)、[本轮精确改善](configuration-bounds/README.md)；有限家族最优不等于全局最优|
|严格真子集约束及有效网格|已证；有显式正有理余量与有限全局误差界|[子集证明](finite-n/SUBSET_PROOF.md)、[混合网格](source-review/HYBRID_GRID.md)；网格上界要求完整穷尽规定网格|
|十点中心正九边形及旧缺件条目|四阶严格局部最大、完整网格、保存整数下界和特定 Hermite 障碍现均已核验|[已更新原审核](earlier-results/general_n_proof.md)；十点该局部最大仍全局次优，旧证明确认归用户原有成果|
|确定的线性渐近项|已证 $`\log M_n=cn^2+\frac n2\log n+\gamma n+o(n)`$，$`c=\frac12\log2-\frac14`$|[外部定理、尺度与常数](asymptotics/README.md)；精确 $`\gamma`$ 及真实有限 n 边界修正仍未知|

## 本轮构型与上下界的联合进展

[最新专题](configuration-bounds/README.md)把完整旧多环、格点、矩平均和保存整数构型合并为基线，再加入两类精确改善。三十点的人口6、12、12及平方半径1/5、4/5、8/5，在正负相位家族中给出相对旧两方案的增益

```math
\left(\frac{90155}{84693}\right)^{12}\gt2.
```

全部相位家族在2≤n≤100已核对，相对完整旧包络提高13项；有理一步上升在n≤30又提高相位之外的13项，n=26至少提高98.30%。小残差只控制驻点条件，不能认证全局最优。

无条件有限圆盘上界保持为

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

每个专题首页均提供证明、Python标准库接受入口、字节清单和离线包。最新[构型与上下界包](configuration-bounds/mordell-configuration-bounds-2026-10-28.zip)在解压根目录运行：

```bash
python -B -S -O verify_all.py
```

其他专题离线包：[完整原稿与联合方法](source-review/mordell-source-review-2026-10-28.zip)、[有限 n 与子集余量](finite-n/mordell-finite-n-2026-10-27.zip)、[全局网格更正包](global-n/mordell-global-n-corrected-2026-10-28.zip)、[五六七点历史专题](earlier-results/mordell-earlier-proofs-2026-10-27.zip)、[渐近专题 Release](https://github.com/harrykane8816627373-cell/mordell-research/releases/tag/general-n-linear-2026-10-27)。历史ZIP保持原字节；当前网页状态及清单另行同步。

完整八点附件为 **mordell-research-complete-2026-10-27.zip**，大小103026838字节，SHA-256：

```text
00970d35ab15ababf310eb53f2e534ca31a2094b9b1c6f3462983956bdaa1654
```

该附件含19个历史及新增证明包、1516个文件和3个必需大型精确对象。2026-10-27实际重放14个套件、24个程序入口全部PASS，另5个由更强包覆盖而跳过。GitHub自动生成的Source code.zip/tar.gz只包含Git目录，不能替代该完整附件。接受程序及覆盖范围见附件内 `docs/REPRODUCING.md`；来源见附件内 `docs/PROVENANCE.md`，引用建议见附件内 `CITATION.md`。

程序验证代数、严格区间及有限覆盖，解析与几何桥梁仍须阅读证明。AI协作交叉审阅不等于外部同行评审或整篇形式化证明；浮点探索和失败搜索不承担不存在性结论。未指定许可证，不作全球优先权声明。

**English:** Jiang Fengcan's ongoing research on finite-point Mordell discriminant maxima. Proved partial theorems, constructive lower bounds, conditional upper bounds and reproducible exact certificates are indexed separately. The complete eight-point sharp bound and a simple general-n sharp formula remain open.
