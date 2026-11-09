# 仓库整理与状态约定

更新于 **2026-11-09**。 [证明状态](../configuration-bounds/PROOF_STATUS.md) · [研究索引](RESEARCH_INDEX.md)

## 十二点论文归档

英文论文及 LaTeX 源文件放在十二点锐界专题的 [论文子目录](../twelve-point-sharp-global/manuscript/README.md)。论文引用固定的证书提交，阅读与复现入口列在仓库首页。当前结果集中在[证明状态](../configuration-bounds/PROOF_STATUS.md)，较早的上界和局部证书按阶段收录在[研究索引](RESEARCH_INDEX.md)。

GitHub About 标注十二点全局锐界、八九点证书复现及一般构型研究。历史 Releases 对应各自快照，证明范围以当时的材料为准。 [历史发布](https://github.com/harrykane8816627373-cell/mordell-research/releases)

## 2026-11-04 整理记录

### 远程核对

整理基于 `main` 的 `5ad6b78d38a25cafe42a26541b19fabd4879194b`。核对时远程只有 `main`，没有其他分支、PR 或 Issue，因此没有可据此合并或删除的远程对象。新的证明与研究通过一次可审阅的更新接入。

### 本次整理

- 首页改为点数、定理和复现入口，撤下“八点全局锐界仍开放”的当前状态表述。整理前完整叙述保存在[首页历史快照](HISTORY_2026-11-04.md)，原字节仍可按原提交查阅。
- 八点全域核的已验收完整覆盖，与历史坐标树的未决盒分开说明。新的全局定理不改写旧树的执行范围。
- 九点文献定理、上游授权代码、解析必要域、实际完整重放分别记录于[九点专题](../nine-point-global/README.md)。本项目自己的成果与外部作者的定理保留准确归属。
- 一般点数的连续相位相容性、环对缺额及可复算研究见[结构专题](../general-structure/README.md)。有限族中的最优不自动升级为全部构型的最优。

### 固定材料与当前网页

历史 ZIP、原始主稿、精确输入、冻结报告及其清单保持原字节。`configuration-bounds/source_manifest.json` 中只有当前网页 `PROOF_STATUS.md` 的绑定随状态更新，历史 `status_updates.zip` 与 `status_update_manifest.json` 不重写。

一些受固定输入清单约束的旧证明页，仍在结尾使用“全八点锐界开放”等当时措辞。应按其日期和局部方法理解；当前结论统一由首页和状态表覆盖。这些记录包括：

|历史主题|当时的范围|
|---|---|
|[八点局部能量](../eight-point/ENERGY_BASIN.md)、[矩管域](../eight-point/MOMENT_TUBE.md)|只证明候选附近的目标值排除|
|[局部力方程续算](../joint-refinement/NOV4_PROGRESS.md)、[力方程接合](../joint-refinement/FORCE_CLOSURES.md)|保存未完成的坐标覆盖森林|
|[有限构型结果](../configuration-bounds/FINITE_N_RESULTS.md)、[原稿复核](../source-review/REVIEW.md)|阶段性下界、来源核对与当时的开放问题|
|[径向层数](../radial-angular/RADIAL_LAYERS.md)、[局部能量](../collective-stability/LOCAL_ENERGY.md)、[角向上界](../joint-bounds/ANGULAR_PROOF.md)|一般 n 必要条件或上界，各自不提供八点全域证明|
|[有效层数](../population-transport/EFFECTIVE_LAYERS.md)、[径向运输](../population-transport/RADIAL_TRANSPORT.md)、[微观覆盖](../distribution-geometry/MICROSCOPIC_COVERING.md)|结构定理与尺度条件，其数学结论继续有效|

不应为了消除历史措辞而改写已冻结证书的哈希。后续若修订某项数学证明，需新建版本、更新绑定并重新运行相应验证。

### 整理后的核验

本次实际运行 `python -B -S -O configuration-bounds/verify_all.py` 通过：34 份当前输入绑定、6 份历史状态快照、6 个构型专题入口，以及 n=2…30 的 29 项完整旧下界基线。程序生成的计时报告在检查后恢复历史原字节，避免把既有报告的日期与本次运行混在一起；当前网页绑定随后单独复核。九点和结构专题接入后，当前中央导航检查通过：8页、249个本地链接全部存在且位于仓库内。

### 后续研究的接受标准

n≥10 的改进构型可用整数或有理坐标直接认证下界；全局锐界需要覆盖所有构型的解析不等式或完整穷尽。规则环族的相位最优、力残差为零或 Hessian 必要条件只承担各自作用，不替代全局证明。精确一般 n 锐值、所有最优构型、精确渐近常数与微观晶格结晶仍是开放研究目标。
