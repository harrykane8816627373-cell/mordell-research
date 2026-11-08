# 复现指南

更新：2026-11-08。

所有命令均在仓库根目录运行，除非另行注明。默认接受只需要 Python 标准库；解析证明与外部论文输入还需阅读正文。

## 候选、必要域与正定核

以下五组正式接受均只需Python标准库，关键判断在`-O`下保留：

```bash
python -B -S -O candidate-global-screen/verify_all.py
python -B -S -O maximizer-domain/verify_all.py
python -B -S -O candidate-basin/verify_all.py
python -B -S -O candidate-kernel/verify_all.py
python -B -S -O twelve-point-kernel-upper/verify_all.py
```

[候选专题](../candidate-global-screen/README.md)核对1000次已保存数值记录的来源与汇总，并重放真实根、精确对称与点对类型证书；它不重新运行搜索，也不把频率变成全局证明。[必要域](../maximizer-domain/README.md)核对已验收的八九点证据绑定，精确重算十点整数下界、Schwarz端点及非线性排除；完整八九点覆盖仍由下述独立重放命令执行。[局部上界](../candidate-basin/README.md)重放真实根与整个扩大盒的Hessian负性，而非只检查候选点处特征值。

[正定核专题](../candidate-kernel/README.md)核验精确矩恒等式及原始数值核的有理失败见证；无余量majorant仍不成立。[十二点全局上界](../twelve-point-kernel-upper/README.md)另行精确接受增加21/10000余量后的连续majorant，重算整数可达值与比值比较、真实候选值下端，以及37053节点的完整树。该入口约30秒，纯标准库；统一默认入口也实际重放完整覆盖，不用已有零未决报告跳过计算。它证明M₁₂<23V₁₂/20，不证明候选锐界或最大构型分类。

数值搜索另外需要NumPy；SDP发现需要NumPy、SciPy、CVXPY及相应求解器。发现依赖不参与正式接受，复现浮点搜索或求解器不能替代连续域证明。各次发现记录及具体运行参数随专题保存。

## 11月8日四环、非线性相位与高点数认证

以下正式接受只需要Python标准库：

```bash
python -B -S -O four-ring-rigidity/verify_all.py
python -B -S -O phase-deficit-sharp/verify_all.py
python -B -S -O family-exclusion-nov8/verify_all.py
python -B -S -O validated-maxima-nov8/verify_all.py
python -B -S -O validated-maxima-nov8/check_integrity.py
```

四组正式检查已接入统一入口。它们固定引用前轮代码、整数种子或角向数据；复制单个新目录而遗漏所列依赖会被拒绝。统一入口在临时副本运行，避免相位与家族程序生成的报告改写冻结文件。新的局部证书仍使用384位整数区间，所有接受判断在`-O`下执行；高精度发现不属于正式接受依赖。

## 11月5日一般点数与局部认证

以下均只需标准库，接受判断在 `-O` 下保留：

```bash
python -B -S -O ring-rigidity/verify_all.py
python -B -S -O phase-packing/verify_all.py
python -B -S -O finite-search-nov5/verify_all.py
python -B -S -O validated-maxima/verify_all.py
python -B -S -O validated-maxima/check_integrity.py
```

统一入口 `scripts/verify_current.py` 也执行上述四组正式核验。相位程序会生成报告，统一入口将它复制到临时目录运行；其余三组默认不改写冻结报告。局部认证拒绝伪造预条件矩阵、Hessian合同、偏移中点和固定数据篡改，记录另见[拒绝检查](../validated-maxima/integrity_report.json)。

[搜索复现](../finite-search-nov5/README.md)另需NumPy，[高精度驻点发现](../validated-maxima/discover.py)另需mpmath；两者只生成候选。正式下界与驻点存在、稳定性接受不依赖数值优化器。实际执行范围及所有作用域见[本轮研究记录](NOV5_RESEARCH.md)。

## 九点与一般点数结构

九点固定上游与完整解析推导见[九点专题](../nine-point-global/README.md)和[证明正文](../nine-point-global/PROOF.md)。首次复核在独立环境安装所列依赖：

~~~powershell
python -m venv .venv-n9
.venv-n9\Scripts\python.exe -m pip install -r nine-point-global/vendor/n9/requirements.txt
.venv-n9\Scripts\python.exe -B nine-point-global/verify_global.py
~~~

该入口默认从完整初始域重新运行，不使用 `--resume` 复用旧零未决结果；`--check-inputs` 仅核对固定字节，不执行证明。只有实际完整运行成功后，才可用 `--record-report` 保存本地报告。本仓库的实际重放与文献原始证书分别记录；不能把下载的零未决报告当成本地重新运行。

[九点实际重放记录](../nine-point-global/reports/replay.json)保存完整初始域的31613个处理盒、15807个终端盒和0未决结果，全部8个上游检查实际执行。记录绑定具体输入与软件版本；默认轻量接受与重新运行完整覆盖应分开理解。

一般点数的连续相位相容性、环对缺额及有限有理见证见[结构研究](../general-structure/README.md)，检查只需 Python 标准库：

```bash
python -B -S -O general-structure/verify_all.py
```

该检查认证环族内的相位结构与明确缺额；不把受限环族中的改善称为全部十二点构型的最优性证明。

十点严格比较只需标准库：

```bash
python -B -S -O scripts/verify_ten_point_transition.py
```

它用精确整数复核 Δ/C₁₀>1，不证明这个见证达到 M₁₀；转换推导见[十点证明](TEN_POINT_TRANSITION.md)。

## 八点完整域证书与桥接

Ji–Xu–Yu固定核证书的首次完整重放使用独立环境：

~~~powershell
python -m venv .venv-n8
.venv-n8\Scripts\python.exe -m pip install -r eight-point-global/vendor/n8/requirements.txt
.venv-n8\Scripts\python.exe -B eight-point-global/verify_global.py
~~~

这会从完整初始盒重新运行六个原始脚本及准确桥接，原始断言保持启用；不读取已保存的零未决结果来跳过覆盖。固定输入及MIT/CC BY署名见[全局目录](../eight-point-global/README.md)。基础统一入口只审核已验收记录、输入绑定并重算桥接；需要统一入口重新运行完整全域时，在上述环境运行 `scripts/verify_current.py --global-eight`。

旧局部树的未决盒作为历史计算记录保留；其 `full_eight_point_exclusion_proved=false` 表示那棵坐标树未展开完成，不撤回新全域核的全局结论。

## 当前幂和与历史八点局部路线

```bash
python -B -S -O scripts/verify_current.py
```

此入口核对当前导航的本地链接，将两组固定输入复制到临时目录，再运行：

- `power-moment/verify_all.py`：四个旧精确入口，包含代数、结构、参考谱及11个目标盒。
- `radial-gap/verify_all.py`：包含已有 `eight-point/verify_all.py` 的两个精确入口、上传统一包的14个标准库入口，以及新的径向—角向常数和实际点见证。新门槛是 E≤1/1200、对数缺口 E/32；上传原包的49条哈希也会核对。可选 SymPy 根配对测试未执行。

仓库文件不因重放而更新计时。这两个历史程序保留 `full_eight_point_exclusion_proved=false`，表示各自局部证书未完成全覆盖；当前全局结论见 `eight-point-global/reports/replay.json`，其完整重放标志为true。

仅重放本轮径向专题及必需输入：

```bash
python -B -S -O radial-gap/verify_all.py
```

[新证明页](../radial-gap/README.md)给出全角度能量比较、一般真实点假设与七点锐界依赖。原统一 ZIP 固定保存，所有脚本在临时副本运行。新包只需标准库，且不依赖32MB首三次图整数对象。

## 联合物理系统与部分覆盖

当前完整续算入口先重放完整父链，再核对上传证书、精确力余项/Newton、相关点差/复余项/二列Gram、共同子集预算、历史119/80半径、二阶Schwarz半径、旧椭圆–Schur半径、19dbbd3前缀防伪、平移子集及力方程总树，共10个接受入口：[11月4日证明与执行范围](../joint-refinement/NOV4_PROGRESS.md)。

~~~bash
python -B -S -O joint-refinement/verify_force.py
~~~

也可在当前基础导航、幂和与径向检查之后一并运行：

~~~bash
python -B -S -O scripts/verify_current.py --refinement
~~~

该命令包含旧联合覆盖，只重放一次旧覆盖树。新专题[证明与报告](../joint-refinement/README.md)明确保留未排空状态；历史离线包仍是当时的冻结快照。

[全局二次证明](../joint-refinement/GLOBAL_QUADRATIC.md)给出150组标量二次极值核对、1152个独立 Cartesian 角点核对及实际非空大盒基准。该阶段节点明确标记算法2或3，较早节点仍按算法1重放。仅核对输入哈希的 `--check-fixed-inputs` 会明确输出 `proof_replay_performed=false`，不能替代完整重放。

算法4继续加入[中心化切片与质量对偶](../joint-refinement/CENTERED_MASS.md)，核对8组切片、1208个可行顶点记录，以及违反交换凸性时的内部极大反例。各版本节点均按记录版本重放。

父阶段算法5加入[精确质量范围与参数包络](../joint-refinement/MASS_DUAL_ENVELOPE.md)，以独立面枚举、255条 Cartesian 顶点记录、264个包络查询及24组穷举参数最小化核对新算术。森林仍逐节点按原算法版本重放。

父阶段[联合证明](../joint-refinement/JOINT_STABILITY.md)给出内点条件 min|qᵢ|²<3/4、严格四阶径向预算与旧椭圆的两个参数包络。新算术独立穷举267220对整数参数，并按系数核对平方和分解、105个 Bernstein 查询与201个质量参数；六份伪造覆盖必须被拒绝。

当前新增 `--nov2-prefix-already-replayed`，仅当19dbbd3此前已经验收时使用：复用94399个旧处理步骤及祖先链，要求70个旧算术输入逐字节不变，仍逐节点重放全部新增步骤。首次复核使用不带跳过参数的默认完整命令。

~~~bash
python -B -S -O joint-refinement/verify_force.py --nov2-prefix-already-replayed
~~~

[执行记录](../joint-refinement/reports/force_execution.json)明确列出本次实际范围。更早的 `--accepted-prefix-already-replayed` 固定8be9f91的60791个处理步骤与61个旧输入；19dbbd3历史阶段曾使用它重放33608个新增步骤。旧选项保留其原语义。

更早的 `--checkpoint-already-replayed` 仍对应768598a的25790步及50个旧算术输入；`--parent-already-replayed` 对应e846b8e父链及36个旧输入。它们不是本次实际使用的前缀模式。

[力方程联合证明](../joint-refinement/FORCE_CLOSURES.md)核对312个复Jacobian条目、24个精确余项、3个已知目标力的Newton像，以及局部圆盘不确定而终点点距可用的条件外包；254个指示向量模式核对完整Hessian必要谱条件，8份伪造覆盖必须被拒绝。上传门槛729/1000按全部多项式系数复核。新增接受检查不使用 `assert`，`-O` 不会跳过它们；原ZIP完整保留。

[新的独立检查](../joint-refinement/EDGE_REMAINDERS.md)另外核对448条点差记录、256个复余项分量、两种必要像各16个已知解，以及84个复Gram主子式；11种伪造检查点必须被拒绝。另有质量7、旧椭圆内的严格探针加强见证，单独运行 `python -B -S -O joint-refinement/verify_mass7_gram.py` 可复算其有理反射、全部127个旧不等式和负主子式。该演示不作为总树叶的依赖。原编号10在19dbbd3阶段已完整关闭；11月4日仅导入旧叶9与55820的两棵部分树，旧14盒没有新增整盒关闭。这条历史路线的原区域仍0–9，保存28个开放后代盒；全部新增步骤已独立重放，但该树未逐盒排空。当前八点全局结论由完整核证书承担。

[新子集检查](../joint-refinement/SUBSET_BUDGET.md)核对1968个准确子集恒等式及重心包围，并验证真实单叶关闭与原父区域关闭的区别。[119/80半径证明](../joint-refinement/TARGET_RADIUS.md)中的四次非负分解、严格有理比较和边界导数恒等式由脚本检查；全体多项式的Schwarz解析桥梁需阅读正文，七点完整驻点分类没有在本轮重建。

[二阶Schwarz核对](../joint-refinement/verify_schwarz2_radius.py)检查231/160常数与0处一次项消失；[旧椭圆–Schur核对](../joint-refinement/verify_ellipse_schur_radius.py)检查23/16的二次项与严格有理间隙；[平移子集核对](../joint-refinement/verify_translated_subset.py)独立核对1968个实际角点恒等式、相关重心包围、5个准确质量7保留盒、连续严格加强、边界等号与作用域。圆盘映射的解析桥梁仍须阅读证明，七点完整驻点分类未在本轮重建。

`reports/force_cover.json` 是最新力方程总账；`reports/joint_cover.json` 和 `reports/cover.json` 分别保留e846b8e、92fbc4d父快照报告。`summarize_frontier.py` 只核对2f25467到92fbc4d的历史边界，不能替代最新增量的算术重放。

基础入口默认不重复耗时较长的60021节点覆盖重放。加上新专题：

~~~bash
python -B -S -O scripts/verify_current.py --joint
~~~

单独重放新专题及其必需旧输入：

~~~bash
python -B -S -O joint-closure/verify_all.py
~~~

接受入口核对等价多项式模型、真实驻点盒和覆盖树；所有21个开放区域都显式保留。接受通过表示这些有限证书正确重放，完整排空状态仍为 false。[解析等价证明与续算命令](../joint-closure/README.md)

## 保留的全局低阶矩耦合

这项检查需要旧首三次图的32MB整数证书，默认轻量入口不会隐式下载或跳过它。用原2026-10-16证明包重放：

```bash
python -B -S -O algebraic-gap/verify_all.py --source-archive PATH_TO_FIRST3_ZIP
```

也可一次执行全部当前入口：

```bash
python -B -S -O scripts/verify_current.py --refinement --first3-archive PATH_TO_FIRST3_ZIP
```

[本轮证明页](../algebraic-gap/README.md)列出准确源包名和哈希；完整新离线包只携带3个必要旧文件，解压后无需额外参数。该入口重新构造结果式、逐项重放14179个整数系数、验证齐次权和七个乘子范数，再用精确分数核算全局楔形下界。它没有重跑旧零矩构型的全部分类；这个局部方法的完整八点排除标志仍为 false。

## 保留的两阶段局部能量证书

只验证最新专题可运行：

```bash
python -B -S -O eight-point/verify_all.py
```

[最新离线包](../eight-point/mordell-eight-point-energy-2026-11-01.zip)解压后也使用这条命令。固定输入由 `eight-point/source_manifest.json` 约束；包内报告应当与重新计算结果完全相同。可用 `eight-point/build_package.py` 从已固定输入确定性重建最新能量 ZIP；首阶段 ZIP 保持原字节。

## 既有一般 n 的累计检查

```bash
python -B -S -O annular-layer/verify_all.py
```

该入口累计重放28个所列程序，覆盖层数、能量及百点条件证书。它不会执行所有历史巨型网格、重建所有六七点外部分类，或完成八点全域消元。

较早构型专题的标准库检查：

```bash
python -B -S -O configuration-bounds/verify_all.py
```

其他专题以[研究索引](RESEARCH_INDEX.md)所链接的各自首页为准。部分旧入口会重写本专题的计时报告；需保留干净工作区时，在副本或对应离线包中运行。

## 历史包与当前网页

历史 ZIP 是冻结快照；当前 README、证明状态页可能随后更新。`configuration-bounds/status_updates.zip` 中的状态快照由自己的清单检查，不要求与当前网页逐字相同。当前 `PROOF_STATUS.md` 在本地专题清单里的字节数与哈希随网页状态更新；历史 ZIP 与其快照清单保持原字节。

哈希保证字节一致；有理程序保证所列运算与有限检查；数学桥梁与全覆盖量词必须由证明及实际穷尽证据承担。[旧八点方法的开放域记录](../eight-point/OPEN_BRANCHES.md)作为历史保留；当前剩余系统由新的全域核接合排除。
