# 复现指南

所有命令均在仓库根目录运行，除非另行注明。默认接受只需要 Python 标准库；解析证明与外部论文输入还需阅读正文。

## 当前幂和与八点

```bash
python -B -S -O scripts/verify_current.py
```

此入口核对当前导航的本地链接，将两组固定输入复制到临时目录，再运行：

- `power-moment/verify_all.py`：四个旧精确入口，包含代数、结构、参考谱及11个目标盒。
- `radial-gap/verify_all.py`：包含已有 `eight-point/verify_all.py` 的两个精确入口、上传统一包的14个标准库入口，以及新的径向—角向常数和实际点见证。新门槛是 E≤1/1200、对数缺口 E/32；上传原包的49条哈希也会核对。可选 SymPy 根配对测试未执行。

仓库文件不因重放而更新计时。两个程序都明确保留 `full_eight_point_exclusion_proved=false`，这表示完整八点覆盖尚未证明。

仅重放本轮径向专题及必需输入：

```bash
python -B -S -O radial-gap/verify_all.py
```

[新证明页](../radial-gap/README.md)给出全角度能量比较、一般真实点假设与七点锐界依赖。原统一 ZIP 固定保存，所有脚本在临时副本运行。新包只需标准库，且不依赖32MB首三次图整数对象。

## 联合物理系统与部分覆盖

当前完整续算入口先重放完整父链，再独立重放上传证书、精确力余项/Newton算术、相关点差/复余项/二列Gram检查及全部力方程总树，共4个新增入口：

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

固定e846b8e父链已验收且36个旧算术输入未变时，可用 `--parent-already-replayed`。若768598a的25790个处理步骤及其父链也已验收，可用 `--checkpoint-already-replayed`；程序固定原证书、报告、元数据，检查50个旧算术输入未变，要求旧处理步骤全部原样保留，再逐节点重放所有后代。本次采用后者，在固定输入临时副本实际重放35001个新增处理步骤。[执行记录](../joint-refinement/reports/force_execution.json)明确报告复用范围；首次复核使用不带跳过参数的默认完整命令。

[力方程联合证明](../joint-refinement/FORCE_CLOSURES.md)核对312个复Jacobian条目、24个精确余项、3个已知目标力的Newton像，以及局部圆盘不确定而终点点距可用的条件外包；254个指示向量模式核对完整Hessian必要谱条件，8份伪造覆盖必须被拒绝。上传门槛729/1000按全部多项式系数复核。新增接受检查不使用 `assert`，`-O` 不会跳过它们；原ZIP完整保留。

[新的独立检查](../joint-refinement/EDGE_REMAINDERS.md)另外核对448条点差记录、256个复余项分量、两种必要像各16个已知解，以及84个复Gram主子式；11种伪造检查点必须被拒绝。另有质量7、旧椭圆内的严格探针加强见证，单独运行 `python -B -S -O joint-refinement/verify_mass7_gram.py` 可复算其有理反射、全部127个旧不等式和负主子式。该演示不作为总树叶的依赖。当前总账有11个未闭合原域、30个开放后代盒，完整排空仍未完成。

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

[本轮证明页](../algebraic-gap/README.md)列出准确源包名和哈希；完整新离线包只携带3个必要旧文件，解压后无需额外参数。该入口重新构造结果式、逐项重放14179个整数系数、验证齐次权和七个乘子范数，再用精确分数核算全局楔形下界。它没有重跑旧零矩构型的全部分类；完整八点排除状态仍为 false。

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

历史 ZIP 是冻结快照；当前 README、证明状态页可能随后更新。`configuration-bounds/status_updates.zip` 中的状态快照由自己的清单检查，不要求与当前网页逐字相同。本次仅同步当前 `PROOF_STATUS.md` 在本地专题清单里的哈希，不追溯改写历史 ZIP。

哈希保证字节一致；有理程序保证所列运算与有限检查；数学桥梁与全覆盖量词必须由证明及实际穷尽证据承担。[八点未关闭部分](../eight-point/OPEN_BRANCHES.md)继续开放。
