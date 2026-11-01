# 复现指南

所有命令均在仓库根目录运行，除非另行注明。默认接受只需要 Python 标准库；解析证明与外部论文输入还需阅读正文。

## 当前幂和与八点

```bash
python -B -S -O scripts/verify_current.py
```

此入口核对当前导航的本地链接，将两组固定输入复制到临时目录，再运行：

- `power-moment/verify_all.py`：四个旧精确入口，包含代数、结构、参考谱及11个目标盒。
- `eight-point/verify_all.py`：两个精确入口：保留首阶段管域，另验证四倍欧氏半径、五倍矩管域、距离 Laplacian 谱、36个 Parseval 内积、两类新见证及共六组错误参数拒绝。

仓库文件不因重放而更新计时。两个程序都明确保留 `full_eight_point_exclusion_proved=false`，这表示完整八点覆盖尚未证明。

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
