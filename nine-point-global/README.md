# 九点全局锐界：固定外部证书的独立重放

2026-11-04。定理、核及原始计算机辅助证书归于 **Xiang Ji、Tiancheng Xu、Wensheng Yu**，见[论文 arXiv:2611.03088v1](https://arxiv.org/abs/2611.03088v1)和[固定代码版本](https://github.com/Jxxxx-2001/mordell-discriminant/tree/31cc934b7dfa0ef2308566a1bce48cbfb2e01559)。本仓库保存署名、授权、精确输入及独立执行记录，并整理[必要域、求和及全部等号的证明](PROOF.md)。

```math
M_9=9^{36}/8^{28}.
```

等号恰为原点加半径 $\sqrt{9/8}$ 的正八边形，允许旋转与重排。证明不使用九点驻点完备分类或预设对称性。先用完整二阶必要条件导出径向矩和共同点对域，再认证三变量核不等式，最后从非负矩二次型的等号恢复构型。

`vendor/` 逐字节保存30个固定上游文件；`UPSTREAM.json` 记录路径、Git blob SHA、SHA256与长度。代码和机器数据采用上游MIT许可，原始说明采用CC BY 4.0，见[授权说明](vendor/LICENSING.md)。上游旧说明中的“尚未证明”描述其早期阶段，不能覆盖最终证明；上游保存的成功报告也不能代替本地重放。

独立环境与完整重放：

```powershell
python -m venv .venv-n9
.venv-n9\Scripts\python.exe -m pip install -r nine-point-global/vendor/n9/requirements.txt
.venv-n9\Scripts\python.exe -B nine-point-global/verify_global.py
```

该入口在临时副本执行八个上游入口，删除旧定义、局部及覆盖输出，**从完整初始域开始**；不使用 `--resume`。父入口可在 `-S -O` 下运行，但子脚本强制启用断言且移除 `PYTHONOPTIMIZE`。`--seconds` 调整覆盖时间预算；预算耗尽而前沿非空时必须失败。`--check-inputs` 只检查字节，不是证明重放，也不能生成成功报告。

日常复现使用上述命令，保留仓库中已冻结的执行证据。只有准备替换冻结记录时才加入 `--record-report`；它在全部成功后覆盖 `reports/` 中本次生成的输出，随后仍须重新审计并更新 `accepted_replay.json` 对报告的哈希，才能把这次执行登记为新的已接受证据。

本地实际执行结果及版本见[重放报告](reports/replay.json)与[原始输出](reports/fresh_replay.log)。报告绑定本地证明、驱动及全部输入；历史数据与本次执行分开。主要软件前提是SymPy 1.14.0、python-flint 0.9.0、NumPy及IEEE-754向外舍入；没有外部同行评审或整篇形式化证明的声明。
