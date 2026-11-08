# 真实驻点与完整局部极大性的整数区间认证

2026-11-05。蒋丰灿研究项目，OpenAI Codex 协作。

本专题把保存的近似坐标升级为真实驻点存在性及完整二阶结论：

|点数|实际认证|
|---|---|
|10、11、12、13、14、17、18、19|各有一个模整体旋转的严格局部极大轨道|
|15|盒内真实驻点有严格增益二阶方向，因而不稳定|

每个规范盒半径为2⁻⁸⁰，使用整数二进区间的压缩映射证明存在与局部唯一，再以完整实 Hessian 的有理合同和严格 Gershgorin 界证明负性。17、18点采用[本轮改进下界](../finite-search-nov5/README.md)的种子。完整推导见[证明](PROOF.md)。

这些结论没有确定任何 n≥10 的全局最优值。小梯度、正则环族的最优和完整局部极大性属于不同层次；15点的真实不稳定分支给出一个明确区别。

在仓库根目录运行（仅需 Python 标准库）：

```bash
python -B -S -O validated-maxima/verify_all.py
```

普通验证不会覆盖冻结报告。`--record-report` 用于有意更新结果。高精度探索入口 `discover.py` 另需 mpmath；它不是正式接受程序的依赖。

[精确数据](certificates.json) · [正式核验](verify_all.py) · [整数区间实现](dyadic.py) · [实际报告](verification_report.json) · [输入绑定](source_manifest.json)。

另运行 `python -B -S -O validated-maxima/check_integrity.py`，会在临时副本拒绝固定字节篡改，并拒绝零预条件矩阵、伪造 Hessian 合同、偏移中心和零不稳定方向；见[拒绝记录](integrity_report.json)。
