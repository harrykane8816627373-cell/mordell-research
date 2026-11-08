# 11月8日扩展：25、35、42、46点的真实局部极大轨道

蒋丰灿（Fengcan Jiang）研究项目，OpenAI Codex 协作。2026-11-08。

已将[11月5日改进下界](../finite-search-nov5/README.md)的四组整数种子升级为**真实驻点存在性、规范盒内唯一性，以及完整 Hessian 的严格局部极大认证**。四个分支都通过了整数区间检查；这不确定任何点数的全局最优值。

|点数 n|规范维数|已认证结论|相对旧保存整数构型的目标值比值包围|
|---|---|---|---|
|25|49|模整体旋转严格局部极大|[1.023141696685, 1.023141696686]|
|35|69|模整体旋转严格局部极大|[1.036339237884, 1.036339237885]|
|42|83|模整体旋转严格局部极大|[1.045381742429, 1.045381742430]|
|46|91|模整体旋转严格局部极大|[1.049316805422, 1.049316805423]|

四个旋转规范盒的半径均为2⁻⁸⁰。存在性由完整力方程的压缩映射证明；局部极大性验证所有点独立移动时的完整实 Hessian，超出指定对称群或多环族。根的 Δ 下端也都严格超过各自11月5日整数种子的准确值。详细推导见[证明](PROOF.md)，完整有理端点见[实际报告](verification_report.json)。

本专题直接导入固定哈希的[既有正式验证程序](../validated-maxima/verify_all.py)和[整数区间实现](../validated-maxima/dyadic.py)，不修改或替换原程序、原证书或原报告。新中点、预条件矩阵、Hessian 合同矩阵与完整整数种子独立保存在[新证书](certificates.json)中。

在仓库根目录复核，只需 Python 标准库：

```bash
python -B -S -O validated-maxima-nov8/verify_all.py
python -B -S -O validated-maxima-nov8/check_integrity.py
```

普通复核重新计算所有连续盒的区间界，并与冻结结果逐项比较，不覆盖报告。`--record-report` 只用于明确更新记录；它仍执行数学接受检查。反例测试拒绝零压缩预条件、奇异 Hessian 合同、错误中点、伪造正方向、篡改发现种子和固定字节。

可选的高精度生成入口另需 mpmath 1.3.0；本次以180位十进制精度生成证书。示例：

```bash
python -B validated-maxima-nov8/discover.py --n 46 --output work/nov8-cert/cert-n46.json
```

依次对25、35、42、46生成候选数据后，才由正式整数验证接受。mpmath 的收敛或 Cholesky 成功只是发现阶段信息，不被当成证明。重放最终证书不需要安装或导入 mpmath、NumPy。

完整 Lagrangian Hessian 的惯性为(2n−1个负方向、1个旋转零方向、0个正方向)；质量球面切空间上为(2n−2,1,0)。表中“唯一”始终只限于明确旋转规范盒，既不是全部驻点的唯一性，也不是全局最优轨道的唯一性。

[正式入口](verify_all.py) · [反例拒绝报告](integrity_report.json) · [源与依赖绑定](source_manifest.json) · [原通用证明](../validated-maxima/PROOF.md)。
