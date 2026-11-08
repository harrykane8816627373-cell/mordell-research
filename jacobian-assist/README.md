# Jacobian计算辅助环境与真实调用记录

2026-11-08。按蒋丰灿的要求接入[morluto/jacobian](https://github.com/morluto/jacobian)，固定官方发行版0.23.0。已在独立Windows Python环境安装，并注册Codex本地stdio MCP服务器`jacobian`；握手、`math.find`检查契约和`math.run`实际矩阵计算均已执行通过，当前Codex会话也已直接调用这两个工具。

官方wheel的SHA256为`72fe5a9210db17e18e756a25cb94473b2c188bba591ad8ab8d7dfff9531e05f8`。环境与版本见[environment.json](environment.json)、[安装锁](requirements-windows.lock)。配置方式遵循[官方Codex MCP文档](https://learn.chatgpt.com/docs/extend/mcp?surface=cli)。没有更改已有其他MCP服务器配置。

## 已实际用于本研究的运算

- 原核和[新十二点核](../twelve-point-kernel-refinement/README.md)的全部12个Gram块，各自调用Jacobian的`matrix.inertia.compute`，均得到全部特征值为正；与独立标准库有理LDL分解和逐项重构一致。
- 对真实径向特征10q(u)²，调用`polynomial.bernstein.coefficients.compute`，FLINT返回的7个精确系数与标准库公式逐项相同。
- 在先前数值搜索遗漏的低谷小盒上，调用`interval.expression.second_jet_enclosure.compute`，Arb严格包围log距离的值、梯度及Hessian，返回`ENCLOSED`；解析中心导数也在区间内。
- 当前会话直接调用MCP，对新核的实际5×5矩块计算精确惯性，并以`polynomial.sos.gram.check`接受10q²的真实秩一Gram平方和证书。

完整输入输出保存在[旧核原生记录](native_report.json)、[新核原生记录](refined_native_report.json)、[stdio调用记录](mcp_report.json)和[当前会话MCP记录](injected_mcp_report.json)。这些是实际计算结果，不是预期输出。

Jacobian辅助构造和交叉核验。新的M₁₂<(103/100)V₁₂全局界仍由[独立标准库接受器](../twelve-point-kernel-refinement/verify_all.py)重放完整87,385节点证明。数值核由CVXPY/Clarabel及广域反例反馈发现；没有将Jacobian调用成功、有限采样或浮点求解状态代替全域证明。

## Windows后端范围

这是明确限定的原生Windows配置。官方发行版的PPLpy依赖没有对应Windows轮子，因此保留官方Jacobian wheel不变，以`--no-deps`安装它，再精确安装其余可用依赖。SymPy、Python-FLINT、Z3和CyPari已安装；PPL、Singular、QEPCAD未启用。依赖元数据仍如实保留对PPLpy的要求，`pip check`的这一项缺失被记录，未伪装成完整依赖栈。

已检查PPL导入仅在相应操作的执行路径中发生，不妨碍上述MCP、矩阵、多项式及区间操作。PPL相关优化、Singular和QEPCAD操作不属于本次可用范围。上游明确测试的完整二进制安装平台是glibc Linux x86-64；见[上游安装说明](https://github.com/morluto/jacobian#quickstart)。

## 复现与核验

在64位Windows CPython3.12或3.13下，从仓库根目录执行：

```powershell
python jacobian-assist/bootstrap.py --dry-run
python jacobian-assist/bootstrap.py --register-mcp
```

[bootstrap.py](bootstrap.py)固定版本、校验官方wheel、创建用户目录下的独立环境。只有显式`--register-mcp`才写Codex配置；已存在的不同Jacobian配置会被保留并拒绝覆盖，首次写入前备份本地配置。脚本使用系统代理设置，保留TLS证书验证，不输出代理值；私人配置备份不进入仓库。

真实后端复算：

```powershell
& "$env:USERPROFILE/.codex/tools/jacobian/0.23.0/venv/Scripts/python.exe" -B jacobian-assist/bridge.py
& "$env:USERPROFILE/.codex/tools/jacobian/0.23.0/venv/Scripts/python.exe" -B jacobian-assist/bridge.py --data twelve-point-kernel-refinement/kernel.json
& "$env:USERPROFILE/.codex/tools/jacobian/0.23.0/venv/Scripts/python.exe" -B jacobian-assist/mcp_probe.py
```

默认不改写记录，显式`--record`才保存。轻量固定记录接受器只需标准库：

```bash
python -B -S -O jacobian-assist/verify_all.py
```

该命令核对输入、记录和独立有理计算，明确报告本次没有重启Jacobian。需要同时重算真实后端时，加`--jacobian-python`及其解释器路径；统一入口提供对应`--jacobian`选项。固定输入见[source_manifest.json](source_manifest.json)，实际接受结果见[verification_report.json](verification_report.json)。

## 继续辅助搜索与构造

[切平面搜索](search/cut_search.py)、[浮点覆盖发现](search/cover_candidate.py)及其[多项式工具](search/bernstein_discovery.py)已保存；原始遗漏区域的数值种子也一起保留。它们固定原特征，优化Gram块，以广域低深度分区和彼此分离的负谷反馈新约束。输出写入仓库的`work/jacobian-search/`，不覆盖已接受证书。

发现另需NumPy、SciPy、CVXPY和Clarabel，版本见[发现环境](../candidate-kernel/discovery_environment.json)，与精确Jacobian环境分开。示例命令：

```bash
python jacobian-assist/search/cut_search.py --iterations 11
python jacobian-assist/search/cover_candidate.py --candidate twelve-point-kernel-refinement/kernel.json --c 0.0004 --seconds 60
```

这些命令只作数值发现。过去11轮中先后调整过oracle深度和多谷策略，求解器也可能给`optimal_inaccurate`；当前发布脚本保留最后的广域策略，不承诺逐位复现历史优化路径。正式结果须把输出解释为精确系数后，重新执行完整整数覆盖与可达值比较。

完整h有602项，部分Jacobian操作有256项等准入限制，因此采用谐波块、径向特征和小盒组合。下一步可用分块精确梯度、严格区间二阶包围及精确PSD面约化来改进搜索；这些操作的返回范围仍须逐项检查。
