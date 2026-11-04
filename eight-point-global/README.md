# 八点全局锐界：外部证书复现与既有结果衔接

2026-11-04。本目录复现 Xiang Ji、Tiancheng Xu、Wensheng Yu 的[预印本 arXiv:2611.03088v1](https://arxiv.org/abs/2611.03088v1)及[固定代码版本](https://github.com/Jxxxx-2001/mordell-discriminant/tree/31cc934b7dfa0ef2308566a1bce48cbfb2e01559)。三变量全局核及其证书是上述作者的成果；本仓库整理复核，并与蒋丰灿项目已有严格子集定理接合。

[桥接证明](PROOF_BRIDGE.md)说明：已有半径与点距条件足以把全部全局最大点的点对放进核证书的域。核求和及全部等号给出

~~~math
M_8=8^{28}/7^{21},\qquad \Delta_q\le7^7\quad(\sum|q_i|^2=7).
~~~

等号恰为原点加正七边形，允许旋转、重排。因此质量7、Δ≥7⁷且s₂s₃≠0的剩余系统为空。该桥接无需七点完整驻点分类。

`vendor/`逐字节保留上游输入、六个脚本、原始说明与授权。代码及机器数据采用上游MIT许可；原始说明采用CC BY 4.0，详见[授权说明](vendor/LICENSING.md)。[来源清单](UPSTREAM.json)记录全部上游路径和哈希。旧局部搜索树继续保留其原始开放盒计数，不能把全局结论误记为它已逐盒完成。

首次复核需独立环境及固定依赖：

~~~powershell
python -m venv .venv-n8
.venv-n8\Scripts\python.exe -m pip install -r eight-point-global/vendor/n8/requirements.txt
.venv-n8\Scripts\python.exe -B eight-point-global/verify_global.py
~~~

该入口从临时副本重新执行全部原始检查、完整初始域覆盖和桥接代数，原始脚本的断言保持启用。`--check-inputs`只核对字节，不是证明重放。使用`--record-report`仅在全部成功后保存[实际执行记录](reports/replay.json)和原始输出。

计算机辅助证明保留明确的软件与解析推导前提；没有声称外部同行评审或整篇形式化证明。
