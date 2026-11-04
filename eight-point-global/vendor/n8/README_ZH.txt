八点 Mordell 全局最优性：复核包

先阅读 PROOF_ZH.txt。结论为 M8=8^28/7^21；全部等号是原点加
半径 sqrt(8/7) 的正七边形，只允许整体旋转和重标号。
这是通常数学论证加严格区间证书的计算机辅助证明。

文件
PROOF_ZH.txt：完整数学推导、域覆盖方法、等号条件、依赖和来源限制。
data/global_pair_bound_exact_contacts.json：唯一核心输入，含精确矩阵与多项式系数。
verify_all.py：从头复核全部必要计算，不依赖此前抽样和优化结果。
data/FINAL_CERTIFICATE.json：最近一次完整复核的总证书。
其它 data 文件：各步骤的输出；单个步骤中的 global_majorant_proved=false
仅表示该步骤本身不是全局证明，请以完整总证书和证明正文为准。

复核方法
安装 Python 3 以及 numpy、sympy、python-flint，然后在本目录运行：
  python verify_all.py
不要使用 python -O；程序必须启用断言。
正常耗时数分钟，视机器而异。结束时应出现
  ALL REQUIRED CERTIFICATES PASSED; UNRESOLVED BOXES: 0
也可单独阅读每个验证程序；没有联网调用、随机采样或 SDP 优化依赖。

系数文件 SHA-256：
405ebbccbe70b6da7c69f13e1fdbb130625f5aba98898729cb2cdf16153995d9
