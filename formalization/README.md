# Mordell 问题的 Lean 形式化

lean-toolchain 固定 Lean 4.28.0，lake-manifest.json 锁定 Mathlib 及传递依赖。从本目录构建并检查主要定理的公理依赖。

~~~sh
lake exe cache get
lake build
lake env lean AxiomAudit.lean
~~~

Windows 上若 Mathlib 的 ProofWidgets 缓存解包遇到 tar 的非 UTF-8 输出错误，可在当前 PowerShell 会话把 Git for Windows 自带的 GNU tar 放到 PATH 前面，再运行缓存命令。

~~~powershell
$env:Path = "C:\Program Files\Git\usr\bin;$env:Path"
$env:LC_ALL = 'C'
lake exe cache get
~~~

## 已在 Lean 中检查

|范围|形式化结论|
|---|---|
|共同解析框架|质量约束的紧性和极大值存在性，平方距离积的缩放和重排，质量惩罚的无约束比较。极大点无碰撞。对质量惩罚后的对数积沿实仿射路径求导，得到复驻点方程和 Ward 恒等式。在质量为 \(n-1\) 时，驻点乘子为 \(n/2\)。Lean 也检查了两点 Gram 核的对角修正和求和。|
|十点候选|五条代数方程在论文的有理盒中恰有一根。Lean 核检查预条件固定点映射的有理展开、盒内映射和压缩上界。由该根定义的构型可行且无碰撞，并给出可达下界。|
|十点等号|径向罚项锁定三个平方半径，有限整数检查推出 \(2,4,4\) 人口。用内外带的 Chebyshev 特征恢复每层一阶与二阶矩，并验证二阶投影支点非零。两点层固定共同相位，四点层的矩确定两对对径点。Lean 验证所得构型与候选构型只相差旋转与重排。|
|十二点分带|在已知两条径向带及 Ward 预算时，用五个有理子区间的 Bernstein 控制值验证四条仿射下界，推出内带恰有三个点。|
|十二点候选和等号|四个三角形轨道的质量与矩。Vandermonde 恢复 \(3,3,6\) 人口及各层角矩。Newton 恒等式验证六点层的根多项式关于 \(z^3\) 为二次式，无碰撞推出两个不同的三次像。内层固定共同相位，外层的两个同模三次像由其非零和确定，恢复整个旋转和重排轨道。|
|覆盖规则|二分盒覆盖的健全性，张量 Bernstein 系数比较，十点所用的三阶对数上界、十二点所用的切线和二阶对数上界。|

`sharp_ten_of_gram_bounds` 与 `sharp_twelve_of_gram_bounds` 将逐边上界和 Gram 零矩接到最大值与唯一性定理。Lean 验证了角因子在三个逐对坐标中的多项式递推，也验证了正定矩阵分块核的有序和，以及等号蕴含每个复特征的零矩。`BlockBound` 接收各块的正定性、逐边上界、零和修正与二阶变分预算。具体证书还需将所用分块特征识别为 `TenGramMoments` 或 `TwelveGramMoments` 中的径向和角向条件。十点根盒系数可用 `tools/generate_root_data.py` 重新生成，Lean 在构建时直接验证这些系数。

## 全域证书的形式化状态

两案的 Lean 锐界定理目前带有显式证书前提。几何必要域、接触系数校正与矩阵正定性、全部局部接触盒和全域 Bernstein 叶子的逐项验收尚未在 Lean 中证明。十二点还需接入论文的 23 维驻点根盒和对称性证书。现有十点 50,267 节点与十二点 728,050 节点由 Python 精确验收器复现。 [十点全域证明](../ten-point-sharp-contacts/README.md) · [十二点全域证明](../twelve-point-sharp-global/README.md)
