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
|十点二阶变分|质量惩罚后的对数积沿任意复位移的二阶变分非正。Lean 将径向位移及其乘以 \(i\) 后的位移化成存档的两种逐边表达式，证明以任意两个实对称正半定矩阵加权时，两块在质量九极大构型上的总和非负。|
|十点零和修正|在质量九的极大构型上，证书中五种 Ward 矩对应的逐边核求和为零。逐点径向力行和质量行允许任意端点权重。分带人口为八时，人口行也有加权零和。给定八个外点的循环次序，三个短程角类减去两倍对径类的端点权重总和为零。|
|十二点分带|在已知两条径向带及 Ward 预算时，用五个有理子区间的 Bernstein 控制值验证四条仿射下界，推出内带恰有三个点。|
|十二点驻点根盒|Lean 使用存档中的 23 个坐标、中心、半径 \(2^{-80}\) 和预条件矩阵。核检查中心残差和逐对有理割线包络，证明全盒压缩以及唯一驻点。驻点到中心的距离至多为 \(2^{-205}\)，质量为 12，由此得到可达下界。|
|十二点候选|两种存档对称置换经旋转到规范坐标后仍落在原盒中。Lean 用唯一性推出它们固定驻点，并将质量 11 的驻点识别为四个三角形轨道。三个平方半径互异，四个种子的三次像互异，第三谐波投影支点和外层三次像之和均非零。|
|十二点等号|四个三角形轨道的质量与矩。Vandermonde 恢复 \(3,3,6\) 人口及各层角矩。Newton 恒等式验证六点层的根多项式关于 \(z^3\) 为二次式，无碰撞推出两个不同的三次像。内层固定共同相位，外层的两个同模三次像由其非零和确定，恢复整个旋转和重排轨道。|
|矩阵区间|十点 12 个 Gram 块和 2 个二阶变分块、十二点 12 个 Gram 块。Lean 对覆在接触修正系数外的有理区间验收三角合同变换的严格对角优势，并由证明推出每个区间中的全部实对称矩阵正定。|
|覆盖规则|仿射重缩放的幂基到张量 Bernstein 系数的精确转换。中点二分和向下舍入维持逐点下界。有限奇次项有理级数给出对数上包络，并用于切线比较、十二点二次上界和十点清分母后的三次上界。|
|十二点已有回放|核验存档 `1.txt` 的 39 个节点及其前序排列，逐项检查二分、17 个切线叶和 3 个几何剔除叶。Lean 核验了 500 个初始 Bernstein 系数的三次坐标变换和向下舍入，也证明被省略高次项的常数下界。导出的有理幂基多项式在两个内层半径、非负角坐标的有序 Gram 可行域上支配距离平方的对数。|
|十二点 Ward 修正|质量为 11 的极大构型满足前四个谐波的 Ward 矩条件。Lean 验证三个秩至多为二的块和一阶稀疏块的有序和为零，并证明分摊对角项后的逐边修正总和为零。|

`sharp_ten_of_gram_bounds` 与 `sharp_twelve_of_gram_bounds` 将逐边上界和 Gram 零矩接到最大值与唯一性定理。Lean 验证了角因子在三个逐对坐标中的多项式递推，也验证了正定矩阵分块核的有序和，以及等号蕴含每个复特征的零矩。`BlockBound` 接收各块的正定性、逐边上界、零和修正与二阶变分预算。具体证书还需将所用分块特征识别为 `TenGramMoments` 或 `TwelveGramMoments` 中的径向和角向条件。十点根盒系数可用 `tools/generate_root_data.py` 重新生成，Lean 在构建时直接验证这些系数。

## 全域证书的形式化状态

两案的 Lean 锐界定理仍带有显式证书前提。`sharp_ten_of_gram_bounds` 和 `TwelveRootBox.isolated_sharp_of_gram_bounds` 分别使用两案的精确候选。十二点 `InnerCoverData.signed_inner_archive_power_bound` 已将这个根盒的逐点比较接到导出的有理幂基多项式。`signed_inner_coefficients_bound` 把后续衔接转成真实核的逐项系数下界。精确接触修正的系数区间、几何必要域、局部接触盒和其余全域叶子仍待验证。十点还需在 Lean 中证明必要角域给出的角类与外点循环次序一致。现有十点 50,267 节点与十二点 728,050 节点由 Python 精确验收器复现。 [十点全域证明](../ten-point-sharp-contacts/README.md) · [十二点全域证明](../twelve-point-sharp-global/README.md)

`tools/generate_gram_data.py` 和 `tools/generate_twelve_root_data.py` 分别生成 `Mordell/GramData.lean` 和 `Mordell/TwelveRootData.lean`。生成器读取仓库已有的证书数据。Lean 重新检查导出的矩阵区间、逐对割线包络和压缩行和。

`tools/generate_inner_cover.py` 生成 `Mordell/InnerCoverData.lean`，参数 `--checks` 生成 `Mordell/InnerCoverChecks.lean`，参数 `--initial` 生成 `Mordell/InnerInitData.lean` 中的幂基系数和变换见证。每个节点的整数不等式分开交给 Lean 核，以控制重放时的内存开销。
