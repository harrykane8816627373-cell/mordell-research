# 用Jacobian组合十二点核的精确解析梯度

2026-11-08。蒋丰灿研究项目，OpenAI Codex协作。本目录固定引用[此前已接受的十二点3%核](../twelve-point-kernel-refinement/kernel.json)，没有随后续候选自动漂移。它保存29次真实Jacobian 0.23.0原生求导的完整结果，并提供不依赖Jacobian的标准库重放。

本轮结论是梯度的完整系数恒等式及物理边界接触条件；**没有由这些导数推出新的全局上界或最大构型分类**。此前全局界仍由独立的连续覆盖证明承担。

## 1. 实际计算了什么

将质量11尺度的核上界多项式写成

\[
h(u,v,t)=h_{\mathrm{affine}}(u,v)+\sum_{k=0}^{11}h_k(u,v,t).
\]

各谐波块包含交互项及对角修正：

\[
h_k=-H_k(u,v,t)F_k(u)^TQ_kF_k(v)
-\frac{u^kF_k(u)^TQ_kF_k(u)+v^kF_k(v)^TQ_kF_k(v)}{22},
\]

其中H₀=1、H₁=t、Hₖ=2tHₖ₋₁−uvHₖ₋₂。代码固定使用[此前覆盖程序](../twelve-point-kernel-upper/verify_cover.py)的精确特征表示；另外用其完整展开独立核对分块之和。

完整h有602项，三个偏导分别有527、527、500项。直接请求整块梯度遇到两个真实操作准入限制：系数最多128位十进制数字，三个导数的总项数最多256。部分源系数约210位，不能只按谐波分块后直接提交。

本目录使用精确线性组合：

1. 每个块取公共分母D，使全部系数成为整数除以D。
2. 将带符号整数系数按2²⁵⁰进制分成小整数多项式Pⱼ。
3. 每次再取最多64个源单项式；这保证总梯度支持至多192项，输入整数系数至多76位。
4. 实际调用Jacobian的polynomial_field.scalar.gradient.compute原生实现。
5. 以精确比例2²⁵⁰ʲ/D恢复每个结果，再逐项求和。

这是对求导线性的精确使用，未修改Jacobian护栏，也未舍入任何系数。[native_blocks.json](native_blocks.json)保留29次真实原生输出及其源多项式、恢复比例、谐波编号和整数分块编号。[gradient.json](gradient.json)保存组合后的完整h和梯度。

## 2. 标准库接受承担哪些检查

从仓库根目录执行：

~~~bash
python -B -S -O jacobian-gradient/verify_all.py
~~~

仅需Python标准库；不导入Jacobian、NumPy或SciPy。接受器：

- 核对本目录证据及旧核、旧覆盖程序的固定哈希。
- 从旧核和特征重新构造全部12个hₖ及仿射项。
- 将29个原生输出保留的源多项式乘其精确恢复比例，逐系数重建对应谐波块。
- 对每一个原生结果检查完整单项式求导恒等式；组合结果还须与完整h的三个普通单项式偏导逐系数相等。
- 核对4个有理点上的残差梯度及接受器的全部固定支持计数。
- 核对[真实MCP边界梯度](boundary_gradient.json)及下面的正负分支方向。

因此接受不只是“抽四个点相等”，也不以报告中的passed字段代替计算。4个点是额外交叉检查；系数恒等式才承担完整导数结论。Jacobian原生计算的实际记录见[生成报告](report.json)，标准库实际接受结果见[重放报告](verification_report.json)。

## 3. 给新核生成梯度，不覆盖已接受数据

只有明确传入--build时才导入并调用Jacobian。使用安装了Jacobian 0.23.0的Python运行：

~~~bash
python -B jacobian-gradient/analytic.py --build
python -B jacobian-gradient/analytic.py --build --kernel work/new-kernel.json --output work/new-gradient
~~~

默认输入仍是上面的已接受3%核；默认输出是work/jacobian-gradient-generation。生成命令禁止向本目录及其子目录写入，接受过的数据不会被后续生成覆盖。新输出的Residual实例明确从新输出目录读取，不会悄悄继续使用旧梯度。

新核生成结果只是新的辅助计算材料，仍需检查对应完整系数恒等式、Gram条件及连续域。把新梯度文件放进数值搜索不会自动产生全局证书。

## 4. 在数值搜索中传入解析jac

Residual默认读取本目录接受的gradient.json。浮点求值使用这些精确系数转换后的双精度数；舍入及高幂抵消仍存在，它是发现工具，不是区间证明。

~~~python
import importlib.util
from scipy.optimize import minimize

spec = importlib.util.spec_from_file_location("kernel_gradient", "jacobian-gradient/analytic.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
residual = module.Residual()
U = 121 / 72
constraints = [
    {"type": "ineq", "fun": lambda p: p[0]*p[1] - p[2]**2,
     "jac": lambda p: (p[1], p[0], -2*p[2])},
    {"type": "ineq", "fun": lambda p: p[0] + p[1] - 2*p[2] - 1/3,
     "jac": lambda p: (1, 1, -2)},
]
result = minimize(
    residual.value, (0.5264, U, 0.9033), jac=residual.jac,
    method="SLSQP", bounds=((0, U), (0, U), (-U, U)),
    constraints=constraints, options={"ftol": 1e-12, "maxiter": 200},
)
~~~

目标为G=h−log d，d=u+v−2t；其解析梯度为∇h−(1,1,−2)/d。[有限数值诊断](numeric_crosscheck.json)在83个点与既有发现程序的梯度公式比对，最大双精度差约5.48×10⁻¹⁰，未观察到公式错误；旧违反点搜索没有向SLSQP传入解析jac。4个Jacobian生成梯度的有理点评估与精确值的最大绝对差约1.15×10⁻¹⁰。这些数值差只针对所测点，不是连续域舍入误差上界，更不是全域极小值证明。

可选的check_discover.py需既有SDP环境中的NumPy、SciPy、CVXPY；它重新运行这项有限诊断及四个SLSQP例子。正式标准库接受不执行它，也不把其收敛状态当作数学结论。

## 5. 物理边界的正确接触方向

设u,v>0、t²=uv，且其余域边界均不活跃。若majorant残差G在该点接触0并在可行域内非负，则KKT条件为

\[
G_u+\frac{t}{2u}G_t=0,\qquad
G_v+\frac{t}{2v}G_t=0,\qquad tG_t\le0.
\]

理由是物理域为g=uv−t²≥0，实际MCP计算给∇g=(v,u,−2t)。局部最小点的法向关系为∇G=λ∇g、λ≥0，所以tGₜ=−2λt²≤0；两个切向等式亦随即得到。

因此t>0时需要Gₜ≤0，t<0时需要Gₜ≥0。统一写Gₜ≥0只覆盖负t分支。这里保留全部[旧数值发现代码](../candidate-kernel/discover_sdp.py)；新搜索示例通过正确的物理约束及其解析梯度处理两个分支。若另有半径端点或点距边界活跃，必须加入相应乘子；若u或v为零，也不能使用含1/u、1/v的这些切向公式。

本目录没有把KKT必要条件误作充分条件。它们用于构造和排查候选核；最终全局上界仍须完整连续证明及正确的等号衔接。
