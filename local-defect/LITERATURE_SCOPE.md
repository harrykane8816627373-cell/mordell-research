# 来源、近期文献与本地认证范围

研究资料与项目发起：蒋丰灿（Fengcan Jiang）；核查及整理包含 OpenAI Codex 协作。2026-10-30。

本轮的有限n主证明自含，不使用外部三角晶格全域最小性或线性渐近余项。它使用用户已证严格子集几何、仓库已有的圆盘精确能量分解，以及本文新的局部圆弧与匹配计数。Gram、Schur、圆盘矩和Onsager涂抹都是标准工具，不作框架发明或全球优先权声明。

经典来源包括：

- Lauritsen, *Floating Wigner Crystal and Periodic Jellium Configurations*, J. Math. Phys. 62, 083305 (2021)，[arXiv:2103.07975v2](https://arxiv.org/abs/2103.07975v2)。本文及仓库旧稿保留其jellium、重整化能量与圆盘下界的来源。本轮只比较项目已列显式上界，没有声称常数优于全部文献。
- Sandier–Serfaty, *2D Coulomb Gases and the Renormalized Energy*，[arXiv:1201.3503v2](https://arxiv.org/abs/1201.3503v2)。γ的本轮条件推论沿用仓库[渐近专题](../asymptotics/README.md)已经说明的这一展开接口；有限n线性缺额本身不依赖它。

检索还发现OpenAI于2026-09-23发布的 *Triangular minimality for planar Coulomb renormalized energy* [预印本目录](https://github.com/openai/math/tree/main/preprints/Triangular-minimality-for-planar-Coulomb-renormalized-energy-September-23-2026)。其[原始README](https://raw.githubusercontent.com/openai/math/main/preprints/Triangular-minimality-for-planar-Coulomb-renormalized-energy-September-23-2026/README.md)列出有限区间证书；[验收说明](https://raw.githubusercontent.com/openai/math/main/preprints/Triangular-minimality-for-planar-Coulomb-renormalized-energy-September-23-2026/verification/README.md)也明确区分计算验收与正文解析桥梁。

本轮没有全篇复核该预印本或执行其C++证书，因此**不将它作为本仓库已认证的定理输入**，也不对其正确性、同行评审状态作判断。“本项目尚未认证精确γ或三角全域最小性”只描述本地证明覆盖范围，不据此断言近期文献不存在相关结果。一般有限n最优构型的完整分类也不由一个无限系统的能量常数自动给出。

本轮同一研究代理完成推导和核验，用独立积分、不同代数表达式、已保存实际输入及完整旧实现作对照；没有独立代理或外部同行评审，也没有整篇形式化验证。
