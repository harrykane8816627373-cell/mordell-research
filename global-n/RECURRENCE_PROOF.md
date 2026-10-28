# A strict cross-cardinality upper bound from subset-moment variance

Research initiated by Fengcan Jiang (蒋丰灿); proof development and verification include OpenAI Codex collaboration. 2026-10-27. The project studies

```math
M_n=\max_{\sum_i|z_i|^2=n}\prod_{i\lt j}|z_i-z_j|^2,
\qquad N_n=\frac{n(n-1)}2.
```

Source correction (2026-10-28): the supplied complete revised manuscript already proves the ordinary recurrence and strict adjacent monotonicity, with the same force-equation argument below. These are restored source results, not first obtained in this note. The general planar subset-variance identity and its quantitative loss are retained as subsequent supplements. See ../source-review/SOURCE_REVIEW.md. This note supplies a structural recurrence and a quantitative strictness factor. It does **not** determine every exact value of M_n. The ordinary deletion/subset comparison is elementary and may already be classical; no priority claim is made. The new use here is to retain its variance loss explicitly in terms of planar fourth moments. In the comparisons actually made below, this recurrence does not improve the previously published disk-smearing upper envelope.

## 1. Exact variance identity

Let n >= 4 and let z_1,...,z_n be any centered planar configuration with sum |z_i|^2=n. Collisions are allowed for the identity. For a uniformly random m-element subset S, 2 <= m < n, put

```math
Q_S=\sum_{i\in S}|z_i-\bar z_S|^2,
\quad R_4=\sum_i|z_i|^4,
\quad P_2=\left|\sum_i z_i^2\right|^2,
\quad h=n-m.
```

Then

```math
\mu:=\mathbb E Q_S=\frac{n(m-1)}{n-1},
```

```math
\boxed{\begin{aligned}
\operatorname{Var}(Q_S)
={}&\frac{n(m-1)h(h-1)}{m(n-1)^2(n-2)}\\
&+\frac{(m-1)h\big((n-1)(m-1)-2\big)}
 {m(n-1)(n-2)(n-3)}(R_4-n)\\
&+\frac{(m-1)h(h-1)}
 {mn(n-1)(n-2)(n-3)}P_2.
\end{aligned}}
```

In particular, since R_4 >= n by Cauchy--Schwarz and P_2 >= 0,

```math
\operatorname{Var}(Q_S)\ge
v_{n,m}:=\frac{n(m-1)(n-m)(n-m-1)}{m(n-1)^2(n-2)}.
```

This bound is positive for 2 <= m <= n-2. Equality holds exactly when all radii are one and sum z_i^2=0 in that range, because both additional coefficients are positive. Such configurations exist, for example the unit regular n-gon. This equality statement concerns the variance bound and does not assert that those polygons maximize the distance product.

For m=n-1, the formula reduces to

```math
\operatorname{Var}(Q_S)=\frac{n}{(n-1)^2}(R_4-n).
```

Indeed deleting i leaves Q_S=n-n|z_i|^2/(n-1). Equal radii therefore make every deletion moment equal. A positive universal variance loss is impossible at this endpoint. The endpoint formula remains valid for n=3, without the n-3 denominators in the full displayed identity.

### Derivation

Write w_{ij}=|z_i-z_j|^2 and W_S=sum_{i<j in S} w_{ij}; then Q_S=W_S/m. Centering gives

```math
\sum_{i\lt j}w_{ij}=n^2,
\qquad s_i:=\sum_{j\ne i}w_{ij}=n(1+|z_i|^2).
```

The exact fourth-distance identity is

```math
A:=\sum_{i\lt j}w_{ij}^2=nR_4+2n^2+P_2.
```

To check it, write z_i=x_i+iy_i and let C=sum_i (x_i,y_i)(x_i,y_i)^T. Expanding the fourth powers yields nR_4+n^2+2 tr(C^2), and 2 tr(C^2)=n^2+P_2. Every term involving sum_i z_i disappears by centering.

Let B be the sum of products w_e w_f over unordered distinct edge pairs sharing a vertex, and let C_0 be the corresponding sum over disjoint edge pairs. Then

```math
2B=\sum_i s_i^2-2A=n^2(R_4+3n)-2A,
\qquad 2C_0=n^4-A-2B.
```

Put p_j=(m)_j/(n)_j, with p_j=0 if j>m. Since a pair of edges uses either two, three or four vertices,

```math
\mathbb E W_S^2
=p_2 A+2p_3 B+2p_4 C_0
=(p_2-2p_3+p_4)A+(p_3-p_4)n^2(R_4+3n)+p_4n^4.
```

Substitute A, subtract (E W_S)^2, and divide by m^2. For explicit algebra checking,

```math
p_2-2p_3+p_4=\frac{m(m-1)h(h-1)}{n(n-1)(n-2)(n-3)},
```

```math
p_3-p_4=\frac{m(m-1)(m-2)h}{n(n-1)(n-2)(n-3)}.
```

This gives the boxed identity, including the vanishing h-1 endpoint.

## 2. Ordinary cross-cardinality comparison

Take a global maximizing n-point configuration, which is centered and collision-free. Every m-subset has

```math
\Delta_m(z_S)\le M_m\left(\frac{Q_S}{m}\right)^{N_m}.
```

This follows by translating the subset to its centroid and scaling it to total moment m. Each pair is contained in exactly binomial(n-2,m-2) subsets, so multiplying all these inequalities and taking logarithms gives

```math
\log M_n\le\frac{N_n}{N_m}\log M_m
 +N_n\big(\mathbb E\log Q_S-\log m\big).
```

Concavity of log and the exact mean produce

```math
M_n\le M_m^{N_n/N_m}
\left(\frac{n(m-1)}{m(n-1)}\right)^{N_n}.
```

Equivalently, the quantities ((n-1)/n) M_n^(1/N_n) are nonincreasing with n. For n=3,m=2,M_2=4 this reproduces M_3<=27, attained by the equilateral triangle.

### Strict adjacent monotonicity

In fact define D_n=((n-1)/n)M_n^(1/N_n). Then

```math
\boxed{D_2=D_3=2,\qquad D_n\lt D_{n-1}\quad(n\ge4).}
```

Suppose equality held in the deletion comparison m=n-1 for some n>=4. Both preceding inequalities must then be equalities. Strict concavity of log implies that all deletion moments Q_i are equal, and the identity

```math
Q_i=n-\frac{n}{n-1}|z_i|^2
```

forces every |z_i|^2=1 and Q_i=n(n-2)/(n-1). Moreover every deletion subset, after the indicated translation and scaling, must itself attain M_(n-1).

For any collision-free maximizing configuration with cardinality k and centered total moment Q, differentiating log Delta under the moment constraint gives

```math
\sum_{\ell\ne j}\frac1{z_j-z_\ell}
=\frac{k(k-1)}{2Q}(\bar z_j-\bar a),
```

where a is its centroid. Thus the n-point stationary equation is (n-1)bar z_j/2. The subset deleting i has centroid -z_i/(n-1) and multiplier (n-1)^2/(2n). Subtracting its stationary equation at j from the full equation gives

```math
\frac1{z_j-z_i}
=\frac{n-1}{2n}(\bar z_j-\bar z_i),
\qquad |z_j-z_i|^2=\frac{2n}{n-1}\quad(i\ne j).
```

This would be an equidistant planar set of n points. Centering its real Gram matrix gives diagonal entries one and off-diagonal entries -1/(n-1), a matrix of rank n-1. A planar real Gram matrix has rank at most two, contradicting n>=4. Hence adjacent monotonicity is strict. Chaining these inequalities also makes the ordinary m-subset comparison strict for every 2<=m<n once n>=4. The argument gives no explicit adjacent loss; the next section supplies one whenever m<=n-2.

## 3. Keeping an explicit strictness factor

Suppose 0<Q_S<=L for every subset under consideration and L>=mu. An elementary strengthened logarithm inequality is

```math
\log x\le x-1-\frac{(x-1)^2}{2\max\{1,x\}},\qquad x\gt0.
```

For x>=1, subtract the rightmost expression from x-1-log x; its derivative is (x-1)^2/(2x^2)>=0 and its value at one is zero. For 0<x<=1, the derivative of x-1-log x-(x-1)^2/2 is -(x-1)^2/x<=0 and its value at one is zero. Applying the inequality to Q_S/mu and replacing max{1,Q_S/mu} by L/mu yields

```math
\mathbb E\log Q_S\le\log\mu-
\frac{\operatorname{Var}(Q_S)}{2\mu L}.
```

Hence, for every n>=4 and 2<=m<=n-2,

```math
\boxed{\displaystyle
M_n\le M_m^{N_n/N_m}
\left(\frac{n(m-1)}{m(n-1)}\right)^{N_n}
\exp\!\left[-\frac{n(n-m)(n-m-1)}{4m(n-2)L}\right].}
```

The simplest completely general choice is L=n, since Q_S<=sum_i|z_i|^2=n. For a global maximizer the previously proved strict subset theorem gives the stronger valid choice

```math
L=\min\left\{n,\frac{2m(n-1)}n,
n-\frac{h(h-1)}{n-1}\right\}.
```

The second term uses the established pointwise radius bound |z_i|^2<2(n-1)/n. The third term uses the centroid decomposition n=Q_S+Q_{S^c}+(mh/n)|bar z_S-bar z_{S^c}|^2 and the established complement moment bound. Taking their closed upper bounds is legitimate.

The published [explicit subset slack](../finite-n/SUBSET_PROOF.md) can be combined directly with this recurrence. With exactly the same positive rational epsilon constants defined there, one may use the smaller rational quantity

```math
L^+_{n,m}=\min\left\{n,\frac{2m(n-1-\varepsilon_{n,n-1})}{n},
n-\frac{h(h-1)+2\varepsilon_{n,h}}{n-1}\right\}
\qquad(h=n-m\ge2).
```

Every global maximizing subset still obeys Q_S<=L^+_{n,m}; in particular its average mu does too. The same proof therefore gives the displayed exponential discount with L replaced by L^+_{n,m}. This is a direct combination of the previous strict-subset method with the fourth-moment variance identity.

Thus the ordinary cross-cardinality upper bound is strictly reduced by a specified factor for every 2<=m<=n-2. This prevents equality in that ordinary comparison even when the unknown extremizer has equal radii and isotropic second moment; the randomness of proper subsets retains a strictly positive variance.

This is not a solution of the finite-n maximization problem. Using the known sharp M_2,...,M_7 in this formula for n=4,...,15 gives bounds above the already available disk-smearing upper bound U_n. The theorem gives an additional exact compatibility condition between different cardinalities and supplies a loss that can be combined with future sharper subset information. It does not justify announcing an improved overall upper envelope without a separate comparison.

## 4. Independent exact checks

`verify_recursion_variance.py` uses Fraction arithmetic. It enumerates every subset directly for several integer/rational planar families, compares the empirical variance with the formula involving R_4 and P_2, separately verifies the fourth-distance identity, checks m=n-1 directly, and checks equality for isotropic equal-radius squares and unions of rotated squares. These finite checks exercise the algebra; the proof above establishes the theorem for all n.
