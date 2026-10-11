# Optimal separable-measurement benchmark for chirality estimation

For $n$ independent copies of an unknown three-spin state, every universally
unbiased estimator based on a fully separable POVM on each copy obeys

$$
n\operatorname{Var}_{I_8/8}(\widehat\chi)\ge\frac9{16}.
$$

The bound permits classical adaptation between copies, arbitrary processing of
the complete measurement record, and local operations with classical
communication (LOCC) within each copy. Six uniformly randomized Pauli-product
settings attain it and have worst-case variance coefficient $9/16$. Thus $9/16$
is the exact minimax coefficient within this measurement class, not merely the
cost of one chosen Pauli decomposition.

Throughout, $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$ and
$\mathbf S=\boldsymbol\sigma/2$, as in [the model](../THEORY.md). The restriction
is on separability across the **three spins within each input copy**. The input
state itself is unrestricted.

## 1. Separable effects have a bounded chirality expectation

For a product state with single-qubit Bloch vectors $\mathbf r_j$,

$$
\operatorname{Tr}\!\left[\chi(\rho_1\otimes\rho_2\otimes\rho_3)\right]
=\frac18\mathbf r_1\cdot(\mathbf r_2\times\mathbf r_3).
$$

Since $\|\mathbf r_j\|\le1$, its absolute value is at most $1/8$. Orthogonal
unit Bloch vectors attain the bound. Convexity extends it to every fully
separable state $\sigma$:

$$
|\operatorname{Tr}(\chi\sigma)|\le\frac18.
$$

In particular, for any nonzero fully separable positive effect $E_y$, define

$$
q_y=\frac{\operatorname{Tr}E_y}{8},\qquad
a_y=\frac{\operatorname{Tr}(\chi E_y)}{\operatorname{Tr}E_y}.
$$

Then $|a_y|\le1/8$. For a POVM, $\sum_yq_y=1$. Zero-trace positive effects
are zero and can be omitted. The same argument holds with integrals for
continuous-outcome POVMs.

The product-state formula and direct Pauli measurements are established
chirality tools; see Reascos et al., [R23](../REFERENCES.md). The bound used here
is the fully separable value $1/8$, distinct from the biseparable-state value
$1/4$ used for genuine tripartite entanglement detection.

## 2. Direct one-copy proof

Suppose a real score $w_y$ is unbiased for chirality on every input state. This
is equivalent to the operator identity

$$
\sum_yw_yE_y=\chi.
$$

The Pauli decomposition gives $\operatorname{Tr}\chi=0$ and
$\operatorname{Tr}\chi^2=3/4$. Hence, on $I_8/8$,

$$
\sum_yq_yw_y=0,\qquad
\sum_yq_yw_ya_y=\frac{\operatorname{Tr}\chi^2}{8}=\frac3{32}.
$$

Weighted Cauchy–Schwarz now gives

$$
\left(\frac3{32}\right)^2
\le\left(\sum_yq_yw_y^2\right)\left(\sum_yq_ya_y^2\right)
\le\frac1{64}\sum_yq_yw_y^2.
$$

The score has zero mean at this input, so

$$
\boxed{\operatorname{Var}_{I_8/8}(w)\ge\frac9{16}.}
$$

Random setting choices are included in $y$ and their probabilities in $E_y$.
Random postprocessing cannot improve the bound: conditioning the score on the
observed record preserves its mean and can only reduce its variance.

This is an elementary witness-bound-to-variance argument. More generally, for
a traceless observable $A$ in dimension $d$, if every allowed normalized
effect has $|\operatorname{Tr}(AE)/\operatorname{Tr}E|\le b$, the same proof
gives $\operatorname{Var}_{I_d/d}(w)\ge[\operatorname{Tr}(A^2)/d]^2/b^2$.
No new general statistical principle is claimed.

## 3. Adaptive measurements and arbitrary estimators

The one-copy proof alone would not establish a bound for all estimators on a
multi-copy record. A one-parameter submodel supplies that extension:

$$
\rho_\theta=\frac{I_8}{8}+\frac43\theta\chi,
\qquad |\theta|<\frac{\sqrt3}{8}.
$$

The eigenvalues of $\chi$ are $0,\pm\sqrt3/4$, so this interval gives
strictly positive states, with $\operatorname{Tr}(\rho_\theta\chi)=\theta$.
For any allowed POVM,

$$
p_y(\theta)=q_y\left(1+\frac{32}{3}\theta a_y\right).
$$

Its classical score and Fisher information at zero satisfy

$$
s_y=\left.\partial_\theta\log p_y(\theta)\right|_0
=\frac{32}{3}a_y,
\qquad
F(0)=\frac{1024}{9}\sum_yq_ya_y^2\le\frac{16}{9}.
$$

Now measure $n$ copies sequentially. Conditional on any preceding classical
record $h$, the next measurement may be any fully separable POVM
$\{E_{y|h}\}$. The choice of measurement is a fixed rule applied to $h$, not
a rule with separate access to the unknown $\theta$. Each conditional score
has conditional mean zero and second moment at most $16/9$. Consequently,
the full-record score $S=\sum_{i=1}^ns_i$ obeys

$$
\mathbb E_0[S^2]=\sum_{i=1}^n\mathbb E_0[s_i^2]
\le\frac{16n}{9}.
$$

Cross terms vanish by conditional expectation. This covers deterministic
allocations, randomized settings, and classical adaptation.

Let $T$ be any finite-variance estimator that is locally unbiased on this
submodel: $\mathbb E_0T=0$ and
$\left.\partial_\theta\mathbb E_\theta T\right|_0=1$. Differentiating the
expectation yields $\mathbb E_0[TS]=1$. Cauchy–Schwarz proves

$$
1\le\operatorname{Var}_0(T)\,\mathbb E_0[S^2]
\le\operatorname{Var}_0(T)\frac{16n}{9},
$$

and hence $n\operatorname{Var}_0(T)\ge9/16$. The conditional scores are
bounded, so at fixed $n$ differentiation is justified for finite-variance
estimators, including continuous outcomes. Universal unbiasedness implies
the required local unbiasedness. No per-copy score representation for $T$
is assumed.

Every LOCC measurement within one copy has fully separable effects. The
restriction therefore includes local ancillas, local filtering, classical
communication among the three spins, and their coarse-grained outcomes,
provided all original input copies are counted and no entanglement resource
is supplied to change the effective measurement class. For the operational
inclusion of LOCC in separable operations, see Chitambar et al.,
[*Everything You Always Wanted to Know About LOCC (But Were Afraid to Ask)*,
Communications in Mathematical Physics **328**, 303–326 (2014)](https://doi.org/10.1007/s00220-014-1953-9),
Section 2.3; [arXiv:1210.4583](https://arxiv.org/abs/1210.4583).

## 4. Attainment and exact minimax value

Write

$$
\chi=\frac18\sum_{k=1}^6\epsilon_k P_k,
\qquad
(P_1,\ldots,P_6)=(XYZ,YZX,ZXY,XZY,YXZ,ZYX),
$$

with signs $(+,+,+,-,-,-)$. Independently choose $k$ uniformly, measure the
three corresponding local Pauli axes, and let $z=\pm1$ be the product of
their outcomes. Assign the score $w=3\epsilon_kz/4$. Then for every state,

$$
\mathbb E_\rho w=\operatorname{Tr}(\rho\chi),\qquad
\mathbb E_\rho w^2=\frac9{16}.
$$

The sample mean from $n$ independent trials is universally unbiased and obeys

$$
n\operatorname{Var}_\rho(\widehat\chi)
=\frac9{16}-\langle\chi\rangle_\rho^2\le\frac9{16}.
$$

Combining this upper bound with the mixed-input lower bound gives

$$
\boxed{\inf_{\mathcal P}\sup_\rho
n\operatorname{Var}_\rho(\widehat\chi_{\mathcal P})=\frac9{16},}
$$

where $\mathcal P$ ranges over fixed-$n$, universally unbiased protocols with
classically adaptive, fully separable measurements on each independent copy.
Uniform randomization proves attainment for every integer $n\ge1$. Equal
deterministic allocation among the six settings also attains the mixed-state
value when the allocation is exact; its variance elsewhere may be smaller.

## 5. Interpreting an advantage over this benchmark

A calibrated native parity estimator with mixed-input variance coefficient
below $9/16$ outperforms every estimator in the specified class on that same
input. An all-input upper bound below $9/16$ additionally proves a strict
minimax sampling advantage. It need not outperform every separable protocol
pointwise on every state.

The comparison concerns unbiased estimation from an unknown input state. If
the state is promised to be $I_8/8$, reporting zero without measuring solves
a different task. Biased or prior-dependent estimators, collective quantum
operations across distinct copies, and retained quantum memory between copies
are outside this theorem. Known entangling measurements, including the
controlled-cycle Hadamard test, are also outside the separable-effects class.

The result compares original input-copy counts under exact calibration and
the declared detector model. Coherent hopping, calibration effort, apparatus
complexity, and bias from model errors remain separate resource costs.
