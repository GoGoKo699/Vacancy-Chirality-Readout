# Finite calibration with unconditional mean-square error

Three independent batches of product references give a finite polynomial
calibration whose risk averages over both reference and target data. There is
no discarded calibration record or unbudgeted confidence-failure event. The
construction has a large but explicit reference cost; reuse can amortize that
cost against the passive-parity lower bound. The counts below are sufficient
mathematical bounds, not measured instrument performance or necessary costs.

## 1. Coordinates, reference preparations, and certified inputs

Use the spin and reflection covariance assumptions of
[direct channel calibration](CHANNEL_CALIBRATION.md). Include every preparation
and recorded outcome. Let

```math
x_\rho=(\langle\chi\rangle_\rho,\langle D\rangle_\rho/4)^T,
\qquad D=P_{12}-P_{23},\qquad h=\sqrt3/4.
```

The identities $\{\chi,D\}=0$, $\chi^2=3(I-\Pi_{\rm s})/16$, and
$D^2=3(I-\Pi_{\rm s})$ imply $\|x_\rho\|\le h$ for every input.
From each randomized two-flux image record the fixed two-feature vector
$Z=(a,v/4)^T$: $a$ is the improved certified chirality score and $v$ the
certified auxiliary $D$ score. Features include their flux signs and have
zero score on the prescribed zero-score events. The complete channel obeys

```math
\mathbb E_\rho Z=Mx_\rho.
```

The known product references
$\rho_x=|+x,+y,+z\rangle\langle+x,+y,+z|$ and
$\rho_d=|+z,+z,-z\rangle\langle+z,+z,-z|$ have coordinates
$(1/8,0)^T$ and $(0,1/4)^T$. Therefore
$M=[8\mathbb E_xZ,4\mathbb E_dZ]$. No separately known loss or detector-error
rate enters this operational identity.

Throughout equal spin-independent loss $0\le\gamma t\le0.01$ and independent
equal parity-bit errors $0\le q\le0.003$, the
[integer/rational certificate](CERTIFICATE.md) supplies the following outward
bounds for these fixed features. Its report derives them from the certified
finite-dimensional effects and the analytic noise inequalities.

```math
B=I-M,\qquad \|B\|\le b,\quad
\|e_1^TB\|\le r,\quad \|M\|\le\mu,
```

```math
\sup_\rho\mathbb E_\rho\|Z\|^2\le H.
```

| Bound | Outward rational value |
|---|---:|
| $b$ | $0.103947447477219$ |
| $r$ | $0.043400571947183$ |
| $\mu$ | $1.065494885645340$ |
| $H$ | $44.219368969990164$ |

For the deterministic rows $s_1=e_1^T$, $s_{j+1}=e_1^T+s_jB$, the same
certificate gives

```math
\sup_\rho\mathbb E_\rho(s_jZ)^2\le V_j.
```

| Row | Second-moment bound |
|---|---:|
| $V_1$ | $0.344784300185690$ |
| $V_2$ | $0.395119522625782$ |
| $V_3$ | $0.370837757092059$ |
| $V_4$ | $0.366919338686148$ |

These are second moments, so they also bound variances. The checker below
requires the produced certificate and verifies all these dependencies; it does
not treat the decimal constants as independent assumptions about the model.

### Simultaneous parameter and noise uncertainty

The same tables also hold throughout the smaller joint region proved in
[the parameter-robustness account](PARAMETER_ROBUSTNESS.md):

```math
\eta\le0.001,\qquad 0\le\gamma t\le0.005,\qquad 0\le q\le0.0015.
```

Here $\eta$ bounds the lossless measurement-dilation distance from the
nominal instrument. A sufficient physical parameter condition is

```math
1.643|u-16|+3.163|\tau-\pi|
+4.086|\Phi_++2\pi/3|\le0.001,
```

with the same $u,\tau$ at both opposite fluxes. Spin and paired-reflection
covariance remain required. The joint certificate derives the raw response
and feature moments after combining parameter variation with the specified
loss and bit errors. Each of its bounds for $b,r,\mu,H,V_1,\ldots,V_4$
is at most the corresponding table entry above.

Every scalar recurrence and risk bound below is nondecreasing in these
nonnegative inputs. Consequently **both finite preparation budgets below
apply unchanged in this joint region**, including the additional preparation
and drift allowances of Section 4. This conclusion uses the joint moment
certificate, rather than combining independently stated tolerances. The
stable example uses one fixed instrument across reference and target trials;
Section 4 accounts separately for bounded differences between them.

The checker verifies the joint certificate's source window and input
dominance, and requires its nominal-certificate hash to match the parity
certificate supplied to the same run.

## 2. Independent polynomial calibration

In stage $j$, collect independent reference batches of sizes $N_{x,j}$ and
$N_{d,j}$, and form

```math
\widehat M_j=[8\overline Z_{x,j},4\overline Z_{d,j}],\qquad
\widehat B_j=I-\widehat M_j.
```

The reference states, instrument, and features are fixed in advance. All three
stages and subsequent target samples are independent. Define random row
weights

```math
w_1=e_1^T,\qquad w_{j+1}=e_1^T+w_j\widehat B_j\quad(j=1,2,3),
\qquad \widehat\chi=w_4\overline Z_\rho.
```

No random inverse is used. Independence and $\mathbb E\widehat B_j=B$ give
$\mathbb Ew_j=s_j$ and

```math
s_4M=e_1^T(I-B^4),\qquad
|\mathbb E\widehat\chi-\langle\chi\rangle|
\le h r b^3.
```

The identity follows by multiplying the degree-three geometric polynomial by
$I-B$. It is an exact expectation statement for finite reference batches.

### Exact covariance recursion

Write $K_j=\mathrm{Cov}(w_j)$, interpreting row covariance as a
$2\times2$ matrix, $Q_j=K_j+s_j^Ts_j$, and
$\Sigma_x=\mathrm{Cov}_xZ$, $\Sigma_d=\mathrm{Cov}_dZ$.
The two columns of $\widehat M_j-M$ are independent and centered. Conditional
on $w_j$, their covariance contribution is diagonal. Consequently

```math
\boxed{
K_{j+1}=B^TK_jB+
\mathrm{diag}\!\left(
\frac{64}{N_{x,j}}\mathrm{Tr}(\Sigma_xQ_j),
\frac{16}{N_{d,j}}\mathrm{Tr}(\Sigma_dQ_j)
\right),\quad K_1=0.
}
```

The cross terms vanish by stage independence and centering. The off-diagonal
covariance between the two features from the same image remains inside each
$\Sigma$; it must not be replaced with two independent feature variances.

For target covariance $\Sigma_\rho$ and mean $m_\rho=Mx_\rho$, conditioning on
$w_4$ gives the exact decomposition

```math
\mathbb E(\widehat\chi-\langle\chi\rangle)^2
=\frac{\mathrm{Tr}(\Sigma_\rho Q_4)}n
+m_\rho^TK_4m_\rho
+(s_4m_\rho-\langle\chi\rangle)^2.
```

This averages over calibration and target sampling. It also shows why shared
calibration errors do not disappear when the same weights are reused.

### Scalar enclosure and finite example

Allocate each stage $4:1$, so $N_{x,j}=4N_j/5$, $N_{d,j}=N_j/5$.
The equal coefficients $64/N_x=16/N_d=80/N_j$ and the previous moment bounds
give $\|K_j\|\le c_j$, where

```math
c_1=0,\qquad
\boxed{c_{j+1}=b^2c_j+\frac{80}{N_j}(V_j+Hc_j).}
```

Indeed $\mathrm{Tr}(\Sigma Q_j)\le V_j+Hc_j$, using positivity and
$\mathrm{Tr}\Sigma\le H$. The uniform risk enclosure is

```math
\boxed{R\le\frac{V_4+Hc_4}{n}+h^2\mu^2c_4+h^2r^2b^6.}
```

For $N_1=1{,}800{,}000$, $N_2=18{,}000{,}000$,
$N_3=180{,}000{,}000$ and $n=382{,}210$, exact rational arithmetic proves

```math
R<9.999796322\times10^{-7}<10^{-6}.
```

All 199.8 million references are counted. No pilot, validation rejection, or
unrecorded failure is omitted: the fixed features come from the mathematical
certificate and the polynomial estimator always returns a value. Additional
experimental characterization or monitoring would carry additional cost.

## 3. A comparison for many distinct unknown states

For the four-mode passive-parity class in
[the passive bound](PASSIVE_PARITY_BOUND.md), the submodel
$\rho_\theta=\Pi_{\rm m}/4+(4/3)\theta\chi$ has per-copy Fisher information

```math
F(\theta)\le I(t)=\frac{16(6+t)}{3(1-t^2)(13+t)},\qquad t=|\theta|/h.
```

Use independent cosine-squared priors on $[-3h/40,3h/40]$ for $K$ distinct
target parameters. Each prior contributes $25600\pi^2/27$ information, and
$I(3/40)=2073600/832093$. A policy may choose which task to sample and which
allowed measurement to use from its whole preceding record. A conditional
likelihood factor for task $k$ depends directly only on $\theta_k$; the
policy itself is a fixed function of the observed history. The conditional
score increments have zero mean and add in Fisher information. Thus the
joint-prior Fisher contribution for $\theta_k$ is at most
$I(3/40)\mathbb E N_k$. This uses the uniform pointwise information bound
on the prior interval; a prior-averaged per-copy bound cannot be substituted
under adaptive allocation. Known references add no parameter information.

Applying the scalar integration-by-parts and Cauchy--Schwarz argument to
coordinate $k$, with all other parameters and data included in the record,
gives

```math
\mathbb E(\widehat\theta_k-\theta_k)^2
\ge\left[I(3/40)\mathbb E N_k+25600\pi^2/27\right]^{-1}.
```

If each task has mean-square error at most $\varepsilon^2$ uniformly in all
target states, summing the necessary bounds yields

```math
\boxed{\mathbb E\sum_{k=1}^K N_k\ge K L_\varepsilon,\qquad
L_\varepsilon=\frac{832093}{2073600}
\left(\varepsilon^{-2}-\frac{25600\pi^2}{27}\right).}
```

The statement covers randomized allocation and stopping with finite expected
sample count under the joint prior. Here every expectation in the
expected-count lower bound averages over both the product prior and the
measurement record. It therefore bounds the worst-case expected total count
from below, without asserting a pointwise count bound at every target tuple.
For a deterministic or almost-sure total budget, that integer budget is at
least $\lceil K L_\varepsilon\rceil$. A worst-case expected budget alone need
not be an integer and should not be rounded up.

For completeness, let $N_{k,m}$ count task-$k$ samples among the first $m$
preparations, and let $S_{k,m}$ be their likelihood score sum, with zero
increments after stopping. Conditional on the parameter vector, these
scores are martingales. Their increments are bounded on the prior interval
by $1/(h-3h/40)$, because $q_z\ge h|c_z|$. Orthogonality and the uniform
information bound give an $L^2$ limit $S_k$ with

```math
\mathbb E|S_k-S_{k,m}|^2
\le I(3/40)\,\mathbb E(N_k-N_{k,m})\longrightarrow0.
```

Clip each estimator to $[-h,h]$, which cannot increase its squared error.
Its conditional expectations given the first $m$ observed records, taken
under the fixed joint prior, are bounded functions of those records and
converge in $L^2$ to the clipped estimator. They are parameter-independent
functions, so finite-record integration by parts applies to each one.
Passing to the limit by Cauchy--Schwarz gives
$\mathbb E[(\widehat\theta_k-\theta_k)
(S_k+\partial_{\theta_k}\log\pi_a(\theta_k))]=1$.
Here $a=3h/40$. The prior score has variance $25600\pi^2/27$ and zero cross
term with $S_k$;
this proves the stated bound for arbitrary finite-risk estimators and
finite-expected-count stopping rules. Removing the clipping only increases
the risk.

At $\varepsilon=10^{-3}$, $L_\varepsilon=397524.305991122\ldots$.
For $K=13{,}100$ distinct tasks the native construction costs
5,206,751,000 preparations, while the passive lower bound on any fixed total
budget is 5,207,568,409. The comparison counts equal-cost preparations and
gives the comparator free reference calibration. It does not establish a
wall-time advantage. The targets must be fixed independently of calibration
randomness. Repeating one identical state is not an arbitrary number of
distinct tasks, and a promised common dynamical model changes the comparison.

## 4. Preparation error and drift

Suppose the reference means instead estimate $\widetilde M=M+E$, while the
target channel has $M_t=M+F$, with $\|E\|\le e$, $\|F\|\le d$. Reference
preparation errors are fixed across the three independent stages and lie in
the same three-spin orbital embedding. Set $C=B-E$, $q_{\ast}=b+e$. The expected
weights are now $s'_4=e_1^T(I+C+C^2+C^3)$, and

```math
s'_4M_t-e_1^T=-e_1^TC^4+s'_4(F-E).
```

Define

```math
S'_4=1+(r+e)(1+q_{\ast}+q_{\ast}^2),\qquad
\beta_{\ast}=h\left[(r+e)q_{\ast}^3+S'_4(e+d)\right].
```

Then $\|s'_4\|\le S'_4$ and the bias is at most $\beta_{\ast}$.
To propagate moment bounds, use

```math
S_1^0=1,\quad S_j^0=1+r\sum_{k=0}^{j-2}b^k\ (j\ge2),\qquad
\delta_1=0,\quad \delta_{j+1}=q_{\ast}\delta_j+eS_j^0,
```

```math
V'_j=(\sqrt{V_j}+\sqrt H\,\delta_j)^2,\qquad
c'_1=0,\quad c'_{j+1}=q_{\ast}^2c'_j+\frac{80}{N_j}(V'_j+Hc'_j).
```

Here $\|s'_j-s_j\|\le\delta_j$ follows by induction; the $L^2$ triangle
inequality gives $V'_j$. The covariance proof above applies with $C$.
Actual imperfect spin references satisfy the same all-input moment bound.

Mean-response drift alone does not control second moments. Also assume

```math
\sup_\rho\left\|\mathbb E_{\rho,t}ZZ^T-
\mathbb E_{\rho,\mathrm{cal}}ZZ^T\right\|\le t_{\ast}.
```

The exact risk decomposition and $\mathrm{Tr}K'_4\le2c'_4$ give

```math
\boxed{R\le
\frac{V'_4+Hc'_4+t_{\ast}[(S'_4)^2+2c'_4]}n
+h^2(\mu+d)^2c'_4+\beta_{\ast}^2.}
```

With $e=0.000081$, $d=0.00002$, $t_{\ast}=0.001$, stages of 2 million,
20 million, and 200 million references and $n=384{,}220$ give
$R<9.999754581\times10^{-7}$.
For $K=16{,}800$, the total is 6,676,896,000 preparations, below the fixed
passive budget lower bound 6,678,408,341. The drift bounds are explicit
instrument assumptions, not estimates obtained from the two references.
Both channels must retain the stated two-coordinate structure; arbitrary
symmetry-breaking errors require separate operator bounds.

### Angular errors with balanced product probes

There is a useful quadratic preparation tolerance for pure spin states.
Prepare the same actual triad of Bloch vectors $(a,b,c)$ in its three cyclic
permutations, chosen independently and uniformly. The averaged $D$ vanishes
exactly, while $\chi=\det(a,b,c)/8$ is unchanged. If each vector is within
angle $\delta$ of its corresponding intended orthogonal axis, then

```math
t_x=\det(a,b,c)\ge\cos^3\delta-\tfrac32\sin^2\delta-\sin^3\delta
\ge1-3\delta^2-\delta^3.
```

For the first inequality expand the determinant. The three quadratic
off-diagonal terms are bounded by half the total off-diagonal squared norm;
the sum of the two cubic terms is at most the product of the three
off-diagonal vector norms. The second uses $\sin\delta\le\delta$ and
$\cos\delta\ge1-\delta^2/2$ in the small-angle range used here.

Prepare the second probe as $(a_d,a_d,c_d)$ with $c_d$ within angle $\delta$
of $-a_d$. Its chirality is exactly zero and its $D$ gain is
$u=(1-a_d\cdot c_d)/2\ge1-\delta^2/4$. Thus the reference coordinate map is
$P=\mathrm{diag}(t_x,u)$ and $E=M(P-I)$ obeys

```math
\|E\|\le\mu(3\delta^2+\delta^3).
```

At $\delta=0.005$ radians, about $0.286$ degrees, this is less than
0.000080045305, within the robust profile's $e=0.000081$. It requires
permutation of the same actual triad or an equivalent site-matched
preparation. Cycling programmed axis labels under unrelated site errors
does not guarantee cancellation. Purity defects are not covered by this
quadratic angular allowance.

For arbitrary spin-reference trace distances $\eta_x,\eta_d$, the general
bound is instead

```math
\boxed{\|E\|\le8h\|M\|\sqrt{4\eta_x^2+\eta_d^2}.}
```

Indeed $\|\delta x\|\le2h\eta$ by the exact norm of every linear
combination of $\chi,D/4$, and the response-column multipliers are 8 and 4.
This spin-space argument does not cover wrong particle number or orbital
loading. Those require a full-preparation error model.

## 5. Product-reference geometry and scope of the allocation

Every separable three-spin reference obeys

```math
64\langle\chi\rangle^2+16\langle D/4\rangle^2\le1.
```

For a pure product, align the middle Bloch vector with $z$ and write
$d=\langle D\rangle=(a_z-c_z)/2$. The transverse determinant satisfies
$|8\langle\chi\rangle|\le\sqrt{(1-a_z^2)(1-c_z^2)}\le1-d^2$.
Therefore its square plus $d^2$ is at most one. Convexity extends the
quadratic bound to all separable mixtures.

In a homoscedastic linear-calibration model let $G$ be the average reference
coordinate outer product. Then $64G_{11}+16G_{22}\le1$; hence
$\lambda_{\min}(G)\le1/80$. The two axis probes with probabilities $4/5,1/5$
attain $G=I/80$. This explains the $4:1$ allocation for the common scalar
variance bound. It is not a global statistical optimum for heterogeneous
reference variances or arbitrary calibration procedures.

The risk guarantees concern each fixed target's unconditional MSE, not
simultaneous confidence coverage of a scan. Reused estimates are correlated
through the common calibration. Independent trials, known preparation
accuracy, and a justified stability interval remain essential. The billions
of preparations are an explicit finite existence construction; they are not
an experimental economy claim.

## 6. Reproduction

After producing the exact parity and parameter-robustness certificates, run

```bash
python research/check_finite_calibration.py \
  --certificate /absolute/output/complete-parity-certificate.json \
  --robustness-certificate /absolute/output/parameter-robustness.json \
  --output /absolute/output/finite-calibration.json
```

The standard-library checker uses rational arithmetic, rational enclosures
for square roots and $\pi$, verifies both certificate dependencies, and
records their hashes. Generated reports must remain outside the repository.
The root verifier runs these dependencies in sequence.
