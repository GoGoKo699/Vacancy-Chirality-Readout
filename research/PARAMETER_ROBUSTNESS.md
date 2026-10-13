# Robustness of the primary complete-parity score

The [primary score certificate](CERTIFICATE.md) extends to a finite parameter
region with exact response recalibration. A separate fixed-score statement
allows general symmetry-breaking measurement errors at finite sample size.
These statements concern the same four-mode, known-three-particle input and
ideal total-site parity detector as [the complete-parity theory](COMPLETE_PARITY_THEORY.md).
Section 5 gives a separate joint parameter/loss/bit-error guarantee and
certifies the inputs needed to extend the
[finite-reference calibration budget](FINITE_CALIBRATION.md) to its smaller
joint window. Section 6 gives microscopic sufficient conditions and a separate
loading-distance allowance. The larger ideal-detector regions are not allowances that can
be added to that joint window.

The exact nominal point is $u=16$, $\tau=\pi$, and
$\Phi_+=-2\pi/3$, with $\Phi_-=-\Phi_+$. The two settings are randomized
independently with equal probabilities. Let $X_a$ be the exactly corrected
primary score and $X_D$ the exactly corrected auxiliary score for
$D=P_{12}-P_{23}$. Their nominal first moments are $\chi$ and $D$, and the
integer/rational certificate supplies

$$
\|T_a\|<v=0.3287,\quad |X_a|<A=1.623,\qquad
\|T_D\|<v_D=646,\quad |X_D|<A_D=194.
\tag{1}
$$

Neither score is asserted to be globally optimal. Their exact definitions,
rather than the displayed rounded score decimals, determine the identities.

## 1. A weighted dilation bound

Let $\mathcal V$ be an isometry for the nominal measurement, including the
equal-probability flux register, and let $S$ be its fixed Hermitian score
operator. Thus the first moment is $\mathcal V^\dagger S\mathcal V$ and
the second moment is $\mathcal V^\dagger S^2\mathcal V$. Suppose

$$
\|S\|\le A,\qquad \|S\mathcal V\|\le\sqrt v,
\qquad \|\mathcal V'-\mathcal V\|\le\eta.
$$

The two isometries must use a common output space and the same score
operator; this is an explicit dilation-distance assumption. Write
$E=\mathcal V'-\mathcal V$. Expanding the first moment gives

$$
\mathcal V'^\dagger S\mathcal V'-\mathcal V^\dagger S\mathcal V
=E^\dagger S\mathcal V+\mathcal V^\dagger SE+E^\dagger SE.
$$

The first two terms each have norm at most $\eta\sqrt v$ and the third
at most $A\eta^2$. The triangle inequality also gives
$\|S\mathcal V'\|\le\sqrt v+A\eta$. Consequently

$$
\boxed{\epsilon_a(\eta)=2\sqrt v\,\eta+A\eta^2,\qquad
Q_a(\eta)=(\sqrt v+A\eta)^2,}
\tag{2}
$$

where the first quantity bounds the first-moment operator error and the
second bounds the perturbed second-moment norm. Substituting $(v_D,A_D)$
gives $\epsilon_D,Q_D$. The estimate uses the nominal score's second
moment as well as its largest absolute value.

## 2. Exact response recalibration

Assume the perturbed measurement remains spin independent and that its two
settings are exact spatial reflections. As in the complete-parity theory,
the odd first moments then have the form

$$
A_a'=(1+x)\chi+yD,\qquad A_D'=z\chi+(1+w)D.
$$

For $h=\sqrt3/4$, the identities
$\chi^2=h^2\Pi_{\rm m}$, $D^2=16h^2\Pi_{\rm m}$, and
$\{\chi,D\}=0$ imply
$\|r\chi+sD\|=h\sqrt{r^2+16s^2}$. Equation (2) therefore implies

$$
|x|\le X=\epsilon_a/h,\quad |y|\le Y=\epsilon_a/(4h),\quad
|z|\le Z=\epsilon_D/h,\quad |w|\le W=\epsilon_D/(4h).
$$

If $W<1$ and

$$
\kappa=1-X-\frac{YZ}{1-W}>0,
$$

the response matrix is invertible. Using its exact true coefficients, report

$$
X_{\rm corrected}=
\frac{X_a-\dfrac{y}{1+w}X_D}
{1+x-\dfrac{yz}{1+w}}.
\tag{3}
$$

Its first moment is exactly $\chi$. For every input state, the classical
$L^2$ triangle inequality bounds the second moment of its numerator;
the denominator has absolute value at least $\kappa$. Thus

$$
\boxed{\|T_{\rm corrected}\|\le
\left[\frac{\sqrt{Q_a}+\dfrac{Y}{1-W}\sqrt{Q_D}}
{\kappa}\right]^2.}
\tag{4}
$$

The [exact-arithmetic checker](check_parameter_robustness.py) encloses all
square roots and verifies these consequences of (1)–(4):

| Dilation distance $\eta$ | Corrected second-moment bound | Passive comparison |
|---:|---:|---|
| $0.0025$ | $<0.394824826$ | Below the unbiased coefficient $13/32$ |
| $0.0015$ | $<0.366080282<0.367$ | Also separates finite all-estimator minimax risk for every $n\ge12000$ |

The finite comparison uses the [passive bound](PASSIVE_PARITY_BOUND.md),
with $I(1/5)=775/297$ and prior information $400\pi^2/3<4000/3$.
Since $0.367(775/297)<1$, its sufficient sample threshold is
$0.367(4000/3)/[1-0.367(775/297)]<12000$.
These are worst-case comparisons over arbitrary input states, not
pointwise dominance for each state.

Exact response coefficients in (3) are a mathematical calibration
assumption. The response is fixed during its characterization and the target
trials; alternatively, its exact coefficients must be known for each trial.
Substitution of finite empirical estimates requires its own
random-inverse or other finite-calibration risk analysis.

## 3. A physical parameter region

Use dimensionless $H_{u,\Phi}=T_\Phi+uW$, where
$W=\tfrac12\sum_jn_j(n_j-1)$. Let $J$ embed the input at occupations
$(0,1,1,1)$; this $J$ is an embedding, not the hopping-energy unit.
The two settings use the same $u,\tau$ and opposite fluxes. Write
$\phi=\Phi_+$ and $\phi_0=-2\pi/3$.

There is no double occupancy initially, so $WJ=0$. The two-step hopping
paths returning to the occupied-input subspace give

$$
J^\dagger H^2J=J^\dagger T_\phi^2J
=6I+2(P_{12}+P_{23}),\qquad \|HJ\|=\sqrt{10}.
\tag{5}
$$

The six directed occupied-input hops contribute the identity terms; the
two occupied nearest-neighbor bonds each supply the two exchange paths.
Forward/backward phases cancel. The sum of the two swaps has norm at
most two and equals two on the symmetric spin subspace, establishing the
last equality. For a duration change at fixed Hamiltonian, integrating
the derivative of $e^{-itH}J$ therefore gives

$$
\|[e^{-i(\tau+\delta\tau)H}-e^{-i\tau H}]J\|
\le\sqrt{10}|\delta\tau|.
\tag{6}
$$

For $u>0$, energy conservation on the same embedded input gives

$$
\|W e^{-isH_{u,\phi}}J\|
\le\frac{\sqrt{10}+\|T_\phi\|}{u}.
$$

Duhamel's formula for the derivative in $u$, integrated along the
repulsion path at duration $\pi$, consequently yields

$$
\|[e^{-i\pi H_{u,\phi}}-e^{-i\pi H_{16,\phi}}]J\|
\le\pi(\sqrt{10}+\|T_\phi\|)|\log(u/16)|.
\tag{7}
$$

The four-mode single-particle hopping eigenvalues are
$-2\cos(k+\phi/4)$, $k\in\{0,\pi/2,\pi,3\pi/2\}$.
On the three-boson space their sums give

$$
\|T_\phi\|=6\max(|\cos(\phi/4)|,|\sin(\phi/4)|),\quad
\|\partial_\phi T_\phi\|=
\frac32\max(|\cos(\phi/4)|,|\sin(\phi/4)|).
$$

When $|\phi-\phi_0|\le0.001$, the maximum is at most
$\sqrt3/2+1/8000$: concavity bounds the cosine, while the sine is at
most $1/2+1/4000$. Hence the two operator norms are respectively below
$5.197$ and $1.300$. Change flux first, then repulsion, then duration.
For $|u-16|\le0.002$, the logarithm in (7) is at most
$|u-16|/15.998$. Equations (6)–(7), the flux derivative bound, and
$\pi<22/7$, $\sqrt{10}<3.163$ prove

$$
\boxed{\|\mathcal V'-\mathcal V\|\le\Delta
=1.643|u-16|+3.163|\tau-\pi|
+4.086|\phi+2\pi/3|.}
\tag{8}
$$

The ideal parity measurement is the same in both dilations, and appending
its outcome register preserves the bound. The two randomized flux branches
obey the same estimate by reflection. For every point satisfying
$\Delta\le0.0025$, the prerequisites $|u-16|<0.002$ and
$|\phi-\phi_0|<0.001$ already follow from (8)'s positive coefficients.

Thus $\Delta\le0.0025$ gives the first row of the table, and
$\Delta\le0.0015$ gives the second. Sufficient equal-coordinate boxes are

```math
\begin{array}{c|c|c}
\max(|u-16|,|\tau-\pi|,|\phi-\phi_0|)&\Delta\text{ upper bound}&
\|T_{\rm corrected}\|\text{ upper bound}\\ \hline
0.00028&0.00248976&0.394824826\\
0.00016&0.00142272&0.366080282
\end{array}
```

These are sufficient dimensionless parameter regions from a conservative
operator bound. They are not measured apparatus tolerances. The mixed-input-score neighborhoods in the nominal certificate remain distinct.

## 4. Fixed-score errors without the symmetry assumption

Keep the nominal primary score instead of recalibrating it. The weighted
dilation lemma still applies if the actual instrument breaks spin or
reflection symmetry. Assume the same unknown input state is freshly
prepared on each trial and every actual instrument obeys the stated
dilation-distance bound. For a common input embedding and parity detector,
a sufficient Hamiltonian condition at fixed duration is
$\int\|H'(t)-H(t)\|\,dt\le\eta$, by Duhamel's formula. Define

$$
\beta_\eta=2\sqrt{0.3287}\,\eta+1.623\eta^2,\qquad
v_\eta=(\sqrt{0.3287}+1.623\eta)^2.
$$

For independent trials, including a deterministic variation of the
instrument between trials, the sample mean obeys

$$
\boxed{\sup_\rho\mathbb E_\rho
(\widehat\chi-\langle\chi\rangle_\rho)^2
\le\frac{v_\eta}{n}+\beta_\eta^2.}
\tag{9}
$$

Indeed each score variance is at most $v_\eta$, the variances add, and
the average bias has absolute value at most $\beta_\eta$. At
$\eta=10^{-4}$, one has $v_\eta<0.328887$ and
$\beta_\eta<0.000114682$. For $n=100000$ this gives

$$
R_{\rm fixed}<3.302022\times10^{-6},\qquad
R_{\rm passive}>3.903384\times10^{-6}.
\tag{10}
$$

The passive lower bound uses $I(1/10)=97600/38907$ and prior information
$1600\pi^2/3<16000/3$. An independent-trial separation also holds for
every integer $4000\le n\le4000000$. To check this whole interval, put
$b=0.000114682$, $v_*=0.328887$, $I_*=775/297$, and $J_*=4000/3$.
The desired strict inequality follows from negativity of

$$
q(n)=b^2I_*n^2+(b^2J_*+v_*I_*-1)n+v_*J_*.
$$

This quadratic is convex. The exact checker proves $q(4000)<-128$ and
$q(4000000)<-17558$, so convexity makes it negative throughout the interval.
This statement is a finite biased-estimator MSE comparison; its nonzero
bias floor does not establish an asymptotic unbiased separation.

If the instrument instead depends on the preceding classical history,
require the dilation bound conditional on **every** history. Write the
centered score differences as martingale differences. Their average has
$L^2$ norm at most $\sqrt{v_\eta/n}$, while the average conditional bias
has absolute value at most $\beta_\eta$. The $L^2$ triangle inequality gives

$$
R_{\rm history}\le
\left(\sqrt{v_\eta/n}+\beta_\eta\right)^2.
\tag{11}
$$

At the same $\eta=10^{-4}$ and $n=100000$, this is below
$3.717980\times10^{-6}$ and still below the passive lower bound in (10).
The independent-trial interval above is not asserted for (11).

## 5. Joint parameter variation, particle loss, and detector errors

Retain the fixed nominal score pair and suppose the perturbed **clean**
measurement obeys the dilation bound $\eta$. During its actual evolution,
allow equal, spin-independent Markov losses
$L_{j\sigma}=\sqrt\gamma\,b_{j\sigma}$, with no gains, and let
$\lambda=\gamma t$ use the **actual** duration after its perturbation.
After evolution, flip each of the four parity bits independently with the
same probability $q$. The rates and detector law are identical for the two
reflection-paired settings, and the Hamiltonian remains number conserving,
spin independent, and paired by reflection. These assumptions preserve the
two-dimensional odd response used in Section 2.

The [noise proof](NOISE_BUDGET.md) applies to the perturbed Hamiltonian itself.
The traces in its number sectors obey the same pure-death equations, so the
number probabilities remain
$p_k={3\choose k}e^{-k\lambda}(1-e^{-\lambda})^{3-k}$.
In particular, its three-particle output is exactly
$sU'\rho U'^\dagger$, where $s=e^{-3\lambda}$: the no-jump imaginary
Hamiltonian is the scalar $-3i\gamma I/2$ in that sector. Every contamination
bound in the noise proof is uniform over the conditional state after losses.
There is therefore no approximation that moves the losses to the end of
the interacting evolution, and no additional closeness assumption on its
one- or two-particle dynamics.

Set

$$
z=e^\lambda-1,\quad \alpha=(1-q)^4,\quad r_q=\frac q{1-q},\qquad
R_b=6r_q^2+r_q^4,\quad R_c=4r_q+4r_q^3.
$$

For $i=a,D$, let $\epsilon_i,Q_i$ be the clean perturbation bounds in (2)
and use the nominal singleton half-ranges $C_a<1.123$, $C_D<107$ from
the exact score certificate. Analyze the normalized features
$Y_i=X_i/(s\alpha)$. Their operator mean errors relative to $\chi,D$,
respectively, and their all-input second moments satisfy

$$
e_i=\epsilon_i+3C_i z^2+A_iR_b(1+3z^2)+3A_iR_cz,
\tag{12}
$$

$$
N_i=\frac{e^{3\lambda}}\alpha
\left[Q_i+R_bA_i^2+3z^2(1+R_b)A_i^2
+(3z+z^3)R_cA_i^2\right].
\tag{13}
$$

Indeed the noise mean error is measured relative to the perturbed clean
first moment; the triangle inequality adds $\epsilon_i$ to it. The clean
three-particle second moment is at most $Q_i$, which replaces the nominal
one in the noise bound. All remaining terms depend only on the implemented
score, the number probabilities, and the detector kernel. This proves the
composition without dropping the loss/detector cross term.

Apply (3)–(4) to the normalized response, substituting
$(e_a,e_D,N_a,N_D)$ for $(\epsilon_a,\epsilon_D,Q_a,Q_D)$.
Every expression in (12)–(13) is nondecreasing in $\eta,\lambda,q$ for
nonnegative rates and $q<1$. The inverse bound is also nondecreasing while
$W<1$ and $\kappa>0$, as verified at the endpoints below. Exact arithmetic
therefore gives uniform guarantees over these joint boxes:

| Clean dilation distance $\eta$ | Actual loss exposure $\lambda$ | Bit-error probability $q$ | Exactly recalibrated second moment |
|---:|---:|---:|---:|
| $\le0.001$ | $\le0.01$ | $\le0.003$ | $<0.395873255<13/32$ |
| $\le0.001$ | $\le0.005$ | $\le0.0015$ | $<0.366684486<0.367$ |

Both rows permit the physical budget $\Delta\le0.001$ in (8), with
matched opposite fluxes and common $u,\tau$. An equal-coordinate radius
$0.00011$ gives $\Delta\le0.00097812<0.001$. The smaller joint window
also separates finite all-estimator passive minimax risk for every
$n\ge12000$, using the same argument as Section 2.

The factor $s\alpha$ is an analytical normalization. Exact inversion of the
**raw** two-feature response absorbs it, so the separately known values of
$\lambda,q$ are not operational inputs. The coefficient statements above
assume exact response coefficients and a fixed response, or exact
coefficients known separately on each trial. Finite references are accounted
for by the following separate argument.

### Finite-reference inputs in the smaller joint window

For raw features $Z=(X_a,X_D/4)^T$ and target coordinates
$x=(\langle\chi\rangle,\langle D\rangle/4)^T$, write
$\mathbb EZ=Mx$ and $M=r_0(I+E)$, with $r_0=s\alpha$.
The mean-error norm identity in Section 2 bounds the two rows of $E$ by

$$
\ell_1=e_a/h,\qquad \ell_2=e_D/(4h),\qquad
\|E\|\le c=\sqrt{\ell_1^2+\ell_2^2}.
$$

Let $r_{\min}=e^{-3\lambda_{\max}}(1-q_{\max})^4$.
Throughout the smaller joint box $r_0\in[r_{\min},1]$ and $c,\ell_1<1$.
For $B=I-M$, triangle inequalities give

$$
\|B\|\le b=1-r_{\min}+r_{\min}c,\qquad
\|e_1^TB\|\le r=1-r_{\min}+r_{\min}\ell_1,
\tag{14}
$$

$$
\|M\|\le\mu=1+c,\qquad
\sup_\rho\mathbb E_\rho\|Z\|^2\le H=N_a+N_D/16.
\tag{15}
$$

The first two bounds use the maximum of $1-r_0+r_0c$ and
$1-r_0+r_0\ell_1$ over the interval for $r_0$; the raw second moments
are at most the normalized ones because $r_0\le1$.

Define the deterministic polynomial rows
$s_j=e_1^T\sum_{k=0}^{j-1}B^k$ as in the finite-calibration theorem.
Their first two second-moment bounds are

$$
V_1=N_a,\qquad
V_2=\left[(1+r)\sqrt{N_a}+\ell_1\sqrt{N_D/16}\right]^2.
\tag{16}
$$

For $s_2=e_1^T(I+B)$, the first entry has absolute value at most $1+r$.
Its second entry obeys $|B_{12}|=r_0|E_{12}|\le\ell_1$; applying the
classical $L^2$ triangle inequality proves (16) with the feature covariance
retained. The exact inverse row $w_*=e_1^TM^{-1}$ has second moment at most
the joint bound $T_{\rm joint}$. Since
$s_j-w_*=-e_1^TB^j(I-B)^{-1}$ and $b<1$,

$$
\|s_j-w_*\|\le\frac{rb^{j-1}}{1-b},\qquad
V_j=\left[\sqrt{T_{\rm joint}}
+\sqrt H\frac{rb^{j-1}}{1-b}\right]^2\quad(j=3,4)
\tag{17}
$$

bounds the remaining moments. The inverse row is used only to prove these
deterministic bounds; the finite-reference estimator uses independent
polynomial calibration and does not compute a random inverse.

| Quantity | Outward bound in the smaller joint window |
|---|---:|
| $b$ | $0.065701062837834$ |
| $r$ | $0.023954914327520$ |
| $\mu$ | $1.045869024670356$ |
| $H$ | $42.632586752575791$ |
| $V_1$ | $0.338062671759970$ |
| $V_2$ | $0.379960845522265$ |
| $V_3$ | $0.367560194372654$ |
| $V_4$ | $0.366741988495767$ |

Every value is strictly below its corresponding nominal-window bound in
[FINITE_CALIBRATION.md](FINITE_CALIBRATION.md). The checker verifies this
domination against the supplied exact nominal certificate. Consequently both
of that theorem's finite-copy profiles hold in this smaller joint window,
with their original independent reference stages, predetermined target tasks,
and reference-preparation and drift assumptions. The base response is fixed
across reference and target sampling; its separately bounded preparation
errors and target drift are treated by that theorem. The larger joint window
above has no such finite-reference count claim here.

## 6. Microscopic control and loading errors

The fixed-score result can be stated directly in terms of microscopic
perturbations, without spin or reflection symmetry of the actual instrument.
Keep the closed four-orbital, two-internal-level model and duration $\pi$.
For either nominal flux branch, write

$$
H'(s)=H_0+d\Gamma(\delta h(s))+R(s),\qquad H_0=T+16W.
$$

Here $\delta h(s)$ is a Hermitian one-particle perturbation, including
spin-dependent hopping or local fields. The second-quantized operator has
norm at most $3\|\delta h(s)\|$ in the three-boson sector.
Let $R(s)$ be a Hermitian sum of on-site pair perturbations, preserving each
spatial occupation block, with each two-spin pair operator bounded in norm
by $\varepsilon(s)$. A block contains exactly $W$ on-site pairs, so the
triangle inequality gives $\|R(s)\|_{\rm block}\le\varepsilon(s)W$.
Since $W$ is scalar on each such block,
$R(s)^2\le\varepsilon(s)^2W^2$. This assumption excludes perturbations
that move a pair between occupation blocks.

Equation (5), energy conservation, and $\|T\|=3\sqrt3$ at the nominal
flux give

$$
\|W e^{-isH_0}J\|\le\frac{\sqrt{10}+3\sqrt3}{16}.
$$

In Duhamel's formula, the perturbation acts on this nominal propagation
and the remaining perturbed propagator is unitary. Consequently

$$
\eta\le\int_0^\pi
\left[3\|\delta h(s)\|+
\frac{\sqrt{10}+3\sqrt3}{16}\varepsilon(s)\right]ds.
\tag{18}
$$

For uniform bounds $\|\delta h(s)\|\le d$ and
$\varepsilon(s)\le e$, this gives the sufficient budget

$$
\boxed{\eta\le9.429d+1.643e.}
\tag{19}
$$

Both flux branches must obey the bound. Equations (18)–(19) feed the
fixed-nominal-score risks (9) and (11); they do not justify two-feature
recalibration when spin or reflection symmetry is broken. They also do not
bound dynamical escape into additional orbitals or internal levels.

Independently, a preparation error can be quantified without discarding
any trial. Suppose the actual input $\sigma$ obeys

$$
\tfrac12\|\sigma-J\rho J^\dagger\|_1\le\delta
\tag{20}
$$

uniformly over the intended unknown spin state $\rho$. The common input
space may include additional particle-number sectors; the actual measurement
channel must act on that same space. If a spectator reference is included,
use the embedding tensored with its identity and the corresponding joint
trace-distance bound. Extend the fixed score to all sixteen parity patterns,
with zero score on even patterns, so $|X_a|\le A=1.623$ on every trial.
Contractivity of trace distance under the actual measurement makes its
outcome distributions differ in total variation by at most $\delta$.
Their score means therefore differ by at most $2A\delta$ and their second
moments by at most $A^2\delta$, using the respective ranges
$[-A,A]$ and $[0,A^2]$.

Combining this with the dilation bound on the intended input gives, for
independent preparations and trials,

$$
\boxed{R_n\le\frac{v_\eta+A^2\delta}{n}
+(\beta_\eta+2A\delta)^2.}
\tag{21}
$$

For $\eta=\delta=10^{-4}$ and $n=100000$, exact arithmetic verifies

$$
R_n<3.484473\times10^{-6}
<3.903384\times10^{-6}<R_{\rm passive}.
\tag{22}
$$

This is a worst-case finite-MSE comparison with the same passive lower bound
as (10). Under history dependence, require both (20) and the dilation bound
conditional on every preceding history. The martingale argument instead gives

$$
R_n\le\left[\sqrt{\frac{v_\eta+A^2\delta}{n}}
+\beta_\eta+2A\delta\right]^2.
\tag{23}
$$

The numerical separation in (22) is asserted only for independent trials.
The metric assumption (20) is not implied by a loading success probability
alone: coherent leakage can have trace distance of order the square root
of its probability. Every preparation is counted and every outcome scored;
there is no heralding or postselection in (21)–(23).

## 7. Reproduction and evidence scope

[`check_parameter_robustness.py`](check_parameter_robustness.py) consumes a
passing exact parity certificate, checks the nominal input bounds in
(1) and both singleton half-ranges, records that report's SHA-256, and runs
five exact-arithmetic groups: response recalibration, physical parameter
budgets, fixed-score risk, and joint parameter/noise bounds with
finite-calibration inputs, and microscopic-control/loading bounds.
For example, from the repository root, choose a new output path outside it:

```bash
python -B research/check_parameter_robustness.py \
  --certificate verification/complete-parity-certificate.json \
  --output /tmp/vacancy-parameter-robustness.json
```

The operator identities, perturbation inequalities, and statistical
arguments are proved above. The checker verifies their numerical
consequences with rational enclosures; arithmetic alone does not prove
those analytical premises or validate an experimental error bound. Run
the [full repository verifier](../REPRODUCE.md) for the integrated check.
