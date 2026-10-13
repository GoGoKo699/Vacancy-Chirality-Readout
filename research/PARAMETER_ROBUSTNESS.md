# Robustness of the primary complete-parity score

The [primary score certificate](CERTIFICATE.md) extends to a finite parameter
region with exact response recalibration. A separate fixed-score statement
allows general symmetry-breaking measurement errors at finite sample size.
These statements concern the same four-mode, known-three-particle input and
ideal total-site parity detector as [the complete-parity theory](COMPLETE_PARITY_THEORY.md).
They do not combine automatically with the loss/bit-error window or the
[finite-reference calibration budget](FINITE_CALIBRATION.md).

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

$$
\begin{array}{c|c|c}
\max(|u-16|,|\tau-\pi|,|\phi-\phi_0|)&\Delta\text{ upper bound}&
\|T_{\rm corrected}\|\text{ upper bound}\\ \hline
0.00028&0.00248976&0.394824826\\
0.00016&0.00142272&0.366080282
\end{array}
$$

These are sufficient dimensionless parameter regions from a conservative
operator bound. They are not measured apparatus tolerances. The older
mixed-input-score neighborhoods in the nominal certificate remain distinct.

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

## 5. Reproduction and evidence scope

[`check_parameter_robustness.py`](check_parameter_robustness.py) consumes a
passing exact parity certificate, checks the four nominal input bounds in
(1), records that report's SHA-256, and runs three exact-arithmetic groups:
response recalibration, physical parameter budgets, and fixed-score risk.
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
