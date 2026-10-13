# Calibration cost and reuse

Complete-image readout has a certified advantage in unknown-input sampling
under exact calibration. Finite reference preparation and validation have
their own cost. This account keeps those resources explicit and distinguishes
finite confidence certificates from asymptotic planning calculations.

For a finite **unconditional** construction within the certified noise
window, see [FINITE_CALIBRATION.md](FINITE_CALIBRATION.md). It uses independent
reference stages and a bounded polynomial correction, includes the reference
preparations in a passive-comparison cost theorem, and treats bounded
preparation error and drift. The present account describes the distinct
held-out-validation and regular inverse-response methods; its numerical
planning profiles use the mixed-input score.

The implemented-channel assumptions, product probes, and empirical Bernstein
certificate are in [CHANNEL_CALIBRATION.md](CHANNEL_CALIBRATION.md). Its
variance-sensitive radii use Maurer and Pontil (2009), Theorem 4
([MP2009](../REFERENCES.md)); they include the finite-sample terms proportional
to $B/(N-1)$ and $B^2/(N-1)$ for a score bounded by $B$.
Preparation error and drift must be included in the final bounds used below.
For the specified loss and bit-error window, the exactly recalibrated mixed-input score
has $V<0.446$ and absolute score below $2.626$
([certificate](CERTIFICATE.md)). These are useful planning bounds for that
score. A score obtained from a random pilot inverse must have its own range
checked before validation; the exact-score range does not automatically bound
an unstable fitted inverse.

## 1. Count the actual preparations

Let $P$ count pilot reference preparations and let $N_x,N_d,N_0$ count the
three held-out validation batches. The reference cost is

$$
M=P+N_x+N_d+N_0.
$$

Add any optional fourth-probe batch and all preparations from earlier failed
attempts. Counting a failed image as a measurement outcome does not remove
its preparation from this budget. The two statistics from one image, its
score and squared score, share that preparation and are not charged twice.

After freezing and validating the score, write its final bias allowance as
$\beta$ and its final second-moment bound as $U$. These include justified
preparation and drift allowances. On the validation confidence event,
independent unknown preparations obey

$$
\mathrm{MSE}\le U/n+\beta^2.
$$

For a desired MSE $\varepsilon^2$ and $\beta<\varepsilon$, the sufficient
integer target-copy budget is

$$
\boxed{n=\left\lceil\frac{U}{\varepsilon^2-\beta^2}\right\rceil.}
$$

If $\beta\ge\varepsilon$, collecting more target copies does not make this
certificate meet the requested accuracy. Improve the calibration or its
error allowances before using this budget.

The resource report is therefore $(n,M)$: unknown-state preparations and
known product-reference preparations. Their costs can differ substantially.
If one reference preparation costs $w$ times one target preparation and the
same valid calibration is used for $K$ target batches, the total weighted
preparation cost is

$$
K n+wM.
$$

Each batch is a separate estimation task on an arbitrary unknown state,
with the same per-batch accuracy requirement and no promised relation
between states. This is not a calculation for pooling all $Kn$ copies of
one common state into a single estimate.

Equal preparation costs correspond to $w=1$. Reuse is justified only while
the calibrated channel and reference preparation remain stable. Shared
calibration bias and drift do not decrease as $1/K$; only their preparation
cost is amortized.

## 2. Two distinct separable comparisons

The ideal uniformly randomized six-Pauli estimator is unbiased and has
worst-case MSE $9/(16n)$. Its sufficient and attained worst-case budget is

$$
n_0=\left\lceil\frac9{16\varepsilon^2}\right\rceil.
$$

For a fixed accepted calibration record, the exact preparation-count
comparison is $Kn+wM<Kn_0$. If $\Delta_0=n_0-n>0$, its smallest strict
integer reuse threshold is

$$
\boxed{K_{\min,0}=
\left\lfloor\frac{wM}{\Delta_0}\right\rfloor+1.}
$$

If $\Delta_0\le0$, that record has no strict advantage under this accounting,
regardless of reuse. This calculation compares the native certificate with
an ideal Pauli implementation; it does not assign a zero physical calibration
cost to an actual Pauli apparatus.

The [separable benchmark](SEPARABLE_BENCHMARK.md) also gives a distinct
necessary target-copy count for arbitrary estimators, including biased ones,
whose worst-case **unconditional** MSE is at most $\varepsilon^2$:

$$
N_{\min}=\left\lceil
\left(\frac3{4\varepsilon}-\pi\right)^2
\right\rceil,
\qquad 0<\varepsilon<\frac1{16\pi}.
$$

The measurement restriction remains fully separable effects within each
three-spin copy, with classical adaptation allowed. This lower bound need
not be attained at finite sample size. Replacing $n_0$ by $N_{\min}$ in
the preceding integer comparison gives the stricter count condition

$$
K_{\min,\mathrm{lb}}=
\left\lfloor\frac{wM}{N_{\min}-n}\right\rfloor+1,
\qquad N_{\min}>n.
$$

Such a count inequality alone does not establish an unconditional
finite-calibration advantage: the validation confidence and rejection
probabilities must also be included in the risk comparison.

## 3. Confidence, rejection, and unconditional risk

Let $G$ be the joint validation confidence event, with
$\Pr(G)\ge1-\delta$. The MSE certificate holds for every fixed pilot and
validation record in $G$. Thus the probability of issuing an incorrect
certificate is at most $\delta$.

This is not coverage $1-\delta$ conditional on passing an acceptance gate
$A$. In general,

$$
\Pr(G^c\mid A)\le\frac{\delta}{\Pr(A)},
$$

which can be much larger when passing is rare. Rejecting a failed validation
and returning no estimate is abstention; it is not an unconditional estimator
with the displayed MSE. A nominal or simulated profile substituted for
observed statistics does not establish that validation will pass.

For example, define a complete procedure that returns zero on rejection and
clips every accepted target estimate to the physical chirality interval
$[-\sqrt3/4,\sqrt3/4]$. If its accepted good-event MSE is at most
$\varepsilon^2$, an unconditional bound is

$$
\mathrm{MSE}_{\rm total}
\le\varepsilon^2+\frac34\delta
+\frac3{16}\Pr(\text{reject}).
$$

Clipping does not increase squared error for a physical chirality value.
The additional terms must be budgeted or improved before comparing this
procedure with an unconditional separable lower bound. In particular, a
fixed $95\%$ calibration confidence does not by itself establish an
arbitrarily small unconditional MSE.

## 4. Leading empirical Bernstein planning

Let $v_x,v_d$ be the variances of the fixed signed score on the two mean
reference states and set $L=\log(20/\delta)$. Ignoring only the lower-order
finite-$N$ terms for planning, the statistical squared bias radius is

$$
r_\beta^2\simeq
\frac{24v_xL}{N_x}+\frac{6v_dL}{N_d}.
$$

For $M_{\rm mean}=N_x+N_d$, this expression is minimized by

$$
N_x:N_d=2\sqrt{v_x}:\sqrt{v_d},\qquad
r_\beta^2\simeq\frac{C_E}{M_{\rm mean}},\qquad
C_E=6L(2\sqrt{v_x}+\sqrt{v_d})^2.
$$

This is an allocation for a confidence radius. Observed residuals, the
empirical Bernstein finite-$N$ terms, the second-moment certificate, pilot
cost, preparation error, and drift remain part of the actual calculation.
In particular, on the confidence event the reported bias bound can be as
large as $\beta_{\rm true}+2r_\beta$, where $\beta_{\rm true}$ is the
fixed score's actual bias norm. A prospective passing guarantee with a known
margin requires

$$
r_\beta\le(\beta_{\rm target}-\beta_{\rm true})/2.
$$

Similarly, the reported second-moment upper bound can exceed the true value
by twice its linear confidence radius. Assuming zero observed residuals
avoids neither fluctuation. The following reuse calculation is therefore a
leading planning model, not a finite-sample acceptance theorem.

Suppose the reference budget contributing to the bias allowance is modeled
as $C/\beta^2$, and set $\beta^2=b\varepsilon^2$. With variance coefficient
$V$, the leading per-batch weighted cost is

$$
\frac1{\varepsilon^2}
\left[\frac{V}{1-b}+\frac{wC}{Kb}\right].
$$

Writing $D=wC/K$, direct minimization gives

$$
\boxed{b_{\ast}=
\frac{\sqrt D}{\sqrt V+\sqrt D},\qquad
\mathrm{cost}_{\min}\simeq
\frac{(\sqrt V+\sqrt D)^2}{\varepsilon^2}.}
$$

Consequently its strict leading break-even condition relative to $9/16$ is

$$
\boxed{K>\frac{wC}{(3/4-\sqrt V)^2},\qquad V<9/16.}
$$

Use $C=C_E$ only with the confidence-planning qualifications above.
Reference-only overhead not included in $C/\beta^2$ must be added. With a
nonzero systematic allowance, the bias cannot be allocated entirely to
finite reference statistics, and this zero-floor formula no longer applies.

## 5. A different regular plug-in asymptotic calculation

A regular inverse-response estimator has a different leading risk account
from a simultaneous confidence bound. Freeze a two-component feature $Z$
independently of the subsequent data, write its response as
$\mathbb E_\rho Z=M(\theta,\zeta)^T$, and set
$w^T=e_1^TM^{-1}$, $X=w^TZ$, where
$\theta=\langle\chi\rangle_\rho$ and $\zeta=\langle D\rangle_\rho$.
The two reference means estimate the columns of $M$ as
$\widehat M=[8\overline Z_x,\overline Z_d]$.
For the plug-in estimate $e_1^T\widehat M^{-1}\overline Z$, differentiating
the inverse gives the leading fluctuation

$$
w^T(\overline Z-\mathbb E_\rho Z)
-8\theta w^T(\overline Z_x-\mathbb E_xZ)
-\zeta w^T(\overline Z_d-\mathbb E_dZ).
$$

Thus independent reference mean errors give the first-order calibration
contribution on a target state $\rho$,

$$
64\langle\chi\rangle_\rho^2\frac{v_x}{N_x}
+\langle D\rangle_\rho^2\frac{v_d}{N_d}.
$$

The anticommuting observables satisfy
$16\langle\chi\rangle_\rho^2+\langle D\rangle_\rho^2\le3$.
The worst value of this contribution is therefore

$$
\max\left\{\frac{12v_x}{N_x},\frac{3v_d}{N_d}\right\}.
$$

For the actual two-probe total $M_{\rm mean}=N_x+N_d$, equalizing the two
terms gives

$$
N_x:N_d=4v_x:v_d,\qquad
C_R=12v_x+3v_d,\qquad
\mathrm{MSE}\le\frac Vn+\frac{C_R}{M_{\rm mean}}
+o(n^{-1}+M_{\rm mean}^{-1}).
$$

This expansion requires a stable inverse-response matrix, independently
prepared reference and target batches, a fixed channel, controlled reference
preparation, and a bounded treatment of singular empirical inverses. It is
an asymptotic statement as all sample sizes grow with controlled ratios.
References must be exact or corrected with no fixed residual systematic
bias; more generally that residual must be
$o(\sqrt{n^{-1}+M_{\rm mean}^{-1}})$. Fixed preparation or drift errors
require an additional bias allowance and invalidate the zero-floor expansion.
For example, keep the inverse only in a fixed well-conditioned neighborhood
of its nonsingular limit and specify a bounded output outside it; bounded
features make departures exponentially rare. Taylor remainders are then
uniform over the compact input-state space. The fitted inverse can have
finite-sample bias; it is not exactly unbiased. Variances $v_j$ refer to the
full combined score $X$, so covariance between its two features is included.
The third probe and pilot overhead must be added;
they do not enter the two mean-variance terms displayed here.

For the nominal profile in [the certificate](CERTIFICATE.md),
$v_x\simeq0.269270494203$, $v_d\simeq0.319649084105$, and
$V<0.389156263913$. The two planning calculations are distinct:

| Account | Mean-reference allocation $N_x/N_d$ | Leading coefficient | Strict equal-cost leading reuse |
|---|---:|---:|---:|
| Empirical Bernstein radius, $\delta=0.05$ | 1.83564 | $C_E\simeq92.39747$ | $K\gtrsim5804$ |
| Regular plug-in risk | 3.36958 | $C_R\simeq4.19019318$ | $K\gtrsim264$ |

The second row maximizes the expected calibration error for each fixed target
state; the first certifies an operator error simultaneously for all inputs.
They have different confidence and risk meanings. Neither row is a necessary
reuse threshold or a proven finite preparation count. The leading thresholds
also omit any non-negligible third-probe and pilot cost.

## 6. One finite planning profile and practical interpretation

The companion [calibration-cost report](../verification/calibration-cost.json)
records a finite planning profile, retaining empirical Bernstein remainder
terms, explicit pilot and validation counts, target-copy ceilings, and integer
reuse thresholds. Model moments used as hypothetical observed statistics
must be labeled as such; they are not acquired validation data or a
guaranteed passing batch.

Take target RMSE $\varepsilon=0.01$, calibration failure budget
$\delta=0.05$, and frozen-score bound $B=2.337$. As an explicit forecast,
substitute the certified nominal first, second, and fourth moments for
hypothetical empirical moments. The unbiased sample variances include
their $N/(N-1)$ factors. Charge an illustrative $P=10{,}000$ pilot trials;
this charge does not establish that such a pilot will select this score.
No preparation error or channel drift is assumed in this example.

| Quantity | Finite planning value |
|---|---:|
| Validation counts $(N_x,N_d,N_0)$ | $(5{,}500{,}000,3{,}000{,}000,100{,}000)$ |
| All pilot and validation preparations $M$ | $8{,}610{,}000$ |
| Reported bias allowance $\beta$ | $0.00335254$ (rounded up) |
| Second-moment upper bound $U$ | $0.39416880$ (rounded up) |
| Target copies from unrounded bounds | $4{,}441$ |
| Ideal six-Pauli copies $n_0$ | $5{,}625$ |
| All-biased separable necessary copies $N_{\min}$ | $5{,}164$ |
| Strict equal-cost reuse versus Pauli | $7{,}272$ batches |
| Strict equal-cost reuse versus the lower bound | $11{,}909$ batches |

The count calculations include the finite Bernstein terms and integer
ceilings. They quantify what this hypothetical certificate would cost;
they are neither observed calibration results nor guaranteed-pass counts.
The comparison retains the conditional-risk qualification in Section 3.

The exact-channel result establishes a sampling improvement using the native
measurement. The present accounting does not establish a one-use total-cost
improvement when all reference and target preparations are charged equally.
Cheaper reference preparation, repeated use of a stable apparatus, and
tighter measured validation statistics can change that comparison. Their
values belong in the resource report. No universal impossibility claim is
made about calibrated native measurements.

For a concrete use case, [repeated chirality scans](REPEATED_STATE_SCANS.md)
translate these copy counts into full-cycle costs and a finite stability
interval. That account also derives when heterogeneous parallel plaquettes
can share a calibration and how unmatched sampling changes the bias budget.
