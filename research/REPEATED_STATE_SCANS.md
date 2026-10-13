# Repeated-state scans with a fixed chirality readout

A candidate use is a scan of chirality against a state-preparation parameter,
such as the duration of a spin quench. Each point prepares an unknown
three-spin state and then applies the same calibrated vacancy-hopping and
parity measurement. This account states when calibration reuse yields a
valid resource comparison; it does not establish a measured time advantage.

The [finite polynomial calibration](FINITE_CALIBRATION.md) now supplies a
separate unconditional comparison with the passive four-mode class, including
all reference copies and bounded preparation error and drift. The conditional
validation method and separable-comparator examples below remain distinct:
their numerical budgets must not be transferred between scores or methods.

## 1. What must remain fixed across the scan

At the start of readout every target must have the orbital embedding
$(0,1,1,1)$ used in [the model](../THEORY.md). The preparation parameter may
change the spin density operator, including its entanglement with an
unmeasured reference. It must not leave additional charge sectors or
spin–motion correlations outside that embedding. The readout dynamics must
also be decoupled from any unmeasured reference as assumed by the model.

The scan can be implemented conceptually by preparing spins while their
occupancies are pinned, ending the preparation dynamics, and applying a
fixed flux-hopping interval and terminal detector. The [channel
calibration](CHANNEL_CALIBRATION.md) is reusable if the resulting measurement
channel is the same for every target and product reference. Equal control
settings alone do not prove this: a longer preparation can change atom
survival, motional excitation, or the noise present during readout.

The comparison treats each scan point as an arbitrary unknown state and
requires the stated accuracy separately at each point. A known dynamical
model or promised smooth curve can allow a comparator to infer unmeasured
points; the independent-task copy bound does not describe that different
estimation problem. Repeating one identical state produces more copies for
one estimate, not an arbitrary number of separate calibration reuses.

Collective spin rotations do not determine the unknown response coefficients:
both $\chi$ and $D$ are invariant under them. Trusted rotational averaging
can enforce spin covariance, but does not calibrate the chirality gain or
the paired reflection relation. The nonparallel product reference axes
still require the relative spin preparations specified in the calibration
account.

Three primary experiments motivate this workflow, with distinct demonstrated
capabilities:

| Source | Evidence relevant to the candidate scan |
|---|---|
| Nascimbène et al. [N12](../REFERENCES.md), Figs. 2–3 and Appendix | Spinful bosonic four-site plaquettes, parallel spin manipulation, and a 5 ms RVB preparation ramp. This is a spin-preparation precedent with four occupied sites and ensemble detection. |
| Aidelsburger et al. [A13](../REFERENCES.md), abstract and Section VI | Internal-structure-independent orbital flux control and single-particle plaquette dynamics. It does not demonstrate the spinful interacting detector considered here. |
| Impertro et al. [I25](../REFERENCES.md), Fig. 2 and Supplement I, IV | Interacting cesium plaquettes, flux control, and parity microscopy. The reported sequence has 220 ms Mott loading and additional 30 ms and 15 ms plaquette ramps; some images are rejected because of phase-lock fluctuations or drift. |

The last example illustrates shared preparation cost and a need to characterize
stability. Its orbital pseudospins are distinct from internal spins; its timing
does not provide a cost ratio for this protocol. These sources establish
components, rather than one apparatus satisfying the full readout assumptions
and reuse inequalities. Single-shot coherence or current lifetimes do not
measure the calibration's stability over repeated experimental cycles.

## 2. Count complete experimental cycles

Let $M$ include all pilot and validation reference preparations, including
unsuccessful calibration runs. Let $n$ be the target-copy count per scan point obtained
from its final calibration certificate. Write $c_R$ for the actual average
cost per counted reference preparation, $c_N$ for a native target trial,
and $c_S$ for a comparator target trial. A cost can mean elapsed time or
another explicitly stated resource. Include loading attempts, state
preparation, readout, imaging, and reset in its definition. If reference
types have different costs, use $M c_R=\sum_j M_j c_{R,j}$.

For $K$ separate target tasks and comparator count budget $m$, the native
budget is strictly smaller precisely when

$$
\boxed{M c_R<K(m c_S-n c_N).}
$$

Unshared setup or switching costs must also be charged on the appropriate
side. The displayed criterion already gives the comparator free calibration;
an actual comparator's measured overhead can instead be included explicitly.
If $m c_S\le n c_N$, reuse cannot make this budget comparison favorable. Otherwise
the smallest strict integer threshold is

$$
K_{\min}=\left\lfloor\frac{M c_R}{m c_S-n c_N}\right\rfloor+1.
$$

For the ideal six-Pauli implementation,
$m=\lceil9/(16\varepsilon^2)\rceil$ is an attained worst-case copy count.
Assigning it a measured $c_S$ compares two specified implementations.
The [all-estimator separable lower bound](SEPARABLE_BENCHMARK.md#6-finite-sample-bound-allowing-biased-estimators)
instead gives the necessary count
$m=\lceil(3/(4\varepsilon)-\pi)^2\rceil$ for
$0<\varepsilon<1/(16\pi)$. At $\varepsilon=0.01$ this is $5164$.
It is not an attained implementation. Turning that count into a cost bound
for the entire separable class requires either equal costs for all target
copies or an independently justified lower bound $c_S$ on their cost.
Copy counts alone do not establish a bound on every apparatus's wall time.

For the passive four-mode class, the [product-prior argument](FINITE_CALIBRATION.md)
allows adaptive allocation across independently specified target tasks and
gives a necessary expected total of
$K(38907/97600)(\varepsilon^{-2}-1600\pi^2/3)$ copies when this is positive.
At RMSE $0.001$, the per-task real bound is greater than $396538.952$.
Round the total once for a fixed integer budget; rounding each task separately
would not be justified for random adaptive allocations. The unconditional
native examples count all three reference stages, but still require enough
distinct tasks within one epoch satisfying the specified stability bounds.

For the conditional validation method below, beating a valid separable-class
cost lower bound is sufficient on the calibration confidence event; failing
to beat that bound does not establish that every actual separable
implementation is faster. The polynomial method instead uses unconditional
risk and the passive-class comparison just described.

Product references share much of the native experimental cycle. If

$$
c_R=c_{\rm shared}+c_{\rm ref\,prep},\qquad
c_N=c_{\rm shared}+c_{\rm target\,prep},
$$

then $c_R/c_N\ge c_{\rm shared}/(c_{\rm shared}+c_{\rm target\,prep})$.
This common cost includes the same loading, native readout, imaging, and
reset, with matched success accounting. A simple spin state does not make
these stages disappear. Different loading yields must be measured and
included rather than absorbed into an assumed small reference weight.

For the existing hypothetical RMSE-$0.01$ profile, $M=8{,}610{,}000$ and
$n=4441$. If native and comparator target cycles cost the same, write
$w=c_R/c_N$. The strict condition is $w<K(m-n)/M$:

| Distinct target tasks $K$ | Required $w$, ideal Pauli $m=5625$ | Required $w$, equal-cost separable lower bound $m=5164$ |
|---:|---:|---:|
| 1 | $<0.0001375145$ | $<0.0000839721$ |
| 100 | $<0.01375145$ | $<0.00839721$ |
| 1000 | $<0.1375145$ | $<0.0839721$ |

Entries round the exact thresholds down. Thus even a thousand required
tasks need reference cycles more than about 7.27 times cheaper than target
cycles for the Pauli comparison, or 11.91 times cheaper for the equal-cost
lower-bound comparison. Under the common-cycle model and zero reference-only
overhead, the latter requires target-only preparation to exceed the shared
cycle by more than 10.91 times. These are requirements, not measured ratios.

Without equal target-cycle costs, beating the chosen budget even with free
references requires
$c_N/c_S<m/n$: approximately $1.266607$ for Pauli or $1.162801$ for the
assumed class cost bound. The full-cycle formula, not the rounded ratios,
determines a strict comparison.

The forecasts and conditional-confidence qualifications in
[CALIBRATION_COST.md](CALIBRATION_COST.md) apply unchanged. In particular,
a hypothetical passing validation record is not acquired evidence, and a
certificate with failure probability $\delta$ is not an unconditional
mean-square-error guarantee for a random calibration procedure.

## 3. Fit the reuse count inside a stable calibration epoch

Let $H$ be a justified interval over which the required channel and
reference-preparation bounds remain valid, including validation acquisition.
Let $M_{\rm epoch}$ references be acquired within it, $S$ be calibration
setup time, $t_{\rm setup}$ be setup time per target task, and $R_R,R_N$
be justified sample throughputs. For sequential calibration followed by the
target scan, the maximum number of complete target batches that fit is

$$
K_{\max}=\max\left\{0,
\left\lfloor\frac{H-S-M_{\rm epoch}/R_R}
{n/R_N+t_{\rm setup}}\right\rfloor\right\}.
$$

Set $M_{\rm epoch}=M$ unless an earlier independent pilot can legitimately
sit outside this epoch; its cost still belongs in $M c_R$. Serial cycle
times are recovered with $R_R=1/c_R$ and $R_N=1/c_N$. The rates must count
statistically justified samples, with the parallel-data conditions below.
The no-drift, equal-cost forecast needs $61{,}497{,}869$ such preparations
to cross the separable count bound at $K=11909$. At a common throughput
$R$, its minimum sequential acquisition time is that number divided by $R$,
before setup or monitoring. No value of $R$ or $H$ is inferred from the
component experiments.

A feasible advantage requires $K_{\min}\le K_{\max}$ and at least
$K_{\min}$ scientifically required target tasks. Count additional calibration
or monitoring within each epoch. A large planned scan does not justify
carrying a calibration past its demonstrated stability interval.

Drift must also fit the accuracy budget. Suppose the validated score has
mean operator $A$, second-moment operator $T$, bias bound $\beta$, and
second-moment bound $U$. During use suppose
$\|A'-A\|\le d_A$ and $\|T'-T\|\le d_T$. Then

$$
\operatorname{MSE}\le\frac{U+d_T}{n}+(\beta+d_A)^2.
$$

Writing $g=\varepsilon^2-U/n-\beta^2$, a sufficient remaining drift budget is

$$
\boxed{\frac{d_T}{n}+2\beta d_A+d_A^2\le g.}
$$

Choosing the smallest integer $n$ allowed by a no-drift certificate leaves
only its rounding slack $g$. The current finite forecast uses nearly all
of its accuracy budget this way. A practical design must reserve drift
margin by tightening calibration, increasing $n$, or shortening the epoch;
it cannot add a drift allowance without revising the count calculation.

In the worked profile, $g\simeq3.75864\times10^{-9}$ at $n=4441$.
With $d_T=0$, this allows only $d_A\simeq5.6052\times10^{-7}$.
Reserving the illustrative allowance $d_A=10^{-4}$ before selecting the
target count instead gives $n=4476$ and $K_{\min}=12515$ against the
equal-cost separable lower bound. This is an operator allowance, not a
measured drift rate or a tolerance on an experimental control parameter.

## 4. When parallel plaquettes can share calibration

**Matched-mixture lemma.** Let $E_{i,s,y}$ be the input POVM effect for
plaquette or channel label $i$, flux $s$, and reported outcome $y$. Suppose
each reference and target trial draws its label independently of the spin
state from the same known distribution $p_i$. The resulting measurement is

$$
\overline E_{s,y}=\sum_i p_i E_{i,s,y}.
$$

Indeed, for any input $\rho$ the outcome probability is
$\sum_i p_i\operatorname{Tr}(E_{i,s,y}\rho)
=\operatorname{Tr}(\overline E_{s,y}\rho)$.
Thus calibration of this mixture is valid for every input. Individual
plaquettes need not be identical. The averaged effects must satisfy the
spin-rotation and paired-reflection assumptions; a common distribution
over individually symmetric paired channels is sufficient.

This is not the replacement
$\sum_i p_i\operatorname{Tr}(E_i\rho_i)
\mapsto\operatorname{Tr}(\overline E\,\overline\rho)$ for arbitrary
different states tied to different channels. If scan points occupy fixed
plaquette subsets, they generally have different mixtures. Use matched
randomized assignments, separate calibrations, or a mismatch bound.

**Mismatch bound.** For one frozen score with $|X|\le B$, let $A_i,T_i$
be its channel-specific mean and second-moment operators. They satisfy
$-BI\le A_i\le BI$ and $0\le T_i\le B^2I$. If reference labels have
weights $p$ and target labels have weights $q$, define
$\operatorname{TV}(p,q)=\frac12\sum_i|p_i-q_i|$. Then

$$
\boxed{\|A_q-A_p\|\le2B\operatorname{TV}(p,q),\qquad
\|T_q-T_p\|\le B^2\operatorname{TV}(p,q).}
$$

For the first inequality use the triangle inequality and $\|A_i\|\le B$.
For the second, $\sum_i(q_i-p_i)=0$ allows replacing each $T_i$ by
$T_i-B^2I/2$, whose norm is at most $B^2/2$. Both follow immediately.
Insert these bounds as $d_A,d_T$ above. They quantify unmatched mixture
weights; physical changes in the component channels need an additional
drift allowance.

## 5. Parallel throughput and the operational decision

The mixture lemma concerns means and moment operators. It does not prove
independence among plaquettes in one camera frame. Common flux noise or
detector fluctuations can correlate their scores and invalidate an
independent-sample confidence calculation. For example, $m$ equal-variance
scores per frame, variance $\sigma^2$, pair correlation $r$, and $F$
independent frames give

$$
\operatorname{Var}(\overline X)
=\frac{\sigma^2[1+(m-1)r]}{mF}.
$$

Justify independence or analyze independent frame blocks with their measured
covariance and bounded scores. Parallelization can reduce elapsed time,
but substituting a variance-based effective sample size into the existing
i.i.d. confidence formula is not justified by this variance identity alone.
Even a shared randomized flux setting can correlate scores within a frame.
The statistical certificate must match the actual setting and noise design,
and the same available capacity must be allowed to the Pauli comparator.
The reuse count $K$ counts separate estimation tasks, not plaquettes.

A concrete assessment therefore needs the preparation sequence and orbital
state at readout, actual cycle times and loading yields, a validated
calibration record, its stability interval and drift bounds, and the
plaquette assignment and noise correlations. Specify whether the output
requires per-point RMSE or simultaneous confidence for an entire curve.
The operator calibration certificate is uniform over inputs, but simultaneous
target-statistics confidence and correlations from shared calibration still
require their own treatment.

The candidate scan is supported by the channel and estimation theorems.
A target-matched full-cycle and stability advantage is not yet established;
the measured inputs above determine whether the operational criterion holds.
