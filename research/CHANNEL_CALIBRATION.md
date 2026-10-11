# Direct channel calibration with product reference states

Complete-parity reconstruction can be calibrated directly on the implemented
measurement channel. Spin-independent loss and reflection-symmetric detector
errors change the calibration coefficients but preserve the two-observable
structure. Their rates need not be inferred separately. Two product reference
states determine the mean response; three product reference states also
determine the complete second-moment operator of any fixed score.

This result complements the parameter-dependent certificates in
[CERTIFICATE.md](CERTIFICATE.md) and the explicit error channels in
[NOISE_BUDGET.md](NOISE_BUDGET.md). It replaces a narrow Hamiltonian-parameter
window with a directly measured response and an explicit statistical
calibration budget. Symmetry, stable reference preparation, and stability of
the channel between calibration and use remain assumptions.

## 1. The complete noisy channel retains the odd two-dimensional space

Include every final event, including lower-particle-number outcomes, and let
$E_{s,y}$ be the resulting POVM effect on the original eight-dimensional
three-spin input. The setting $s=+$ denotes the first flux, nominally
$-2\pi/3$, and $s=-$ denotes its reflected partner. For the sixteen observed
parity strings, write $r(y_0,y_1,y_2,y_3)=(y_0,y_3,y_2,y_1)$.
Assume

$$
[E_{s,y},U^{\otimes3}]=0\quad\text{for every spin rotation }U,
\qquad E_{-,r(y)}=R E_{+,y}R,\quad R=P_{13}.
$$

The effects sum to the identity at each setting. If an acquisition can fail
to return an image, add a separate reflection-invariant failure label, normally
with score zero, beyond the sixteen recorded parity strings. Missing events
must be included as outcomes; conditioning on successful trials generally changes the
linear measurement model.

**Channel theorem.** Choose any real parity feature $f(y)$. Randomize the two
fluxes equally and assign $X=f(y)$ at the first flux and
$X=-f(r(y))$ at the second. Its mean operator and second-moment operator obey

$$
A_f=\frac12\sum_y f(y)(E_{+,y}-R E_{+,y}R)=c\chi+dD,
\qquad
T_f=\frac12\sum_y f(y)^2(E_{+,y}+R E_{+,y}R),
$$

where $D=P_{12}-P_{23}$. In particular, $T_f$ commutes with collective
rotations and reflection and has the three scalar blocks described below.

**Proof.** Collective-spin invariance places every effect in the commutant
spanned by $I,P_{12},P_{23},P_{13},\chi$. The first displayed combination is
reflection-odd, whose Hermitian subspace is exactly
$\operatorname{span}\{\chi,D\}$. The squared-score combination is
reflection-even. This argument acts on the original input effects and does
not require the final state to contain three particles.

Uniform spin-independent hopping and interactions, equal spin-independent
loss, and a spin-blind final detector satisfy the first assumption. The
reflected dynamics and a classical parity-error kernel satisfying
$K_-(r(y)\mid r(p))=K_+(y\mid p)$ satisfy the second. Independent equal-rate
bit flips are one example. The symmetry theorem holds at arbitrary rates;
useful reconstruction additionally requires a nonvanishing, well-conditioned
response. Complete erasure of input information cannot be calibrated away.

## 2. Two product probes determine the complete mean response

Use the same orbital loading, dynamics, and detector as for the unknown input,
changing only the prepared spin state. Define

$$
\rho_x=|{+x},{+y},{+z}\rangle\langle{+x},{+y},{+z}|,
\quad
\rho_d=|{+z},{+z},{-z}\rangle\langle{+z},{+z},{-z}|,
\quad
\rho_0=|{+z},{+z},{+z}\rangle\langle{+z},{+z},{+z}|.
$$

Here $|{+a}\rangle$ is the positive eigenstate of Pauli axis $a$.
The first two probes obey

$$
(\langle\chi\rangle_x,\langle D\rangle_x)=(1/8,0),
\qquad
(\langle\chi\rangle_d,\langle D\rangle_d)=(0,1).
$$

Indeed a product state with Bloch vectors $\mathbf r_j$ has
$\langle\chi\rangle=\mathbf r_1\cdot(\mathbf r_2\times\mathbf r_3)/8$
and $\langle D\rangle=\mathbf r_2\cdot(\mathbf r_1-\mathbf r_3)/2$.
Both probe responses attain their respective maximum absolute values over
fully separable states. Consequently the factors $8$ and $1$ below cannot
be reduced by choosing other separable probes that isolate the same
coordinates. This is a sensitivity statement, not a minimum-shot theorem.

For two simultaneously recorded signed features $Z=(Z_1,Z_2)^T$, let
$M$ have rows $(c_j,d_j)$. Then

$$
\mathbb E_\rho Z=M
\begin{pmatrix}\langle\chi\rangle_\rho\\\langle D\rangle_\rho\end{pmatrix},
\qquad
\boxed{M=\bigl[\,8\mathbb E_xZ\,,\,\mathbb E_dZ\,\bigr].}
$$

If $\det M\ne0$, $w^T=(1,0)M^{-1}$ gives an exactly unbiased score
$X=w^TZ$ for the complete noisy channel. No loss-rate or detector-error-rate
estimate enters this identity. The same two probe batches can calibrate any
number of predetermined features extracted from each image.

For exact $M$ and independent unknown inputs, the variance coefficient is
$w^T\operatorname{Cov}_\rho(Z)w$. The features come from the same image:
their covariance must be retained rather than adding two independent-noise
estimates. Two independent probe response vectors are sufficient and, for
an unrestricted two-coordinate response, necessary. Preparing the product
axes is an additional calibration resource, even though the unknown-input
measurement uses only spin-independent dynamics and parity readout.

For the nominal finite-repulsion point, the [integer/rational certificate](CERTIFICATE.md) supplies a concrete nonzero-noise example. Two exactly defined features, calibrated to $\chi$ and $D$ in the ideal model, retain an invertible full-channel response throughout $0\le\gamma t\le0.01$, $0\le q\le0.003$ for equal loss and independent bit errors. Recalibrating their raw responses gives an exactly unbiased score with all-input second moment below $0.444512137719$. No separately known rate enters the operational inversion. The assertion uses exact channel coefficients; finite reference data require the validation below.

## 3. Freeze the score before validating finite calibration

A practical procedure has three stages:

1. Use pilot reference data to choose features and a bounded score $f$, for
   example by inverting an empirical two-feature response matrix.
2. Freeze that score and collect independent validation data on the reference
   states below. Certify its actual mean and second moment.
3. Apply the fixed score to independent unknown preparations, counting every
   preparation in the estimator's sample size.

An empirical inverse is not an exactly known inverse. In general
$\mathbb E[\widehat M^{-1}]\ne M^{-1}$, and an inverse can become unstable
near a singular empirical response. Independence of pilot and unknown data
does not remove the resulting calibration bias. The validation step bounds
that bias without treating the inverse as unbiased or differentiating a
random self-normalization formula. If the score is changed after inspecting
validation results, use fresh validation data or a confidence statement that
is uniform over the entire score-selection procedure.

For the fixed score define $m_j=\mathbb E_jX$ and $Q_j=\mathbb E_jX^2$.
The exact first-moment coefficients are simply

$$
c=8m_x,\qquad d=m_d.
$$

Let $\Pi_{\rm s}$ be the symmetric spin-$3/2$ projector and
$\Pi_{\rm e},\Pi_{\rm o}$ the mixed-spin doublet projectors with reflection
eigenvalues $+1,-1$. The three probe weights in these blocks are

| Probe | $\Pi_{\rm s}$ | $\Pi_{\rm e}$ | $\Pi_{\rm o}$ |
|---|---:|---:|---:|
| $\rho_x$ | $1/2$ | $1/4$ | $1/4$ |
| $\rho_d$ | $1/3$ | $1/6$ | $1/2$ |
| $\rho_0$ | $1$ | $0$ | $0$ |

Since $T_f=t_{\rm s}\Pi_{\rm s}+t_{\rm e}\Pi_{\rm e}+t_{\rm o}\Pi_{\rm o}$,
these three second moments determine it exactly:

$$
\boxed{t_{\rm s}=Q_0,\qquad
t_{\rm o}=3Q_d-2Q_x,\qquad
t_{\rm e}=6Q_x-3Q_d-2Q_0.}
$$

Thus three product probes suffice for both response and variance; three
independent second-moment equations are necessary to reconstruct an
unrestricted three-block operator. These equalities use actual noisy outcome
frequencies, including losses and mislabeled images.

Applying the same construction to every outcome-indicator feature reconstructs
each individual effect: its odd part has the two calibrated coefficients and
its even part has the three reconstructed block values. The three preparation
types therefore suffice for complete measurement calibration within this
five-dimensional symmetry class. Pilot data can support full-image score
optimization; an independently validated fixed score avoids having to certify
that optimization or every inverted calibration coefficient separately.

There is also an exact all-input bias formula. The identities
$\{\chi,D\}=0$, $\chi^2=3(I-\Pi_{\rm s})/16$, and
$D^2=3(I-\Pi_{\rm s})$ imply

$$
\boxed{\sup_\rho|\mathbb E_\rho X-\langle\chi\rangle_\rho|
=\|A_f-\chi\|
=\sqrt{\frac3{16}(c-1)^2+3d^2}.}
$$

The randomized worst-case variance still equals $\|T_f\|$ even if the mean
is slightly miscalibrated: a state in the largest reflection block has zero
expectation of every reflection-odd $A_f$.

## 4. Explicit finite-sample confidence certificate

Suppose $|f(y)|\le B$ and collect independent validation batches of sizes
$N_x,N_d,N_0\ge2$, fixed before observing those batches. Let
$\widehat m_x,\widehat m_d$ be the sample means and
$\widehat Q_x,\widehat Q_d,\widehat Q_0$ the sample second moments. The
default variance-sensitive certificate also uses the unbiased sample variances

$$
\widehat v_j=\frac1{N_j-1}\sum_i(X_{j,i}-\widehat m_j)^2,\qquad
\widehat w_j=\frac1{N_j-1}\sum_i(X_{j,i}^2-\widehat Q_j)^2.
$$

For a failure probability $0<\delta<1$, take $L_E=\log(20/\delta)$ and

$$
\boxed{\epsilon_j=
\sqrt{\frac{2\widehat v_jL_E}{N_j}}
+\frac{14BL_E}{3(N_j-1)},\qquad
\nu_j=\sqrt{\frac{2\widehat w_jL_E}{N_j}}
+\frac{7B^2L_E}{3(N_j-1)}.}
$$

These are affine rescalings of
[Maurer and Pontil (2009), Theorem 4](https://arxiv.org/abs/0907.3740),
applied to $X\in[-B,B]$ and $X^2\in[0,B^2]$. Applying the one-sided
theorem to both tails and allocating failure probability $\delta/5$ to each
of the five statistics gives the displayed logarithm. With probability at
least $1-\delta$, all five mean and second-moment errors are within their
respective radii. The theorem permits the mean and sample variance to use
the same observations; no independence among statistics from the same
probe batch is assumed. It does require the score to have been frozen
independently of these validation data. Repeatedly checking the same batch
until a bound passes requires a separate confidence allocation over the
looks or a valid confidence sequence.

The simpler range-only alternative remains

$$
L_H=\log(10/\delta),\qquad
\epsilon_j=B\sqrt{2L_H/N_j},\qquad
\nu_j=B^2\sqrt{L_H/(2N_j)}.
$$

It follows from Hoeffding's inequality. Choose the method before inspecting
validation results, or divide the failure budget between the two methods
before intersecting their confidence intervals. Either set of radii can be
inserted into the same operator certificate. On its confidence event define

$$
b_c=|8\widehat m_x-1|+8\epsilon_x,\qquad
b_d=|\widehat m_d|+\epsilon_d,\qquad
\beta=\sqrt{3b_c^2/16+3b_d^2},
$$

and

$$
\begin{aligned}
U_{\rm s}&=\widehat Q_0+\nu_0,\\
U_{\rm o}&=3\widehat Q_d-2\widehat Q_x+3\nu_d+2\nu_x,\\
U_{\rm e}&=6\widehat Q_x-3\widehat Q_d-2\widehat Q_0
             +6\nu_x+3\nu_d+2\nu_0,\\
U_T&=\min\{B^2,\max\{0,U_{\rm s},U_{\rm o},U_{\rm e}\}\}.
\end{aligned}
$$

Then the actual, fixed score obeys $\|A_f-\chi\|\le\beta$ and
$\|T_f\|\le U_T$ for **every** input state. For the mean of $n$ independent
unknown-input scores this gives the conditional mean-square-error guarantee

$$
\boxed{\mathbb E_\rho[(\widehat\chi-\langle\chi\rangle_\rho)^2
\mid\text{fixed pilot and validation data}]
\le U_T/n+\beta^2.}
$$

The inequality holds for all validation data in the stated confidence event;
it is not an unconditional assertion that the random calibration is exact.
If a certificate is issued only when its observed bounds pass a threshold,
the probability of issuing a false certificate is at most $\delta$.
This does not imply $1-\delta$ coverage conditional on passing: that
conditional error probability can be as large as
$\delta/\Pr(\mathrm{pass})$. Rejecting a validation run without reporting
an estimate is abstention, not an unconditional estimation procedure with
the displayed mean-square error.
Independently, for any $0<\eta<1$, the probability of reporting an estimate
whose error $|\widehat\chi-\langle\chi\rangle_\rho|$
exceeds $\beta+B\sqrt{2\log(2/\eta)/n}$ is at most $\delta+\eta$.
The operator certificate is uniform in the input and also applies when its
three spins are entangled with an unmeasured reference.

## 5. Calibration cost, conditioning, and reuse

The known probes consume preparations and characterization time. They are
additional resources to the $n$ unknown input copies. Finite preparation
counts and the number of uses needed to amortize them are treated in
[the calibration-cost account](CALIBRATION_COST.md).

Let $r_\beta=\sqrt{12\epsilon_x^2+3\epsilon_d^2}$. For the Hoeffding
radii its squared value is

$$
r_\beta^2=B^2L_H\left(\frac{24}{N_x}+\frac6{N_d}\right).
$$

At fixed $N_x+N_d=N_{\rm mean}$, this sufficient bound is minimized by
$N_x:N_d=2:1$, giving $54B^2L_H/N_{\rm mean}$. This allocation optimizes
the Hoeffding radius, not the information-theoretic calibration cost. With
empirical Bernstein radii, the observed variances and the finite-$N$
linear terms must also be included.

Crucially, $r_\beta$ equals the reported bias bound only when the observed
mean residuals are zero. A forecast that substitutes their nominal
expectations is not a guarantee that a validation batch will pass. On the
confidence event, if the actual fixed score has bias norm $\beta_{\rm true}$,
the reported bound satisfies

$$
\beta_{\rm true}\le\beta\le\beta_{\rm true}+2r_\beta.
$$

The second inequality follows by applying the triangle inequality to the two
weighted mean errors. Similarly, writing
$r_T=\max\{\nu_0,3\nu_d+2\nu_x,6\nu_x+3\nu_d+2\nu_0\}$,
one has $U_T\le\|T_f\|+2r_T$ on that event. A prospective passing
guarantee must budget these fluctuations as well as the confidence radii.

For a desired mean-square error $\varepsilon^2$, the displayed certificate
requires $\beta<\varepsilon$ and permits
$n\ge U_T/(\varepsilon^2-\beta^2)$ unknown-input preparations, rounded up.
Tight systematic-error tolerances can therefore require a large reference
sample even when the subsequent unknown-input variance is small.

The three-probe second-moment inversion is economical in preparation types
but can amplify statistical errors. An optional fourth probe
$\rho_v=|{+z},{-z},{+z}\rangle\langle{+z},{-z},{+z}|$ gives

$$
t_{\rm e}=(3Q_v-Q_0)/2,\qquad
t_{\rm o}=2Q_d-(Q_0+Q_v)/2.
$$

For equal second-moment error radii $\nu$, this bounds the largest eigenvalue
with additive uncertainty at most $3\nu$, instead of $11\nu$ in the
three-probe inversion. It can therefore reduce the required number of
reference preparations despite adding one preparation type. The two mean
probes remain $\rho_x,\rho_d$.

The same fixed-score certificate can be reused for different unknown input
states and repeated batches while the physical channel remains stable.
This can amortize calibration effort; it does not make calibration free.
If all preparations are charged equally, compare an actual total resource
budget containing both known and unknown preparations. A finite-validation
estimator can have residual bias, so the exact unbiased separable benchmark
must be compared using the displayed mean-square-error budget rather than
variance alone.

## 6. What the reference data do not establish

These few probes do not themselves establish collective-spin covariance or
the paired reflection relation. A symmetry-breaking operator can vanish on
all tested probes and still bias an untested input. Those structural
assumptions require physical justification or additional characterization.

Reference states must have the stated spin expectations and use the same
orbital preparation and channel as the unknown inputs. If a reference state
has trace-distance error at most $\eta_j$, adding $2B\eta_j$ to its
mean-error radius and $B^2\eta_j$ to its second-moment-error radius is
always sufficient. The following sharper correction can be much smaller.
For that sharper correction, the preparation error must be a spin-state
error within the same three-spin orbital embedding: the actual reference
and its target are states of the same eight-dimensional input space.
Wrong particle number, different orbital loading, or leakage outside that
embedding do not satisfy this assumption. When such errors have an
independently established full-preparation trace-distance bound, use the
coarse $2B\eta_j$ and $B^2\eta_j$ allowances instead.

Let $\beta_{\rm stat},U_{\rm stat}$ denote the full observed certificate
above, including the observed mean residuals but before preparation and
drift corrections. Set

$$
\kappa_A=\sqrt{48\eta_x^2+12\eta_d^2},\qquad
\kappa_T=\max\{\eta_0,3\eta_d+2\eta_x,
                         6\eta_x+3\eta_d+2\eta_0\}.
$$

If both are less than one, a sufficient corrected certificate is

$$
\boxed{\beta_{\rm use}\le
\frac{\beta_{\rm stat}+\kappa_A\sqrt3/4}{1-\kappa_A}+d_A,
\qquad
U_{\rm use}\le
\min\left\{B^2,\frac{U_{\rm stat}}{1-\kappa_T}\right\}+d_T.}
$$

Here $d_A$ and $d_T$ bound the operator-norm drift of the mean and
second-moment operators between validation and use. To prove the preparation
correction, a reference perturbation changes its mean by at most
$2\eta_j\|A_f\|$ and its second moment by at most $\eta_j\|T_f\|$.
Propagating the first pair of errors through the exact bias norm gives
$\beta_{\rm cal}\le\beta_{\rm stat}+\kappa_A\|A_f\|$.
Use $\|A_f\|\le\sqrt3/4+\beta_{\rm cal}$ and rearrange.
The three second-moment equations similarly give
$\|T_f\|\le U_{\rm stat}+\kappa_T\|T_f\|$; positivity of $T_f$
is what permits the smaller second-moment perturbation bound. Finally apply
the triangle inequality to the channel drift. The coarse radius additions
remain available if these denominators are not useful.

Use $\beta_{\rm use},U_{\rm use}$ in the mean-square-error budget. These
are independently justified error allowances, not quantities established by
the three probe data themselves. In particular, spin-dependent loss,
unpaired flux errors, or spin-selective detector responses are not
automatically absorbed by the two-coordinate calibration theorem.

The exact statement is therefore a calibration theorem for a symmetry class
of implemented measurement channels, with finite-data certification of the
resulting estimator. It is not a claim of exact finite-sample inversion,
calibration-free readout, or coherent retained-state measurement.
