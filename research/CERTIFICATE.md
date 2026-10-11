# Certified complete-parity reconstruction

The standard-library checker [`parity_certificate.py`](parity_certificate.py)
certifies a feasible, exactly calibrated complete-parity estimator at

$$
u=16,\qquad \tau=\pi,\qquad \Phi_+=-2\pi/3,\qquad \Phi_-=-\Phi_+.
$$

All enclosures use integer or rational arithmetic. The certificate includes the
transcendental constants, propagator truncation, arithmetic rounding, calibration
rank, and second-moment bound. It does not promote floating-point eigensolver
residuals to error bounds. The scores were motivated by a floating-point
optimization; feasibility and the following bounds are independently certified.

## Estimator and exact calibration

Let $E_p$ be the plus-flux effect for parity pattern $p$, and let
$R=P_{13}$. Reflection changes both the flux and the parity pattern. On a
plus-flux trial, report $b_p$; on a minus-flux trial with observed pattern
$p$, report $-b_{\bar p}$, where
$\bar p=(p_0,p_3,p_2,p_1)$. Randomize the flux with equal probabilities, or use
equal independent sample allocations. Write

$$
O_p=\frac{E_p-RE_pR}{2}=c_p\chi+d_pD,
\qquad F_p=\frac{E_p+RE_pR}{2}.
$$

The $F_p$ form a POVM. Exact unbiasedness for every input state is equivalent to
$\sum_p b_pc_p=1$ and $\sum_p b_pd_p=0$. Its second-moment operator is

$$
T=\sum_p b_p^2F_p.
$$

For randomized fluxes, $n\operatorname{Var}(\widehat\chi)
=\operatorname{Tr}(\rho T)-\operatorname{Tr}(\rho\chi)^2\le\|T\|$.
Fixed equal flux allocations have variance no larger than this expression.

The following finite decimals define exact rational reference scores $b_p^0$:

| Pattern | Rational reference score |
|---|---:|
| 0001 | 0.736714019013 |
| 0010 | 0.252900613464 |
| 0100 | 2.336375499791 |
| 0111 | 0.730697335072 |
| 1000 | −0.660422527924 |
| 1011 | −0.288874449056 |
| 1101 | −0.305946882115 |
| 1110 | −0.540240033331 |

They are not themselves asserted to be exactly unbiased. Let $p=0111$,
$q=0100$, and define

$$
\Delta=c_pd_q-c_qd_p,\qquad
r_c=1-\sum_j b_j^0c_j,\qquad r_d=-\sum_j b_j^0d_j.
$$

The exact scores are $b_p=b_p^0+\alpha$, $b_q=b_q^0+\beta$, and
$b_j=b_j^0$ otherwise, where

$$
\alpha=\frac{r_cd_q-c_qr_d}{\Delta},\qquad
\beta=\frac{c_pr_d-r_cd_p}{\Delta}.
$$

Here the coefficients denote exact finite-model quantities. The checker encloses

$$
-0.003664877543543<\Delta<-0.003664877543542,
$$

$$
3.76658\times10^{-13}<\alpha<3.76659\times10^{-13},\qquad
-4.83984\times10^{-13}<\beta<-4.83983\times10^{-13}.
$$

Consequently the two-row calibration has rank two and these corrected scores
exist. The rounded coefficients in a report are not substitutes for exact
calibration or experimental parameter characterization.

## Certified performance

The full three-spin operator obeys

$$
\boxed{\|T\|<0.389156263913<0.4.}
$$

For the maximally mixed input, the randomized-flux variance coefficient obeys

$$
0.284895494200<n\operatorname{Var}(\widehat\chi)
<0.284895494206.
$$

These bounds include the exact score correction. They assume three input
particles, ideal total parity readout, and the calibrated uniform model.

The certificate also gives two conservative neighborhoods. In both, the scores
are recalibrated by the same two-outcome correction of $b^0$.

* If every plus-flux effect changes by operator norm at most $10^{-5}$, while
  the effects remain a complete spin-rotation-invariant POVM and the two fluxes
  remain related by reflection, the calibration minor stays nonzero and
  $\|T\|<0.45$.
* In the uniform spin-independent ring, if each of
  $|u-16|$, $|\tau-\pi|$, and $|\Phi_++2\pi/3|$ is less than $10^{-7}$,
  with $\Phi_-=-\Phi_+$, then
  $\|T\|<0.395388567053<0.4$.

These are sufficient continuity bounds, not estimates of experimental tolerance.
Uncharacterized errors with nominal scores, detector misclassification, and loss
require the separate bias budgets in the complete-parity analysis. In particular,
using single-odd parity patterns does not inherit return-event loss rejection.

## Stronger stability with all-outcome recalibration

Using all eight score corrections gives a larger certified neighborhood than
the two-outcome correction above. If each effect obeys
$\|E'_p-E_p\|\le0.0015$, and the perturbed effects remain a complete
spin-rotation-invariant POVM with reflection-related opposite fluxes, there are
exactly unbiased recalibrated scores satisfying

$$
\boxed{\|T'\|<0.536350926403<0.54<9/16.}
$$

This is a statement about existence and construction of calibrated scores;
using nominal scores without recalibration is not covered. In the uniform ring,
it applies throughout the sufficient box

$$
|u-16|,\ |\tau-\pi|,\ |\Phi_++2\pi/3|<0.00009,
\qquad \Phi_-=-\Phi_+,
$$

because the effect-norm bound is less than $0.001388254316470$ there. This
sufficient box is a mathematical continuity guarantee, not an experimentally
measured tolerance. The derivation below uses the known common eigenspaces of
all reflection-even effects to avoid a loose maximum-score-square estimate.

## Certified loss and bit-error window

The same rational score corrections certify

$$
A<2.337,\qquad C_1<1.499,\qquad
Q_0(I_8/8)<0.285,\qquad \|T\|<0.390,
$$

where $A$ is the maximum absolute score and $C_1$ is half the range of
singleton-parity scores. Thus the sufficient window in
[the noise budget](NOISE_BUDGET.md) applies to these exactly calibrated scores.
For every $0\le\lambda\le0.01$ and $0\le q\le0.003$, known
survival-and-no-flip correction gives

$$
|\mathbb EY-\langle\chi\rangle|<0.001430,\qquad
\mathbb EY^2<0.411.
$$

For maximally mixed input the second-moment bound improves to $0.302$.
These statements assume equal spin-independent loss and independent symmetric
parity-bit errors. They are conditional channel guarantees, not apparatus
characterization. Calibration uncertainty and uncertainty in the correction
rates require the separate budgets given in the noise analysis.

The checker evaluates the analytical endpoint formulas with rational upper
bounds $A=2.337$, $C_1=1.499$, $Q=0.390$ (or $0.285$). It
bounds $e^{0.01}$ and $e^{0.03}$ without floating arithmetic: for a positive
argument $x$, the Taylor sum through order 20 is a lower bound, and the
remaining positive tail is at most

$$
\frac{x^{21}}{21!}\frac{1}{1-x/22}.
$$

The error and second-moment budgets increase throughout the displayed
rectangle, as proved in the noise analysis. Exact rational interval arithmetic
therefore extends the endpoint certificate to the entire rectangle.

## Exact channel recalibration throughout the noise window

A second fixed image feature removes the modeled noise bias through direct
channel calibration. Let $v$ be an auxiliary score that reconstructs $D$ in
the ideal nominal model. The checker starts from four times
`RIGHT_INVERSE_SECOND` and applies the same two-outcome correction, now with
target coefficients $(0,1)$. It certifies

$$
\max_p|v_p|<194,\qquad C_{1,v}<107,\qquad \|T_v\|<646.
$$

For the full channel with $0\le\lambda\le0.01$ and $0\le q\le0.003$, extract
both signed features from every image, including zero-score events. Calibrating
their actual two-dimensional mean response gives an exactly unbiased chirality
score with

$$
\boxed{n\operatorname{Var}(\widehat\chi)
\le\|T_{\rm calibrated}\|
<0.444512137719<0.446<9/16.}
$$

The final per-trial score has absolute value less than $2.625215826513$.
This theorem assumes the exact full-channel response coefficients. It requires
neither separately known loss rates nor separately known bit-error rates:
the response can be obtained from the product probes in
[CHANNEL_CALIBRATION.md](CHANNEL_CALIBRATION.md). Finite empirical response
estimates must instead use that document's independent validation and bias
bounds. The $0.446$ guarantee here concerns the stated noise channel at the
nominal Hamiltonian parameters; it is not a joint assertion over the preceding
Hamiltonian perturbation box.

**Certificate.** For the proof only, divide both raw features by the common
factor $r=e^{-3\lambda}(1-q)^4$. Let their normalized mean responses be

$$
A_a=(1+x)\chi+yD,\qquad A_v=z\chi+(1+w)D.
$$

The same analytical noise budget, now also applied to the certified auxiliary
score, gives

$$
\|A_a-\chi\|<0.001429313207303=:\epsilon_a,\qquad
\|A_v-D\|<0.113367439057571=:\epsilon_v,
$$

with normalized second moments bounded by
$V_a<0.410838798136739$ and $V_v<701.993354716872$.
For coefficient bounds set

$$
X=4\epsilon_a/\sqrt3,\quad Y=\epsilon_a/\sqrt3,
\quad Z=4\epsilon_v/\sqrt3,\quad W=\epsilon_v/\sqrt3.
$$

Anticommutation of $\chi$ and $D$ implies
$|x|\le X$, $|y|\le Y$, $|z|\le Z$, and $|w|\le W$.
The checker evaluates the unrounded rational error bounds and encloses the
Schur complement using

$$
\kappa=1-X-\frac{YZ}{1-W}>0.
$$

The normalized response determinant is at least $(1-W)\kappa$, certified
larger than $0.931246421082804$. The unnormalized, directly measured response
has determinant larger than $0.856186181626231$, so it too is invertible.
For the first row $(w_a,w_v)$ of the normalized inverse,

$$
|w_a|\le\kappa^{-1}<1.003544558545665,\qquad
|w_v|\le\frac{Y}{(1-W)\kappa}<0.000886139636417.
$$

Keeping the Schur complement before bounding avoids an unnecessarily loose
ratio of independent numerator and determinant bounds. Cauchy–Schwarz then
gives the uniform second-moment bound

$$
\mathbb E(w_aY_a+w_vY_v)^2
\le(|w_a|\sqrt{V_a}+|w_v|\sqrt{V_v})^2
<0.444512137719.
$$

The two features share an image, and this inequality includes their covariance.
Finally, the raw response matrix is $r$ times the normalized response matrix;
its inverse absorbs $1/r$. Therefore an operational reconstruction using raw
reference responses implements the same score without knowing $r$.

## Nominal calibration-moment profiles

For variance-sensitive planning, the checker also encloses the second and fourth
moments of the exactly corrected nominal score on the product reference states.
These are ideal-model profiles, not measured variances of an implemented channel
or of a score fitted from pilot data. They use independent, equally randomized
flux settings.

| Reference | Exact mean | Certified upper bound on $\operatorname{Var}(X)$ | Certified upper bound on $\operatorname{Var}(X^2)$ |
|---|---:|---:|---:|
| $xyz$ | $1/8$ | 0.269271 | 0.255520 |
| $001$ | 0 | 0.319650 | 0.324633 |
| $000$ | 0 | 0.200631 | 0.170496 |
| $010$ (optional) | 0 | 0.299654 | 0.192583 |

The fourth-moment operator is $\sum_p b_p^4F_p$. It has the same three known
symmetry sectors as the second-moment operator. Multiplying its certified sector
eigenvalues by the product-probe weights in
[CHANNEL_CALIBRATION.md](CHANNEL_CALIBRATION.md) gives $\mathbb E X^4$.
The reported variance intervals then follow from
$\operatorname{Var}(X)=\mathbb E X^2-(\mathbb E X)^2$ and
$\operatorname{Var}(X^2)=\mathbb E X^4-(\mathbb E X^2)^2$.
The `calibration_moments` report section supplies outward-rounded intervals for
all these quantities and the score range, including the optional fourth probe.
Actual calibration guarantees must use a justified variance bound or
independently validated statistics for the score and channel being used.

## Arithmetic enclosure

The computation uses the 40-state sector with two spin-up and one spin-down
boson and evolves its three input columns $001,010,100$. This restriction
retains one component of each total-spin representation: the spin-3/2 scalar and
the full spin-1/2 multiplicity space. Therefore it determines every
spin-rotation-invariant input effect and its full eight-dimensional spectrum.

Set $s=2^{1024}$. Integer square roots enclose each needed $\sqrt m$ between
$\lfloor s\sqrt m\rfloor/s$ and that rational plus $1/s$.
Machin's identity

$$
\pi=16\arctan(1/5)-4\arctan(1/239)
$$

is enclosed using alternating rational series. Each fixed-point division has
error less than $1/s$; the first omitted term bounds the remaining series.
There are 219 and 64 retained terms respectively. The flux phase is algebraic,
$e^{-i\pi/6}=(\sqrt3-i)/2$, so no trigonometric evaluation is needed.

Let $A=-i\pi H$, and let $\widetilde A$ be its complex fixed-point midpoint
matrix. In the submultiplicative row norm using entry magnitude
$|\operatorname{Re}z|+|\operatorname{Im}z|$, the checker verifies

$$
\|A\|_\infty,\|\widetilde A\|_\infty\le B=256,
\qquad
\|A-\widetilde A\|_\infty\le206388/s.
$$

It computes the Taylor action through order 1280 on each input column (1281 terms including the constant term). Each real and
imaginary component is divided with integer floor rounding, giving complex
rounding error at most $2/s$ per component and recurrence step. With
$N=1280$, the uniform entry error of the final propagated columns is at most

$$
\varepsilon=
3^B\left[\frac{206388+2N}{s}
+\frac{B^{N+1}}{(N+1)!}\right]<10^{-170}.
$$

The first term combines the exponential's perturbation bound with all recurrence
rounding; the second bounds the exponential-series tail. The inequality
$e^B<3^B$ avoids transcendental arithmetic. For the recurrence error, each term
is bounded by $(2/s)e^B$, and summing $N$ terms gives the stated factor.
Integer summation itself introduces no error.

True propagated amplitudes have modulus at most one. Each effect matrix entry is
a sum of at most 40 amplitude products, so its error is at most

$$
\eta=40(2\varepsilon+\varepsilon^2)+2/s.
$$

The coefficient formulas on the three-input sector are
$c_p=\operatorname{Tr}(\chi E_p)/(3/8)$ and
$d_p=\operatorname{Tr}(DE_p)/6$. Their elementary entrywise enclosures have
radii $4\eta$ and $\eta$, respectively. Rational interval arithmetic then
certifies the determinant and exact score-correction bounds above.

## Spectrum and perturbation bounds

For any spin-rotation-invariant, reflection-even Hermitian operator on the
three-input sector, write its first row as $(a,b,c)$. Commutation with reflection
and with the symmetric-spin projector fixes the orthogonal eigenvectors to

$$
(1,1,1),\quad(1,0,-1),\quad(1,-2,1),
$$

with eigenvalues $a+b+c$, $a-c$, and $a+c-2b$. In particular, $b,c$ are
real. The checker bounds these linear expressions directly. No computed
eigenvalue is used in the proof. The maximally mixed second moment is
$(2\lambda_{\rm sym}+\lambda_{\rm odd}+\lambda_{\rm even})/4$, because the
symmetric eigenvalue has full-space multiplicity four and the other two each
have multiplicity two.

Since the $F_p$ form a POVM, changing scores changes the second-moment operator
by at most $\max_p|b_p^2-(b_p^0)^2|$. The checker adds the bound
$2|b_p^0|\,|\delta b_p|+|\delta b_p|^2$ for the two corrected outcomes.

For an effect perturbation bounded by $b_*$, duality between trace norm and
operator norm gives

$$
|\delta c_p|\le4b_*/\sqrt3,\qquad
|\delta d_p|\le b_*/\sqrt3.
$$

The reference-score second moment changes by at most
$b_*\sum_p(b_p^0)^2$. Applying the same interval determinant and correction
formulas proves the effect-neighborhood result. For uniform-ring parameter
changes, the sufficient effect bound used in the parameter box is

$$
b_*\le6|\delta\tau|
+\frac32\tau_{\max}(|\delta u|+|\delta\Phi|).
$$

Here $\|T_{\rm hop}\|\le6$, $\|W-3I/2\|=3/2$, and
$\|\partial_\Phi T_{\rm hop}\|\le3/2$. The commutator of a Hermitian operator
with a projector is bounded by the norm of the operator after an arbitrary
scalar shift. A Duhamel interpolation therefore gives the displayed constants.

## Proof of the all-outcome stability bound

Scale the calibration columns as $g_p=(c_p,4d_p)^T$ and let $G$ be the
resulting $2\times8$ matrix. On the total-spin-1/2 multiplicity space,
$\chi$ and $D/4$ are anticommuting Pauli directions with common norm
$\sqrt3/4$. Thus an effect perturbation of norm at most $\epsilon$ implies

$$
\|\delta g_p\|_2\le r:=4\epsilon/\sqrt3.
$$

Let $b$ be the exactly calibrated nominal score above. Define a rational
$8\times2$ matrix $B_0$: its first column is $b^0$, and its second column is
the eight exact finite decimals in `RIGHT_INVERSE_SECOND` in the checker. Set

$$
B=B_0(GB_0)^{-1},\qquad GB=I_2.
$$

The checker encloses the determinant of $GB_0$ between
$0.999999999999640$ and $0.999999999999641$, so this definition is nonsingular.
All entries of $B$ are obtained by exact rational interval formulas for the
$2\times2$ inverse. No floating-point inverse supplies a certified value.

For a perturbed calibration $G'=G+\delta G$, choose

$$
b'=b+Bz,\qquad
z=-(I_2+\delta GB)^{-1}\delta Gb.
$$

This satisfies $G'b'=(1,0)^T$ exactly. Write $B_{p,:}$ for its rows and
$L=\sum_p\|B_{p,:}\|_2$. If $rL<1$, the Neumann bound gives

$$
\|z\|_2\le\zeta:=\frac{r\|b\|_1}{1-rL}.
$$

At $\epsilon=0.0015$, the checker proves
$rL<0.270626586743$ and $\zeta<0.027794427232$.

Let $f_{pj}$ be the eigenvalue of the nominal $F_p$ in one of the three
common sectors $j$. They are nonnegative and obey $\sum_p f_{pj}=1$.
Define

$$
t_j=\sum_p f_{pj}b_p^2,\qquad
\ell_j=\sum_p f_{pj}b_p B_{p,:}^{T},\qquad
K_j=\sum_p f_{pj}B_{p,:}^{T}B_{p,:}.
$$

The nominal effects at the corrected score have sector second moment
$t_j+2\ell_j^Tz+z^TK_jz$. Since $K_j$ is positive semidefinite,
$\|K_j\|\le\operatorname{Tr}K_j$. Changing each effect by at most
$\epsilon$ adds at most $\epsilon\|b+Bz\|_2^2$. Consequently,

$$
\lambda_j(T')\le
 t_j+2\zeta\|\ell_j\|_2
 +\zeta^2\operatorname{Tr}K_j
 +\epsilon(\|b\|_2+\|B\|_F\zeta)^2.
$$

Every constant on the right is enclosed from the original rational effect
certificate, using integer square roots for the Euclidean and Frobenius norms.
The certified bounds for the symmetric, reflection-odd, and remaining
reflection-even sectors are respectively
$0.294482501767$, $0.536350926403$, and $0.423599386487$.
The largest gives the stated all-input second-moment bound.

## Reproduction

From the repository root, with ordinary (non-optimized) Python:

```bash
python research/parity_certificate.py --output /absolute/new/path/certificate.json
```

The script has no third-party dependencies. Its assertions check the rational
bounds, and all displayed decimals are rounded outward. Generated reports must
remain outside the repository; report creation is exclusive and refuses to
overwrite existing evidence.
