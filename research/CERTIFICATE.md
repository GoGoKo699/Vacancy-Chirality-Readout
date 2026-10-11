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

## Reproduction

From the repository root, with ordinary (non-optimized) Python:

```bash
python research/parity_certificate.py --output /absolute/new/path/certificate.json
```

The script has no third-party dependencies. Its assertions check the rational
bounds, and all displayed decimals are rounded outward. Generated reports must
remain outside the repository; report creation is exclusive and refuses to
overwrite existing evidence.
