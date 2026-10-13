# Loss, parity errors, and calibration for complete-image readout

This account concerns the one-duration complete-parity estimator in
[COMPLETE_PARITY_THEORY.md](COMPLETE_PARITY_THEORY.md). It uses the same four-site,
spin-independent Hamiltonian as [the model](../THEORY.md). The inequalities
below are analytical. The [exact certificate](CERTIFICATE.md) establishes the
score and operator bounds used in the uniform noise window in Section 5.
The finer numerical table evaluates the bounds with the floating-point scores
in [complete_parity.py](complete_parity.py). The loss and detector channels
are diagnostic models; their parameters are not measured apparatus data.
[Direct channel calibration](CHANNEL_CALIBRATION.md) includes the whole
noisy response without separately known rates. For the worst-case-oriented
score, the [exact certificate](CERTIFICATE.md) gives an unbiased second
moment below $0.366461$ throughout $\lambda\le0.01$, $q\le0.003$.
The mixed-input score has a separate bound below $0.446$ in the same window.
Both assume exact response coefficients. The [finite-calibration theorem](FINITE_CALIBRATION.md)
accounts for finite reference data and preparation/drift errors. All generic
budgets below must be evaluated with the particular score being implemented.
The [parameter-robustness theorem](PARAMETER_ROBUSTNESS.md) combines
Hamiltonian variation with these uniform number-sector budgets. It proves
joint parameter/loss/detector regions and freshly certifies finite-reference
moment bounds in the smaller joint region. Its generic fixed-score and
loading-error guarantees are a distinct bias-aware treatment.

## 1. Score convention and input-trial accounting

Let $\mathcal O$ be the eight four-bit strings of odd Hamming weight. Choose
the two fluxes independently with probability $1/2$ on each input preparation.
At flux $-2\pi/3$, report score $a_p$. At flux $+2\pi/3$, reflect the
reported string, $r(p)=(p_0,p_3,p_2,p_1)$, and report $-a_{r(p)}$.
Every reported even-weight string has score zero. The ideal calibration is

$$
\frac12\sum_{p\in\mathcal O}a_p
   (E_{+,p}-E_{-,r(p)})=\chi.
$$

Write $X$ for the resulting per-trial score. All preparations, including
zero-score trials, count toward the sample size $n$. Define

$$
A=\max_{p\in\mathcal O}|a_p|,\qquad
A_1=\max_{|p|=1}|a_p|,\qquad
C_1=\frac{\max_{|p|=1}a_p-\min_{|p|=1}a_p}{2}.
$$

The ideal second-moment operator is

$$
T=\frac12\sum_{p\in\mathcal O}a_p^2
        (E_{+,p}+E_{-,r(p)}),\qquad
Q_0(\rho)=\operatorname{Tr}(T\rho).
$$

Thus $n\operatorname{Var}(\widehat\chi)=Q_0(\rho)-\langle\chi\rangle^2$
for randomized flux and ideal calibration. Fixed equal allocations have no
larger variance. Every operator bound below remains valid for inputs entangled
with an external reference: tensoring an effect error with the reference
identity preserves its operator norm.

## 2. Equal-rate loss: quadratic contamination after survival correction

Assume independent Markovian loss operators
$L_{j\sigma}=\sqrt\gamma\,b_{j\sigma}$, the same rate for all eight modes,
no gains, and a number-conserving Hamiltonian. Set

$$
\lambda=\gamma t,\quad \eta=e^{-\lambda},\quad
p_k={3\choose k}\eta^k(1-\eta)^{3-k},\quad s=p_3=\eta^3.
$$

These number probabilities are independent of the spin input and flux.
Indeed, taking the trace in each number block of the Lindblad equation gives
$\dot p_k=\gamma[(k+1)p_{k+1}-kp_k]$, with $p_3(0)=1$.
The solution is the displayed binomial law. In the three-particle block the
no-jump Hamiltonian is $H-i3\gamma I/2$, so its entire output equals
$sU\rho U^\dagger$, including any external reference.

The parity of the sum of reported ideal bits equals the parity of particle
number. Consequently one and three losses have zero score. Two losses leave
one particle, whose singleton parity can coincide with a three-particle
collision outcome. The complete-image estimator therefore does not have the
return event's exact rejection of every loss.

For known $\lambda$, use $Y=X/s$. Then

$$
\boxed{\left|\mathbb E Y-\langle\chi\rangle\right|
 \le \frac{p_1}{s}C_1
 =3C_1(e^\lambda-1)^2.}
$$

**Proof.** The three-particle contribution is exactly
$s\langle\chi\rangle$. For either flux the one-particle sector has trace
$p_1$. Its conditional mean score lies between the smallest and largest
singleton scores. Half the difference between the two flux means therefore
has absolute value at most $p_1C_1$. This proves the inequality without
assuming anything about the interaction history before the losses.

In particular, if all four singleton scores are equal, this contamination
cancels exactly for every input. Generic calibrated scores need not have that
property. Taking $C_1\le A_1\le A$ gives simpler, weaker bounds.

The second moment satisfies

$$
\boxed{\mathbb E Y^2\le
\frac{Q_0(\rho)+3A_1^2(e^\lambda-1)^2}{s}.}
$$

This follows by bounding the one-particle squared score by $A_1^2$.
The upper bound divided by $n$ also bounds the variance of the empirical
mean. Loss correction raises sampling cost even when its bias is small.

## 3. Independent parity-bit errors and loss–detector cross terms

Assume each of the four parity bits is flipped independently with known
probability $q<1$, independent of the state, loss history, and flux. This
classical channel acts on the total site-parity bits at the end of evolution.
It is invariant under the reflection used to relabel the second flux.
For strings $p,r$, write

$$
K_q(p\mid r)=q^{d(p,r)}(1-q)^{4-d(p,r)},
$$

where $d$ is Hamming distance. Useful probabilities are

$$
\alpha=(1-q)^4,\qquad
b=6q^2(1-q)^2+q^4,\qquad
c=4q(1-q)^3+4q^3(1-q).
$$

Here $\alpha$ is the no-flip probability, $b$ the probability of two or
four flips, and $c$ the probability of an odd number of flips. Odd true
parity survives the total-parity rejection test with probability
$\alpha+b$; even true parity can become a scored string with probability
$c$.

Use the simple survival-and-no-flip corrected score

$$
Y=\frac{X}{s\alpha}.
$$

The no-flip correction is necessary for the stated quadratic dependence on
$q$: rejecting a single flipped bit otherwise attenuates the clean signal
at first order. This correction does not remove every detector error.

Let $z=e^\lambda-1$. A bound valid for every input is

$$
\boxed{B_{\rm noise}:=
3C_1z^2+\frac{Ab}{\alpha}(1+3z^2)
 +\frac{3Ac}{\alpha}z,
\qquad
|\mathbb EY-\langle\chi\rangle|\le B_{\rm noise}.}
$$

Its terms are, respectively, two-loss aliases, wrong odd-to-odd parity labels,
and one-loss events that an odd number of detector errors makes look odd.
The leading dependences are $O(\lambda^2)$, $O(q^2)$, and
$O(q\lambda)$. A loss and a detector error cannot be treated as two
separate independent bias budgets: their product opens an additional accepted
event. Three losses give the vacuum; its mean contribution cancels exactly
between the fluxes because the vacuum and its detector-error law are the same
for both settings. Vacuum events still contribute to the second moment.

### A sharper finite-outcome bound

Define two explicit functions on the sixteen strings,

$$
f(r)=\sum_{p\in\mathcal O}a_pK_q(p\mid r),\qquad
g(r)=\sum_{p\in\mathcal O}a_p^2K_q(p\mid r).
$$

For a real function $v$ on a finite set $S$, set
$\operatorname{osc}_S(v)=(\max_S v-\min_S v)/2$. Use
$S_1=\{r:|r|=1\}$ and $S_2=\{r:|r|=0\text{ or }2\}$. Then

$$
k_3=\operatorname{osc}_{\mathcal O}(f-\alpha a),\qquad
k_1=\operatorname{osc}_{S_1}(f),\qquad
k_2=\operatorname{osc}_{S_2}(f)
$$

give the tighter bound

$$
\boxed{|\mathbb EY-\langle\chi\rangle|
\le\frac{k_3+3z^2k_1+3zk_2}{\alpha}.}
$$

**Proof.** Condition on particle number. Both fluxes have the same sector
trace $p_k$. The difference of their conditional expectations of any
function lies within its range, so half that difference is bounded by its
oscillation. The three-particle no-flip part is calibrated exactly. The
remaining three-particle function is $f-\alpha a$, the one-particle
function is $f|_{S_1}$, and the two-particle function is $f|_{S_2}$.
The vacuum function is a constant and hence cancels. Finally,
$p_1/s=3z^2$ and $p_2/s=3z$.

The simpler boxed bound follows from
$k_3\le bA$, $k_1\le\alpha C_1+bA$, and $k_2\le cA$.
All of these inequalities are uniform in the post-loss spin and spatial
states; the actual loss trajectory need not be reconstructed.

### Sampling and mean-square error

Set $G_1=\max_{S_1}g$, $G_2=\max_{S_2}g$, and
$G_0=g(0000)$. A usable second-moment bound is

$$
\boxed{\mathbb EY^2\le
\frac{\alpha Q_0(\rho)+bA^2+3z^2G_1+3zG_2+z^3G_0}
     {s\alpha^2}=:V_{\rm budget}(\rho).}
$$

For the three-particle term use
$g(r)\le\alpha a_r^2+bA^2$; use the displayed maxima for the other
number blocks. Replacing $Q_0(\rho)$ by $\|T\|$ gives an all-input
bound. If the bias is bounded by $B$, then

$$
\mathbb E[(\widehat\chi-\langle\chi\rangle)^2]
\le V_{\rm budget}(\rho)/n+B^2.
$$

Thus a small variance alone does not prove lower estimation error at
arbitrarily large $n$: an unknown systematic bias sets an error floor.
Independence of the preparations also gives the distribution-free statement

$$
\Pr\!\left[
|\widehat\chi-\langle\chi\rangle|
>B+\frac{A}{s\alpha}\sqrt{\frac{2\log(2/\delta)}n}
\right]\le\delta.
$$

### Optional exact detector recalibration

For known $0\le q<1/2$, one may replace the odd-image scores by
$\widetilde a=K_{\mathcal O,\mathcal O}^{-1}a$, keep even-image scores
zero, and divide only by $s$. The retained eight-by-eight kernel is
invertible. With $r=1-2q$, its eigenvalues are
$(1+r^4)/2$, $(r+r^3)/2$, and $r^2$, with multiplicities
$1,4,3$. One obtains them by restricting Walsh characters to one parity
coset; complementary characters become the same character up to sign.

This exactly recovers the target score on every odd true pattern, so the
three-particle detector bias vanishes and the one-particle alias bound is
again $3C_1z^2$. The two-particle term remains, bounded by
$3z\operatorname{osc}_{S_2}(K\widetilde a)$. The vacuum mean still cancels.
Second moments must be recomputed using $\widetilde a^2$, and diverging
inverse weights near $q=1/2$ prevent a uniform sampling-cost claim.

## 4. Hamiltonian and score calibration

No reflection symmetry of the parameter error is needed for the following
norm bound. Work in dimensionless units with
$H(u,\phi)=T_{\rm hop}(\phi)+uW$ on the complete three-particle sector.
Here $0\le W\le3I$, $\|T_{\rm hop}\|\le6$, and

$$
\|T_{\rm hop}(\phi')-T_{\rm hop}(\phi)\|
\le12\left|\sin\frac{\phi'-\phi}{8}\right|
\le\frac32|\phi'-\phi|.
$$

The hopping estimate follows first on the single-particle ring and then by
second quantization to three bosons. Subtracting the physically irrelevant
scalar $3uI/2$ centers the interaction spectrum. Duhamel's formula then
gives, up to a scalar phase,

$$
\|U'-U\|\le\epsilon:=
\frac32|\tau'u'-\tau u|+6|\tau'-\tau|
 +12|\tau'|\left|\sin\frac{\phi'-\phi}{8}\right|.
$$

Let $C=(\max_{\mathcal O}a-\min_{\mathcal O}a)/2$. Centering the
three-particle score observable gives, for potentially different errors in
the two flux settings,

$$
\boxed{B_H\le
C\,[\min\{1,\epsilon_+\}+\min\{1,\epsilon_-\}].}
$$

Indeed each setting's score observable has centered norm $C$, its
conjugation changes by at most $2C\min\{1,\epsilon_\pm\}$, and each
setting enters with coefficient $1/2$. For common error bounds and
$u=16,\tau=\pi$, its small-error form is

$$
B_H\le2C\left[
\frac{3\pi}{2}|\delta u|+30|\delta\tau|
+\frac{3\pi}{2}|\delta\phi|
+\frac32|\delta\tau|(|\delta u|+|\delta\phi|)
\right].
$$

This is a conservative sufficient budget, not a measured apparatus
sensitivity. All angles are in radians.

If the nominal finite-matrix reconstruction has operator residual
$\beta_{\rm rec}$, add it to $B_H$. For score rounding or calibration
error $\delta a$, the ideal-setting mean changes by at most
$(\max\delta a-\min\delta a)/2$, and in particular by at most
$\|\delta a\|_\infty$. Evaluate the noise constants using the scores
actually implemented. Exact calibration can instead absorb these changes
into the reconstruction residual.

For a second-moment perturbation, use
$C_2=(\max a_p^2-\min a_p^2)/2$ in place of $C$.
Then $Q_0(\rho)$ changes by at most
$C_2[\min\{1,\epsilon_+\}+\min\{1,\epsilon_-\}]$.
Insert this increase into the second-moment budget above.

Finally, suppose the implemented normalization is
$\widehat r=e^{-3\widehat\lambda}(1-\widehat q)^4$, while the true
clean-event factor is $r=e^{-3\lambda}(1-q)^4$. Then

$$
|\mathbb E(X/\widehat r)-\langle\chi\rangle|
\le |r/\widehat r-1|\,\frac{\sqrt3}{4}
 +(r/\widehat r)(B_{\rm noise}+B_H+\beta_{\rm rec}).
$$

If $|\lambda-\widehat\lambda|\le d_\lambda$,
$|q-\widehat q|\le d_q$, and both bit-error rates are at most
$q_{\max}<1$, one may bound

$$
|r/\widehat r-1|\le
\exp\!\left(3d_\lambda+\frac{4d_q}{1-q_{\max}}\right)-1.
$$

These formulas require equal, spin-independent loss and the stated bit-error
channel. Unequal loss, particle gains, correlated readout errors, and a
flux-dependent error law require their own calibrated channel.

## 5. A certified uniform window and a numerical diagnostic

For the worst-case-oriented score at $u=16$, $\tau=\pi$, the
[certificate](CERTIFICATE.md) proves $A<1.623$, $C_1<1.123$, and
$\|T\|<0.3287$. The formulas above give, uniformly over
$0\le\lambda\le0.01$, $0\le q\le0.003$,

$$
|\mathbb EY-\langle\chi\rangle|<0.001018,\qquad
\mathbb EY^2<0.344785.
$$

These are bounds for survival-and-no-flip correction; the bias remains.
The two-feature full-channel inverse removes that bias under exact
calibration and has second moment below $0.366461<13/32$. Without particle
loss, inversion of the odd-parity bit-error kernel instead gives an exactly
unbiased coefficient below $0.382273$ for known $q\le0.03$. The certificate
derives and checks these two different corrections separately.

### Mixed-input score and diagnostic table

For the exactly calibrated mixed-input estimator at $u=16$, $\tau=\pi$,
the [rational-arithmetic certificate](CERTIFICATE.md#mixed-input-score-loss-and-bit-error-window)
proves $A<2.337$, $C_1<1.499$, and $\|T\|<0.390$.
For every $0\le\lambda\le0.01$, $0\le q\le0.003$, and every input,
known survival-and-no-flip correction then gives

$$
|\mathbb EY-\langle\chi\rangle|<0.001430,
\qquad \mathbb EY^2<0.411.
$$

These are certified guarantees for the stated equal-loss and independent
bit-error channels with known rates. To verify uniformity, use
$G_1\le(\alpha+b)A^2$ and
$G_2,G_0\le cA^2$ in the second-moment formula. Both budgets increase
throughout this rectangle: $b/\alpha=6r_q^2+r_q^4$ and
$c/\alpha=4r_q+4r_q^3$, where $r_q=q/(1-q)$, are increasing,
as are $1/\alpha$, $z$, and $1/s$. The certificate evaluates the endpoint
budgets with rational exponential enclosures and includes the exact score
correction. Its additional bound $Q_0(I_8/8)<0.285$ gives a second moment
below $0.302$ for maximally mixed input. Calibration and normalization
uncertainty remain separate additions as specified in Section 4. The window
does not assert that these rates or calibration accuracies have been achieved
in an apparatus.

For the nominal $u=16,\tau=\pi$ mixed-input optimal scores,
$A=A_1=2.3363754998$, $C_1=1.4983990139$,
$Q_0(I_8/8)\simeq0.284895494$, and $\|T\|\simeq0.389156264$.
At the illustrative diagnostic values $\lambda=0.01$, $q=0.003$,
using known survival and no-flip correction gives:

| Quantity | Upper bound, evaluated with nominal floating-point scores |
|---|---:|
| Two-loss contribution to simple bias bound | 0.000454042 |
| Wrong odd-to-odd labels | 0.000126964 |
| Loss–detector cross term | 0.000847866 |
| Total simple bias bound | 0.001428871 |
| Sharper finite-outcome bias bound | 0.000678311 |
| Second moment for maximally mixed input | 0.299782851 |
| Second moment for arbitrary input | 0.408517793 |

The finite-outcome values use
$k_3=9.3958552\times10^{-6}$, $k_1=1.480490134$, and
$k_2=0.0070379270$. Parameter uncertainty and normalization uncertainty
must be added through Section 4. These sampling bounds remain below the
ideal separately sampled six-Pauli coefficient $9/16$ in this diagnostic,
but a comparison of total error must include the displayed squared bias and
the comparator's own noise and calibration budget.
