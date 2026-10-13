# Vacancy chirality readout: model, results, and proofs

A localized vacancy in a four-site ring converts three-spin chirality into position statistics. Complete parity images permit exact reconstruction at one duration and two opposite fluxes. A certified interacting readout beats a universal precision bound for passive four-mode interferometers with the same terminal parity detector. The passive comparison permits arbitrary estimators and adaptation across fresh copies; its finite-risk formulation is given in the [passive theorem](research/PASSIVE_PARITY_BOUND.md). The separate separable-effects comparison, return-only reconstruction, and retained-state filter are complementary results. For a known total particle number at most three, total site parity recognizes the return event exactly. Retaining an entangled output additionally requires spin-coherent measurement back-action.

## 1. Model and measurement

Three spin-½ bosons initially occupy sites 1, 2, 3; site 0 is empty. Their input spin state $\rho$ is arbitrary and may be entangled with an external reference. With $\hbar=1$, the spin-independent Hamiltonian is

$$
H_{U,\Phi}=-J\sum_{j,\sigma}\left(e^{i\Phi/4}b^\dagger_{j\sigma}b_{j+1,\sigma}+\mathrm{h.c.}\right)+UW,
\qquad W=\frac12\sum_j n_j(n_j-1).
$$

Indices are modulo four and $n_j$ counts both spin components. The four bonds are equal; interactions are independent of spin. Write $u=U/J$ and $\tau=Jt$. The fixed-three-particle Hilbert space has dimension 120. The no-double-occupancy projector $P$ has rank 32.

Let $V$ embed the eight spin states into occupancy $(0,1,1,1)$ and $S=VV^\dagger$. The return operator and effect are

$$
K_{U,\Phi}(t)=V^\dagger e^{-itH_{U,\Phi}}V,
\qquad E_{U,\Phi}(t)=K_{U,\Phi}(t)^\dagger K_{U,\Phi}(t).
$$

The probability is $\operatorname{Tr}(E\rho)$. An ideal spin-preserving accepted instrument produces the unnormalized state $K\rho K^\dagger$. Other outcomes terminate the trial. At finite $U$, detecting only that site 0 is empty does not identify the return event.

Define the cyclic spin permutation and scalar chirality by

$$
C|s_1s_2s_3\rangle=|s_2s_3s_1\rangle,
\qquad \chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3),\qquad \mathbf S=\boldsymbol\sigma/2.
$$

With this convention $C-C^\dagger=4i\chi$. The eigenvalues of $\chi$ are $0,\pm\sqrt3/4$, with projectors $\Pi_0,\Pi_+,\Pi_-$ of ranks 4, 2, 2. For the displayed hopping convention, flux $-2\pi/3$ selects positive chirality in the hard-core filter below. Labels $+$ and $-$ refer to selected chirality, not the sign of the flux.

## 2. Parity identifies the accepted occupancy exactly

For $s=(0,1,1,1)$ define the projector onto total site parity $s$:

$$
Q_s=\prod_{j=0}^3\frac{I+(-1)^{s_j}(-1)^{n_j}}2.
$$

**Parity identity.** On the entire particle-number subspace $N\le3$, $Q_s=S$, extending $S$ by zero on $N<3$.

**Proof.** Odd occupation at each of sites 1, 2, 3 requires at least three particles. With at most three particles, each occupation is exactly one and site 0 is empty. The converse is immediate. This proves a projector equality, including coherent superpositions and reference-entangled inputs. Consequently,

$$
V^\dagger e^{itH}Q_s e^{-itH}V=K^\dagger K.
$$

Ideal total-parity images therefore implement the required effect even when intermediate doublon and triplon amplitudes are appreciable. Hopping must be stopped during the modeled readout, and the reported bits must represent the pre-imaging total site parities. Additional particles invalidate the identity: $(0,3,1,1)$ has the same parity at $N=5$. Loss alone cannot produce three odd occupied sites from fewer than three particles.

Parity imaging is an established technique [B24](REFERENCES.md); the identity specifies when its event equals this protocol's return event. It says nothing by itself about the retained spin state.

## 3. Exact finite-repulsion reconstruction

Set $R=P_{13}$ and $D=P_{12}-P_{23}$, where $P_{jk}$ swaps spins. Then $\|D\|=\sqrt3$ and $\|\chi\|=\sqrt3/4$. Define $E_+(t)=E_{U,-2\pi/3}(t)$ and $E_-(t)=E_{U,+2\pi/3}(t)$.

**Contrast theorem.** For every $U$ and time in the uniform model,

$$
E_+(t)-E_-(t)=g_u(\tau)\chi+h_u(\tau)D.
$$

The real coefficients are the Hilbert–Schmidt projections onto the orthogonal operators $\chi$ and $D$.

**Proof.** Collective spin rotations commute with the dynamics and occupancy projection. Their three-qubit commutant is the five-dimensional span of $I,P_{12},P_{23},P_{13},\chi$ [EW01](REFERENCES.md). Spatial reflection fixes the vacancy and return pattern, exchanges sites 1 and 3, and reverses the flux. Hence $K_-=RK_+R$, and the effect difference is odd under $R$. The odd Hermitian subspace is precisely the span of $\chi$ and $D$.

Reciprocity also gives $K_-=K_+^T$ in the real occupation basis, but it does not imply $E_-=E_+^*$ for a nonnormal return operator. The $D$ term generally survives. Reflection-symmetric scalar imperfections preserve the symmetry argument with recalibrated coefficients; arbitrary finite-$U$ scalar disorder need not do so.

Measure independent input copies at durations $t_1,t_2$ and both flux signs. Let $d_i=p_{i,+}-p_{i,-}$ and

$$
M=\begin{pmatrix}g_1&h_1\\g_2&h_2\end{pmatrix},\qquad \Delta=g_1h_2-g_2h_1.
$$

Whenever $\Delta\ne0$,

$$
\boxed{\langle\chi\rangle=\frac{h_2d_1-h_1d_2}{\Delta}=w_1d_1+w_2d_2.}
$$

This is an arbitrary-input expectation identity. It neither assigns a chirality eigenvalue to a single return nor corrects the accepted state's back-action.

**Existence of suitable durations.** At every fixed finite $u>0$, exact short-time coefficients give

$$
g_u(\tau)=4\sqrt3\,\tau^4+O_u(\tau^6),\qquad
h_u(\tau)=\sqrt3\left(\frac{7u}{72}+\frac{u^3}{2520}\right)\tau^9+O_u(\tau^{11}).
$$

Reciprocity makes $g$ even and $h$ odd in time. Taking $\tau_1=s$, $\tau_2=2s$ yields

$$
\Delta=5952\left(\frac{7u}{72}+\frac{u^3}{2520}\right)s^{13}+O_u(s^{15})>0
$$

for sufficiently small positive $s$. The coefficients are certified by exact integer/rational arithmetic in $\mathbb Q(\sqrt2,\sqrt3,i)[u]$ on the complete 40-state fixed-magnetization sector; see the [certificate](source/resource-review-2026-10-10/prior/exact_short_time.py) and [derivation, Section 4](source/resource-review-2026-10-10/prior/FOLLOWUP.md). This fixed-$u$ existence result does not imply efficient estimation at very short times or a uniform large-$u$ limit.

An order-one numerical example uses $u=16$:

| $\tau$ | $g$ | $h$ |
|---|---:|---:|
| $\pi$ | 1.76499071140 | 0.026588956448 |
| $3\pi/4$ | −0.39157486049 | 0.020164887920 |

Here $\Delta\simeq0.04600240679$ and $\langle\chi\rangle\simeq0.43834419389d_1-0.57799055104d_2$. Exact coefficients are defined by the finite-matrix functions, not these rounded decimals. The values are floating-point diagnostics rather than certified numerical intervals or device calibrations. A single duration at $\tau=\pi$ leaves the bias $(h_1/g_1)\langle D\rangle$; even $|\uparrow\uparrow\downarrow\rangle$ has zero chirality but contrast $h_1$.

### Complete parity images: one duration and two fluxes

At exactly three particles, the eight possible parity strings have odd Hamming weight. Let $E_p$ be the effect for pattern $p$ at flux $-2\pi/3$, and let $r(p_0,p_1,p_2,p_3)=(p_0,p_3,p_2,p_1)$. At the opposite flux, the **reflected physical pattern** $r(p)$ has effect $RE_pR$. Each contrast obeys

$$
O_p=\frac{E_p-RE_pR}{2}=c_p\chi+d_pD.
$$

This follows from the same collective-spin and reflection symmetries as the return contrast. Choose either flux independently with probability $1/2$. Assign score $a_p$ to pattern $p$ at the first flux and $-a_p$ to pattern $r(p)$ at the second. Then

$$
\boxed{\sum_pc_pa_p=1,\qquad\sum_pd_pa_p=0}
$$

is exactly the condition for an unbiased chirality estimate on every input. Rank two of the matrix with columns $(c_p,d_p)^T$ suffices. All original preparations count, and even-weight parity strings receive score zero. With $F_p=(E_p+RE_pR)/2$ and $T=\sum_pa_p^2F_p$, the randomized sample mean obeys

$$
n\operatorname{Var}_\rho(\widehat\chi)=\operatorname{Tr}(T\rho)-\langle\chi\rangle_\rho^2,
\qquad \sup_\rho n\operatorname{Var}_\rho(\widehat\chi)=\|T\|.
$$

The last equality follows because $T$ has three scalar total-spin/reflection blocks, and a state in its largest-eigenvalue block has zero chirality. Fixed equal allocations between fluxes have no larger variance. The [complete proof](research/COMPLETE_PARITY_THEORY.md) also derives the weighted Gram formula for the minimum mixed-input variance.

At $u=16$, $\tau=\pi$, an [integer/rational certificate](research/CERTIFICATE.md) proves a nonzero calibration minor and constructs exactly feasible scores by correcting two rational reference coordinates. The score designed for worst-case precision satisfies

$$
\boxed{\sup_\rho n\operatorname{Var}_\rho(\widehat\chi)<0.328697776905<13/32.}
$$

This is a certified feasible score, without a claim of global optimality. The [passive theorem](research/PASSIVE_PARITY_BOUND.md) derives $13/32$ as a lower bound for arbitrary spin-independent passive $U(4)$ followed by complete site parity, including product spin unitaries before mixing. It also proves three-mode impossibility, gives an explicit passive four-mode reconstruction, and extends the comparison to biased estimators through a finite minimax risk bound. Number resolution, extra spatial modes, intermediate detection, and spin controls interleaved with spatial mixing are different measurement resources.

Exact recalibration of the modeled channel retains a worst-case coefficient below $0.366461$ for equal-rate loss $\gamma t\le0.01$ and independent bit errors $q\le0.003$. The [finite-calibration theorem](research/FINITE_CALIBRATION.md) replaces exact response coefficients by independent product-reference batches, with unconditional error and preparation-cost bounds, including explicit preparation and drift allowances.

A separate score minimizes the variance at the maximally mixed input. Its nearby exactly corrected version satisfies

$$
\|T_{\rm mixed}\|<0.389156263913,\qquad
0.284895494200<n\operatorname{Var}_{I_8/8}(\widehat\chi_{\rm mixed})<0.284895494206.
$$

For this mixed-input score, the certificate also proves $\|T\|<0.536350926403<9/16$ after all-outcome recalibration whenever each parity effect changes by operator norm at most $0.0015$, preserving the stated symmetries and completeness. In the uniform ring this includes the sufficient box $|u-16|,|\tau-\pi|,|\Phi_++2\pi/3|<0.00009$ with $\Phi_-=-\Phi_+$. These are sufficient mathematical neighborhoods, not measured control tolerances. They concern that score and the separable benchmark; they do not extend the passive separation to the same parameter box. The [noise budget](research/NOISE_BUDGET.md) treats loss aliases, detector errors, calibration, and normalization separately. The return-only existence theorem remains established for every fixed finite $u>0$.

## 4. Hard-core return instrument and retained entanglement

In the strict hard-core limit, transporting spin labels with the vacancy reduces each cycle eigenspace $C=e^{i\theta}$ to a scalar four-site ring with flux $\Phi+\theta$ [K22](REFERENCES.md). Its return amplitude is

$$
A(\tau;f)=\frac12\left[\cos\left(2\tau\cos\frac f4\right)+\cos\left(2\tau\sin\frac f4\right)\right],
\qquad K_{\infty,\Phi}=\sum_{\theta=0,\pm2\pi/3}A(\tau;\Phi+\theta)\Pi_\theta.
$$

At $\tau=\pi$ and flux $\mp2\pi/3$,

$$
K_{\infty,\pm}=\Pi_\pm+a(I-\Pi_\pm),\quad
a=\frac{\cos(\sqrt3\pi)-1}{2},\quad r=a^2=0.02786714004\ldots.
$$

Thus $E_{\infty,\pm}=rI+(1-r)\Pi_\pm$ and $\langle\chi\rangle=\sqrt3(p_+-p_-)/[4(1-r)]$. A single accepted return is a coherent filter, with nonzero acceptance outside the selected sector.

From input $I_8/8$, the positive filter yields

$$
p_1=\frac{1+3r}{4}=0.2709003550\ldots,\qquad
\rho_1=\frac{\Pi_++r(I-\Pi_+)}{2+6r},\qquad
\langle\chi\rangle_{\rho_1}=\frac{\sqrt3}{4}\frac{1-r}{1+3r}=0.3884692918\ldots.
$$

Every biseparable three-qubit state obeys $|\langle\chi\rangle|\le1/4$ [R23](REFERENCES.md). For a state product across $1|23$, rotate the first spin's polarization to $z$: its magnitude is at most $1/2$, while $\|S_{2x}S_{3y}-S_{2y}S_{3x}\|=1/2$. The bound follows, then extends to the other cuts and their convex mixtures. Hence $\rho_1$ is genuinely tripartite entangled. Its one-spin marginals are maximally mixed and its pair marginals are separable.

With ideal coherent readout after every accepted round,

$$
K_{\infty,+}^m=\Pi_++a^m(I-\Pi_+),\qquad
p_m=q+(1-q)r^m,\quad q=\operatorname{Tr}(\Pi_+\rho).
$$

For $I_8/8$, $p_m\to1/4$ and the conditional state tends to the mixed doublet $\Pi_+/2$. The probability counts the entire successful record per original preparation. Spin coherence and atom survival between rounds are necessary assumptions. The [hard-core derivation](source/resource-review-2026-10-10/prior/prior/PILOT.md) gives the cycle reduction, pair marginals, and static scalar-disorder extension in detail.

## 5. Uniform finite-$U$ operation bound

For every $u>0$, $\tau\ge0$, and flux, the accepted operators satisfy

$$
\boxed{\|K_{U,\Phi}(t)-K_{\infty,\Phi}(t)\|\le\epsilon_u(\tau),\qquad
\epsilon_u(\tau)=\min\{2,20\sqrt3\,\tau/u\}.}
$$

**Proof.** Write $H=T+UW$, $Q=I-P$, $A=PTP$, and $B=QTP$. The one-particle hopping norm gives $\|T\|\le6J$; the hard-core vacancy graph gives $\|A\|\le2J$. Counting two-hop paths gives $\|B\|^2\le12J^2$: at zero flux the nonnegative Gram matrix has every row sum equal to $12J^2$, and flux cannot increase absolute path sums. Thus $\|TP\|\le4J$.

For a unit vector $\psi\in P$, energy conservation and $W\ge Q$ give

$$
U\|Qe^{-itH}\psi\|\le U\|We^{-itH}\psi\|
\le\|H\psi\|+\|T\|\le10J.
$$

Therefore $\|Qe^{-itH}P\|\le\min\{1,10J/U\}$. Duhamel's formula for the projected equation then bounds its deviation from $e^{-itA}P$ by $t\|B\|(10J/U)=20\sqrt3 J^2t/U$. Compression with $V$ proves the claim. The [full proof and path-count verification](source/resource-review-2026-10-10/prior/FOLLOWUP.md) distinguish leakage probability from coherent spin-operation error.

The success maps obey $\|\mathcal K_U-\mathcal K_\infty\|_\diamond\le2\epsilon_u$, including arbitrary external references. With a terminal failure flag and discarded failed spins, the complete flagged instruments differ by at most $\min\{2,4\epsilon_u\}$ in diamond norm.

For a specified input with ideal success probability $p_0$, finite-$U$ probability $p_U$, and $p_U>0$,

$$
|\sqrt{p_U}-\sqrt{p_0}|\le\epsilon_u,\qquad
D(\rho_U^{\mathrm{succ}},\rho_0^{\mathrm{succ}})\le\min\{1,\epsilon_u/\sqrt{p_0}\}.
$$

Here $D$ is trace distance. Purification reduces the proof to two unnormalized accepted vectors separated by at most $\epsilon_u$; distance to the line of the second vector bounds their normalized pure-state angle. Partial trace preserves the bound. For $m$ rounds, replace $\epsilon_u$ by $\min\{2,m\epsilon_u\}$ and use the whole-record ideal success probability.

At $u=2048$, $\tau=\pi$, these conservative bounds prove $p_U>0.2184$ and $\langle\chi\rangle>0.3000>1/4$ for the maximally mixed input. This is a sufficient analytical example, not a necessary repulsion threshold. At $u=16$, full-matrix diagnostics give four-round success $0.1141865390$ and chirality $0.4324656707$. Fixed-round convergence as $u\to\infty$ does not establish an infinite-round survival plateau at finite $u$.

## 6. Detector back-action, losses, and calibration

Correct return probabilities do not fix the retained state. A trace-preserving spin channel $\mathcal N$ after ideal acceptance gives

$$
\mathcal S_{\mathcal N}(\rho)=\mathcal N(K\rho K^\dagger),\qquad
\operatorname{Tr}\mathcal S_{\mathcal N}(\rho)=\operatorname{Tr}(K\rho K^\dagger).
$$

Complete local dephasing produces a fully separable spin output with exactly the ideal success probabilities. Destructive terminal imaging can therefore serve expectation estimation without realizing the retained-state instrument.

For the explicit diagnostic of independent phase damping with single-spin coherence factor $\eta$, $(\mathcal D_\eta^{\otimes3})^\dagger(\chi)=\eta^2\chi$. The one-round hard-core output remains certified by the chirality witness precisely when $\eta>0.8022166239\ldots$. This is a witness threshold for that channel. Under the same disturbance after every accepted round, the exact recurrence has spectral radius below one for every $\eta<1$, so total success tends to zero; see [the recurrence and proof, Section 3](source/resource-review-2026-10-10/READOUT_AND_RESOURCES.md).

For equal-rate Markovian loss $L_{j,\sigma}=\sqrt\gamma\,b_{j,\sigma}$ with no gains, every loss trajectory is rejected by the ideal **return parity event** $0111$. Since $\sum L^\dagger L=\gamma N$, the accepted return map is exactly

$$
\mathcal S_\gamma(\rho)=e^{-3\gamma t}K\rho K^\dagger.
$$

This particular loss model reduces success without changing the conditional state. Unequal rates need not have that property. The damping and loss channels are stated diagnostics, not fitted detector parameters.

For the complete-image estimator, one loss changes the parity weight to even and receives zero score, while two losses can imitate a singleton collision pattern. If $\lambda=\gamma t$ and each parity bit independently flips with known probability $q$, dividing scores by $e^{-3\lambda}(1-q)^4$ leaves bounded contributions of orders $\lambda^2$, $q^2$, and $\lambda q$. For the exactly calibrated nominal score, the [certificate](research/CERTIFICATE.md) and [analytical noise budget](research/NOISE_BUDGET.md) give the uniform sufficient window

$$
0\le\lambda\le0.01,\quad0\le q\le0.003:
\qquad |\operatorname{bias}|<0.001430,\qquad \mathbb EY^2<0.411.
$$

On the maximally mixed input the second-moment bound is $0.302$. Thus the mean-square error is bounded by $0.411/n+(0.001430)^2$ for arbitrary input, before adding calibration and rate uncertainty. These guarantees concern the stated loss and detector channels, and the residual bias must be included when comparing total estimation error.

If the implemented occupation effect obeys $\|F_{i,\pm}-S\|\le b_{i,\pm}$, compression after evolution gives reconstructed systematic bias at most $\sum_i|w_i|(b_{i,+}+b_{i,-})$. Hamiltonian parameters, timing, symmetry errors, and coefficient uncertainty require their own calibration budgets. Calibration copies and characterization time are additional resources.

### Direct calibration of the implemented channel

Spin-rotation invariance and paired reflection constrain the input effects even when atoms are lost or parity labels are mistaken. Include all outcomes, with no successful-trial postselection. Any signed two-flux image score then has mean operator $A=c\chi+dD$ and a reflection-even second-moment operator $T$ with three scalar blocks. The channel's noise rates need not be separately known.

Use three product references with the same orbital loading and measurement channel:

$$
\rho_x=|{+x},{+y},{+z}\rangle\langle{+x},{+y},{+z}|,
\quad \rho_d=|{+z},{+z},{-z}\rangle\langle{+z},{+z},{-z}|,
\quad \rho_0=|{+z},{+z},{+z}\rangle\langle{+z},{+z},{+z}|.
$$

Their score means determine $c=8m_x$ and $d=m_d$. Two image features with a nonsingular two-coordinate response can therefore be combined into an exactly unbiased estimator when that response is known. Their second moments $Q_j=\mathbb E_jX^2$ determine all three blocks:

$$
t_{\rm s}=Q_0,\qquad t_{\rm o}=3Q_d-2Q_x,\qquad
 t_{\rm e}=6Q_x-3Q_d-2Q_0.
$$

The same identities reconstruct every outcome effect within the five-dimensional symmetry class. They apply at arbitrary noise strength when the response remains informative; complete erasure cannot be inverted. At the nominal operating point, the [exact certificate](research/CERTIFICATE.md) constructs a second score for $D$ and proves that their full noisy response stays nonsingular throughout $0\le\gamma t\le0.01$, $0\le q\le0.003$ under the stated equal-loss and bit-error channels. Exact channel recalibration removes the mean bias and gives $\|T\|<0.444512137719<0.446<9/16$. The rates are used to establish this sufficient window; operational recalibration uses the raw reference responses, without those rates.

With finite calibration data, freeze a score chosen from pilot data and use independent reference data to validate it. Its exact all-input bias is

$$
\|A-\chi\|=\sqrt{3(c-1)^2/16+3d^2}.
$$

The [channel-calibration proof](research/CHANNEL_CALIBRATION.md) converts the two means and three second moments into simultaneous confidence bounds $\beta$ and $U_T$, giving mean-square error at most $U_T/n+\beta^2$ for every unknown input. It includes reference-preparation errors and the cost of known reference copies separately from unknown input copies. Three probes do not establish the symmetry assumptions themselves, and a random inverse of estimated calibration coefficients is not automatically unbiased.

## 7. Matched measurement resources

For $n$ independent input preparations and exact calibration, with ideal detection except in the stated noisy row, compare the same normalization of $\chi$. All failed return trials remain included in $n$.

| Method | Primitive and allocation | Worst-case coefficient or upper bound on $n\operatorname{Var}(\widehat\chi)$ |
|---|---|---:|
| Optimal fully separable measurements per input copy | Classical adaptation allowed; six uniformly randomized Pauli settings attain the minimax value | $9/16=0.5625$ exactly |
| Controlled-cycle Hadamard test | Ancilla quadrature after controlled cycle | $1/4=0.25$ |
| Spectral chirality measurement | Outcomes $0,\pm\sqrt3/4$; Fourier/count implementation follows [N26](REFERENCES.md) | $3/16=0.1875$ |
| Hard-core vacancy return | Two fluxes, equal allocation | $3/[16(1-r)^2]=0.1984038180\ldots$ |
| Finite-$U$ complete parity, $u=16$, $\tau=\pi$ | One duration, two equally randomized fluxes; certified corrected score | $<0.389156263913$ |
| Finite-$U$ channel-calibrated readout | Same nominal point; equal loss $\gamma t\le0.01$, independent bit errors $q\le0.003$ | $<0.444512137719$ |
| Finite-$U$ return-only example above | Four settings, allocation proportional to $|w_i|$ per sign | $(|w_1|+|w_2|)^2=1.0329363137\ldots$ |

The Pauli decomposition is $\chi=(XYZ+YZX+ZXY-XZY-YXZ-ZYX)/8$; it needs no entangling measurement gate. The Hadamard test uses $\operatorname{Im}\langle C\rangle=2\langle\chi\rangle$. These established alternatives and projective sector preparation are described in [R23](REFERENCES.md).

The separable benchmark is a lower bound for the whole declared measurement class. A normalized fully separable effect has chirality expectation at most $1/8$ in absolute value. With $\operatorname{Tr}\chi^2=3/4$, Cauchy–Schwarz yields $n\operatorname{Var}_{I_8/8}(\widehat\chi)\ge9/16$ for unbiased scoring. A classical score/Fisher-information argument extends the bound to arbitrary estimators locally unbiased around that input, even with classical adaptation between copies. Randomized Pauli scores $\pm3/4$ attain worst-case coefficient $9/16$. The [self-contained proof](research/SEPARABLE_BENCHMARK.md) establishes the exact minimax value without claiming a new general statistical principle; LOCC measurements within one copy are included [CLMOW14](REFERENCES.md).

The certified complete-parity bound therefore proves a strict minimax sampling advantage over universally unbiased protocols using fully separable effects across the three spins of each copy. It also improves the mixed-input coefficient from $9/16$ to the certified interval above. This does not assert pointwise superiority on every state, superiority to the controlled-cycle test, or a bound on collective measurements across copies or estimators exploiting a state promise.

For arbitrary biased estimators in the same per-copy measurement class, a uniform mean-square error at most $\varepsilon^2$ requires

$$
N\ge\left(\frac{3}{4\varepsilon}-\pi\right)^2,
\qquad 0<\varepsilon<\frac1{16\pi}.
$$

The [finite-sample proof](research/SEPARABLE_BENCHMARK.md#6-finite-sample-bound-allowing-biased-estimators) combines the submodel's separable Fisher-information bound with the scalar van Trees inequality. This is a necessary count, not an attained finite-sample minimax formula. It provides a bias-aware comparison for a fixed validated native estimator; the confidence of calibration must still be distinguished from unconditional estimation risk. The [calibration-cost account](research/CALIBRATION_COST.md) charges pilot and reference preparations explicitly and finds no established one-use total-preparation advantage. Its reuse calculations distinguish finite confidence certificates from regular plug-in asymptotics.

Signed internal-state information in spatial counts has optical antecedents [M17](REFERENCES.md). Fourier interferometry already reads cyclic observables with entangled inputs [N26](REFERENCES.md), and phase-sensitive scattering measurements isolate static Hubbard spin chirality [NGH25](REFERENCES.md). The implementation studied here combines interacting four-site vacancy dynamics, parity-only detection, and finite-repulsion calibration. Correlated inputs alone do not distinguish it from these general cycle-readout methods.

For the four-setting finite-$u=16$ return-only allocation, the coefficient on input $I_8/8$ is $0.8759834953$. Retaining complete parity images reduces this coefficient to the certified interval near $0.2848954942$, with one evolution duration. Both methods count all original preparations. The original return-only derivations remain in [the matched resource comparison, Section 6](source/resource-review-2026-10-10/READOUT_AND_RESOURCES.md); the complete-image and separable-benchmark proofs above specify the stronger comparison.

The vacancy implementation uses coherent spin-independent flux hopping, a localized motional ancilla, calibrated durations, and total site parity. Retained entanglement adds atom-preserving spin-coherent detection; repetition requires those properties after every cycle. Component experiments in [the references](REFERENCES.md) establish relevant imaging and flux capabilities separately. The mathematical result is the specified native instrument and its calibrated input readout, with these distinct operational requirements.
