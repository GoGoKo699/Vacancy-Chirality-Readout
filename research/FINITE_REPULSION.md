# Finite-repulsion vacancy filtering and calibrated chirality readout

## 1. Fixed model, actual accepted event, and two distinct tasks

Use three spin-1/2 bosons on four sites with

$$
H_U=T_\Phi+UW,\qquad
T_\Phi=-J\sum_{j,\sigma}(e^{i\Phi/4}b^\dagger_{j\sigma}b_{j+1,\sigma}+\mathrm{h.c.}),\qquad
W=\frac12\sum_jn_j(n_j-1).
$$

The two spin components have identical hopping and interactions. All four bonds are equal, and there is no spin-dependent potential. The full fixed-three-particle space has dimension 120; W has eigenvalues 0,1,3. Let P project onto all 32 no-double-occupancy configurations. Let V embed the eight spin states into the particular occupancy pattern (0,1,1,1), and S=VV†. The accepted Kraus operator is

$$K_{U,\Phi}(t)=V^\dagger e^{-itH_U}V. \qquad\text{(1)}$$

This definition presumes an ideal spin-preserving projection onto that occupancy pattern. Checking only that site 0 is empty is insufficient at finite U. The returned spins remain at sites 1,2,3. All other outcomes terminate the trial; no state-preservation claim is made for rejected configurations. An unread detector channel that resolves or dephases internal spins would change (1).

The t*=pi/J and Phi=−2pi/3 filter is paired with duration 3pi/(4J) on separate input copies for expectation-value readout. No spin rotations, controlled permutations, interaction-sign reversal, auxiliary spin gates, feedback during motion, or new couplings are added. Doubling the readout settings and calibrating their coefficients are explicit costs.

Task A is the quantum operation on the accepted state. Task B is estimating the *input* chirality expectation from independent-copy return probabilities. Correcting the classical estimate in Task B does not reverse the measurement back-action in Task A.

## 2. Uniform finite-U bound on the accepted quantum operation

Write tau=Jt and u=U/J. Let K∞ denote the hard-core return operator at the same time and flux. For all u>0 and tau>=0,

$$
\boxed{\|K_{U,\Phi}(t)-K_{\infty,\Phi}(t)\|\le
\epsilon_u(\tau):=\min\{2,20\sqrt3\,\tau/u\}.} \qquad\text{(2)}
$$

This deliberately conservative estimate is uniform over input spin states, including arbitrary external reference systems. It holds at every time, though it becomes uninformative when tau/u is large. It is not an optimal required U/J.

### Elementary proof without assuming that rare double occupancy means small spin error

On the full three-boson space, ||T||<=6J: the one-body four-ring hopping has norm at most 2J, and there are three bosons. In P the vacancy graph has weighted degree two, so ||PTP||<=2J. Put B=QTP, Q=I−P. Direct counting of the two-hop paths gives

$$\|B\|^2\le12J^2,\qquad \|TP\|^2=\|(PTP)^2+B^\dagger B\|\le16J^2. \qquad\text{(3)}$$

For clarity, the counting bound can be established at zero flux. The nonnegative matrix B†B/J² has row sum exactly 12 for every hard-core basis configuration. For spins around the three occupied sites, the diagonal contribution is 8 for three equal spins, 6 for an equal neighboring pair and an opposite third spin, and 4 for alternating spins along the occupied segment. The off-diagonal contributions are respectively 4,6,8. A flux can only reduce the absolute value of each entry relative to its zero-flux path sum. The maximum absolute row sum is consequently at most 12. The [finite-repulsion checker](../source/check_finite_repulsion.py) reproduces this with exact integer products of the squared one-hop occupation factors; no numerical spectral maximum is needed for (3).

For any normalized initial vector psi in P, energy conservation gives

$$\|H_Ue^{-itH_U}\psi\|=\|T\psi\|\le4J.$$

Since W>=Q and ||T||<=6J,

$$U\|Qe^{-itH_U}\psi\|\le U\|We^{-itH_U}\psi\|
\le\|H_Ue^{-itH_U}\psi\|+\|Te^{-itH_U}\psi\|\le10J.$$

This yields the all-time leakage bound

$$\boxed{\|Qe^{-itH_U}P\|\le\min\{1,10J/U\}.} \qquad\text{(4)}$$

The low-space equation differs from evolution under A=PTP by the forcing B†Qe^{-itH_U}P. Duhamel and (3)–(4) give

$$\|Pe^{-itH_U}P-e^{-itA}P\|
\le t\sqrt{12}J\,10J/U=20\sqrt3\,J^2t/U.$$

Compressing further with V proves (2). All spaces are finite, so no domain or particle-number truncation issue enters this proof. The probability of leaving P is O((J/U)²), whereas coherent accepted-operation error can be O(J²t/U): the former alone would not justify the latter.

### Completely bounded operation and postselection estimates

Define the trace-nonincreasing success map Kcal_U(rho)=K_U rho K_U†. Then

$$\boxed{\|\mathcal K_U-\mathcal K_\infty\|_\diamond\le2\epsilon_u.} \qquad\text{(5)}$$

This follows by adding and subtracting K∞rho K_U† and using that both Kraus operators are contractions, also after tensoring an identity reference map. If failure is represented by a terminal classical flag and its spin state is discarded, the complete trace-preserving flagged instrument satisfies

$$\boxed{\|\mathfrak I_U-\mathfrak I_\infty\|_\diamond\le\min\{2,4\epsilon_u\}.} \qquad\text{(6)}$$

The failure functional's difference is bounded by ||K_U†K_U−K∞†K∞||<=2epsilon_u. Equation (6) does not compare the unretained physical states of failed bosons.

For any specified possibly reference-entangled input, let p0 and pU be its ideal and finite-U success probabilities. Purify the input, and call the two unnormalized successful vectors x,y. Since ||x−y||<=epsilon_u,

$$|\sqrt{p_U}-\sqrt{p_0}|\le\epsilon_u,\qquad
\boxed{D(\rho_U^{\mathrm{succ}},\rho_0^{\mathrm{succ}})
\le\min\{1,\epsilon_u/\sqrt{p_0}\}.} \qquad\text{(7)}$$

The conditional statement requires pU>0. In particular epsilon_u<sqrt(p0) guarantees that. To prove the sharper distance bound, the distance from x to the complex line spanned by y is at most ||x−y||. Dividing by ||x|| bounds the sine of the pure-state angle by epsilon_u/sqrt(p0); tracing the purification cannot increase trace distance. This avoids treating a normalized conditional map as a linear channel or silently assuming a favorable acceptance probability.

For m retained returns, contraction gives ||K_U^m−K∞^m||<=min{2,m epsilon_u}; use this error in (5)–(7) with the *whole-record* ideal probability. Fixed m is controlled as U/J→infinity. No infinite-round finite-U guarantee follows.

### A non-asymptotic sufficient entanglement example

At tau=pi, the ideal maximally mixed input gives p0=(1+3r)/4=0.2709003550… and chi0=0.3884692918…. At u=2048, (2) gives epsilon_u=0.05313865325…. Thus

$$p_U\ge(\sqrt{p_0}-\epsilon_u)^2>0.2184,$$
$$\langle\chi\rangle_{\rho_U}\ge\chi_0-2\|\chi\|\epsilon_u/\sqrt{p_0}>0.3000>1/4.$$

This proves that genuine-tripartite-entanglement generation is not confined to infinite repulsion. The sufficient value 2048 is a conservative analytical example, not a necessary experimental ratio. Lower-U values discussed below are full finite-matrix diagnostics, not consequences of this loose bound. The scalar chirality witness itself is inherited from the spin-chirality literature and proved in the [hard-core account](HARD_CORE_FILTER.md).

## 3. The finite-U readout error has exactly one contaminating observable

Let R=P13 exchange spins at sites 1 and 3, and define

$$D=P_{12}-P_{23}=2(\mathbf S_1\cdot\mathbf S_2-\mathbf S_2\cdot\mathbf S_3),\quad
\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3).$$

The swap convention makes ||D||=sqrt3 and ||chi||=sqrt3/4. D is an imbalance of two pair correlations, not an additional three-spin handedness.

For one fixed time, define positive-chirality-labeled K+=K_{U,−2pi/3} and K−=K_{U,+2pi/3}, E±=K±†K±. At **every U and time** in the declared uniform model,

$$\boxed{E_+(t)-E_-(t)=g_U(t)\chi+h_U(t)D.} \qquad\text{(8)}$$

No other spin operator can occur in this difference. The real coefficients are exactly defined by the Hilbert–Schmidt projections onto chi and D, which are orthogonal.

Proof: the Hamiltonian and occupancy projector commute with collective spin rotations. On three qubits the commutant is the five-dimensional span of I, P12, P23, P13, chi. This is the established permutation/collective-unitary algebra [EW01](../REFERENCES.md). Spatial reflection fixes the source vacancy and the accepted pattern, interchanges sites 1 and 3, and reverses the flux. Consequently K−=R K+ R and the effect difference is odd under R. The only odd Hermitian combinations in the stated commutant are chi and P12−P23. That proves (8) without a 1/U expansion.

In the real occupation basis there is also reciprocity K−=K+^T. This does **not** imply E−=E+* unless K+ is normal, which is not generally true. Treating return probabilities as though this extra normality held would incorrectly discard the D term.

At U/J=16 and Jt=pi, complete-matrix evaluation gives

$$g_1=1.7649907114\ldots,\qquad h_1=0.02658895645\ldots.$$

Thus even after calibrating the chi coefficient, the single-time estimate Delta p/g1 has systematic error (h1/g1)<D>. Its worst possible magnitude is sqrt3 |h1/g1|≈0.02609. A real D eigenstate has zero chirality and attains this bias. Even the fully separable input |up,up,down> has <chi>=0, <D>=1 and produces a nonzero contrast h1. These are operator-level counterexamples, not merely changes in one generated state's fidelity.

Equation (8) needs reflection symmetry as well as spin-independent dynamics. Reflection-symmetric static scalar imperfections can be included with new coefficients; arbitrary scalar disorder does not generally preserve this two-dimensional form at finite U. The [hard-core disorder statement](HARD_CORE_FILTER.md#6-exact-structure-under-static-spin-independent-disorder) applies in its own limit; the finite-repulsion contrast requires the reflection symmetry stated here.

## 4. Exact expectation reconstruction with four native settings

For two durations t1,t2 on independent preparations, let

```math
d_i=p(-2\pi/3,t_i)-p(+2\pi/3,t_i),\quad
M=\begin{pmatrix}g_1&h_1\\g_2&h_2\end{pmatrix}.
```

Whenever DetM=g1 h2−g2 h1 is nonzero, the exact identity is

$$\boxed{\langle\chi\rangle=
\frac{h_2d_1-h_1d_2}{g_1h_2-g_2h_1}=w_1d_1+w_2d_2.} \qquad\text{(9)}$$

It holds for every unknown spin state, including a subsystem of a larger entangled state. It uses two times times two flux signs—four occupancy measurement configurations. No measurement of D, spin rotation, or controlled-SWAP is needed during these runs. One observes frequencies on independent copies; this is not a single-shot chirality assignment.

### Two suitable durations exist at every fixed positive finite repulsion

The nonsingularity condition is not only a numerical possibility. At the fixed flux magnitude 2pi/3, the exact short-time series at any fixed u>0 is

$$
g_u(\tau)=4\sqrt3\,\tau^4+O_u(\tau^6),\qquad
h_u(\tau)=\sqrt3\left(\frac{7u}{72}+\frac{u^3}{2520}\right)\tau^9+O_u(\tau^{11}). \qquad\text{(11)}
$$

Here g and h are functions of dimensionless time tau=Jt. Reciprocity implies g is even and h is odd in time. The displayed coefficients are exact finite path sums, not floating-point fits. For two durations tau=s and tau=2s, the determinant in (9) consequently obeys

$$
\det M=5952\left(\frac{7u}{72}+\frac{u^3}{2520}\right)s^{13}+O_u(s^{15})>0
$$

for all sufficiently small positive s. Thus **for every fixed finite u>0 there are two durations whose flux-reversed contrasts identify both input observables**. This statement grants no spin controls. It also does not assert that the particular order-one durations below work at every U, or that very small durations are statistically efficient: their chirality contrast vanishes as s^4. The expansion holds at fixed u as time tends to zero, not uniformly in a joint large-u/time limit.

A reproducible algebraic certificate for (11) is [exact_short_time.py](../source/exact_short_time.py). In the complete N_up=2,N_down=1 sector, use the 40 occupation basis states. All entries of 2H/J are integer polynomials in u with coefficients in Q(sqrt2,sqrt3,i). Set G_0=V and G_(n+1)=(2H/J)G_n, retaining powers of u explicitly. Then the coefficient of tau^ell in E is exactly

$$
\frac1{2^\ell\ell!}\sum_{n=0}^{\ell}\binom\ell n
 i^n(-i)^{\ell-n}(V^\dagger G_n)^\dagger(V^\dagger G_{\ell-n}).
$$

Subtract its reflection conjugate. In the ordered input basis (up,up,down), (up,down,up), (down,up,up), D_00=1 and chi_01=i/4, so the diagonal (0,0) entry gives h and four times the (0,1) imaginary part gives g. The recurrence gives zero h coefficients through order eight and precisely the two u and u^3 coefficients at order nine in (11). It gives zero g coefficients below order four and 4sqrt3 at order four. Exact integer and rational arithmetic, including the multiplication rules (sqrt2)^2=2, (sqrt3)^2=3, i^2=−1, verifies every entry. Finite-dimensional analyticity and the even/odd symmetry supply the stated remainders. This is an exact coefficient certificate, not interval validation of the longer-time example.

For the same finite ratio u=16, choose tau1=pi and tau2=3pi/4. The complete 120-state calculation gives

| tau | g | h |
|---|---:|---:|
| pi | 1.76499071140 | 0.026588956448 |
| 3pi/4 | −0.39157486049 | 0.020164887920 |

The determinant is approximately 0.04600240679, and

$$\langle\chi\rangle\simeq0.43834419389\,d_1-0.57799055104\,d_2.$$

The unrounded coefficients are defined by the exact finite-matrix functions in (8). The displayed rounded values give an approximate numerical formula, not an exact decimal identity. Independent Taylor scaling/squaring and sparse exponential evolution agree with the spectral computation at the 10^-14 level. These are not interval-certified parameter values or a calibration of a real device. No uniform determinant lower bound for all U is asserted. A singular or poorly conditioned choice of two times must be rejected; the helper does so rather than inverting numerical zero. Zero flux provides no such odd-observable calibration.

### Statistical and calibration costs

For n independent input trials divided among the four settings, use counts in proportions |w_i| for each sign of time i. The Bernoulli variance bound yields

$$\mathrm{Var}(\widehat\chi)\le (|w_1|+|w_2|)^2/n.$$

The numerical u=16 example gives a numerator 1.0329363…. This is a worst-case bound and an allocation statement, not a measured variance or a sample-complexity advantage over the circuit predecessor. The ideal hard-core two-setting numerator would be [sqrt3/(4(1−r))]^2≈0.1984; the four-setting calibration has a statistical cost.

Calibration is not assumed cost-free. If the implemented contrasts have operator errors ||Delta E_i−(g_i chi+h_i D)||<=b_i, then the reconstructed expectation has systematic error at most

$$\boxed{|w_1|b_1+|w_2|b_2.} \qquad\text{(10)}$$

Parameter uncertainty, asymmetric bonds, or imperfect readout must be included in those budgets. Knowing that a device is SU(2)-invariant is insufficient to know g_i,h_i. The finite Hamiltonian, timing, and measurement must be known or independently characterized; no hidden spin-resolved calibration apparatus is granted here. If g_i,h_i estimates have errors, an additional budget |Delta g_i|sqrt3/4+|Delta h_i|sqrt3 can be used in (10).

### Why this is not an exact quantum-nondemolition instrument

Equation (9) is a linear combination of ensemble frequencies, with a negative weight in the example. It is not one physical CP map and cannot alter K_U's back-action. The finite-U return operation generally does not commute with chi. Thus calibrated input-expectation readout can be exact while the accepted state is disturbed and while repeated filtering fails to retain the hard-core success plateau. These are distinct operational claims.

## 5. Repeated filtering at finite U: high chirality has a survival cost

The exact successful record remains K_U^m rho (K_U†)^m. For the maximally mixed input at u=16, tau=pi:

| Successful rounds m | Probability per original trial | Target-sector weight | Output chirality |
|---|---:|---:|---:|
| 1 | 0.2716923706 | 0.7571463353 | 0.3045132972 |
| 4 | 0.1141865390 | 0.9991250391 | 0.4324656707 |
| 16 | 0.01081184194 | 0.9996427368 | 0.4327033029 |
| 64 | 8.70860e−7 | 0.9996427368 | 0.4327033029 |

The maximal chirality eigenvalue is sqrt3/4≈0.4330127. Four returns therefore give a high-chirality genuinely tripartite-entangled **mixed** state with about 11.42% success, not the ideal 25% plateau. SU(2) symmetry leaves the internal doublet mixed for this input. The numerical spectral radius of K_U at u=16 is 0.9064726…, and its operator norm is 0.9074959…. For any exact contraction with norm below one, survival is bounded by ||K_U||^(2m); the finite-matrix diagnostics demonstrate the corresponding decay here. A theorem excluding accidental perfect-return parameters at every finite U is not claimed.

At u=128 four returns give success 0.24733838… and target weight 0.9999955863…: the four-round ideal limit is approached when U increases at fixed round count. At fixed finite U one cannot freely send the number of rounds to infinity and retain the ideal success yield. The all-input bound in Section 2 explicitly accounts for the round count.

These numbers use the complete three-particle sector. They do not model detector disturbance, losses, or parameter uncertainty, and require spin-preserving nondestructive occupancy readout between successful rounds. Destructive number detection can estimate independent-copy probabilities in (9), but cannot retain or repeatedly filter the entangled output.

## 6. First-order interpretation

Let A=PTP and B=QTP. Standard finite-dimensional block elimination [BDL11](../REFERENCES.md), on the isolated W=0 spectral group as J/U→0, gives

$$H_{\mathrm{eff}}=A-U^{-1}B^\dagger W_Q^{-1}B+O(J^3/U^2).$$

Here B reaches only a single doublon, so W_Q^{-1}B=B. Projected endpoint dressing is O((J/U)²). At fixed tau, this yields

$$K_U=K_\infty+\frac{i}{U}V^\dagger\int_0^t e^{-iA(t-s)}B^\dagger B e^{-iAs}ds\,V
+O_\tau((J/U)^2).$$

The high-space phases may oscillate rapidly with U, but their accepted-block amplitude is of second order. This is a finite-time expansion, not a uniform-in-round-count assertion. The helper computes the first correction as a finite spectral sum; independent tests check the residual scaling. These tests do not supply an explicit second-order remainder constant. The nonperturbative bound (2) is the quantitative theorem used for arbitrary inputs.

The correction represents ordinary virtual particle crossings and exchange. It explains why small doublon occupation does not by itself make the accepted spin map ideal. The Schrieffer-Wolff method, the spin-permutation algebra, and linear calibration are inherited mathematical tools; none is a new general principle here.

## 7. Attribution and verification

The [reference account](../REFERENCES.md) attributes the spin-cycle flux to Kesharpu and Kochetov [K22](../REFERENCES.md), chirality measurement and sector preparation to Reascos et al. [R23](../REFERENCES.md), collective-unitary permutation algebra to Eggeling and Werner [EW01](../REFERENCES.md), and block elimination to Bravyi et al. [BDL11](../REFERENCES.md). The bounds here specify the accepted operation and calibrated expectation reconstruction for the declared spin-independent, reflection-symmetric ring. Retained outputs also require spin-preserving occupancy readout.

The [finite-repulsion model](../source/finite_repulsion.py) implements the complete 120-state bosonic sector. The [independent checker](../source/check_finite_repulsion.py) runs five groups. Its 512-dimensional first-quantized construction is projected onto that sector as a separate Hamiltonian check; spectral evolution is compared with scaling/squaring Taylor evolution and sparse exponential action. The Taylor remainder is an exact-arithmetic truncation estimate, with floating-point roundoff checked separately rather than included in an interval certificate.

The [short-time coefficient checker](../source/exact_short_time.py) uses integer and rational arithmetic for the coefficients in (11). See [Reproduce](../REPRODUCE.md) for the full suite, report identity checks, and the distinction between analytical certificates and numerical diagnostics.
