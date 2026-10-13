# From many-particle interference to chirality readout

**One external tutorial, followed by a bridge to this repository.** The goal
is to understand how spin-independent motion and spin-insensitive images
can estimate a three-spin observable, and why interactions improve precision
when spatial and detector resources are fixed.

The external foundation is Mattia Walschaers,
[*Signatures of Many-Particle Interference*](https://arxiv.org/abs/1908.08370).
Section and equation numbers below refer to
[arXiv:1908.08370v1](https://arxiv.org/html/1908.08370v1).
Basic density matrices, tensor products, and the Born rule are sufficient
starting knowledge. The full project statements remain in [Theory](THEORY.md);
this guide explains how to approach their proofs.

## 1. Read the tutorial with a destination

Read **§§2.1–2.3 → §§3.1–3.3 → §4.1**, focusing on Eqs. (138)–(142).
Their worked input has product internal states and distinct input and output
ports; §4.1 assumes fully indistinguishable particles and omits internal
degrees of freedom.

Our input spins may be correlated, and parity groups together outcomes with
collisions. The operator treatment below supplies these extensions and leads
to the repository's precision comparison. The guiding questions are: which
spin information survives the image, how can its mean be reconstructed,
and how many preparations are needed?

## 2. Translate the notation

| Object | Meaning here |
|---|---|
| Spatial space $`\mathcal H_E`$ | Four modes, labeled $`0,1,2,3`$ |
| Internal space $`\mathcal H_I`$ | One spin-½ space $`\mathbb C^2`$ |
| $`a^\dagger(e_j\otimes f_\sigma)=b^\dagger_{j\sigma}`$ | Creation of spin $`\sigma`$ in spatial mode $`j`$ |
| Spatial unitary $`U`$ | Passive comparator $`U\in U(4)`$ |
| Detector output | Total occupation parity at each site |

Here $`U\in U(4)`$ must be distinguished from the interaction energy $`U`$ in
the Bose–Hubbard Hamiltonian. The ratio $`u=U/J`$ always denotes interaction
strength relative to hopping.

Initially site 0 is empty and each other site holds one boson. The embedding

```math
V|s_1s_2s_3\rangle
=b^\dagger_{1s_1}b^\dagger_{2s_2}b^\dagger_{3s_3}|\mathrm{vac}\rangle
```

identifies the input with three site-labeled qubits. Bosonic symmetry is
already enforced by the creation operators; it does not require the
three-qubit density matrix $`\rho`$ to be permutation-symmetric.

## 3. Replace overlap products by spin operators

A passive device transforms creation operators with the same spatial
matrix for both spin components. Matching three input positions to three
output positions produces permutations of the spin labels. For a product
spin vector $`|\psi_1\psi_2\psi_3\rangle`$, their expectations are products
of overlaps. For an arbitrary correlated input, the quantity to retain is
instead $`\operatorname{Tr}(\rho P_\pi)`$, where

```math
P_\pi|s_1s_2s_3\rangle
=|s_{\pi(1)}s_{\pi(2)}s_{\pi(3)}\rangle.
```

For example, compressing a product of output parities on a set $`\mathcal S`$
back to the occupied input modes gives

```math
\begin{aligned}
W_{\mathcal S}&=U^\dagger\operatorname{diag}
\bigl((-1)^{1_{j\in\mathcal S}}\bigr)U,\\
B_{\mathcal S}&=\sum_{\pi\in S_3}
\left[\prod_{i=1}^3(K_{\mathcal S})_{i,\pi(i)}\right]P_\pi,
\end{aligned}
```

where $`K_{\mathcal S}`$ is the occupied-input $`3\times3`$ submatrix of
$`W_{\mathcal S}`$. This is an operator identity before any input is chosen.
Consequently it covers mixtures, entangled spins, and inputs entangled with
an external reference. The creation-operator derivation is in
[the passive proof, §1](research/PASSIVE_PARITY_BOUND.md#1-compressing-a-spatial-parity-observable).

Chirality is the part sensitive to the orientation of a three-cycle. Our
conventions are

```math
C|s_1s_2s_3\rangle=|s_2s_3s_1\rangle,\qquad
\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3),\qquad
\mathbf S=\boldsymbol\sigma/2,
```

```math
C-C^\dagger=4i\chi.
```

Thus the two cycle terms in $`B_{\mathcal S}`$ contribute a chirality
coefficient $`-4\operatorname{Im}(K_{12}K_{23}K_{31})`$. Pair exchanges alone
cannot supply this orientation-sensitive component. The sign of the
chirality coefficient depends on the stated cycle convention; $`\chi`$ has
eigenvalues $`0,\pm\sqrt3/4`$.

## 4. Understand what the parity image retains

Number detection reports an occupation vector $`\mathbf n`$. Parity reports
only $`p_j=n_j\bmod2`$. If $`E_{\mathbf n}`$ is the count effect, then

```math
E_p=\sum_{\mathbf n:\,\mathbf n\bmod2=p}E_{\mathbf n}.
```

This sum can erase interference information. With exactly three particles,
four-mode images have one or three odd bits. The pattern $`(0,1,1,1)`$ fixes
the occupation exactly: three odd sites already require all three particles.
A singleton pattern instead includes collision possibilities.

For three passive modes, every parity effect has zero chirality component.
The two opposite chiral-sector states therefore remain indistinguishable
regardless of the number of fresh copies. A fourth, initially empty mode
allows informative two-output parity correlations. The
[passive proof, §§2–3 and §7](research/PASSIVE_PARITY_BOUND.md)
establishes the obstruction and an explicit four-mode reconstruction.
The vacancy enables access; access alone does not determine precision.

## 5. Convert an image into an estimator

For the interacting ring, use $`u=16`$, $`\tau=Jt=\pi`$, and fluxes
$`\mp2\pi/3`$. Define $`E_p`$ at the negative flux, $`R=P_{13}`$, and
$`r(p_0,p_1,p_2,p_3)=(p_0,p_3,p_2,p_1)`$. Reflection identifies the
positive-flux effect for the **physical pattern $`r(p)`$** with $`RE_pR`$.
Spin symmetry then gives

```math
O_p=\frac{E_p-RE_pR}{2}=c_p\chi+d_pD,
\qquad D=P_{12}-P_{23}.
```

Choose the flux independently with equal probabilities. Report $`a_p`$ for
the negative-flux pattern $`p`$, and $`-a_p`$ for the positive-flux pattern
$`r(p)`$. The two calibration conditions

```math
\sum_pa_pc_p=1,\qquad \sum_pa_pd_p=0
```

make the mean score equal to $`\langle\chi\rangle`$ for every input.
The second condition removes a pair-correlation nuisance. All preparations
count, including zero-score outcomes.

With $`F_p=(E_p+RE_pR)/2`$ and $`T=\sum_pa_p^2F_p`$, the sample mean obeys

```math
n\operatorname{Var}_\rho(\widehat\chi)
=\operatorname{Tr}(T\rho)-\langle\chi\rangle_\rho^2.
```

The [complete-parity proof](research/COMPLETE_PARITY_THEORY.md) shows that
the worst case equals $`\|T\|`$: a state in a largest-eigenvalue reflection
block has zero chirality. The [exact certificate](research/CERTIFICATE.md)
constructs feasible scores with $`\|T\|<0.328697776905`$. Their objective is
worst-case precision; the separately documented mixed-input score optimizes
a different quantity. Neither is claimed globally optimal over devices.

## 6. See what makes the comparison universal

Testing one passive interferometer would only compare two devices. The
proof instead bounds every permitted $`U(4)`$. It uses the physical spin
family

```math
\rho_\theta=\frac{\Pi_{\rm m}}4+\frac43\theta\chi,
\qquad |\theta|<h=\frac{\sqrt3}{4},
```

where $`\Pi_{\rm m}`$ projects onto the two chiral sectors. Its chirality is
$`\theta`$, and each outcome probability has the form $`q_p+\theta c_p`$.
Unitarity constrains the chirality coefficients; positivity constrains
$`q_p\ge h|c_p|`$; normalization fixes $`\sum_pq_p=1`$. Together these imply
the per-copy Fisher-information bound $`F(0)\le32/13`$.

For universally unbiased estimators, the resulting worst-case variance
coefficient is at least $`13/32=0.40625`$. Classical adaptation and independent
spin prerotations are included. The inequality

```math
0.328697776905<13/32
```

therefore separates an interacting construction from a measurement class.
For biased estimators the appropriate statement uses mean-square error;
[§6 of the passive proof](research/PASSIVE_PARITY_BOUND.md#6-finite-mean-square-error-including-biased-estimators)
provides the finite bound and the asymptotic minimax bound. The unbiased
variance constant alone does not establish a finite-sample biased comparison.

The resource boundary matters: three bosons, four spatial modes, one vacuum
input, terminal spin-insensitive parity, and fresh copies. Extra modes,
intermediate detection, collective-copy measurements, and interleaved spin
controls change the class. An explicit eight-mode passive construction has
coefficient $`27/88`$, below the present interacting score.

## 7. Keep implementation guarantees attached to their assumptions

The [parameter theorem](research/PARAMETER_ROBUSTNESS.md) distinguishes
exact recalibration at varied parameters from keeping nominal scores and
budgeting bias. Its joint theorem combines specified parameter variation,
equal-rate loss, and independent detector-bit errors. General
symmetry-breaking errors use the fixed-score bound.

[Finite calibration](research/FINITE_CALIBRATION.md) includes random
reference data in unconditional MSE and counts reference preparations.
Its sufficient examples require substantial reuse of a stable instrument;
they establish a theoretical resource guarantee, not measured throughput.
Expectation estimation also does not require the spin-coherent back-action
needed to retain a useful postmeasurement state. For attribution and
reproduction, continue to [References](REFERENCES.md) and
[Reproduce](REPRODUCE.md).

## Worked self-checks

<details>
<summary>What chirality corresponds to a cycle eigenvalue exp(2πi/3)?</summary>

Using $`\chi=(C-C^\dagger)/(4i)`$ gives
$`\chi=\sin(2\pi/3)/2=\sqrt3/4`$. Reversing the cycle convention reverses
the chirality sign.

</details>

<details>
<summary>Why can the pattern (0,1,1,1) certify a return, while (1,0,0,0) cannot?</summary>

With at most three particles, the first pattern requires occupation
$`(0,1,1,1)`$. At exactly three particles, the second is compatible with
$`(3,0,0,0)`$ and $`(1,2,0,0)`$, among others. The particle-number assumption
is essential: five particles could produce $`(0,3,1,1)`$ with the return
parity pattern.

</details>

<details>
<summary>If calibration leaves mean χ + εD, does increasing the target sample size remove the error?</summary>

No. On input $`\rho`$, the bias is $`\epsilon\langle D\rangle_\rho`$.
Its worst-case magnitude is $`\sqrt3|\epsilon|`$ because $`\|D\|=\sqrt3`$.
More independent target copies reduce variance, while this fixed bias
remains. This is why noise and finite-calibration comparisons use MSE.

</details>
