# Vacancy chirality: parity readout, retained-state requirements, and matched resources

## 1. Model and measurement tasks

Three spin-1/2 bosons start at sites 1,2,3 of a four-site ring, with site 0 empty. The initial spin state rho is arbitrary, including a subsystem of a larger entangled state. The spin-independent Bose-Hubbard Hamiltonian is

$$
H_U=-J\sum_{j,\sigma}(e^{i\Phi/4}b^\dagger_{j\sigma}b_{j+1,\sigma}+\mathrm{h.c.})
+\frac U2\sum_j n_j(n_j-1).
$$

Let V embed the eight spin states into occupancy (0,1,1,1), and S=VV†. The ideal accepted spin operator is K=V† exp(-itH_U)V. At finite U, the exact uniform-ring effect difference between fluxes -2pi/3 and +2pi/3 is g(t) chi+h(t) D, with D=P12-P23 and chi=S1·(S2×S3). Two calibrated durations cancel D when their 2x2 coefficient matrix is invertible. The independent-copy input-expectation estimate is not a correction to the accepted quantum back-action. The [finite-repulsion proof](FINITE_REPULSION.md) establishes these identities.

There are three distinct tasks: estimate input chirality; herald and retain an entangled spin state; repeat the filter on that same surviving state. Their detector requirements differ. The resource comparison below uses the detector requirements appropriate to each task.

## 2. Exact parity readout suffices for the accepted event

Write the total site occupation, summed over both spin states, as n_j. Define the four-bit parity projector

$$
Q_{0111}=\prod_{j=0}^{3}\frac{I+(-1)^{s_j}(-1)^{n_j}}2,
\qquad (s_0,s_1,s_2,s_3)=(0,1,1,1).
$$

**Parity lemma.** On the entire particle-number subspace N<=3,

$$
\boxed{Q_{0111}=S,}
$$

where S is zero on N<3. In particular the equality holds on the full 120-dimensional N=3 space, including every doublon and triplon configuration.

Proof: an odd occupation at each of sites 1,2,3 requires at least three particles. With N<=3, all three occupations must equal one and n_0=0. Conversely, (0,1,1,1) has precisely the required parity. This proves an equality of projectors, hence also of outcome effects for arbitrary coherent superpositions and reference-entangled inputs. No assumption of negligible doublon probability enters.

Consequently,

$$
V^\dagger e^{itH_U}Q_{0111}e^{-itH_U}V=K^\dagger K.
$$

The four finite-U calibrated measurements can therefore use ideal total-occupation parity images, rather than full atom-number resolution. Tunneling must be frozen or otherwise stopped during a readout modeled this way. Light-assisted pair removal can implement parity projection without conserving particle number during the image; the measured bits must reliably represent the pre-imaging site parities. Neither an empty-site report alone nor a spin-selective fluorescence image is the same effect.

The number restriction is essential. For example, (0,3,1,1) at N=5 has the accepted parity pattern but is not the accepted occupancy. Preparation errors with additional atoms or gains can therefore cause false acceptance. Loss alone cannot do so: N<3 cannot produce three odd occupied sites. Detection errors are a separate issue, treated below.

Parity projection in fluorescence quantum-gas imaging is established [B24](../REFERENCES.md). The contribution here is this elementary use of the known total particle number to implement the protocol's exact accepted effect. It is not a new parity-imaging method or an assertion that a particular microscope realizes the complete spinful flux ring.

## 3. Same probabilities, different retained states

Suppose a readout has the correct accepted effect but disturbs the retained spins by a trace-preserving channel N. Its success map is

$$
\mathcal S_N(\rho)=\mathcal N(K\rho K^\dagger).
$$

For every input, Tr S_N(rho)=Tr K rho K†. This remains true with an external reference, since N†(I)=I. Therefore all independent-copy return frequencies and the input-chirality estimator can be unchanged while the accepted output state is very different.

For a concrete counterexample, take complete local dephasing in the computational spin basis:

$$
\mathcal D_0^{\otimes3}(X)=\sum_{s\in\{0,1\}^3}|s\rangle\langle s|X|s\rangle\langle s|.
$$

The output is a positive mixture of local product basis states and hence fully separable, for every accepted spin input. Its return probability is exactly the ideal one. This does not change the unknown input being estimated; the dephasing occurs during or after the terminal detection, after the return effect has acted.

Thus input-expectation readout can use a destructive final image. Retained entanglement requires a sufficiently coherent accepted instrument, and repeated filtering additionally requires atom survival, stable positions, and coherent state retention between rounds. A high atom-retention probability or a small spin-flip probability does not alone establish preserved spin superpositions. A readout that distinguishes spins can dephase them even without flipping them.

### 3.1 A quantified coherence requirement for the ideal filter

As a controlled diagnostic, not an apparatus noise fit, let D_eta be independent phase damping on each spin: diagonal entries are unchanged, and single-qubit off-diagonal entries are multiplied by eta in [0,1]. Every Pauli term in chi contains exactly two transverse factors, so

$$
(\mathcal D_\eta^{\otimes3})^\dagger(\chi)=\eta^2\chi.
$$

For the hard-core one-return maximally mixed input, the chirality chi_0=0.388469291814954... becomes eta² chi_0. The witness |<chi>|>1/4 certifies genuine tripartite entanglement exactly when

$$
\boxed{\eta>\sqrt{1/(4\chi_0)}=0.8022166239\ldots.}
$$

This is a threshold for this witness under this stated phase-damping channel, not a necessary and sufficient entanglement threshold for all detector noise. The success probability remains 0.2709003550... at every eta. These equal success probabilities alone do not distinguish the coherent entangler from the fully dephased separable output.

### 3.2 Repetition: exact recurrence for readout dephasing

Keep the ideal hard-core K=Pi_+ +a(I-Pi_+), r=a²=0.02786714004..., and apply D_eta independently after every accepted return. This intentionally isolates detector disturbance; it does not replace the separate finite-U dynamics.

In each three-dimensional one-up or two-up magnetization block, cyclic symmetry leaves an unnormalized target-cycle weight x_m and equal weights y_m on each of the other two cycle eigenstates. Initially x_0=y_0=1/8. Off-diagonal computational-basis entries within a block differ on two spins, so D_eta is multiplication by lambda=eta² on those entries. It follows exactly that

$$
\begin{pmatrix}x_{m+1}\\y_{m+1}\end{pmatrix}
=\begin{pmatrix}
(1+2\lambda)/3&2r(1-\lambda)/3\\
(1-\lambda)/3&r(2+\lambda)/3
\end{pmatrix}
\begin{pmatrix}x_m\\y_m\end{pmatrix}.
$$

The two fully polarized states contribute r^m/4 in total. Hence

$$
p_m=2(x_m+2y_m)+r^m/4,\qquad
\langle\chi\rangle_m=\frac{\sqrt3(x_m-y_m)}{2p_m}.
$$

For every eta<1, the positive 2x2 recurrence has spectral radius strictly below one. To prove this without a numerical scan, set X=x and Y=2y. The resulting nonnegative matrix has column sums 1 and r<1 and a strictly positive transition from X to Y. Its square has both column sums strictly below one. Its induced 1-norm is therefore strictly below one after two steps, so its spectral radius is below one. Consequently p_m tends to zero. At eta=1 the target column never leaks, giving the ideal 1/4 plateau.

For eta=.9, four successful rounds give p_4=0.1680156... and <chi>_4=0.3485571..., compared with the ideal p_4=0.25000045... and <chi>_4=0.4330117.... These are evaluations of the exact recurrence, not predicted parameters of a detector. The calculation identifies an independent reason that atom retention alone does not justify arbitrarily repeated ideal filtering.

## 4. Equal-rate particle loss: rejection without conditional spin bias

A second explicit diagnostic concerns particle loss during the hopping stage. Let the Markovian jump operators be L_(j,sigma)=sqrt(gamma) b_(j,sigma), with the same rate for every site and spin, no atom gains, and no other error. Then sum L†L=gamma N and H_U preserves N.

Every trajectory with at least one loss ends in N<3 and is rejected by ideal Q_0111. In the surviving N=3 trajectory the non-Hermitian term is the scalar -3i gamma/2. Therefore the accepted map is exactly

$$
\boxed{\mathcal S_\gamma(\rho)=e^{-3\gamma t}K\rho K^\dagger.}
$$

This is valid for all inputs and references. Conditional accepted states are unchanged by this specific loss model, but success per original trial is reduced. It is not a claim that loss is harmless in general. Unequal site/spin rates introduce a nonscalar no-jump evolution, and gains or detector errors can invalidate the rejection argument.

For independent-copy chirality measurements, the two-time contrasts become d_i^obs=s_i d_i with s_i=exp(-3 gamma t_i). Known survival factors may be divided out, replacing w_i by w_i/s_i. A conservative variance bound is then (sum |w_i|/s_i)²/n. This costs samples and survival calibration. It does not turn a retained-state measurement with spin dephasing into a coherent one.

The check integrates the complete number-block Lindblad equation in N=0,1,2,3, rather than substituting the no-jump answer. The Hilbert-sector dimensions are 1,8,36,120, and the block-diagonal density-matrix evolution has 15,761 components. These finite diagnostics test the stated channel, not an experimental loss rate.

## 5. Detector errors and calibration costs

If F_s is the implemented accepted effect on occupations and ||F_s-S||<=b_s, compression after unitary evolution gives

$$
\|V^\dagger U_t^\dagger(F_s-S)U_tV\|\le b_s.
$$

For a reconstructed estimate sum_i w_i(p_(i,+)-p_(i,-)), the unknown-input systematic bias is consequently bounded by

$$
\boxed{\sum_i |w_i|(b_{i,+}+b_{i,-}).}
$$

This is ordinary operator-error propagation, not a calibration-free promise. For example, independent final parity-bit flips at probability .003 give a finite-U16 chirality operator bias about .00517666 in the displayed four-setting implementation, below the conservative .02428249 bound. This is a diagnostic of that classical confusion channel only.

Hamiltonian parameters and timing still need independent characterization or valid uncertainty budgets. The reconstruction below assumes independently characterized response coefficients; its four unknown-input frequencies are not used as calibration data. Calibration copies and characterization time are excluded from the ideal sampling comparison below; they are not silently free resources in a hardware claim.

## 6. Compare the correct tasks and input-trial normalization

The strongest baseline is not only a controlled-SWAP circuit. The chirality operator already has the direct local-spin decomposition

$$
\chi=\tfrac18(XYZ+YZX+ZXY-XZY-YXZ-ZYX).
$$

This measures six three-spin Pauli products using separately chosen local axes. It requires no entangling measurement gate. Reascos et al. [R23](../REFERENCES.md) explicitly compare this strategy with their controlled-cycle Hadamard test and also describe a qutrit phase-estimation measurement of chirality. The vacancy protocol trades those spin-addressing/circuit resources for orbital flux, coherent hopping, a localized motional ancilla, and site-occupation readout. It is not ancilla-free merely because it adds no spin qubit.

### 6.1 Matched ideal sampling accounting

Use the same chi=S1·(S2×S3), n independent preparations, exact calibration, ideal detection, and each method's stated shot allocation. No failed return trial is removed from n.

| Method | Input settings/primitive | Worst-case bound on n Var(estimated chi) |
|---|---|---:|
| Six Pauli products, equal allocation | Six local-axis configurations and spin-resolved readout | 9/16 = .5625 |
| Controlled-cycle Hadamard test | One quadrature of an ancilla after the controlled cycle | 1/4 = .25 |
| Ideal spectral chirality measurement | Outcomes 0,+sqrt3/4,-sqrt3/4 | 3/16 = .1875 |
| Hard-core vacancy return | Two opposite fluxes at one duration | 3/[16(1-r)²] = .1984038180... |
| Finite-U16 vacancy return, two-duration design | Four parity settings and calibrated weights | $`(\lvert w_1\rvert+\lvert w_2\rvert)^2=1.0329363137\ldots`$ |

Derivation: each Pauli estimator has outcomes +/-1 and coefficient +/-1/8. Giving n/6 trials to each yields n Var<=36/64=9/16, saturated by I8/8. In the Hadamard test Im<C>=2<chi>; measuring its ancilla quadrature therefore uses an outcome +/-1/2 (up to orientation), giving the 1/4 bound. A spectral chi measurement has bounded outcomes +/-sqrt3/4, giving 3/16. These rescale the normalization chi_hat=4chi/sqrt3 in [R23](../REFERENCES.md); they are not new measurement theorems.

For the hard-core protocol set c=sqrt3/[4(1-r)]. Equal halves of the trials on each sign give n Var<=c². The bound is attainable: choose equal sector probabilities q_+=q_-=(1/2-r)/(1-r) and q_0=r/(1-r), making both Bernoulli return probabilities 1/2. These q values are nonnegative and sum to one.

For the finite-U protocol, allocate n_(i,+)=n_(i,-)=n |w_i|/[2W], W=|w1|+|w2|. Then

$$
n\operatorname{Var}(\widehat\chi)=2W\sum_i|w_i|[p_{i,+}(1-p_{i,+})+p_{i,-}(1-p_{i,-})]\le W^2.
$$

This last upper bound alone would not prove poorer performance. On the explicit input I8/8, however, the full finite-U16 matrices give p_(1,+)=p_(1,-)=.271692370645... and p_(2,+)=p_(2,-)=.334877621834.... The actual variance coefficient for this allocation is .875983495266..., larger than .5625 for direct Pauli products and .25 for the Hadamard test on the same input. Thus this finite-U parameter/estimator choice has no sample advantage in that matched example. No optimality claim across all vacancy times, fluxes, full-image statistics, or adaptive designs is made.

At fixed ideal resources the hard-core pair of flux measurements has a smaller worst-case coefficient than six independently sampled Pauli products, but it is not a universal physical advantage: the primitives, calibration, cycle duration, and attainable error rates differ. The projective comparator also outputs an actual chirality eigensector and accepts positive chirality from I8/8 with probability 1/4. Our ideal repeated-filter plateau is therefore not a new fundamental yield advantage.

### 6.2 Minimal operational comparison

| Task | Vacancy implementation needs | What cannot be inferred |
|---|---|---|
| Estimate input chi on independent copies | Known three-particle spinful preparation; calibrated spin-independent/reflection-symmetric flux hopping; two durations/two flux signs; final total site parity | Coherent postmeasurement state or exact single-shot eigenvalue |
| Herald and retain an entangled output | The above generation dynamics plus atom-preserving, sufficiently spin-coherent accepted readout | Atom survival alone does not validate the Kraus operator |
| Repeated filtering on the same trial | Reliable coherent accepted readout after every cycle; no added atoms; stable ring and timing | The hard-core 1/4 plateau is not automatic at finite U or with readout dephasing |

Neither direct local-spin readout nor the known circuit implementation is automatically available on a hardware platform offering only the vacancy resources. Conversely, a processor already offering high-quality local rotations/readout gains no demonstrated advantage merely by adding a flux plaquette. The assessment must name a platform and its available operations before comparing overall hardware cost.

## 7. Attribution and verification

The [reference account](../REFERENCES.md) credits cyclic spin readout, parity imaging, synthetic flux, and spin-coherence measurements, and distinguishes component capabilities from a complete spinful-ring apparatus. The native return estimator and retained-state filter have different operational requirements; terminal outcome probabilities alone do not validate coherent measurement back-action.

The [readout model](../source/readout.py) and [resource checker](../source/check_resources.py) cover five groups: the exact parity projector over all particle numbers up to three and a five-particle alias; complete number-block Lindblad evolution; dephasing and the exact repetition recurrence; matched sampling coefficients; and classical parity-bit error budgets. These finite-matrix diagnostics supplement the analytical identities and do not certify apparatus parameters or floating-point intervals.

Run these checks with the complete suite in [Reproduce](../REPRODUCE.md). The complete-parity estimator, which uses the entire image rather than only the return event, has its separate [proof and precision bounds](COMPLETE_PARITY_THEORY.md).
