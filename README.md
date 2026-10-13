# Vacancy Chirality Readout

**An interacting chirality readout that beats passive four-mode parity measurements.**

Three spin-½ bosons occupy three sites of a four-site ring. Spin-independent hopping and interactions in an imposed flux convert their scalar spin chirality into site-parity statistics. At $U/J=16$, one duration $Jt=\pi$ and two opposite fluxes reconstruct the chirality of an arbitrary input spin state. A universal bound on passive four-mode interferometers, followed by the same detector, makes the improvement from this interacting readout quantitative. Every prepared copy contributes to the estimate.

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

## The result

Let $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$, with $\mathbf S=\boldsymbol\sigma/2$. For each of the eight odd-weight parity patterns, pair the effect $E_p$ at flux $-2\pi/3$ with the reflected outcome at flux $+2\pi/3$. Spin and reflection symmetries give

$$
\frac{E_p-RE_pR}{2}=c_p\chi+d_p(P_{12}-P_{23}),\qquad R=P_{13}.
$$

Calibrated scores $a_p$ satisfying $\sum_pa_pc_p=1$ and $\sum_pa_pd_p=0$ cancel the pair-correlation nuisance term. With equally randomized fluxes, score $a_p$ at the first flux and $-a_p$ for its reflected outcome at the other. The sample mean is unbiased for arbitrary input spin states, including inputs entangled with an external reference.

An [integer/rational certificate](research/CERTIFICATE.md) establishes an exactly calibrated score at the stated operating point with

$$
\boxed{\sup_\rho n\operatorname{Var}(\widehat\chi)<0.328697776905<\frac{13}{32}.}
$$

The [passive-measurement theorem](research/PASSIVE_PARITY_BOUND.md) proves the lower bound $13/32=0.40625$ for **every spin-independent passive four-mode interferometer followed by complete total-site parity detection**. It permits independent local spin rotations before mixing, randomized and adaptive settings across fresh copies, and all universally unbiased estimators. The native upper bound is over 19.1% below this lower bound. The passive bound is not claimed to be tight. Extra spatial modes, intermediate measurements, interleaved spin controls, and stronger detectors define different resource classes.

The same proof gives finite worst-case mean-square-error bounds for **arbitrary, possibly biased estimators**. With exact channel calibration, the interacting readout beats the passive minimax risk for every $n\ge12{,}000$, even throughout the modeled equal-loss and independent bit-error window below.

| Interacting readout | Certified worst-case variance coefficient |
|---|---:|
| Ideal model, exactly calibrated score | $<0.328698$ |
| Independent bit errors $q\le0.03$, no loss, known $q$ | $<0.382273$ |
| Equal-rate loss $\gamma t\le0.01$, independent bit errors $q\le0.003$, exact full-channel calibration | $<0.366461$ |

The loss and detector rates specify analytical models, not measured apparatus performance. Full-channel calibration uses two known product preparations to determine the mean response, without separately known noise rates. A third preparation determines its second moment within the spin and reflection symmetry class; see [direct channel calibration](research/CHANNEL_CALIBRATION.md).

The [finite-calibration theorem](research/FINITE_CALIBRATION.md) includes reference randomness and preparation costs in an unconditional risk bound. Its independent-stage polynomial correction has a concrete total-preparation advantage after sufficient reuse, with extensions to bounded preparation error and drift. The sufficient reference and reuse counts are large. A [repeated chirality scan](research/REPEATED_STATE_SCANS.md) specifies the stability, cycle-cost, and parallel-channel assumptions needed to interpret such reuse physically.

The mode restriction has an operational meaning: three-mode passive parity measurements cannot identify chirality at all, while one extra vacuum mode makes passive identification possible. Interactions improve precision within that four-mode setting. The separate [separable-effects benchmark](research/SEPARABLE_BENCHMARK.md) has tight unbiased coefficient $9/16$. Number-resolved Fourier readout and ancilla-controlled circuits provide stronger measurement resources; their established results and the interaction-enabled cyclic-parity precedent are credited in [References](REFERENCES.md).

A complementary [return-only reconstruction](THEORY.md#3-exact-finite-repulsion-reconstruction), using two durations and two fluxes, has a general existence proof for every fixed finite $U/J>0$. The two-setting passive separation is certified at the stated Hamiltonian point. Separate parameter neighborhoods are certified for the mixed-input score's comparison with $9/16$; they do not establish passive separation throughout those neighborhoods.

## What the detector must do

| Task | Required readout |
|---|---|
| Estimate input chirality with two settings | Complete total-site-parity image; terminal destructive detection suffices |
| Estimate using the four-setting return method | Identify parity pattern $0111$ |
| Retain the heralded entangled spin state | Accepted occupancy measurement preserving the required spin coherence |
| Repeat the filter on the same state | Coherent accepted readout, atom survival, and stable retained positions |

For known particle number at most three, parity pattern $0111$ is exactly occupancy $(0,1,1,1)$, even at finite repulsion. Complete-parity estimation assumes three input particles: one lost particle produces an even-weight pattern and receives score zero, whereas two losses can mimic singleton patterns. The [noise budget](research/NOISE_BUDGET.md) includes this contamination and detector errors. Correct probabilities or atom survival alone do not establish coherent retained-state back-action.

## Read the research

| Start here | Contents |
|---|---|
| [Theory and proofs](THEORY.md) | Model, both reconstructions, retained-state filter, and physical requirements |
| [Complete-parity derivation](research/COMPLETE_PARITY_THEORY.md) | One-duration calibration, score optimization, and exact worst-case variance |
| [Passive parity bound](research/PASSIVE_PARITY_BOUND.md) | Three-mode impossibility, four-mode Fisher bound, constructions, and finite risk with arbitrary estimators |
| [Separable measurement benchmark](research/SEPARABLE_BENCHMARK.md) | Tight lower bound, adaptive measurements, and scope |
| [Exact certificate](research/CERTIFICATE.md) | Validated finite-repulsion performance and sufficient neighborhoods |
| [Direct channel calibration](research/CHANNEL_CALIBRATION.md) | Three product references, unknown noise rates, and finite-sample validation |
| [Calibration cost and reuse](research/CALIBRATION_COST.md) | Reference counts, finite biased-estimator comparison, confidence, and amortization |
| [Finite calibration theorem](research/FINITE_CALIBRATION.md) | Independent-stage correction, unconditional risk, total preparation counts, and preparation/drift tolerances |
| [Repeated chirality scans](research/REPEATED_STATE_SCANS.md) | Candidate workflow, experimental precedents, full-cycle costs, stability, and parallel-channel pooling |
| [Noise budget](research/NOISE_BUDGET.md) | Loss, parity errors, calibration uncertainty, and mean-square error |
| [References and attribution](REFERENCES.md) | Mathematical predecessors, measurement baselines, and component experiments |
| [Reproduce the results](REPRODUCE.md) | Scientific check groups, exact certificates, and source integrity |

The runtime protocol uses coherent flux hopping, a localized vacancy, calibrated dynamics, and parity detection. It requires no local spin rotations, controlled permutations, or spin-resolved measurements. Direct channel calibration additionally uses known local spin-axis preparations. Calibration, survival, and sampling costs remain part of the comparison with other measurement resources.

## Reproduce

With Python 3.12:

```bash
python -m pip install -r requirements.txt
python -B verify.py --output-dir ../vacancy-verification
```

Choose a new output directory outside this repository. The verifier reports scientific assertions, exact-certificate completion, preserved numerical-report identity, and source integrity separately. See [REPRODUCE.md](REPRODUCE.md) for interpretation.

Code and original research materials are available under the [MIT license](LICENSE).
