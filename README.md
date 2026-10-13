# Vacancy Chirality Readout

**Interactions improve chirality-readout precision at fixed spatial and detector resources.**

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

The [passive-measurement theorem](research/PASSIVE_PARITY_BOUND.md) proves the lower bound $13/32=0.40625$ for **every spin-independent passive four-mode interferometer followed by complete total-site parity detection**. It permits independent local spin rotations before mixing, randomized and adaptive settings across fresh copies, and all universally unbiased estimators. The native upper bound is about 19.1% below this lower bound. The passive bound is not claimed to be tight. Extra spatial modes, intermediate measurements, interleaved spin controls, and stronger detectors define different resource classes.

The same proof gives finite worst-case mean-square-error bounds for **arbitrary, possibly biased estimators**. With exact channel calibration and nominal control parameters, the interacting readout beats the passive minimax risk for every $n\ge12{,}000$ throughout the modeled window $\gamma t\le0.01$, $q\le0.003$.

| Interacting readout | Certified worst-case variance coefficient |
|---|---:|
| Ideal model, exactly calibrated score | $<0.328698$ |
| Parameter budget $\Delta\le0.0025$, ideal parity, exact response calibration | $<0.394825$ |
| Parameter budget $\Delta\le0.0015$, ideal parity, exact response calibration | $<0.366081$ |
| Independent bit errors $q\le0.03$, no loss, known $q$ | $<0.382273$ |
| Equal-rate loss $\gamma t\le0.01$, independent bit errors $q\le0.003$, exact full-channel calibration | $<0.366461$ |
| Joint $\Delta\le0.001$, $\gamma t\le0.01$, $q\le0.003$, exact response calibration | $<0.395874$ |
| Joint $\Delta\le0.001$, $\gamma t\le0.005$, $q\le0.0015$, exact response calibration | $<0.366685$ |

The [parameter-robustness theorem](research/PARAMETER_ROBUSTNESS.md) uses

$$
\Delta=1.643|u-16|+3.163|\tau-\pi|+4.086|\Phi_++2\pi/3|,
\qquad \Phi_-=-\Phi_+.
$$

Both flux settings use the same $u$ and $\tau$. The ideal-parity region $\Delta\le0.0015$ and the smaller joint region in the final row both beat the passive minimax risk for every $n\ge12{,}000$, allowing biased passive estimators. The joint rows include parameter variation and the stated loss and detector errors simultaneously; $t$ is the actual evolution duration. The theorem also bounds fixed-score finite MSE under generic spin-dependent perturbations and imperfect loading.

These control and noise bounds are analytical assumptions, not measured apparatus performance. Full-channel calibration uses two known product preparations to determine the mean response, without separately known noise rates. A third preparation determines its second moment within the spin and reflection symmetry class; see [direct channel calibration](research/CHANNEL_CALIBRATION.md).

The [finite-calibration theorem](research/FINITE_CALIBRATION.md) includes reference randomness and preparation costs in an unconditional risk bound. Its independent-stage polynomial correction has a concrete total-preparation advantage after sufficient reuse, with extensions to bounded preparation error and drift. A fresh moment certificate extends both finite budgets to the final joint region above without exact operational response coefficients. The sufficient reference and reuse counts remain large. A [repeated chirality scan](research/REPEATED_STATE_SCANS.md) specifies the single-orbital mode resource, stability, cycle-cost, and parallel-channel assumptions needed to interpret such reuse physically.

The mode restriction has an operational meaning: three-mode passive parity measurements cannot identify chirality at all, while one extra vacuum mode makes passive identification possible. Interactions improve precision within that four-mode setting. With **eight spatial modes**, a [passive Fourier-and-splitter construction](research/PASSIVE_PARITY_BOUND.md) attains coefficient $27/88\simeq0.306818$, below the current interacting score, using the same type of parity detector on more outputs. This locates the resource boundary; it does not prove that eight modes are necessary. The separate [separable-effects benchmark](research/SEPARABLE_BENCHMARK.md) has tight unbiased coefficient $9/16$. Established cyclic-readout and interaction-enabled parity precedents are credited in [References](REFERENCES.md).

A complementary [return-only reconstruction](THEORY.md#3-exact-finite-repulsion-reconstruction), using two durations and two fluxes, has a general existence proof for every fixed finite $U/J>0$. A separate mixed-input score and its own neighborhoods address the $9/16$ comparison. The primary score's passive-separation region is the $\Delta$ region above.

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
| [Parameter robustness](research/PARAMETER_ROBUSTNESS.md) | Primary-score parameter regions, exact recalibration, and fixed-score finite MSE |
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
