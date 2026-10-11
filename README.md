# Vacancy Chirality Readout

**Measuring three-spin chirality through the motion of a single vacancy.**

Three spin-½ bosons occupy three sites of a four-site ring. Spin-independent hopping in an imposed flux converts their scalar spin chirality into site-parity statistics. At the certified finite-repulsion point $U/J=16$, one duration $Jt=\pi$ and two opposite fluxes reconstruct the input chirality from complete parity images. Every prepared copy contributes to the estimate.

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

## The result

Let $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$, with $\mathbf S=\boldsymbol\sigma/2$. For each of the eight odd-weight parity patterns, pair the effect $E_p$ at flux $-2\pi/3$ with the reflected outcome at flux $+2\pi/3$. Spin and reflection symmetries give

$$
\frac{E_p-RE_pR}{2}=c_p\chi+d_p(P_{12}-P_{23}),\qquad R=P_{13}.
$$

Calibrated scores $a_p$ satisfying $\sum_pa_pc_p=1$ and $\sum_pa_pd_p=0$ cancel the pair-correlation nuisance term. With equally randomized fluxes, score $a_p$ at the first flux and $-a_p$ for its reflected outcome at the other. The sample mean is unbiased for arbitrary input spin states, including inputs entangled with an external reference.

An [integer/rational certificate](research/CERTIFICATE.md) establishes an exactly calibrated score at the stated operating point with

$$
\boxed{\sup_\rho n\operatorname{Var}(\widehat\chi)<0.389156263913<\frac9{16}.}
$$

The value $9/16$ is the exact minimax coefficient for universally unbiased estimation using **fully separable measurement effects on each three-spin copy**, including adaptive local measurements. The [benchmark proof](research/SEPARABLE_BENCHMARK.md) states the precise measurement class. The certified result reduces the worst-case copy requirement by over 30.8% within this comparison. It does not equate apparatus costs or cover collective measurements across copies.

The implemented readout can also be calibrated directly using **three known product-state preparations**. They determine the full response within the spin-rotation and reflection symmetry class, including loss and detector errors, without separately known noise rates. At the nominal operating point, exact channel recalibration retains a certified variance coefficient below $0.446<9/16$ throughout the specified equal-loss and bit-error window $\gamma t\le0.01$, $q\le0.003$. Independent validation bounds finite-calibration bias and sampling cost; see [direct channel calibration](research/CHANNEL_CALIBRATION.md).

**Calibration cost limits the practical claim.** Variance-sensitive validation and a finite-sample separable lower bound that permits biased estimators now support an explicit [preparation-cost audit](research/CALIBRATION_COST.md). It does not establish a one-use total-cost advantage when reference and unknown-state preparations cost the same. Cheap product references or many uses of a stable calibration can change that conclusion. The worked finite profile is a labeled forecast, not experimental data or an unconditional risk guarantee.

A complementary [return-only reconstruction](THEORY.md#3-exact-finite-repulsion-reconstruction), using two durations and two fluxes, has a general existence proof for every fixed finite $U/J>0$. The two-setting result is certified at the stated point and in a sufficient neighborhood; it is not asserted for every repulsion.

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
| [Separable measurement benchmark](research/SEPARABLE_BENCHMARK.md) | Tight lower bound, adaptive measurements, and scope |
| [Exact certificate](research/CERTIFICATE.md) | Validated finite-repulsion performance and sufficient neighborhoods |
| [Direct channel calibration](research/CHANNEL_CALIBRATION.md) | Three product references, unknown noise rates, and finite-sample validation |
| [Calibration cost and reuse](research/CALIBRATION_COST.md) | Reference counts, finite biased-estimator comparison, confidence, and amortization |
| [Noise budget](research/NOISE_BUDGET.md) | Loss, parity errors, calibration uncertainty, and mean-square error |
| [References and attribution](REFERENCES.md) | Mathematical predecessors, measurement baselines, and component experiments |
| [Reproduce the results](REPRODUCE.md) | Twenty-seven check groups, exact certificates, and source integrity |

The runtime protocol uses coherent flux hopping, a localized vacancy, calibrated dynamics, and parity detection. It requires no local spin rotations, controlled permutations, or spin-resolved measurements. Direct channel calibration additionally uses known local spin-axis preparations. Calibration, survival, and sampling costs remain part of the comparison with other measurement resources.

## Reproduce

With Python 3.12:

```bash
python -m pip install -r requirements.txt
python verify.py --output-dir ../vacancy-verification
```

Choose a new output directory outside this repository. The verifier reports scientific assertions, exact-certificate completion, preserved numerical-report identity, and source integrity separately. See [REPRODUCE.md](REPRODUCE.md) for interpretation.

Code and original research materials are available under the [MIT license](LICENSE).
