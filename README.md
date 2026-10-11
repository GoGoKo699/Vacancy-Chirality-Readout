# Vacancy Chirality Readout

**Measuring three-spin chirality through the motion of a single vacancy.**

Three spin-½ bosons occupy three sites of a four-site ring. Spin-independent hopping in an imposed flux converts their scalar spin chirality into a position-return probability. At finite repulsion, two calibrated durations and two flux signs reconstruct the input chirality from four sets of return frequencies. With known total particle number at most three, total on-site parity identifies the accepted occupancy exactly.

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

## The result

Let $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$, with $\mathbf S=\boldsymbol\sigma/2$. In the uniform spin-independent model, the flux-reversed return effects obey

$$
E_+(t)-E_-(t)=g_U(t)\chi+h_U(t)(P_{12}-P_{23}).
$$

The second term is a pair-correlation imbalance. Two durations cancel it: if $d_i=p(-2\pi/3,t_i)-p(+2\pi/3,t_i)$ and $g_1h_2-g_2h_1\ne0$, then

$$
\boxed{\langle\chi\rangle=\frac{h_2d_1-h_1d_2}{g_1h_2-g_2h_1}.}
$$

Suitable durations exist for every fixed finite $U/J>0$. The identity holds for arbitrary input spin states, including inputs entangled with an external reference. It is an expectation-value measurement on independently prepared copies.

## What the detector must do

| Task | Required readout |
|---|---|
| Estimate input chirality | Total site parity; a terminal destructive image suffices |
| Retain the heralded entangled spin state | An accepted occupancy measurement that preserves the required spin coherence |
| Repeat the filter on the same state | Coherent accepted readout, atom survival, and stable retained positions between rounds |

For particle number at most three, parity pattern $0111$ is exactly occupancy $(0,1,1,1)$, even at finite repulsion. Correct return probabilities or atom survival alone do not establish coherent retained-state back-action.

## Read the research

| Start here | Contents |
|---|---|
| [Theory and proofs](THEORY.md) | Model, exact reconstruction, hard-core filter, finite-repulsion bounds, and detector requirements |
| [References and attribution](REFERENCES.md) | Mathematical predecessors, measurement baselines, and component-level experimental precedents |
| [Reproduce the results](REPRODUCE.md) | Sixteen check groups, exact coefficient certificates, and source-integrity checks |

The protocol uses coherent flux hopping, a localized vacancy, calibrated dynamics, and parity detection. Its runtime readout uses no local spin rotations, controlled permutations, or spin-resolved measurements. These are different physical resources from direct Pauli or controlled-cycle measurements; their calibration, survival, and sampling costs are included in the account.

## Reproduce

With Python 3.12:

```bash
python -m pip install -r requirements.txt
python verify.py --output-dir ../vacancy-verification
```

Choose a new output directory outside this repository. The verifier checks sixteen scientific groups and reports assertion success, canonical numerical-report identity, and source preservation separately. See [REPRODUCE.md](REPRODUCE.md) for the result map and interpretation.

Code and original research materials are available under the [MIT license](LICENSE).
