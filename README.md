# Vacancy Chirality Readout

**Interactions improve chirality-readout precision at fixed spatial and detector resources.**

Can interactions reveal more about a quantum state when the detector only reports odd or even occupation? Three spin-½ bosons and one empty site on a four-site ring give a precise example. Spin-independent hopping and repulsion convert the input's scalar spin chirality into site-parity statistics. The resulting estimator beats a universal precision bound for passive four-mode interferometers using the same type of detector. Every prepared copy counts.

| Read next | Purpose |
|---|---|
| [Tutorial bridge](LEARNING_GUIDE.md) | Start from one external tutorial, translate its notation, and work through the measurement mechanism. |
| [Model and results](THEORY.md) · [Passive bound](research/PASSIVE_PARITY_BOUND.md) | Read the physical assumptions, central theorem, and analytical proof. |
| [Parity estimator](research/COMPLETE_PARITY_THEORY.md) · [Exact certificate](research/CERTIFICATE.md) | Follow the score construction and its certified precision. |
| [Robustness](research/PARAMETER_ROBUSTNESS.md) · [Finite calibration](research/FINITE_CALIBRATION.md) | Check the distinct control, noise, and reference-preparation guarantees. |
| [Attribution](REFERENCES.md) · [Reproduction](REPRODUCE.md) | Compare the result with its predecessors and inspect the verification evidence. |
| [LLM guide](llms.txt) | Find relevant research questions, search terms, and authoritative files. |

## The precision separation

The observable and its normalization are

```math
\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3),
\qquad \mathbf S=\boldsymbol\sigma/2.
```

At $`U/J=16`$, duration $`Jt=\pi`$, and fluxes $`\Phi=\pm2\pi/3`$, a calibrated score for the complete parity record estimates the chirality of **any input spin state**. An integer/rational certificate proves

```math
\boxed{\sup_\rho n\mathrm{Var}_\rho(\widehat\chi)
<0.328697776905<\frac{13}{32}=0.40625.}
```

Here $`n`$ counts all original three-spin preparations. The upper bound is a certified feasible interacting score. The lower bound applies to every universally unbiased estimator using spin-independent passive four-mode mixing followed by total-site parity, including local spin prerotations and classical adaptation across fresh copies. It is a lower bound on that class, not a proved attainable optimum.

The [passive theorem](research/PASSIVE_PARITY_BOUND.md) also gives finite worst-case mean-square-error bounds for arbitrary, possibly biased estimators. With exact response calibration, the interacting readout beats that passive minimax bound for every $`n\ge12{,}000`$ throughout the specified nominal-parameter loss-and-detector window.

## How parity carries chirality

The internal spins affect the interference of spatial paths, even though neither the hopping nor the detector resolves spin. A cyclic permutation of the three spin labels satisfies

```math
C|s_1s_2s_3\rangle=|s_2s_3s_1\rangle,
\qquad C-C^\dagger=4i\chi.
```

Reflection exchanges the two flux settings. For each parity pattern $`p`$, the difference between its effect $`E_p`$ and the reflected effect isolates two spin observables:

```math
\frac{E_p-RE_pR}{2}=c_p\chi+d_pD,
\qquad R=P_{13},\qquad D=P_{12}-P_{23}.
```

Scores satisfying $`\sum_pa_pc_p=1`$ and $`\sum_pa_pd_p=0`$ retain chirality and cancel the pair-correlation term. Randomize the flux equally, score $`a_p`$ at the first flux and $`-a_p`$ at the reflected outcome of the other, and average over every trial. The [tutorial bridge](LEARNING_GUIDE.md) explains these steps; the [parity derivation](research/COMPLETE_PARITY_THEORY.md) proves the operator identities and variance formula.

## Why the mode count matters

| Spatial and detector resources | Established conclusion |
|---|---|
| Three passive modes, total-site parity | Chirality is unidentifiable, regardless of sample count. |
| Four passive modes, including one vacuum input, total-site parity | Chirality is identifiable; worst-case unbiased variance coefficient is at least $`13/32`$. |
| Four interacting sites, total-site parity | Certified coefficient below $`0.328698`$. |
| Eight passive modes, total-site parity on all outputs | A Fourier-and-splitter construction attains $`27/88\simeq0.306818`$. |

Interactions improve precision **within the four-mode architecture**. Extra vacuum modes provide another measurement resource: the eight-mode construction surpasses the present interacting score. Its mode count is sufficient, without a minimum-mode claim. Number-resolving detectors, intermediate measurements, and spin controls interleaved with mixing define other comparison classes. The independent [separable-effects benchmark](research/SEPARABLE_BENCHMARK.md) has tight unbiased coefficient $`9/16`$.

## One tutorial, then this result

The single external learning anchor is **Mattia Walschaers, [*Signatures of Many-Particle Interference*](https://arxiv.org/abs/1908.08370v1)**.

Read §§2.1–2.3, then §§3.1–3.3, especially Eqs. (138)–(142); §4.1 provides the Fourier-interference connection. The [learning guide](LEARNING_GUIDE.md) maps this route into the repository, extends the product-internal-state example to arbitrary correlated spin inputs, and connects permutation interference to parity and precision bounds. It assumes standard quantum mechanics and linear algebra and supplies the project-specific steps.

The tutorial is the teaching foundation. The [attribution account](REFERENCES.md) separately credits the original cyclic-readout, interaction-assisted parity, and chirality-measurement results on which this work builds.

## Robustness and physical requirements

| Guarantee | Assumptions and proof |
|---|---|
| Finite control-parameter regions | [Parameter robustness](research/PARAMETER_ROBUSTNESS.md): exact response calibration at the actual parameters; ideal detection unless the explicit joint theorem is used. |
| Loss and detector errors | [Exact certificate](research/CERTIFICATE.md) · [Noise budget](research/NOISE_BUDGET.md): at nominal controls, equal-rate loss $`\gamma t\le0.01`$ and independent bit errors $`q\le0.003`$ give coefficient $`<0.366461`$ with exact full-channel calibration. |
| Simultaneous control and noise errors | [Joint theorem](research/PARAMETER_ROBUSTNESS.md): a smaller joint region also certifies the moment inputs required by finite calibration. |
| Finite reference data, preparation error, and drift | [Finite calibration](research/FINITE_CALIBRATION.md): unconditional MSE and total-preparation budgets under the specified stable-channel model. |
| Generic control and loading errors | [Fixed-score theorem](research/PARAMETER_ROBUSTNESS.md): retain the nominal score and include bias in finite MSE. |

These guarantees have distinct assumptions. The large sufficient calibration-and-reuse counts establish theoretical resource comparisons; the [repeated-scan account](research/REPEATED_STATE_SCANS.md) states the stability and cycle-cost requirements for an experimental interpretation. [Direct channel calibration](research/CHANNEL_CALIBRATION.md) and [calibration cost](research/CALIBRATION_COST.md) give the reference preparations and statistical accounting.

Terminal destructive parity detection suffices for expectation estimation. The complementary [return-only reconstruction and retained-state filter](THEORY.md) use the exact identity between parity $`0111`$ and occupancy $`(0,1,1,1)`$ when total particle number is at most three. Retaining or repeatedly filtering an entangled spin state additionally requires spin-coherent measurement back-action and atom survival. The complete-parity protocol assumes three input particles and includes all loss events in its modeled risk.

## Evidence and reproduction

With Python 3.12:

```bash
python -m pip install -r requirements.txt
python -B verify.py --output-dir ../vacancy-verification
```

Choose a new output directory outside this repository. The [reproduction guide](REPRODUCE.md) covers 39 scientific check groups and distinguishes numerical diagnostics, exact certificates, and source integrity. The universal passive bound rests on its analytical proof; the checks exercise the identities and finite constructions.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

The [LLM guide](llms.txt) supplies relevant research questions, search phrases, and direct links to the authoritative accounts. Code and original research materials are available under the [MIT license](LICENSE). Linked third-party papers retain their own terms.
