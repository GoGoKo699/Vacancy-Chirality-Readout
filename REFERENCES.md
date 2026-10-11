# References and attribution

The result is a specified vacancy-return instrument and its finite-repulsion chirality reconstruction. The spin-cycle mechanism, chirality observables, permutation algebra, and standard measurement tools have established antecedents.

## Mathematical and measurement predecessors

**[K22]** K. K. Kesharpu and E. A. Kochetov, *Spin-flip-induced superfluidity in a ring of spinful hard-core bosons*, Physical Review A **106**, 033323 (2022). [arXiv:2209.06149](https://arxiv.org/abs/2209.06149). The cyclic spin permutation on the closing link and the resulting spin-sector boundary flux underlie the hard-core reduction. The repository specializes this mechanism to an explicit vacancy-return measurement.

**[R23]** L. I. Reascos, B. Murta, E. F. Galvão, and J. Fernández-Rossier, *Quantum circuits to measure scalar spin chirality*, Physical Review Research **5**, 043087 (2023). [Article](https://doi.org/10.1103/PhysRevResearch.5.043087) · [arXiv:2306.14804](https://arxiv.org/abs/2306.14804). This work supplies direct Pauli-string measurements, a controlled-cycle Hadamard test, qutrit phase estimation, and measurement-based chirality-sector preparation. Its normalized observable is $\widehat\chi=4\chi/\sqrt3$. The chirality witness and sector-preparation principle are inherited. The direct Pauli baseline requires no entangling measurement gate.

**[EW01]** T. Eggeling and R. F. Werner, *Separability properties of tripartite states with $U\otimes U\otimes U$ symmetry*, Physical Review A **63**, 042111 (2001). [arXiv:quant-ph/0010096](https://arxiv.org/abs/quant-ph/0010096). The collective-unitary commutant is permutation-generated. For three qubits its Hermitian span is five-dimensional; combining it with spatial reflection gives the two-observable finite-repulsion contrast used here.

**[BDL11]** S. Bravyi, D. P. DiVincenzo, and D. Loss, *Schrieffer–Wolff transformation for quantum many-body systems*, Annals of Physics **326**, 2793–2826 (2011). [arXiv:1105.0675](https://arxiv.org/abs/1105.0675). This is the background for the optional fixed-time block-elimination interpretation. The quantitative operation bound in this repository instead follows directly from energy conservation and Duhamel’s formula.

**[M17]** A. J. Menssen and collaborators, *Distinguishability and many-particle interference*, Physical Review Letters **118**, 153603 (2017). [Article](https://doi.org/10.1103/PhysRevLett.118.153603) · [arXiv:1609.09804](https://arxiv.org/abs/1609.09804). Independently prepared photons in a tritter exhibit spatial-count sensitivity to a collective internal-state triad phase. Their supplementary Eqs. (20)–(21), with $R=r_{12}r_{23}r_{31}$, give

$$
P_{120}=\frac{1-2R\cos(\varphi+\pi/3)}9,
\qquad
P_{021}=\frac{1-2R\cos(\varphi-\pi/3)}9.
$$

Subtracting these expressions yields

$$
P_{120}-P_{021}=\frac{2\sqrt3}{9}R\sin\varphi.
$$

This algebraic consequence explicitly credits a precedent for **signed** internal-state information in spatial interference. The present derivations concern arbitrary correlated three-spin inputs, vacancy-return dynamics, and finite-repulsion calibration. The broad conversion of internal-state information into spatial counts is established prior art. This attribution supersedes the narrower “cosine-type” characterization in the preserved source record.

## Experimental component precedents

| Source | Established component | Scope of its use here |
|---|---|---|
| **[B24]** S. Buob et al., *A Strontium Quantum-Gas Microscope*, PRX Quantum **5**, 020316 (2024). [Article](https://doi.org/10.1103/PRXQuantum.5.020316) | On-site parity projection through light-assisted collisions during cooling | Supports the parity-imaging primitive; the stated experiment uses bosonic strontium |
| **[A13]** M. Aidelsburger et al., *Experimental realization of strong effective magnetic fields in optical superlattice potentials*, Applied Physics B **113**, 1–11 (2013). [Article](https://doi.org/10.1007/s00340-013-5418-1) · [arXiv:1212.2911](https://arxiv.org/abs/1212.2911) | Photon-assisted orbital phases independent of internal atomic structure; four-site single-particle plaquette dynamics | Supports the synthetic-flux control concept |
| **[F25]** A. Muzi Falconi et al., *Microsecond-Scale High-Survival and Number-Resolved Detection of Ytterbium Atom Arrays*, Physical Review Letters **135**, 203402 (2025). [Article](https://doi.org/10.1103/n3bg-7yw7) · [arXiv:2507.01011](https://arxiv.org/abs/2507.01011) | Fast number-resolved imaging and atom survival | Supports those detector capabilities; arbitrary retained internal-spin coherence requires its own evidence |
| **[N24]** Y. Nakamura et al., *Hybrid Atom Tweezer Array of Nuclear Spin and Optical Clock Qubits*, Physical Review X **14**, 041062 (2024). [Article](https://doi.org/10.1103/PhysRevX.14.041062) · [arXiv:2406.12247](https://arxiv.org/abs/2406.12247) | Nuclear-spin coherence tested while imaging a different isotope | Demonstrates a relevant coherence test in a dual-isotope architecture |

These are component-level comparisons. Combining their capabilities into the particular spinful flux ring and coherent retained detector is an additional apparatus question. The analytical result specifies that detector requirement explicitly.

## Sampling comparison

Use the common $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$ normalization and count $n$ original independent input trials. On the maximally mixed input, the declared $U/J=16$ four-setting example has $n\operatorname{Var}(\widehat\chi)\simeq0.8759835$ under its stated allocation. Separately sampled Pauli products give $0.5625$ and the controlled-cycle Hadamard test gives $0.25$. These ideal sampling coefficients compare the named implementations; they do not assign equal hardware cost to their controls.

For derivations and the exact allocation, see [the resource account, §6](source/resource-review-2026-10-10/READOUT_AND_RESOURCES.md#6-compare-the-correct-tasks-and-input-trial-normalization). The [preserved source reading record](source/resource-review-2026-10-10/SOURCES.md) documents the original targeted comparison. The present references refine attribution without asserting exhaustive priority clearance.
