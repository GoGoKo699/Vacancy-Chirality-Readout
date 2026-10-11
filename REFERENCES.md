# References and attribution

The results specify a vacancy-return instrument, finite-repulsion chirality reconstruction from return or complete-parity statistics, and a certified sampling comparison with separable measurements. The spin-cycle mechanism, chirality observables, permutation algebra, and standard measurement tools have established antecedents.

## Mathematical and measurement predecessors

**[K22]** K. K. Kesharpu and E. A. Kochetov, *Spin-flip-induced superfluidity in a ring of spinful hard-core bosons*, Physical Review A **106**, 033323 (2022). [arXiv:2209.06149](https://arxiv.org/abs/2209.06149). The cyclic spin permutation on the closing link and the resulting spin-sector boundary flux underlie the hard-core reduction. The repository specializes this mechanism to an explicit vacancy-return measurement.

**[R23]** L. I. Reascos, B. Murta, E. F. Galvão, and J. Fernández-Rossier, *Quantum circuits to measure scalar spin chirality*, Physical Review Research **5**, 043087 (2023). [Article](https://doi.org/10.1103/PhysRevResearch.5.043087) · [arXiv:2306.14804](https://arxiv.org/abs/2306.14804). This work supplies direct Pauli-string measurements, a controlled-cycle Hadamard test, qutrit phase estimation, and measurement-based chirality-sector preparation. Its normalized observable is $\widehat\chi=4\chi/\sqrt3$. The chirality witness and sector-preparation principle are inherited. The direct Pauli baseline requires no entangling measurement gate.

**[EW01]** T. Eggeling and R. F. Werner, *Separability properties of tripartite states with $U\otimes U\otimes U$ symmetry*, Physical Review A **63**, 042111 (2001). [arXiv:quant-ph/0010096](https://arxiv.org/abs/quant-ph/0010096). The collective-unitary commutant is permutation-generated. For three qubits its Hermitian span is five-dimensional; combining it with spatial reflection gives the two-observable finite-repulsion contrast used here.

**[CLMOW14]** E. Chitambar, D. Leung, L. Mancinska, M. Ozols, and A. Winter, *Everything You Always Wanted to Know About LOCC (But Were Afraid to Ask)*, Communications in Mathematical Physics **328**, 303–326 (2014). [Article](https://doi.org/10.1007/s00220-014-1953-9) · [arXiv:1210.4583](https://arxiv.org/abs/1210.4583). Section 2.3 establishes the inclusion of LOCC in separable operations. This supports the operational scope of the per-copy separable-effects benchmark; it does not supply the chirality-specific constant derived here.

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

This algebraic consequence explicitly credits a precedent for **signed** internal-state information in spatial interference. The broad conversion of internal-state information into spatial counts is established prior art. This attribution supersedes the narrower “cosine-type” characterization in the preserved source record.

**[N26]** L. Novo, M. Robbio, E. F. Galvão, and N. J. Cerf, *Native linear-optical protocol for efficient multivariate trace estimation*, [arXiv:2601.14204](https://arxiv.org/abs/2601.14204) (2026). Equations (3)–(12) diagonalize the cyclic permutation with Fourier interferometry and internal-state-insensitive particle counts. Theorem 1 explicitly permits entangled inputs. Thus arbitrary correlated inputs alone do not distinguish the present protocol. The count statistic is $f(\mathbf S)=\sum_{j=0}^{M-1}jS_j\pmod M$. At $M=3$, occupations $(3,0,0)$ and $(1,2,0)$ have identical parity $(1,0,0)$ but $f=0$ and $f=2$, respectively: total-site parity cannot recover this count statistic. The distinction pursued here is an interacting four-site vacancy implementation with parity-only detection and finite-repulsion calibration.

**[NGH25]** G. Nambiar, A. Grankin, and M. Hafezi, *Diagnosing Electronic Phases of Matter Using Photonic Correlation Functions*, Physical Review X **15**, 041020 (2025). [Article](https://doi.org/10.1103/67zs-hqf3) · [arXiv:2410.24215](https://arxiv.org/abs/2410.24215). Equation (70) relates an outgoing optical quadrature to a material scattering operator. Section VII, particularly Eq. (96), extracts a lattice scalar-chirality sum using phase-sensitive filtered two-photon quadratures and polarization combinations. This is a precedent for indirect static chirality readout and symmetry-based isolation. Its derivation treats a half-filled Fermi-Hubbard Mott insulator in the regime $t\ll|\omega_L-U|\ll U$; the present instrument instead uses three spinful bosons and a vacancy on four sites, with explicitly calibrated finite-repulsion parity effects.

## Experimental component precedents

| Source | Established component | Scope of its use here |
|---|---|---|
| **[B24]** S. Buob et al., *A Strontium Quantum-Gas Microscope*, PRX Quantum **5**, 020316 (2024). [Article](https://doi.org/10.1103/PRXQuantum.5.020316) | On-site parity projection through light-assisted collisions during cooling | Supports the parity-imaging primitive; the stated experiment uses bosonic strontium |
| **[A13]** M. Aidelsburger et al., *Experimental realization of strong effective magnetic fields in optical superlattice potentials*, Applied Physics B **113**, 1–11 (2013). [Article](https://doi.org/10.1007/s00340-013-5418-1) · [arXiv:1212.2911](https://arxiv.org/abs/1212.2911) | Photon-assisted orbital phases independent of internal atomic structure; four-site single-particle plaquette dynamics | Supports the synthetic-flux control concept |
| **[F25]** A. Muzi Falconi et al., *Microsecond-Scale High-Survival and Number-Resolved Detection of Ytterbium Atom Arrays*, Physical Review Letters **135**, 203402 (2025). [Article](https://doi.org/10.1103/n3bg-7yw7) · [arXiv:2507.01011](https://arxiv.org/abs/2507.01011) | Fast number-resolved imaging and atom survival | Supports those detector capabilities; arbitrary retained internal-spin coherence requires its own evidence |
| **[N24]** Y. Nakamura et al., *Hybrid Atom Tweezer Array of Nuclear Spin and Optical Clock Qubits*, Physical Review X **14**, 041062 (2024). [Article](https://doi.org/10.1103/PhysRevX.14.041062) · [arXiv:2406.12247](https://arxiv.org/abs/2406.12247) | Nuclear-spin coherence tested while imaging a different isotope | Demonstrates a relevant coherence test in a dual-isotope architecture |

These are component-level comparisons. Combining their capabilities into the particular spinful flux ring and coherent retained detector is an additional apparatus question. The analytical result specifies that detector requirement explicitly.

## Sampling comparison

Use the common $\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$ normalization and count $n$ original independent input trials. The [separable benchmark](research/SEPARABLE_BENCHMARK.md) proves exact minimax coefficient $9/16=0.5625$ for universally unbiased estimation with fully separable effects on each three-spin copy, including classical adaptation between copies. Six uniformly randomized Pauli-product settings attain it. The attained coefficient is already implicit in [R23], Eq. (11) and note [42], after rescaling its variance bound $3/n$ by $3/16$; the additional result here is the lower bound over the full declared separable-effects class. Its proof applies elementary Cauchy–Schwarz and classical Fisher-information arguments to the product-state chirality bound; no new general statistical principle or exhaustive priority result is claimed.

Applying [N26] to one boson per input with arbitrary correlated qubit spins gives a spectral measurement of $\chi$: the cyclic eigenphases produce scores $0,\pm\sqrt3/4$, hence worst-case variance coefficient $3/16$. This chirality consequence is a deduction here from their theorem, not an explicit claim of that source. It requires access to the count statistic above and supplies a stronger sampling reference outside the separable-effects class.

For complete parity images at $U/J=16$, $Jt=\pi$, the [rational certificate](research/CERTIFICATE.md) proves a worst-case variance coefficient below $0.389156263913$ and a mixed-input coefficient in $(0.284895494200,0.284895494206)$. This gives strict minimax sampling advantage over the declared separable-effects class. The controlled-cycle Hadamard test has coefficient $0.25$ on the mixed input and uses an entangling measurement outside that class. These ideal comparisons do not equate hardware costs or establish pointwise superiority on every input.

The return-only four-setting example has mixed-input coefficient approximately $0.8759835$ under its stated allocation. Its reconstruction remains valid for every fixed finite positive repulsion with a suitable pair of durations. Complete-image scoring extracts additional information from the detector record, with the separate loss and calibration conditions in [the noise budget](research/NOISE_BUDGET.md).

For derivations and the exact allocation, see [the resource account, §6](source/resource-review-2026-10-10/READOUT_AND_RESOURCES.md#6-compare-the-correct-tasks-and-input-trial-normalization). The [preserved source reading record](source/resource-review-2026-10-10/SOURCES.md) documents the original targeted comparison. The present references refine attribution without asserting exhaustive priority clearance.
