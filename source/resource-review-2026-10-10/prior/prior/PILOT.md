# A vacancy as a native three-spin chirality filter

10 October 2026. A fresh, bounded research pilot. No repository or earlier scientific project is accessed or modified. The historical electronic-decoherence attachment is not the scientific basis of this calculation.

## Result in one paragraph

Three spin-1/2 hard-core bosons and one vacancy occupy an isolated four-site ring. Uniform spin-independent hopping with vacancy flux −2π/3, evolution time π/J, and a spin-insensitive projection onto return of the vacancy implement the exact spin Kraus operator K = Π+ + a(I−Π+), where Π+ is the positive scalar-chirality projector and a = [cos(√3π)−1]/2. Thus the false-return probability on the other chirality sectors is r = a² = 0.02786714004…. Starting with three maximally mixed spins, one return heralds a genuinely tripartite-entangled mixed state with probability 0.2709003550…, witnessed by scalar chirality 0.3884692918… above the biseparable bound 1/4. The spin-charge permutation mechanism and chirality measurement through a cycle are established; the candidate additional result is a native vacancy-return implementation and its explicit instrument, not a new entanglement witness or a new general measurement principle.

## 1. Physical resources and conventions

Label the ring clockwise by j = 0,1,2,3. At a site the hard-core basis is |0>, |↑>, |↓>. Exactly three particles are present; the complete Hilbert space has dimension 4×2³ = 32. Initially site 0 is empty and the spin density operator on sites 1,2,3 is arbitrary ρ. The entanglement-generation example uses ρ = I8/8. Set ℏ = 1.

For clarity define the phase for **vacancy** motion, which is opposite to the corresponding atom-motion flux:

\[
H_\Phi=-J\sum_{j=0}^3\sum_{\sigma=\uparrow,\downarrow}
\left[e^{i\Phi/4}|\sigma_j,0_{j+1}\rangle\langle0_j,\sigma_{j+1}|+\mathrm{h.c.}\right].
\tag{1}
\]

The omitted sites are unchanged, and indices are modulo four. The motion does not rotate or resolve the spin. There is no separate spin-exchange term, particle loss, feedback during evolution, or Hamiltonian chirality term. This is the strict hard-core, exchange-frozen model, not an exact finite-repulsion Hubbard model. Uniform hopping defines the explicit calibrated example; Section 6 permits static spin-blind scalar disorder.

The measurement asks whether the original occupancy pattern has returned. In the strict model, detecting the vacancy at 0 is equivalent to occupancies (0,1,1,1). On this outcome the spin state again belongs to the original three occupied sites. The measurement must preserve the spins to retain the prepared entanglement; repeated filtering additionally requires nondestructive position/occupancy readout. Ordinary destructive imaging is not silently granted those properties. The localized vacancy and postselection are resources; maximally mixed initial spins do not make the entire input maximally mixed.

## 2. Spin cycle and scalar chirality

Write

\[
C|s_1s_2s_3\rangle=|s_2s_3s_1\rangle,
\qquad \chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3),
\qquad \mathbf S=\boldsymbol\sigma/2.
\]

A clockwise vacancy circuit makes exactly this cyclic spin permutation. It follows either by tracking the four allowed vacancy hops or by multiplying the corresponding swaps. Since C = P23 P12,

\[
C-C^\dagger=4i\chi,\qquad C^3=I. \tag{2}
\]

The phases θ = 0,+2π/3,−2π/3 of C have projectors

\[
\Pi_\theta=\tfrac13[I+e^{-i\theta}C+e^{-2i\theta}C^2], \tag{3}
\]

of ranks 4,2,2. The scalar chirality eigenvalues are respectively 0,+√3/4,−√3/4. Denote the latter two projectors Π+,Π−. The θ=0 subspace is the symmetric total-spin-3/2 subspace. The two nonzero chirality sectors each retain a spin-1/2 doublet; selecting one chirality does not prepare a particular spin orientation or a pure W state.

Equations (2)–(3) are standard spin-permutation algebra. Reascos et al. [R23] use the normalized operator χhat = 4χ/√3, so their eigenvalues are 0,±1 and their special unitary is C† with the present orientation. Kesharpu–Kochetov [K22] use the reverse cyclic shift. Neither sign convention changes the physics; the flux selecting Π+ here is −2π/3.

## 3. Exact return instrument

Transport the spin labels along the vacancy's tree path 0→1→2→3. In this basis the first three edges act trivially on the spin register, while the closing edge carries C. A gauge transformation moves all Peierls phases to the closing edge. Restricted to C=e^{iθ}, the charge problem is the ordinary four-site ring with flux Φ+θ.

For uniform hopping, its return amplitude is

\[
A(t;f)=\frac14\sum_{n=0}^3
\exp\!\left[2iJt\cos\frac{2\pi n+f}{4}\right]
=\frac12\left[\cos\!\left(2Jt\cos\frac f4\right)+
\cos\!\left(2Jt\sin\frac f4\right)\right]. \tag{4}
\]

This amplitude is real, even in f, and 2π periodic. Therefore the full spin return Kraus operator is

\[
K_\Phi(t)=\langle0_{\rm hole}|e^{-itH_\Phi}|0_{\rm hole}\rangle
=\sum_\theta A(t;\Phi+\theta)\Pi_\theta. \tag{5}
\]

This is an identity for all initial spin states, including mixed states and states entangled with an external spectator. It is not obtained by evaluating a single chosen input.

At τ=Jt=π,

\[
A(t;0)=1,\qquad A(t;\pm2\pi/3)=a:=\frac{\cos(\sqrt3\pi)-1}{2}.
\]

Consequently the two flux settings give

\[
K_\pm:=K_{\mp2\pi/3}(\pi/J)=\Pi_\pm+a(I-\Pi_\pm),\qquad
E_\pm=K_\pm^\dagger K_\pm=rI+(1-r)\Pi_\pm,
\quad r=a^2. \tag{6}
\]

The plus/minus label denotes **selected chirality**, not flux. Numerically,

\[
a=-0.166934538198736\ldots,
\qquad r=0.0278671400436253\ldots.
\]

The desired sector returns with probability one. Its internal spin-1/2 state is unchanged. The other sectors return with probability r, not zero: a single return is a filter, **not an exact projective measurement of chirality**. A unitary-evolution-and-position-readout implementation remains quantum coherent; no individual hopping history has been measured.

### Measuring the chirality expectation

Let p± be the return frequencies obtained on independently prepared copies using the two flux signs. Then

\[
\boxed{\langle\chi\rangle_\rho
=\frac{\sqrt3}{4(1-r)}(p_+-p_-).} \tag{7}
\]

No two-spin correlator or local magnetization appears in this expression. The distinction between the measurement of an unknown input and heralded state generation matters: Equation (7) concerns ensemble frequencies on copies; it does not assign an exact chirality to every one-shot return.

## 4. Entanglement from a maximally mixed spin input

For input I8/8 and selected positive chirality, the unnormalized output is E+/8. Its success probability and normalized state are

\[
p_1=\frac{1+3r}{4},\qquad
\rho_1=\frac{\Pi_++r(I-\Pi_+)}{2+6r}. \tag{8}
\]

The numerical values are p1=0.270900355032719… and

\[
\langle\chi\rangle_{\rho_1}
=\frac{\sqrt3}{4}\frac{1-r}{1+3r}
=0.388469291814954\ldots. \tag{9}
\]

**Every biseparable three-qubit state obeys |<χ>|≤1/4.** To prove this, take a product state across 1|23. Rotational invariance allows the first spin's Bloch vector to lie along z with length at most 1/2. The two-spin operator S2x S3y−S2y S3x has norm 1/2, giving the bound 1/4. The same holds for either other cut; convexity gives the bound for arbitrary mixtures over cuts. A polarized first spin times (|↑↓>+i|↓↑>)/√2 saturates it. Thus (9) establishes genuine tripartite entanglement, not merely entanglement across one bipartition. The witness and its normalized equivalent in [R23] are inherited, not a new criterion.

Interestingly, all individual spin marginals are I2/2 and every two-spin marginal is separable. Indeed,

\[
\rho_1=w\,\Pi_+/2+(1-w)I_8/8,\qquad
w=\frac{1-r}{1+3r}.
\]

SU(2) and cyclic symmetry imply that a pair marginal is the two-qubit Werner state with singlet-mixture parameter w/3≤1/3. Its partial transpose is positive, and it is separable. Thus the generated state has genuine three-spin entanglement even though neither an individual polarization nor pairwise entanglement detects it. This statement concerns this explicit hard-core output, not arbitrary finite-U outputs.

## 5. Repeated return filtering and its true success probability

After a successful return, apply the same evolution again and again, retaining a trial only while every readout reports the original vacancy position. Nothing is reinitialized within a surviving trial. With m successful rounds,

\[
K_+^m=\Pi_++a^m(I-\Pi_+). \tag{10}
\]

For an arbitrary initial state with q=Tr(Π+ρ), the probability of the entire successful record is q+(1−q)r^m and its target-sector weight is q/[q+(1−q)r^m]. The measurement cannot create target-sector weight from q=0. For I8/8,

\[
p_m=\frac{1+3r^m}{4},\qquad
\rho_m=\frac{\Pi_++r^m(I-\Pi_+)}{2+6r^m},\qquad
\operatorname{Tr}(\Pi_+\rho_m)=\frac1{1+3r^m}. \tag{11}
\]

| Successful rounds | Total evolution Jt | Total success per original trial | Selected-sector weight |
|---|---:|---:|---:|
| 1 | π | 0.2709003550 | 0.9228485506 |
| 2 | 2π | 0.2505824331 | 0.9976756826 |
| 4 | 4π | 0.2500004523 | 0.9999981908 |

The success rate tends to 1/4, not to (1/4)^m, because later outcomes act on an already filtered state. The limit is the **mixed state Π+/2**, not a pure state. All four-round headline numbers presuppose ideal, nondestructive vacancy readout and the frozen-spin hard-core dynamics through all rounds.

## 6. Exact structure under static spin-independent disorder

There is a useful structural control that does not require equal bond magnitudes. Give the four bonds arbitrary fixed real magnitudes Jj and the vacancy arbitrary real onsite potentials vj. Keep a single ring, no particle crossings, no spin dynamics, and a known total vacancy flux Φ. The tree-basis reduction still gives

\[
K_\Phi(t)=\sum_\theta A(t;\Phi+\theta)\Pi_\theta,
\]

where A is now the return amplitude of the corresponding disordered scalar four-site ring. Reciprocity yields Hcharge(−f)=Hcharge(f)^T after a gauge choice, and hence

\[
A(t;-f)=A(t;f).
\]

At Φ=∓2π/3 one therefore still has exactly

\[
K_\pm(t)=A(t;0)\Pi_\pm+A(t;2\pi/3)(I-\Pi_\pm). \tag{12}
\]

The return contrast is

\[
p_+-p_-=
\frac4{\sqrt3}\left(|A(t;0)|^2-|A(t;2\pi/3)|^2\right)\langle\chi\rangle.
\tag{13}
\]

More generally arbitrary opposite nonzero flux settings also isolate χ, with their explicitly calibrated scalar difference. Static spin-blind bond or onsite disorder therefore changes the calibration and filter selectivity, but introduces **no other spin observable** into that flux-reversed difference. This is not a disorder-independent contrast or a topological protection result. The contrast may be small or zero; timing or parameter errors can destroy useful filtering. Time-dependent noncommuting controls do not automatically inherit the transpose argument.

Removing one ring bond is a decisive control: all fluxes become gauge-removable and every closed vacancy path restores the spin order. The return operator is scalar on the spin register, so this readout cannot filter chirality or entangle the maximally mixed spins. At zero flux the two handed sectors also have equal return probability. A nonzero external orbital flux supplies the handedness reference.

## 7. Finite-repulsion diagnostic: generation survives before ideal readout does

To test the assumption carrying the interpretation, separately assemble the finite-U spinful Bose-Hubbard Hamiltonian

\[
H_U=-J\sum_{j,\sigma}
[e^{i\Phi/4}b^\dagger_{j\sigma}b_{j+1,\sigma}+\mathrm{h.c.}]
+\frac U2\sum_j n_j(n_j-1). \tag{14}
\]

There are three particles in eight modes: its full Hilbert space has dimension C(10,3)=120. The matrix calculation uses that whole sector, with no occupation truncation inside it. In this control the accepted readout is explicitly **(0,1,1,1)**; allowing a doublon and simply calling site 0 vacant would change the measurement. Evolve for π/J with flux +2π/3 and the same maximally mixed spin preparation.

| U/J | Occupancy-return probability | Magnitude of signed output chirality |
|---:|---:|---:|
| 8 | 0.20243399997 | 0.18519038040 |
| 16 | 0.27169237065 | 0.30451329715 |
| 32 | 0.27507594737 | 0.36216748808 |
| 64 | 0.27228318704 | 0.38141609054 |
| 128 | 0.27126864655 | 0.38667317535 |
| Strict hard core | 0.27090035503 | 0.38846929181 |

These are finite-matrix diagnostics, not a uniform finite-U theorem, a certified numerical interval, or a threshold determination. The witness exceeds 1/4 at the displayed U/J≥16 samples; its failure at U/J=8 proves only that this witness does not certify entanglement there.

The exact arbitrary-input readout identity does **not** survive automatically. The finite-U return block no longer commutes exactly with C. At U/J=16, the difference of the two flux-reversed return POVM elements has a nonchirality remainder of operator norm 0.0460534… after projecting it onto χ in the Hilbert–Schmidt inner product. Consequently a scalar recalibration of Equation (7) is not an exact cure for every unknown input. The distinction between a robust generated-state example and a valid state-independent measurement is essential.

The data at finite repulsion motivate, but do not complete, an instrument-level error analysis. A microscopic apparatus also has to justify spin-independent Peierls phases, isolation of the four sites, negligible additional spin evolution, and spin-preserving occupancy readout. The ideal model does not assert that an existing apparatus simultaneously supplies all of them.

## 8. Strongest predecessor comparison

**[K22] K. K. Kesharpu and E. A. Kochetov**, *Spin-flip-induced superfluidity in a ring of spinful hard-core bosons*, Physical Review A 106, 033323 (2022), arXiv:2209.06149. Their main-text Eq. (2) already places the cyclic spin permutation on the closing hopping link, and Eqs. (6)–(8) turn its eigenphase into a charge-sector flux. This covers the structural spin-charge reduction used here. Setting three particles on four sites and adding an external flux leads directly to the scalar ring spectrum underlying (4). We cannot claim a new spin-generated gauge field or a new spin-charge separation principle.

**[R23] L. I. Reascos, B. Murta, E. F. Galvão and J. Fernández-Rossier**, *Quantum circuits to measure scalar spin chirality*, Physical Review Research 5, 043087 (2023), arXiv:2306.14804. Their Eq. (1) uses χhat=4χ/√3. Eqs. (5)–(6) and the accompanying text use a two-SWAP cyclic unitary in a Hadamard test. Their qutrit-QPE discussion, Eqs. (12)–(13), already provides projective chirality readout and explicitly proposes preparing well-defined chirality sectors by the measurement. Thus indirect chirality measurement, its entanglement-witness interpretation, and chirality-sector state preparation are inherited. Our finite-time return filter is not a superior ideal projective measurement compared with their circuit; it uses different native controls.

**[A13] M. Aidelsburger et al.**, *Experimental realization of strong effective magnetic fields in optical superlattice potentials*, Applied Physics B 113, 1–11 (2013), arXiv:1212.2911. Its primary abstract describes widely tunable photon-assisted gauge fields that do not rely on atomic internal structure. This supports the relevance of spin-independent orbital flux as a control concept, not validation of our complete three-particle filtering apparatus. The earlier PRL 107, 255301 (2011), arXiv:1110.5314, reports local tunneling-phase and single-atom cyclotron observations; those experiments are not demonstrations of the present filter.

Source reading: primary parsed bodies of [K22] and [R23], especially the equations identified above, were inspected. The web screenshot of the [R23] circuit page and an independent PDF download failed. No circuit-layout interpretation, plotted datum, or device performance was inferred from those failed visual retrievals. [A13] and the 2011 experiment were used only at the primary abstract level. Targeted searches for vacancy/hole-return chirality measurement, spin-cycle readout, heralded chirality preparation, and four-site flux plaquettes did not identify this exact native return instrument. This is **not exhaustive priority clearance**. Nearby quantum-dot chirality and generic mobile-impurity interferometry references are not treated as equivalent without a model-matched comparison.

After attribution, the candidate remainder is: a four-site, spin-independent native implementation of a known spin-cycle measurement, with an exact return Kraus operator, calibrated chirality contrast, a nonvanishing-yield repeated filter, and a clear boundary between static scalar disorder and finite-U particle-crossing errors. Whether that implementation is a substantial physical improvement over established measurement options remains unresolved.

Primary URLs:
- [K22] https://arxiv.org/pdf/2209.06149 ; https://doi.org/10.1103/PhysRevA.106.033323
- [R23] https://arxiv.org/pdf/2306.14804 ; https://doi.org/10.1103/PhysRevResearch.5.043087
- [A13] https://arxiv.org/abs/1212.2911 ; https://doi.org/10.1007/s00340-013-5418-1
- 2011 experiment: https://arxiv.org/abs/1110.5314 ; https://doi.org/10.1103/PhysRevLett.107.255301

## 9. Allocation and reproducibility

**Retain for one bounded continuation, not repository promotion.** The ideal protocol gives a concrete positive operation rather than another singular-velocity or approximation-failure example. Its readout needs neither controlled-SWAP gates nor spin-resolved measurement in the declared native model. But cyclic spin flux and chirality projections are known, and the ideal resources are substantive. The next discriminating question is an instrument-level finite-repulsion bound, keeping the same native controls and occupancy readout, and assessing whether the claimed hardware simplification survives that comparison. The generated-state example alone cannot establish an accurate measurement on arbitrary unknown inputs.

`check_pilot.py --output /path/to/a/new/report.json` runs five groups: spin algebra/witness, full 32-dimensional return instrument, repeated heralding and reduced states, arbitrary static scalar disorder and broken-ring controls, and a separately built full 120-dimensional finite-U control. It refuses to overwrite an existing report. All formal groups passed on the first run and on repetition. No scientific assertion, tolerance, or equation was changed to obtain a pass. The two reports are compared in RECORD.json. Exploratory scripts are preserved under development; their optional long-time recurrence calculations are not part of the proposed result. No old scientific checker was imported or rerun, no manuscript was written, and no repository was created or modified.
