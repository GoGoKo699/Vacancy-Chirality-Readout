# Reproduce the results

The root verifier runs **39 scientific check groups**, verifies source membership and content, and compares four numerical outputs with fixed reference reports. It records code hashes before and after execution. Analytical proofs establish the universal statements; numerical cross-checks exercise their identities and finite constructions.

## Run

Use Python 3.12 and the pinned numerical dependencies:

```bash
python -m pip install -r requirements.txt
python -B verify.py --output-dir ../vacancy-verification
```

Choose a new output directory outside the repository. The wrapper supplies single-thread numerical-library settings, disables bytecode writes, and keeps Python assertions enabled in subprocesses. Existing outputs are never overwritten. The numerical reference reports use Python **3.12.14**, NumPy **2.3.5**, and SciPy **1.17.0**.

The complete-parity exact certificate can also run independently using only the Python standard library:

```bash
python research/parity_certificate.py --output ../parity-certificate/certificate.json
```

## Thirty-nine scientific check groups

| Groups | Implementation | Coverage |
|---:|---|---|
| 5 | [Passive-parity cross-checks](research/check_passive_parity.py) | Three-mode blindness; four-mode Fisher information and local prerotations; constructive identifiability; finite-risk constants and the matched passive ring; six- through nine-mode Fourier-and-splitter construction |
| 2 | [Finite-calibration checks](research/check_finite_calibration.py) | Exact dependency-bound unconditional risk and total counts; preparation and drift extension, including the certified joint window |
| 5 | [Parameter-robustness checks](research/check_parameter_robustness.py) | Exact-response correction; physical parameter budgets; fixed-score finite MSE; joint parameter/noise and finite-reference inputs; microscopic control and loading, with nominal-certificate dependency checks |
| 3 | [Calibration-cost cross-checks](research/check_calibration_cost.py) | Exact moments against the numerical model; confidence-box propagation; biased separable information/prior bound; finite counts, full-cycle cost, stability, and matched-channel pooling |
| 3 | [Channel-calibration cross-checks](research/check_channel_calibration.py) | Product-reference identities; full noisy-effect reconstruction; calibrated mean and variance; a symmetry-breaking counterexample |
| 4 | [Complete-parity cross-checks](research/check_complete_parity.py) | Physical opposite-flux effects, exact reconstruction identities checked numerically, symmetry-block variance, separable benchmark attainment, and loss/detector probability budgets |
| 1 | [Complete-parity exact certificate](research/parity_certificate.py) | Rational/integer enclosure of propagator, calibration rank, both score families, moments, mixed-input neighborhoods, passive-separation noise windows, and finite-calibration inputs |
| 5 | [Readout resource checker](source/check_resources.py) | Exact return-parity projector; complete number-block loss evolution; detector dephasing and recurrence; matched variance; detection-error budgets |
| 5 | [Finite-repulsion checker](source/check_finite_repulsion.py) | Complete bosonic Hamiltonian and independent constructions; accepted-operation bounds; two-observable contrast and calibration; finite-round filtering; first-order diagnostic |
| 1 | [Short-time certificate](source/exact_short_time.py) | Integer/rational coefficient arithmetic proving the leading chirality and nuisance terms |
| 5 | [Hard-core checker](source/check_hard_core.py) | Cycle algebra, return filter, entanglement witness, repeated filtering, disorder and finite-repulsion controls |


## Models and dependencies

The [complete-parity model](research/complete_parity.py) uses the [readout](source/readout.py) and [finite-repulsion](source/finite_repulsion.py) implementations. The fixed-three-particle calculation includes the complete 120-dimensional bosonic sector, including doublon and triplon states. The hard-core checker independently assembles its Hamiltonians.

The exact checker constructs a separate Hamiltonian in the full 40-state fixed-magnetization sector. Collective spin symmetry makes its three input columns sufficient to certify the eight-dimensional input effects. It uses integer/rational enclosures without a floating matrix exponential or eigensolver; see [the certificate account](research/CERTIFICATE.md).

The channel-calibration diagnostic uses terminal uniform attenuation and parity-bit errors. The readout resource checker uses continuous-time Lindblad loss. These are distinct diagnostic models. The [noise budget](research/NOISE_BUDGET.md) and [parameter theorem](research/PARAMETER_ROBUSTNESS.md) specify the certified channels and joint regions.

The calibration-cost checker reads the exact report produced in the same run and compares its moment bounds with the numerical model. Its planning profiles use hypothetical empirical moments; [calibration cost](research/CALIBRATION_COST.md) distinguishes these forecasts from finite guarantees. The parallel-channel calculations check the operator algebra and sampling assumptions in [repeated scans](research/REPEATED_STATE_SCANS.md).

The parameter-robustness checker reads and hashes the nominal exact report. The finite-calibration checker then reads both reports, verifies their common dependency, and checks the joint-window moment inputs before evaluating unconditional risk and total preparation counts. To run them separately, supply `--certificate`, and for finite calibration also `--robustness-certificate`, along with a new external `--output` path. Their proofs are in [Parameter robustness](research/PARAMETER_ROBUSTNESS.md) and [Finite calibration](research/FINITE_CALIBRATION.md).

## Interpret the output

The root `RECEIPT.json` records each suite's exit status, group count, and report path. Overall `status: PASS` requires all 39 groups to pass, source membership and hashes to match [source/MANIFEST.json](source/MANIFEST.json), code and reference files to remain unchanged during the run, and the four foundational reports to match [verification/reference/](verification/reference/) byte for byte.

Reports and per-suite logs are written into the chosen output directory. Scientific assertion status, numerical-report identity, and input integrity are separate checks. A different numerical-library environment can change report bytes while leaving mathematical assertions true; the receipt records both facts. Reference reports must not be replaced merely to obtain a passing comparison.

Numerical agreement supplements the analytical proofs. Only the specifically enclosed quantities in the exact certificates have validated numerical intervals. Device parameters and calibration data must satisfy the assumptions stated in the relevant theorem.

## Verification evidence

The [complete verification receipt](verification/SCIENTIFIC_RECEIPT.json) records the full run and its code hashes. The associated reports are:

| Report | Contents |
|---|---|
| [Exact certificate](verification/complete-parity-certificate.json) | Both score families, noise windows, and finite-calibration moment inputs |
| [Complete-parity diagnostics](verification/complete-parity-diagnostics.json) | Independent full-120-state comparison |
| [Passive parity](verification/passive-parity.json) | Measurement-class identities, finite constants, and extra-mode constructions |
| [Parameter robustness](verification/parameter-robustness.json) | Control, joint-noise, fixed-score, and loading guarantees |
| [Finite calibration](verification/finite-calibration.json) | Unconditional MSE and total-preparation budgets |
| [Channel calibration](verification/channel-calibration.json) | Reference-state identities, response reconstruction, and symmetry test |
| [Calibration cost](verification/calibration-cost.json) | Conditional validation, resource accounting, and matched-channel pooling |

The four foundational reference reports cover [hard-core filtering](verification/reference/hard-core.json), [finite repulsion](verification/reference/finite-repulsion.json), [short-time coefficients](verification/reference/short-time.json), and [readout resources](verification/reference/readout-resources.json). They are also the fixed byte-comparison targets. The complete receipt records the reproduced reports' hashes.
