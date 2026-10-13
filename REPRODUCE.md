# Reproduce the results

The root verifier runs all 16 preserved scientific groups, four complete-parity numerical cross-checks, one independent integer/rational certificate, three channel-calibration cross-checks, three calibration-cost cross-checks, five passive-parity cross-checks, two finite-calibration checks, and five parameter-robustness checks. It verifies the bytes and membership of the preserved research package before and after execution, and records hashes of the current verification code.

## Run

Use Python 3.12 and the pinned numerical dependencies:

```bash
python -m pip install -r requirements.txt
python -B verify.py --output-dir ../vacancy-verification
```

The output directory must be new and outside the repository. The wrapper supplies single-thread numerical-library settings, disables bytecode writes, and keeps Python assertions enabled in subprocesses. Existing outputs are never overwritten.

The new exact certificate can also run independently using only the Python standard library:

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
| 5 | [Preserved resource checker](source/resource-review-2026-10-10/check_resources.py) | Exact return-parity projector; complete number-block loss evolution; detector dephasing and recurrence; matched variance; detection-error budgets |
| 5 | [Preserved finite-repulsion checker](source/resource-review-2026-10-10/prior/check_followup.py) | Complete bosonic Hamiltonian and independent constructions; accepted-operation bounds; two-observable contrast and calibration; finite-round filtering; first-order diagnostic |
| 1 | [Preserved short-time certificate](source/resource-review-2026-10-10/prior/exact_short_time.py) | Integer/rational coefficient arithmetic proving the leading chirality and nuisance terms |
| 5 | [Preserved hard-core checker](source/resource-review-2026-10-10/prior/prior/check_pilot.py) | Cycle algebra, return filter, entanglement witness, repeated filtering, disorder and finite-repulsion controls |

The floating-point [complete-parity model](research/complete_parity.py) uses the preserved [readout](source/resource-review-2026-10-10/readout.py) and [finite-repulsion](source/resource-review-2026-10-10/prior/finite_repulsion.py) implementations. Their fixed-three-particle calculation uses the complete 120-dimensional bosonic sector, without truncating doublon or triplon states.

The channel-calibration diagnostic includes uniform particle attenuation after the hopping interval and parity-bit errors. This terminal-loss model is distinct from the continuous-time Lindblad model in the preserved checks. Its displayed numerical examples are not interval-certified operating regions. The channel theorem itself follows analytically from the stated symmetries.

The calibration-cost checker reads the exact-certificate report produced earlier in the same run. It independently compares its nominal moments with the numerical model and exercises the confidence and counting algebra. Its finite planning profile substitutes model moments for hypothetical empirical moments, including the unbiased sample-variance factors. It is not an experimental record or a guaranteed-pass sample budget. The [cost account](research/CALIBRATION_COST.md) separates confidence-qualified conditional risk, unconditional risk, and leading asymptotic forecasts. The calculator is [calibration_cost.py](research/calibration_cost.py); its floating-point results are not interval certificates.

The same checker tests full-cycle reuse and stability boundaries, adds a specified operator-drift allowance, and constructs valid heterogeneous POVMs to check matched-mixture calibration and its failure under state-channel correlation. These are analytical cross-checks for [repeated scans](research/REPEATED_STATE_SCANS.md), not measurements of cycle time, spatial independence, or laboratory stability.

The new exact checker constructs its own Hamiltonian in the complete 40-state fixed-magnetization sector. Collective spin symmetry makes its three input columns sufficient to certify the full eight-dimensional input effects. Its numerical bounds use no floating matrix exponential or eigensolver. The [certificate account](research/CERTIFICATE.md) explains the enclosure, while the [separable benchmark](research/SEPARABLE_BENCHMARK.md) is an analytical theorem; the numerical checker verifies its saturating local protocol.

The [passive bound](research/PASSIVE_PARITY_BOUND.md) is an analytical theorem;
the numerical checker exercises its operator identities, information bounds,
constructions, and finite constants. Random settings are cross-checks, not
a proof of the universal bound.

The [parameter-robustness checker](research/check_parameter_robustness.py)
reads the nominal exact certificate and records its SHA-256 hash. Its five
exact-arithmetic groups verify moment inputs, recalibration bounds, weighted
parameter constants, fixed-score finite risk, joint parameter/noise bounds,
and microscopic-control/loading consequences. The joint group certifies
norms and polynomial moments throughout its smaller window and checks that
they lie below the existing finite-calibration inputs. The operator and
statistical proofs are in [PARAMETER_ROBUSTNESS.md](research/PARAMETER_ROBUSTNESS.md).
To run it separately, supply `--certificate` and `--output` with the nominal
report and a new external output path.

The [finite-calibration checker](research/check_finite_calibration.py) runs
after that robustness check. It reads both newly produced reports, checks
their common nominal dependency, validates the joint source window and all
imported bounds, and records both report hashes. Its rational arithmetic
certifies the sufficient finite examples in [FINITE_CALIBRATION.md](research/FINITE_CALIBRATION.md),
including their extension to the smaller joint window. The underlying
statistical, adaptive-stopping, and direct-sum statements are proved in that
document. To run it separately, supply `--certificate`,
`--robustness-certificate`, and `--output` with the corresponding report paths.
The separate ideal-detector and fixed-score allowances cannot be added to
the explicit joint window.

## Read the output

The root `RECEIPT.json` separates scientific assertions, certificate completion, preserved-report byte identity, and input integrity. Overall `status: PASS` requires all 39 groups to succeed, the preserved reports to remain byte-identical, and both preserved-source and current-code integrity checks to pass.

The preserved numerical evidence is written under `source-run/`, with wrapper output in `source-verifier.log`. Diagnostic, channel-calibration, calibration-cost, exact-certificate, passive-parity, finite-calibration, and parameter-robustness reports and logs are written alongside it. The receipt gives their paths, return codes, and current-code hashes.

A different numerical-library environment can change floating-point report bytes while leaving scientific assertions true. Inspect the separate fields and recorded differences in that case; preserve the original evidence. Floating-point agreement supplements the proofs. Only the specifically enclosed quantities in the exact certificates have validated numerical intervals. No verifier certifies laboratory parameters.

## Recorded verification

The [scientific verification receipt](verification/SCIENTIFIC_RECEIPT.json) records **PASS** for all 39 groups. The [exact-certificate report](verification/complete-parity-certificate.json) supplies both scores, noise windows, and finite-calibration inputs. The [numerical cross-check report](verification/complete-parity-diagnostics.json) checks both scores independently in the full 120-state model. The [passive-parity report](verification/passive-parity.json), [finite-calibration report](verification/finite-calibration.json), and [parameter-robustness report](verification/parameter-robustness.json) record the measurement-class comparisons and their distinct finite guarantees. The [channel-calibration report](verification/channel-calibration.json) and [calibration-cost report](verification/calibration-cost.json) cover the three-reference and conditional-validation methods. These reports come from the same completed external verification run; the receipt records code hashes, and both dependent exact checkers bind their nominal certificate by hash, and the finite-calibration checker also binds the joint-robustness report.

The [37-group](verification/FIXED_RESOURCE_RECEIPT.json), [33-group](verification/PASSIVE_SEPARATION_RECEIPT.json), [27-group](verification/COMPLETE_PARITY_RECEIPT.json), and [16-group](verification/RECEIPT.json) receipts retain evidence for their recorded code versions. The imported scientific reports remain byte-identical under Python **3.12.14**, NumPy **2.3.5**, and SciPy **1.17.0**. The [import manifest](source/IMPORT_MANIFEST.json) covers all 79 incoming files, including the original root manifest. No preserved report or source file is regenerated in place.
