# Reproduce the results

The root verifier runs all 16 preserved scientific groups, four complete-parity numerical cross-checks, one independent integer/rational certificate, and three channel-calibration cross-checks. It verifies the bytes and membership of the preserved research package before and after execution, and records hashes of the current verification code.

## Run

Use Python 3.12 and the pinned numerical dependencies:

```bash
python -m pip install -r requirements.txt
python verify.py --output-dir ../vacancy-verification
```

The output directory must be new and outside the repository. The wrapper supplies single-thread numerical-library settings, disables bytecode writes, and keeps Python assertions enabled in subprocesses. Existing outputs are never overwritten.

The new exact certificate can also run independently using only the Python standard library:

```bash
python research/parity_certificate.py --output ../parity-certificate/certificate.json
```

## Twenty-four scientific check groups

| Groups | Implementation | Coverage |
|---:|---|---|
| 3 | [Channel-calibration cross-checks](research/check_channel_calibration.py) | Product-reference identities; full noisy-effect reconstruction; calibrated mean and variance; a symmetry-breaking counterexample |
| 4 | [Complete-parity cross-checks](research/check_complete_parity.py) | Physical opposite-flux effects, exact reconstruction identities checked numerically, symmetry-block variance, separable benchmark attainment, and loss/detector probability budgets |
| 1 | [Complete-parity exact certificate](research/parity_certificate.py) | Rational/integer enclosure of propagator, calibration rank, exactly corrected scores, second moments, sufficient neighborhoods, and a uniform noise window |
| 5 | [Preserved resource checker](source/resource-review-2026-10-10/check_resources.py) | Exact return-parity projector; complete number-block loss evolution; detector dephasing and recurrence; matched variance; detection-error budgets |
| 5 | [Preserved finite-repulsion checker](source/resource-review-2026-10-10/prior/check_followup.py) | Complete bosonic Hamiltonian and independent constructions; accepted-operation bounds; two-observable contrast and calibration; finite-round filtering; first-order diagnostic |
| 1 | [Preserved short-time certificate](source/resource-review-2026-10-10/prior/exact_short_time.py) | Integer/rational coefficient arithmetic proving the leading chirality and nuisance terms |
| 5 | [Preserved hard-core checker](source/resource-review-2026-10-10/prior/prior/check_pilot.py) | Cycle algebra, return filter, entanglement witness, repeated filtering, disorder and finite-repulsion controls |

The floating-point [complete-parity model](research/complete_parity.py) uses the preserved [readout](source/resource-review-2026-10-10/readout.py) and [finite-repulsion](source/resource-review-2026-10-10/prior/finite_repulsion.py) implementations. Their fixed-three-particle calculation uses the complete 120-dimensional bosonic sector, without truncating doublon or triplon states.

The channel-calibration diagnostic includes uniform particle attenuation after the hopping interval and parity-bit errors. This terminal-loss model is distinct from the continuous-time Lindblad model in the preserved checks. Its displayed numerical examples are not interval-certified operating regions. The channel theorem itself follows analytically from the stated symmetries.

The new exact checker constructs its own Hamiltonian in the complete 40-state fixed-magnetization sector. Collective spin symmetry makes its three input columns sufficient to certify the full eight-dimensional input effects. Its numerical bounds use no floating matrix exponential or eigensolver. The [certificate account](research/CERTIFICATE.md) explains the enclosure, while the [separable benchmark](research/SEPARABLE_BENCHMARK.md) is an analytical theorem; the numerical checker verifies its saturating local protocol.

## Read the output

The root `RECEIPT.json` separates scientific assertions, new certificate completion, preserved-report byte identity, and input integrity. Overall `status: PASS` requires all 24 groups to succeed, the preserved reports to remain byte-identical, and both preserved-source and current-code integrity checks to pass.

The old numerical evidence remains under `source-run/`, with wrapper output in `source-verifier.log`. New diagnostic, channel-calibration, and exact-certificate reports and logs are written alongside it. The receipt gives their paths, return codes, and current-code hashes.

A different numerical-library environment can change floating-point report bytes while leaving scientific assertions true. Inspect the separate fields and recorded differences in that case; preserve the original evidence. Floating-point agreement supplements the proofs. Only the specifically enclosed quantities in the exact certificates have validated numerical intervals. No verifier certifies laboratory parameters.

## Recorded verification

The [complete-parity integration receipt](verification/COMPLETE_PARITY_RECEIPT.json) records **PASS** for all 24 groups. The compact [exact-certificate report](verification/complete-parity-certificate.json) and [numerical cross-check report](verification/complete-parity-diagnostics.json) record the parity results. The [channel-calibration report](verification/channel-calibration.json) checks the three-reference reconstruction. The receipt binds them to the checked implementation through its input hashes. These curated reports were copied from a completed external verification run.

The earlier [16-group integration receipt](verification/RECEIPT.json) is retained as historical evidence. The imported scientific reports remain byte-identical under Python **3.12.14**, NumPy **2.3.5**, and SciPy **1.17.0**. The [import manifest](source/IMPORT_MANIFEST.json) still covers all 79 incoming files, including the original root manifest. No preserved report or source file is regenerated in place.
