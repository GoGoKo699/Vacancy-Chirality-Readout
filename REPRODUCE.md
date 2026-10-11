# Reproduce the results

The verifier runs the complete finite model, the exact short-time coefficient certificate, and the detector/resource checks. It also verifies the bytes and membership of the preserved research package before and after execution.

## Run

Use Python 3.12 and the pinned numerical dependencies:

```bash
python -m pip install -r requirements.txt
python verify.py --output-dir ../vacancy-verification
```

The output directory must be new and outside the repository. The wrapper supplies single-thread numerical-library settings, disables bytecode writes, and keeps Python assertions enabled in its subprocesses. Existing outputs are never overwritten.

## Sixteen scientific check groups

| Groups | Implementation | Coverage |
|---:|---|---|
| 5 | [Resource checker](source/resource-review-2026-10-10/check_resources.py) | Exact parity projector; complete number-block loss evolution; detector dephasing and recurrence; matched variance; detection-error budgets |
| 5 | [Finite-repulsion checker](source/resource-review-2026-10-10/prior/check_followup.py) | Complete bosonic Hamiltonian and independent constructions; accepted-operation bounds; two-observable contrast and calibration; finite-round filtering; first-order diagnostic |
| 1 | [Exact short-time certificate](source/resource-review-2026-10-10/prior/exact_short_time.py) | Integer/rational coefficient arithmetic proving the leading chirality and nuisance terms |
| 5 | [Hard-core checker](source/resource-review-2026-10-10/prior/prior/check_pilot.py) | Cycle algebra, return filter, entanglement witness, repeated filtering, disorder and finite-repulsion controls |

The helper models are [readout.py](source/resource-review-2026-10-10/readout.py) and [finite_repulsion.py](source/resource-review-2026-10-10/prior/finite_repulsion.py). The fixed-three-particle calculation uses the complete 120-dimensional bosonic sector. It does not truncate doublon or triplon states.

## Read the output

The root `RECEIPT.json` reports three separate checks:

| Field | Meaning |
|---|---|
| `scientific_assertions_pass` | All sixteen check groups completed and their assertions passed |
| `canonical_reports_byte_identical` | Generated scientific JSON reports exactly match the preserved numerical reports |
| `source_integrity_pass` | All 79 imported files match their manifest before and after the run, with no change to the import manifest |

Overall `status: PASS` requires all three checks and successful completion of the original verifier. The archived wrapper log is `source-verifier.log`; full generated reports and nested receipts are under `source-run/`.

A different numerical-library environment can change floating-point report bytes while leaving scientific assertions true. Inspect the separate fields and recorded differences in that case; preserve the original evidence. Numerical agreement supplements the analytical proofs and does not certify laboratory parameters or floating-point intervals.

## Recorded verification

The [integration receipt](verification/RECEIPT.json) records **PASS** for all sixteen groups, canonical-report identity, and source integrity. It was produced with Python **3.12.14**, NumPy **2.3.5**, and SciPy **1.17.0** using the root verifier. All four scientific JSON reports were byte-identical to their preserved references.

The [import manifest](source/IMPORT_MANIFEST.json) covers the complete incoming package, including its original root manifest. New verification outputs remain outside the immutable source snapshot. Only the compact integration receipt is included here; the original detailed numerical evidence is already present with the source records.
