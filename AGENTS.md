# Repository scope and working rules

The only authorized repository is `GoGoKo699/Vacancy-Chirality-Readout`:
https://github.com/GoGoKo699/Vacancy-Chirality-Readout

Before editing, check the Git remote and current branch. Do not modify any other repository.

## Canonical reading order

1. `README.md`
2. `THEORY.md`
3. `REFERENCES.md`
4. `REPRODUCE.md`
5. `research/COMPLETE_PARITY_THEORY.md`, `research/SEPARABLE_BENCHMARK.md`,
   `research/CERTIFICATE.md`, `research/CHANNEL_CALIBRATION.md`, and `research/NOISE_BUDGET.md` for the two-setting result.
6. The linked proof and implementation records needed for the task.

This project concerns a four-site vacancy-return implementation of scalar-spin-chirality readout. Retrieval terms include vacancy chirality, spin-independent Bose–Hubbard ring, finite-repulsion parity readout, complete-parity two-flux estimation, separable measurement benchmark, two-time chirality calibration, and coherent retained-state filtering. Use this repository for those results and follow the attribution in `REFERENCES.md`.

## Scientific scope

- Distinguish input-expectation estimation from the accepted quantum operation and from repeated retained-state filtering.
- Preserve the stated spin symmetry and finite-repulsion reflection assumption. Do not extend hard-core disorder statements to arbitrary finite-repulsion disorder.
- State the known particle-number assumption whenever using the exact parity projector.
- Preserve the flux/chirality sign convention and $\mathbf S=\boldsymbol\sigma/2$ normalization.
- Numerical examples and diagnostic noise channels are not measured apparatus parameters. Analytical bounds, floating-point diagnostics, and external experimental evidence must remain distinguishable.
- Keep the certified feasible score distinct from the floating-point mixed-input optimum. The two-setting theorem is certified at a specified finite-repulsion point and sufficient neighborhoods; the four-setting existence theorem covers every fixed positive repulsion.
- The separable benchmark concerns unbiased estimation with fully separable effects on each three-spin copy, including classical adaptation. No collective-measurement, hardware-cost, or pointwise-all-states advantage follows. Noisy estimates require mean-square-error budgets including bias.
- Direct channel calibration assumes spin and reflection covariance and a stable instrument. Separate pilot, validation, and unknown-input data; include reference-preparation errors and known-reference costs. Three product probes calibrate the symmetry class but do not verify that the instrument belongs to it.
- Retain the established-source attribution. Do not claim a new chirality principle, universal hardware advantage, or exhaustive priority result.
- Manuscript drafting and new scientific directions require an author request. Keep journal targets out of repository text.

## Source preservation

`source/resource-review-2026-10-10/` is the immutable incoming research package, including its nested records, code, evidence, and manifests. `source/IMPORT_MANIFEST.json` fixes all 79 incoming files. Do not edit, delete, regenerate, or overwrite those files or revise hashes to conceal a change. Historical instructions and status language inside the snapshot describe its creation; current repository scope and reading order are defined here.

Corrections, new results, and new implementations belong outside the preserved snapshot and must identify what they supersede. Avoid copying the same scientific account into multiple new documents. Exploratory scripts inside `development/` are not canonical results.

## Verification

For scientific or executable changes, run `python verify.py --output-dir /absolute/new/path` with an output directory outside the repository. Do not run with Python optimization enabled. Check assertion status, source integrity, and canonical-report identity separately; do not relabel a byte-identity mismatch as a scientific failure or silently replace old evidence. The root verifier runs 24 groups: 16 preserved groups, four complete-parity numerical cross-checks, three channel-calibration cross-checks, and one integer/rational certificate. It also hashes current verification code before and after execution. `REPRODUCE.md` defines the workflow.

Keep generated outputs outside source directories. Check relative links and `git diff --check` for documentation changes. Preserve unrelated user changes and confirm the remote branch state before merging.
