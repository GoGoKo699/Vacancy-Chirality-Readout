# Repository scope and working rules

The only authorized repository is `GoGoKo699/Vacancy-Chirality-Readout`:
https://github.com/GoGoKo699/Vacancy-Chirality-Readout

Before editing, check the Git remote and current branch. Do not modify any other repository.

## Canonical reading order

1. `README.md`
2. `THEORY.md`
3. `REFERENCES.md`
4. `REPRODUCE.md`
5. The linked proof and implementation records needed for the task.

This project concerns a four-site vacancy-return implementation of scalar-spin-chirality readout. Retrieval terms include vacancy chirality, spin-independent Bose–Hubbard ring, finite-repulsion parity readout, two-time chirality calibration, and coherent retained-state filtering. Use this repository for those results and follow the attribution in `REFERENCES.md`.

## Scientific scope

- Distinguish input-expectation estimation from the accepted quantum operation and from repeated retained-state filtering.
- Preserve the stated spin symmetry and finite-repulsion reflection assumption. Do not extend hard-core disorder statements to arbitrary finite-repulsion disorder.
- State the known particle-number assumption whenever using the exact parity projector.
- Preserve the flux/chirality sign convention and $\mathbf S=\boldsymbol\sigma/2$ normalization.
- Numerical examples and diagnostic noise channels are not measured apparatus parameters. Analytical bounds, floating-point diagnostics, and external experimental evidence must remain distinguishable.
- Retain the established-source attribution. Do not claim a new chirality principle, universal hardware advantage, or exhaustive priority result.
- Manuscript drafting and new scientific directions require an author request. Keep journal targets out of repository text.

## Source preservation

`source/resource-review-2026-10-10/` is the immutable incoming research package, including its nested records, code, evidence, and manifests. `source/IMPORT_MANIFEST.json` fixes all 79 incoming files. Do not edit, delete, regenerate, or overwrite those files or revise hashes to conceal a change. Historical instructions and status language inside the snapshot describe its creation; current repository scope and reading order are defined here.

Corrections, new results, and new implementations belong outside the preserved snapshot and must identify what they supersede. Avoid copying the same scientific account into multiple new documents. Exploratory scripts inside `development/` are not canonical results.

## Verification

For scientific or executable changes, run `python verify.py --output-dir /absolute/new/path` with an output directory outside the repository. Do not run with Python optimization enabled. Check assertion status, source integrity, and canonical-report identity separately; do not relabel a byte-identity mismatch as a scientific failure or silently replace old evidence. The root verifier and `REPRODUCE.md` define the current workflow.

Keep generated outputs outside source directories. Check relative links and `git diff --check` for documentation changes. Preserve unrelated user changes and confirm the remote branch state before merging.
