# Vacancy chirality readout: resource assessment

The existing four-site vacancy protocol reconstructs spin chirality from calibrated position-return frequencies. This continuation proves that total on-site parity is sufficient for the exact accepted occupancy event with known particle number at most three. It distinguishes final destructive readout from coherent retained-state filtering and compares the physical and statistical resources with established chirality measurements.

Read [the analytical resource account](READOUT_AND_RESOURCES.md), [the allocation review](REVIEW.md), and [primary-source comparisons](SOURCES.md). The original finite-repulsion and hard-core results are preserved under `prior/`.

To rerun five new groups and the eleven unchanged incoming groups, use a new output directory outside this package:

```bash
python verify.py --output-dir /absolute/new/path
```

The helper reuses the previous finite-U Hamiltonian explicitly. New loss and dephasing channels are diagnostic models, not measured device parameters. Numerical checks do not establish an experimental operating window, exact interval bounds, external proof review, or global priority. No repository was accessed or modified.
