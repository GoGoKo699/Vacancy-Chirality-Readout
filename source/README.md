# Bosonic models and foundational checks

These implementations support the complete-parity model and independently check the return instrument, finite-repulsion reconstruction, loss, and retained-state filtering.

| File | Role |
|---|---|
| [finite_repulsion.py](finite_repulsion.py) | Full 120-state bosonic model, spin operators, return amplitudes, and calibration |
| [readout.py](readout.py) | Parity effects, number-sector loss, and retained-state dephasing |
| [check_hard_core.py](check_hard_core.py) | Five groups: cycle reduction, entanglement, repeated filtering, disorder, and finite-repulsion controls |
| [check_finite_repulsion.py](check_finite_repulsion.py) | Five groups: independent model assembly, operation bounds, calibrated contrast, filtering, and perturbative diagnostics |
| [exact_short_time.py](exact_short_time.py) | One group: integer/rational short-time coefficients |
| [check_resources.py](check_resources.py) | Five groups: exact parity event, loss, dephasing, sampling, and detector errors |

The mathematical accounts are [Hard-core filtering](../research/HARD_CORE_FILTER.md), [Finite repulsion](../research/FINITE_REPULSION.md), and [Readout resources](../research/READOUT_RESOURCES.md). Start with [Theory](../THEORY.md) for the complete result and [Reproduce](../REPRODUCE.md) to run all 39 groups.

[MANIFEST.json](MANIFEST.json) records the six implementations. The verifier checks their membership and hashes and compares the four foundational output reports against [numerical reference reports](../verification/reference/).
