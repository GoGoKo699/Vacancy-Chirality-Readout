# Vacancy chirality: first research pilot

Read [PILOT.md](PILOT.md) for the exact hard-core instrument, three-spin entanglement proof, repeat-filter yield, static-disorder statement, finite-U countercheck, and source attribution.

The retained question is whether this native spin-independent vacancy implementation supplies a useful, controllably accurate chirality instrument beyond its ideal limit. This package is a bounded pilot, not a repository handoff or manuscript.

Run the checks with Python, NumPy, and SciPy:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python check_pilot.py --output /absolute/new/report.json
```

The output must not exist. The saved formal reports and logs are in `evidence/`; exact environment and comparison details are in `RECORD.json`. The largest full Hilbert space is 120. Numerical checks are not an external theorem review, priority certificate, or proof that all required experimental controls have been realized together.
