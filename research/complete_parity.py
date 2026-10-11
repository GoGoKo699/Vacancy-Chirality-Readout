"""One-duration chirality readout from complete parity images.

Floating-point diagnostics; the exact enclosure is a separate certificate.
Flux is randomized equally on each independent input preparation.
"""
from __future__ import annotations

import itertools
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
SOURCE = Path(__file__).resolve().parents[1] / "source/resource-review-2026-10-10"
sys.path.insert(0, str(SOURCE))
import readout as rd

fr = rd.fr
PATTERNS = tuple(p for p in itertools.product((0, 1), repeat=4) if sum(p) % 2)


def reflect(pattern):
    return (pattern[0], pattern[3], pattern[2], pattern[1])


def effects(u=16.0, tau=np.pi, phi=-fr.PHI):
    """Eight parity effects on the original three-spin input space."""
    evolved = rd.evolved_inputs(float(u), float(tau), float(phi))
    occupation = fr.model(float(phi))[4]
    parity = [tuple(v % 2 for v in rd.site_occupancy(n)) for n in occupation]
    return np.array([evolved[[i for i, p in enumerate(parity) if p == pattern]].conj().T
                     @ evolved[[i for i, p in enumerate(parity) if p == pattern]]
                     for pattern in PATTERNS])


def calibration(u=16.0, tau=np.pi):
    """Mixed-input optimal antisymmetric score with calibrated unbiasedness."""
    e = effects(u, tau)
    reflected = np.array([fr.R @ x @ fr.R for x in e])
    odd, even = (e - reflected) / 2, (e + reflected) / 2
    a_matrix = np.array([[np.vdot(x, o).real / np.vdot(x, x).real
                          for o in odd] for x in (fr.CHI, fr.D)])
    q = np.trace(e, axis1=1, axis2=2).real / 8
    if np.min(q) <= 0:
        raise ArithmeticError("A zero-probability outcome requires support reduction.")
    gram = (a_matrix / q) @ a_matrix.T
    if np.linalg.eigvalsh(gram)[0] < 1e-12:
        raise ArithmeticError("Numerically unresolved odd-observable calibration.")
    scores = (a_matrix.T @ np.linalg.solve(gram, [1.0, 0.0])) / q
    estimate = np.einsum('p,pij->ij', scores, odd)
    second = np.einsum('p,pij->ij', scores**2, even)
    return dict(effects=e, odd=odd, even=even, coefficients=a_matrix,
                probabilities=q, gram=gram, scores=scores, estimate=estimate,
                second_moment=second,
                mixed_variance=float(q @ scores**2),
                worst_variance=float(np.linalg.eigvalsh(second)[-1]),
                residual=float(np.linalg.norm(estimate - fr.CHI, 2)))


def bit_flip_observables(cal, error):
    """Mean and second moment for fixed ideal scores; even reports score zero."""
    if not 0 <= error <= 1:
        raise ValueError("Use a bit-flip probability in [0,1].")
    score = dict(zip(PATTERNS, cal['scores']))
    mean = np.zeros((8, 8), complex)
    second = np.zeros((8, 8), complex)
    for sign in (1, -1):
        for pattern, e in zip(PATTERNS, cal['effects']):
            physical = pattern if sign == 1 else reflect(pattern)
            event = e if sign == 1 else fr.R @ e @ fr.R
            for reported in itertools.product((0, 1), repeat=4):
                distance = sum(x != y for x, y in zip(physical, reported))
                probability = error**distance * (1 - error)**(4 - distance)
                key = reported if sign == 1 else reflect(reported)
                value = sign * score.get(key, 0.0)
                mean += probability * value * event / 2
                second += probability * value**2 * event / 2
    return mean, second


def summary(cal):
    return {
        'patterns': [''.join(map(str, p)) for p in PATTERNS],
        'scores': cal['scores'].tolist(),
        'odd_coefficients': cal['coefficients'].tolist(),
        'mixed_probabilities': cal['probabilities'].tolist(),
        'mixed_variance': cal['mixed_variance'],
        'worst_variance': cal['worst_variance'],
        'max_abs_score': float(np.max(np.abs(cal['scores']))),
        'reconstruction_residual': cal['residual'],
        'scope': 'Floating-point diagnostics; see the separate exact certificate.',
    }


if __name__ == '__main__':
    import json
    print(json.dumps(summary(calibration()), indent=2))
