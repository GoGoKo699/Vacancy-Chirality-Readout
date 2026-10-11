#!/usr/bin/env python3
"""Three numerical cross-checks of the symmetry-based channel calibration theorem.

The diagnostic loss channel here acts AFTER coherent hopping. It is not the
continuous-time Lindblad loss model used in the preserved source checks.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
import complete_parity as cp

PATTERNS = tuple(itertools.product((0, 1), repeat=4))
FR = cp.fr


def product(*vectors):
    state = np.kron(np.kron(vectors[0], vectors[1]), vectors[2])
    return np.outer(state, state.conj())


UP = np.array([1., 0.], complex)
DOWN = np.array([0., 1.], complex)
PLUS_X = (UP + DOWN)/np.sqrt(2)
PLUS_Y = (UP + 1j*DOWN)/np.sqrt(2)
PROBES = (product(PLUS_X, PLUS_Y, UP), product(UP, UP, DOWN), product(UP, UP, UP))


def noisy_effects(u, tau, phi, survival, bit_error):
    """Exact finite-model effects for terminal uniform attenuation and bit errors.

Each particle survives with probability eta. Conditional on pre-loss site
occupation n, its surviving parity has odd probability [1-(1-2eta)^n]/2.
A classical independent bit channel follows. No survival postselection.
"""
    evolved = cp.rd.evolved_inputs(float(u), float(tau), float(phi))
    occ = FR.model(float(phi))[4]
    site_n = np.array([cp.rd.site_occupancy(n) for n in occ])
    odd_probability = (1-(1-2*survival)**site_n)/2
    reported_odd = bit_error + (1-2*bit_error)*odd_probability
    es = []
    for p in PATTERNS:
        weight = np.prod(np.where(np.array(p), reported_odd, 1-reported_odd), axis=1)
        es.append(evolved.conj().T @ (weight[:, None]*evolved))
    return np.array(es)


def reconstruct_moments(q):
    qc, qd, qs = q
    return np.array([qs, 6*qc-3*qd-2*qs, 3*qd-2*qc])


def run():
    symmetric = (FR.permutation((1, 0, 2)) + FR.permutation((0, 2, 1)) + FR.R)/3
    blocks = (symmetric, (FR.I8-symmetric)@(FR.I8+FR.R)/2,
              (FR.I8-symmetric)@(FR.I8-FR.R)/2)
    coords = np.array([[np.trace(rho @ a).real for a in (FR.CHI, FR.D)] for rho in PROBES])
    weights = np.array([[np.trace(rho @ p).real for p in blocks] for rho in PROBES])
    assert np.allclose(coords, [[1/8, 0], [0, 1], [0, 0]], atol=1e-14)
    assert np.allclose(weights, [[1/2, 1/4, 1/4], [1/3, 1/6, 1/2], [1, 0, 0]], atol=1e-14)
    assert np.linalg.norm(FR.CHI @ FR.D + FR.D @ FR.CHI, 2) < 1e-14
    assert np.linalg.norm(FR.CHI @ FR.CHI - (3/16)*(FR.I8-symmetric), 2) < 1e-14
    assert np.linalg.norm(FR.D @ FR.D - 3*(FR.I8-symmetric), 2) < 1e-14
    group1 = dict(probe_coordinates=coords.tolist(), block_weights=weights.tolist())

    diagnostics = []
    for u, tau, phi, eta, q in ((16., np.pi, -FR.PHI, 1., 0.),
                               (15.7, 1.01*np.pi, -1.02*FR.PHI, .985, .006)):
        plus = noisy_effects(u, tau, phi, eta, q)
        physical_minus = noisy_effects(u, tau, -phi, eta, q)
        minus = physical_minus[[PATTERNS.index(cp.reflect(p)) for p in PATTERNS]]
        reflection_error = max(np.linalg.norm(em-FR.R@ep@FR.R, 2) for ep, em in zip(plus, minus))
        assert reflection_error < 2e-12
        assert np.linalg.norm(plus.sum(axis=0)-FR.I8, 2) < 2e-12
        assert min(np.linalg.eigvalsh(ep)[0] for ep in plus) > -2e-12
        odd, even = (plus-minus)/2, (plus+minus)/2
        # Recover the complete response using only two product references.
        a = np.array([[factor*np.trace(rho@o).real for o in odd]
                      for factor, rho in zip((8., 1.), PROBES[:2])])
        decomp_error = max(np.linalg.norm(o-c*FR.CHI-d*FR.D, 2) for o, c, d in zip(odd, *a))
        assert decomp_error < 2e-12
        # Reconstruct every effect, not only the final scored observable,
        # from the five symmetry coordinates of the same three probes.
        effect_reconstruction_error = 0.
        for ep, op, fp, c, d in zip(plus, odd, even, *a):
            reference_probabilities = [np.trace(rho@fp).real for rho in PROBES]
            sector_values = reconstruct_moments(reference_probabilities)
            reconstructed = c*FR.CHI+d*FR.D+sum(t*p for t, p in zip(sector_values, blocks))
            effect_reconstruction_error = max(effect_reconstruction_error,
                                              np.linalg.norm(ep-reconstructed, 2))
        assert effect_reconstruction_error < 2e-12
        prob = np.trace(even, axis1=1, axis2=2).real/8
        support = prob > 1e-13
        gram = (a[:, support]/prob[support])@a[:, support].T
        assert np.linalg.eigvalsh(gram)[0] > 1e-8
        scores = np.zeros(16)
        scores[support] = (a[:, support].T@np.linalg.solve(gram, [1., 0.]))/prob[support]
        first = np.einsum('p,pij->ij', scores, odd)
        second = np.einsum('p,pij->ij', scores*scores, even)
        assert np.linalg.norm(first-FR.CHI, 2) < 2e-12
        ref_second = np.array([np.trace(rho@second).real for rho in PROBES])
        inferred = reconstruct_moments(ref_second)
        direct = np.array([np.trace(p@second).real/np.trace(p).real for p in blocks])
        assert np.max(np.abs(inferred-direct)) < 2e-12
        assert abs(max(inferred)-np.linalg.eigvalsh(second)[-1]) < 2e-12
        diagnostics.append(dict(u=u, tau=tau, phi=phi, terminal_survival=eta, bit_error=q,
                                calibration_residual=float(np.linalg.norm(first-FR.CHI, 2)),
                                symmetry_residual=float(decomp_error),
                                reflection_residual=float(reflection_error),
                                full_effect_reconstruction_residual=float(effect_reconstruction_error),
                                inferred_moment_eigenvalues=inferred.tolist(),
                                max_score=float(max(abs(scores))),
                                mixed_second_moment=float(np.trace(second).real/8)))
    group2 = dict(cases=diagnostics,
                  scope='Floating-point diagnostics, not interval-certified parameter windows or device data.')

    # Validate the analytic residual norm for a genuinely imperfect frozen score.
    dc, d = .012, -.009
    residual = dc*FR.CHI+d*FR.D
    exact_bias = np.sqrt(3*dc*dc/16+3*d*d)
    assert abs(np.linalg.norm(residual, 2)-exact_bias) < 1e-14
    # A symmetry-breaking input rotation escapes the two-response model;
    # the reference scheme must not be advertised as a test of its assumptions.
    local_z = np.kron(np.diag([1., -1.]), np.eye(4))
    broken = .9*FR.CHI + .1*local_z@FR.CHI@local_z
    c_bad = 8*np.trace(PROBES[0]@broken).real
    d_bad = np.trace(PROBES[1]@broken).real
    escaped = np.linalg.norm(broken-c_bad*FR.CHI-d_bad*FR.D, 2)
    assert escaped > .01
    group3 = dict(residual_norm=float(exact_bias), symmetry_breaking_unseen_residual=float(escaped),
                  scope='The confidence theorem requires independently justified spin and reflection covariance.')
    return dict(status='PASS', scientific_groups=3,
                groups=dict(product_reference_identities=group1, calibrated_noisy_instrument=group2,
                            bias_norm_and_scope=group3))


def main():
    if not __debug__:
        raise RuntimeError('Run with Python assertions enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    target = args.output.resolve()
    if target.is_relative_to(Path(__file__).resolve().parents[1]):
        parser.error('Use a new report outside the repository.')
    report = run()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS: three channel-calibration cross-check groups.')


if __name__ == '__main__':
    main()
