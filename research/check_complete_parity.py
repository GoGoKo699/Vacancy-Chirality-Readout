#!/usr/bin/env python3
"""Four floating-point cross-check groups, separate from the exact certificate."""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import complete_parity as cp
from parity_certificate import PASSIVE_SCORE_DECIMALS, PIVOTS


def run():
    fr = cp.fr
    cal = cp.calibration()
    e, odd, scores, moment = (cal[k] for k in ('effects', 'odd', 'scores', 'second_moment'))
    opposite = cp.effects(phi=fr.PHI)
    reflection_error = max(np.linalg.norm(opposite[cp.PATTERNS.index(cp.reflect(p))]
                                         - fr.R @ ep @ fr.R, 2)
                           for p, ep in zip(cp.PATTERNS, e))
    decomposition_error = max(np.linalg.norm(op - c*fr.CHI - d*fr.D, 2)
                              for op, c, d in zip(odd, *cal['coefficients']))
    assert np.linalg.norm(e.sum(axis=0) - fr.I8, 2) < 1e-12
    assert reflection_error < 1e-12 and decomposition_error < 1e-12
    assert min(np.linalg.eigvalsh(ep)[0] for ep in e) > -1e-12
    assert cal['residual'] < 1e-12
    group1 = dict(reflection_error=reflection_error, decomposition_error=decomposition_error,
                  reconstruction_residual=cal['residual'])

    symmetric = (fr.permutation((1, 0, 2)) + fr.permutation((0, 2, 1)) + fr.R)/3
    projectors = [symmetric, (fr.I8-symmetric)@(fr.I8+fr.R)/2,
                  (fr.I8-symmetric)@(fr.I8-fr.R)/2]
    coefficients = []
    for p in projectors:
        rho = p / np.trace(p).real
        value = np.trace(moment @ rho).real
        assert abs(np.trace(fr.CHI @ rho)) < 1e-14
        assert np.linalg.norm(moment @ p - value*p, 2) < 1e-12
        coefficients.append(float(value))
    assert abs(max(coefficients)-cal['worst_variance']) < 1e-12
    assert cal['worst_variance'] < .39
    assert abs(cal['mixed_variance'] - .284895494202985) < 1e-12
    assert abs(cal['probabilities'] @ scores) < 1e-12
    group2 = dict(block_coefficients=coefficients, **cp.summary(cal))

    # Reconstruct the new reference score independently from the complete
    # 120-state floating model; import only its exact decimal definition.
    primary_scores = np.array(PASSIVE_SCORE_DECIMALS, dtype=float)
    primary_residual = np.array([1., 0.])-cal['coefficients']@primary_scores
    primary_correction = np.linalg.solve(cal['coefficients'][:, PIVOTS], primary_residual)
    primary_scores[list(PIVOTS)] += primary_correction
    primary_mean = np.einsum('p,pij->ij', primary_scores, odd)
    primary_moment = np.einsum('p,pij->ij', primary_scores**2, cal['even'])
    primary_error = float(np.linalg.norm(primary_mean-fr.CHI, 2))
    assert primary_error < 1e-12
    primary_blocks = []
    for projector in projectors:
        value = float(np.trace(primary_moment@projector).real/np.trace(projector).real)
        assert np.linalg.norm(primary_moment@projector-value*projector, 2) < 1e-12
        primary_blocks.append(value)
    primary_worst = float(np.linalg.eigvalsh(primary_moment)[-1])
    assert abs(max(primary_blocks)-primary_worst) < 1e-12
    assert primary_worst < .328697776905 < 13/32
    group2['primary_minimax_score'] = dict(
        score_family='primary minimax-oriented feasible score',
        patterns=[''.join(map(str, p)) for p in cp.PATTERNS],
        scores=primary_scores.tolist(), two_pivot_corrections=primary_correction.tolist(),
        reconstruction_residual=primary_error, block_coefficients=primary_blocks,
        worst_variance=primary_worst,
        mixed_variance=float(cal['probabilities']@primary_scores**2),
        checked_upper_bound=.328697776905, passive_reference=13/32,
        scope='Independent 120-state numerical cross-check; exact enclosure is in parity_certificate.py.')

    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.])
    axes = (x, y, z)
    def tensor3(parts):
        return np.kron(np.kron(parts[0], parts[1]), parts[2])
    first = np.zeros((8, 8), complex)
    second = np.zeros((8, 8), complex)
    fisher = 0.
    for perm in itertools.permutations(range(3)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        for outcome in itertools.product((-1, 1), repeat=3):
            effect = tensor3([(np.eye(2)+s*axes[k])/2 for s, k in zip(outcome, perm)])/6
            score = .75*sign*np.prod(outcome)
            first += score*effect
            second += score**2*effect
            probability = np.trace(effect).real/8
            derivative = (4/3)*np.trace(fr.CHI @ effect).real
            fisher += derivative**2/probability
            normalized_chi = np.trace(fr.CHI @ effect).real/np.trace(effect).real
            assert abs(abs(normalized_chi)-.125) < 1e-14
    assert np.linalg.norm(first-fr.CHI, 2) < 1e-14
    assert np.linalg.norm(second-(9/16)*fr.I8, 2) < 1e-14
    assert abs(fisher-16/9) < 1e-13
    group3 = dict(six_Pauli_mixed_variance=9/16, fisher_at_origin=float(fisher),
                  benchmark='Analytical separable lower bound is proved in SEPARABLE_BENCHMARK.md.')

    q, lam = .003, .01
    all_patterns = tuple(itertools.product((0, 1), repeat=4))
    kernel = np.array([[q**sum(a != b for a, b in zip(p, t))
                        * (1-q)**(4-sum(a != b for a, b in zip(p, t)))
                        for p in cp.PATTERNS] for t in all_patterns])
    f, g = kernel@scores, kernel@(scores**2)
    odd_indices = [all_patterns.index(p) for p in cp.PATTERNS]
    single = [i for i, p in enumerate(all_patterns) if sum(p) == 1]
    two = [i for i, p in enumerate(all_patterns) if sum(p) in (0, 2)]
    alpha = (1-q)**4
    b = 6*q*q*(1-q)**2 + q**4
    c = 4*q*(1-q)**3 + 4*q**3*(1-q)
    eta, zz = math.exp(-lam), math.expm1(lam)
    number = [math.comb(3,k)*eta**k*(1-eta)**(3-k) for k in range(4)]
    assert abs(sum(number)-1) < 1e-14
    s = number[3]
    amax = np.max(np.abs(scores))
    singleton_scores = scores[[cp.PATTERNS.index(all_patterns[i]) for i in single]]
    c1 = np.ptp(singleton_scores)/2
    bias_bound = 3*c1*zz**2 + amax*b/alpha*(1+3*zz**2) + 3*amax*c/alpha*zz
    moment_bound = (alpha*cal['worst_variance']+b*amax**2
                    + 3*zz**2*max(g[single])+3*zz*max(g[two])+zz**3*g[0])/(s*alpha**2)
    # Deliberately adversarial lower-number pattern distributions. They need
    # not arise from this Hamiltonian: the analytical bounds cover all of them.
    max_bias = 0.
    for rho in (fr.I8/8, fr.PPLUS/2, fr.PMINUS/2, projectors[-1]/2):
        plus = np.array([np.trace(ep@rho).real for ep in e])
        minus = np.array([np.trace(fr.R@ep@fr.R@rho).real for ep in e])
        mean = (s*(plus-minus)@f[odd_indices]/2
                + number[1]*np.ptp(f[single])/2
                + number[2]*np.ptp(f[two])/2)/(s*alpha)
        bias = abs(mean-np.trace(fr.CHI@rho).real)
        assert bias <= bias_bound+1e-12
        max_bias = max(max_bias, float(bias))
    clean_mean, clean_second = cp.bit_flip_observables(cal, q)
    assert np.linalg.norm(clean_mean/alpha-fr.CHI, 2) <= amax*b/alpha+1e-12
    assert np.linalg.eigvalsh(clean_second/alpha**2)[-1] <= (
        cal['worst_variance']/alpha+b*amax**2/alpha**2)+1e-12
    # Independent known-detector correction at the certified window endpoint.
    # The direct probability kernel and full effects avoid using the analytic
    # complement-matrix inverse or its operator-norm bound.
    detector_q = .03
    odd_kernel = np.array([[detector_q**sum(x != y for x, y in zip(p, t))
                            *(1-detector_q)**(4-sum(x != y for x, y in zip(p, t)))
                            for p in cp.PATTERNS] for t in cp.PATTERNS])
    detector_scores = np.linalg.solve(odd_kernel, primary_scores)
    detector_mean, detector_moment = cp.bit_flip_observables(
        dict(effects=e, scores=detector_scores), detector_q)
    detector_error = float(np.linalg.norm(detector_mean-fr.CHI, 2))
    detector_worst = float(np.linalg.eigvalsh(detector_moment)[-1])
    assert detector_error < 1e-12
    assert detector_worst < .382273 < 13/32

    group4 = dict(lambda_loss=lam, bit_flip=q, number_probabilities=number,
                  simple_bias_bound=float(bias_bound), second_moment_bound=float(moment_bound),
                  adversarial_examples_max_bias=max_bias,
                  scope='Probability-budget and detector-channel cross-checks, not a device fit.')
    group4['primary_detector_only'] = dict(
        bit_flip=detector_q, atom_loss=False, corrected_scores=detector_scores.tolist(),
        reconstruction_residual=detector_error, worst_variance=detector_worst,
        checked_upper_bound=.382273, passive_reference=13/32,
        scope='Single endpoint numerical cross-check; the uniform window has a separate exact certificate.')
    return dict(status='PASS', scientific_groups=4,
                groups=dict(effects_and_reconstruction=group1, symmetry_variance=group2,
                            separable_attainment=group3, noise_budget=group4))


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
    print('PASS: four complete-parity cross-check groups.')


if __name__ == '__main__':
    main()
