#!/usr/bin/env python3
"""Four cross-check groups for the analytical passive-parity bounds.

Finite numerical examples check identities; they do not certify a global
optimum. Rational arithmetic checks the Fisher and finite-risk constants.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

import complete_parity as cp

fr = cp.fr
PERMUTATIONS = tuple(itertools.permutations(range(3)))
PERMUTATION_MATRICES = np.array([fr.permutation(p) for p in PERMUTATIONS])
PM = fr.I8 - fr.PSYM
H = math.sqrt(3) / 4
B = 4 / (3 * math.sqrt(3))


def unitary(size, seed):
    """A fixed diagnostic example, not a unitary optimization or search."""
    rng = np.random.default_rng(seed)
    q, r = np.linalg.qr(rng.normal(size=(size, size))
                        + 1j * rng.normal(size=(size, size)))
    return q @ np.diag(np.diag(r) / np.abs(np.diag(r)))


def patterns(modes):
    return tuple(p for p in itertools.product((0, 1), repeat=modes)
                 if sum(p) % 2)


def compressed_parity(u, subset):
    """Permutation expansion of a pulled-back product of output parities."""
    modes = len(u)
    occupied = tuple(range(modes - 3, modes))
    signs = np.diag([-1 if j in subset else 1 for j in range(modes)])
    w = u.conj().T @ signs @ u
    k = w[np.ix_(occupied, occupied)]
    return sum(np.prod([k[i, p[i]] for i in range(3)]) * matrix
               for p, matrix in zip(PERMUTATIONS, PERMUTATION_MATRICES))


def passive_effects(u):
    modes = len(u)
    subsets = tuple(tuple(i for i in range(modes) if mask & (1 << i))
                    for mask in range(1 << modes))
    observables = [compressed_parity(u, s) for s in subsets]
    return np.array([sum((-1) ** sum(p[i] for i in s) * obs
                         for s, obs in zip(subsets, observables)) / 2**modes
                     for p in patterns(modes)])


def fock_effects(u):
    """Independent creation-operator expansion, with occupation factorials."""
    modes = len(u)
    occupied = tuple(range(modes - 3, modes))
    occupations = []
    for multiset in itertools.combinations_with_replacement(range(2*modes), 3):
        n = [0] * (2*modes)
        for j in multiset:
            n[j] += 1
        occupations.append(tuple(n))
    lookup = {n: i for i, n in enumerate(occupations)}
    evolved = np.zeros((len(occupations), 8), complex)
    for col, spins in enumerate(fr.SPINS):
        for outputs in itertools.product(range(modes), repeat=3):
            n = [0] * (2*modes)
            for out, spin in zip(outputs, spins):
                n[2*out + spin] += 1
            amplitude = np.prod([u[out, src] for out, src in zip(outputs, occupied)])
            amplitude *= math.sqrt(math.prod(math.factorial(k) for k in n))
            evolved[lookup[tuple(n)], col] += amplitude
    assert np.linalg.norm(evolved.conj().T @ evolved-fr.I8, 2) < 2e-13
    parity = [tuple((n[2*j]+n[2*j+1]) % 2 for j in range(modes))
              for n in occupations]
    return np.array([evolved[[i for i, p in enumerate(parity) if p == target]].conj().T
                     @ evolved[[i for i, p in enumerate(parity) if p == target]]
                     for target in patterns(modes)])


def coefficients(effects):
    return np.einsum('ij,pji->p', fr.CHI, effects).real * (4/3)


def mixed_sector_probabilities(effects):
    return np.einsum('ij,pji->p', PM, effects).real / 4


def information_bound(t):
    if not 0 <= t < 1:
        raise ValueError('The normalized prior radius must lie in [0, 1).')
    return 16*(6+t)/(3*(1-t*t)*(13+t))


def explicit_unitary(x):
    v, r = math.sqrt(x), math.sqrt((1-x)/2)
    rt2 = math.sqrt(2)
    return np.array([[v/rt2, r/rt2, 1/rt2, r/rt2],
                     [r, (-v+1j)/2, 0, (-v-1j)/2],
                     [r, (-v-1j)/2, 0, (-v+1j)/2],
                     [v/rt2, r/rt2, -1/rt2, r/rt2]])


def pi_interval():
    """Machin's identity and alternating-series remainders, all rational."""
    def atan(q):
        count = 32
        partial = sum((F((-1)**k, (2*k+1)*q**(2*k+1))
                       for k in range(count)), F(0))
        # The next term is positive because count is even.
        return partial, partial + F(1, (2*count+1)*q**(2*count+1))
    a5, a239 = atan(5), atan(239)
    return 16*a5[0]-4*a239[1], 16*a5[1]-4*a239[0]


def ceiling(x):
    return -((-x.numerator)//x.denominator)


def run():
    fock_errors, three_mode_slopes = [], []
    for modes, seeds in ((3, (13, 29)), (4, (7, 19))):
        for seed in seeds:
            u = unitary(modes, seed)
            effects = passive_effects(u)
            direct = fock_effects(u)
            fock_errors.append(float(max(np.linalg.norm(a-b, 2)
                                         for a, b in zip(effects, direct))))
            assert fock_errors[-1] < 2e-13
            assert np.linalg.norm(effects.sum(axis=0)-fr.I8, 2) < 2e-13
            assert min(np.linalg.eigvalsh(e)[0] for e in effects) > -2e-13
            if modes == 3:
                three_mode_slopes.append(float(np.max(np.abs(coefficients(effects)))))
                assert three_mode_slopes[-1] < 2e-13
    f3 = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/math.sqrt(3)
    ef = passive_effects(f3)
    for p, e in zip(patterns(3), ef):
        target = fr.PSYM/3 if sum(p) == 3 else fr.I8/3-fr.PSYM/9
        assert np.linalg.norm(e-target, 2) < 2e-13
    chiral_difference = fr.PPLUS/2-fr.PMINUS/2
    assert max(abs(np.trace(e @ chiral_difference)) for e in ef) < 2e-13
    group1 = dict(max_permutation_vs_fock_error=max(fock_errors),
                  max_three_mode_chirality_coefficient=max(three_mode_slopes),
                  three_mode_minimax_mse_exact='3/16',
                  scope='Direct amplitudes independently check the compressed-permutation formula.')

    rotation_errors, fisher_ratios, j_errors = [], [], []
    x, y, z = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])
    local = [expm(-1j*.37*x), expm(1j*.61*y), expm(-1j*.29*z)]
    rotation = np.kron(np.kron(local[0], local[1]), local[2])
    g = (4/3)*np.trace(fr.CHI @ rotation.conj().T @ fr.CHI @ rotation).real
    assert abs(g) <= 1+1e-14
    for seed in (2, 11, 31):
        u = unitary(4, seed)
        e = passive_effects(u)
        c = coefficients(e)
        b = np.array([(4/3)*np.trace(fr.CHI @ compressed_parity(u, (0, j))).real
                      for j in (1, 2, 3)])
        assert np.max(abs(b)) <= B+1e-13
        expected = np.array([sum((-1)**(p[0]+p[j])*b[j-1] for j in (1, 2, 3))/8
                             for p in patterns(4)])
        assert np.max(abs(c-expected)) < 2e-13
        assert H*sum(abs(c)) <= .5+1e-13
        er = np.array([rotation.conj().T @ ep @ rotation for ep in e])
        cr = coefficients(er)
        rotation_errors.append(float(np.max(abs(cr-g*c))))
        assert rotation_errors[-1] < 2e-13
        for effects in (e, er):
            q, slope = mixed_sector_probabilities(effects), coefficients(effects)
            assert np.min(q-H*abs(slope)) > -2e-13
            for normalized_theta in (-.7, 0., .4):
                probability = q+H*normalized_theta*slope
                keep = probability > 1e-14
                fi = float(np.sum(slope[keep]**2/probability[keep]))
                bound = information_bound(abs(normalized_theta))
                assert fi <= bound+2e-12
                fisher_ratios.append(fi/bound)
        # Check both compression eigenvalues and the Vandermonde/J identity.
        w = u.conj().T @ np.diag([-1., -1., 1., 1.]) @ u
        k = w[1:, 1:]
        ev, vec = np.linalg.eigh(k)
        assert abs(ev[0]+1) < 2e-13 and abs(ev[-1]-1) < 2e-13
        j = (vec[0, 0]*vec[1, 1]*vec[0, 1].conjugate()*vec[1, 0].conjugate()).imag
        vandermonde = math.prod(abs(ev[i]-ev[jj]) for i in range(3) for jj in range(i+1, 3))
        j_errors.append(abs(abs((k[0, 1]*k[1, 2]*k[2, 0]).imag)-abs(j)*vandermonde))
        assert j_errors[-1] < 2e-13
        assert abs(j) <= 1/(6*math.sqrt(3))+2e-13
    # At a cube corner use a_z=h*c_z, so all calculations are rational.
    scaled_c = [F(-1, 8)]*2 + [F(1, 24)]*6
    assert sum(scaled_c) == 0 and sum(map(abs, scaled_c)) == F(1, 2)
    vertex_values = {}
    for t in (F(0), F(1, 10), F(1, 5), F(7, 10)):
        values = []
        for sign in (-1, 1):
            for j in range(8):
                probabilities = [abs(v)+t*sign*v for v in scaled_c]
                probabilities[j] += F(1, 2)
                values.append(sum((F(16, 3)*v*v/q for v, q in zip(scaled_c, probabilities)), F(0)))
        assert max(values) == information_bound(t)
        vertex_values[str(t)] = str(max(values))
    assert information_bound(F(0)) == F(32, 13)
    group2 = dict(max_local_rotation_coefficient_error=max(rotation_errors),
                  max_vandermonde_identity_error=float(max(j_errors)),
                  largest_example_fisher_to_bound_ratio=max(fisher_ratios),
                  rational_vertex_fisher_bounds=vertex_values,
                  unbiased_passive_variance_lower_bound='13/32',
                  scope='Numerical examples and exact cube-vertex arithmetic supplement the convexity proof.')

    def check_construction(value):
        u = explicit_unitary(value)
        assert np.linalg.norm(u.conj().T @ u-np.eye(4), 2) < 2e-13
        e = passive_effects(u)
        q, c = mixed_sector_probabilities(e), coefficients(e)
        formula_f = 32*value*(1-value*value)/((value+3)*(3*value+1))
        assert abs(np.sum(c*c/q)-formula_f) < 2e-13
        a = c/(q*formula_f)
        first, second = np.zeros((8, 8), complex), np.zeros((8, 8), complex)
        for p, matrix in zip(PERMUTATIONS, PERMUTATION_MATRICES):
            sign = (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
            rotated = np.array([matrix.conj().T @ ep @ matrix for ep in e])
            first += sign*np.einsum('p,pij->ij', a, rotated)/6
            second += np.einsum('p,pij->ij', a*a, rotated)/6
        ts = (3*value+1)**2/(16*(1-value)*(value+1)**2)
        tm = 1/formula_f
        assert np.linalg.norm(first-fr.CHI, 2) < 2e-13
        assert np.linalg.norm(second-ts*fr.PSYM-tm*PM, 2) < 2e-13
        assert ts <= tm+2e-13
        assert abs(np.trace(second @ PM/4).real-tm) < 2e-13
        return tm, ts
    tm, ts = check_construction(.4)
    assert abs(tm-float(F(935, 1344))) < 2e-13
    optimum = brentq(lambda a: 3*a**4+20*a**3+12*a*a-3, .3, .5, xtol=1e-15)
    improved, _ = check_construction(optimum)
    assert improved < tm and abs(improved-.695135389809721) < 2e-13
    group3 = dict(algebraic_choice_x='2/5', coefficient_exact='935/1344',
                  symmetric_second_moment=ts, mixed_second_moment=tm,
                  family_optimum_x=optimum, family_optimum_coefficient=improved,
                  scope='An explicit unbiased construction; the one-parameter optimum is not a global passive optimum.')

    pi_lo, pi_hi = pi_interval()
    assert pi_lo < F(22, 7) and pi_hi*pi_hi < 10
    eps, radius = F(1, 1000), F(1, 10)
    fi = information_bound(radius)
    counts = [(1/(eps*eps)-16*pi*pi/(3*radius*radius))/fi for pi in (pi_hi, pi_lo)]
    assert ceiling(counts[0]) == ceiling(counts[1]) == 396539
    comparison_fi = information_bound(F(1, 5))
    native = F(367, 1000)
    threshold = native*F(4000, 3)/(1-native*comparison_fi)
    assert comparison_fi == F(775, 297)
    assert threshold == F(5813280, 503) and threshold < 12000
    assert native/12000 < 1/(12000*comparison_fi+F(4000, 3))
    ring_errors, ring_model_errors = [], []
    for tau, phi in ((.39, .71), (.83, -1.2), (math.pi, -2*math.pi/3),
                     (math.pi, 2*math.pi/3)):
        hopping = np.zeros((4, 4), complex)
        for j in range(4):
            k = (j+1) % 4
            hopping[j, k] = -np.exp(1j*phi/4)
            hopping[k, j] = -np.exp(-1j*phi/4)
        u = expm(-1j*tau*hopping)
        xx, yy = 2*tau*math.cos(phi/4), -2*tau*math.sin(phi/4)
        first_row = np.array([(math.cos(xx)+math.cos(yy))/2,
                              (1j*math.sin(xx)+math.sin(yy))/2,
                              (math.cos(xx)-math.cos(yy))/2,
                              (1j*math.sin(xx)-math.sin(yy))/2])
        assert np.max(abs(u[0]-first_row)) < 2e-13
        e = passive_effects(u)
        kappa = math.sin(2*xx)*math.sin(2*yy)*(math.cos(2*xx)-math.cos(2*yy))/8
        expected = np.array([(-1)**(p[0]+p[2])*kappa for p in patterns(4)])
        ring_errors.append(float(np.max(abs(coefficients(e)-expected))))
        assert ring_errors[-1] < 2e-13
        assert abs(kappa) <= 1/(6*math.sqrt(3))+2e-13
        full = cp.effects(u=0., tau=tau, phi=phi)
        ring_model_errors.append(float(max(np.linalg.norm(a-b, 2) for a, b in zip(e, full))))
        assert ring_model_errors[-1] < 3e-12
        if tau == math.pi:
            assert abs(kappa) < 2e-14
    # In scaled units h*kappa=1/24, 1/h^2=16/3.
    kk = F(1, 24)
    ring_fi = F(16, 3)*(7*kk+kk*kk/(1-7*kk))
    assert ring_fi == F(80, 51)
    group4 = dict(epsilon='0.001', prior_radius_over_h='1/10',
                  necessary_passive_count_interval=[float(c) for c in counts],
                  necessary_passive_integer_copies=396539,
                  noisy_native_coarse_coefficient='0.367',
                  all_biased_strict_comparison_sufficient_n=12000,
                  exact_sufficient_threshold=str(threshold),
                  uniform_ring_fisher_bound='80/51',
                  uniform_ring_variance_lower_bound='51/80',
                  max_ring_coefficient_error=max(ring_errors),
                  max_ring_vs_full_bosonic_model_error=max(ring_model_errors),
                  nominal_noninteracting_ring_is_blind=True,
                  scope='Rational finite-risk constants; ring formulas checked against the independent full bosonic implementation.')
    return dict(status='PASS', scientific_groups=4,
                groups=dict(compression_and_three_modes=group1,
                            four_mode_information_and_prerotations=group2,
                            explicit_four_mode_construction=group3,
                            finite_risk_and_uniform_ring=group4))


def main():
    if not __debug__:
        raise RuntimeError('Run with Python assertions enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        help='New report outside the repository; default is a new temporary directory.')
    args = parser.parse_args()
    target = (args.output if args.output is not None else
              Path(tempfile.mkdtemp(prefix='vacancy-passive-parity-'))/'passive-parity.json').resolve()
    if target.is_relative_to(Path(__file__).resolve().parents[1]):
        parser.error('Use a new report outside the repository.')
    if target.exists():
        parser.error('The report already exists; choose a new output path.')
    report = run()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(f'PASS: four passive-parity cross-check groups. Report: {target}')


if __name__ == '__main__':
    main()
