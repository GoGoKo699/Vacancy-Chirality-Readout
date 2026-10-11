#!/usr/bin/env python3
"""Integer/rational certificate for the calibrated complete-parity estimator.

This checker uses only Python's standard library.  All certified quantities are
rational intervals.  No floating eigensolver, matrix exponential, transcendental
library call, or floating tolerance enters the certificate.  See CERTIFICATE.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools
import json
from math import factorial, isqrt
from pathlib import Path

BITS = 1024
SCALE = 1 << BITS
ORDER = 1280
NORM_BOUND = 256
PATTERNS = tuple(p for p in itertools.product((0, 1), repeat=4) if sum(p) % 2)
SCORE_DECIMALS = (
    '0.736714019013', '0.252900613464', '2.336375499791',
    '0.730697335072', '-0.660422527924', '-0.288874449056',
    '-0.305946882115', '-0.540240033331',
)
# The exact correction uses these two columns of the two-row calibration map.
PIVOTS = (PATTERNS.index((0, 1, 1, 1)), PATTERNS.index((0, 1, 0, 0)))
# Rational approximation to a second right-inverse column for the scaled map
# G=(c;4d).  Multiplication by (G B0)^(-1) below makes the inverse exact.
RIGHT_INVERSE_SECOND = (
    '-12.579980610295', '-7.146759029809', '-48.311368416521',
    '1.385838754865', '5.088534362922', '0.402902729639',
    '-0.010440970833', '2.443795589669',
)


def ceil_fraction(x: F) -> int:
    return -((-x.numerator) // x.denominator)


class Interval:
    """Closed rational interval; every arithmetic operation is exact."""
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = as_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_interval(other))

    def __rsub__(self, other):
        return as_interval(other) - self

    def __mul__(self, other):
        other = as_interval(other)
        endpoints = [x*y for x in (self.lo, self.hi)
                     for y in (other.lo, other.hi)]
        return Interval(min(endpoints), max(endpoints))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_interval(other)
        assert not other.lo <= 0 <= other.hi, 'Unresolved denominator.'
        return self * Interval(1/other.hi, 1/other.lo)

    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def centered_widen(self, error):
        return Interval(self.lo - error, self.hi + error)


def as_interval(x):
    return x if isinstance(x, Interval) else Interval(x)


def decimal_bound(x: F, places=15, upper=False):
    """Outward decimal display, still without binary floating point."""
    scale = 10**places
    n = ceil_fraction(x * scale) if upper else (x * scale).__floor__()
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n // scale}.{n % scale:0{places}d}'


def display(x: Interval, places=15):
    return [decimal_bound(x.lo, places), decimal_bound(x.hi, places, True)]


def fixed_interval(x: Interval):
    return (x.lo*SCALE).__floor__(), ceil_fraction(x.hi*SCALE)


def from_fixed(bounds):
    return Interval(F(bounds[0], SCALE), F(bounds[1], SCALE))


def arctan_inverse(q):
    # Each retained division has error < 1/SCALE; the first omitted term is
    # < 1/SCALE, and the alternating-series theorem bounds the remaining tail.
    power = q
    k = 0
    total = 0
    while True:
        quotient = SCALE // ((2*k + 1)*power)
        if quotient == 0:
            break
        total += quotient if k % 2 == 0 else -quotient
        k += 1
        power *= q*q
    return Interval(F(total-k-1, SCALE), F(total+k+1, SCALE)), k


def constants():
    a5, n5 = arctan_inverse(5)
    a239, n239 = arctan_inverse(239)
    # Machin's identity is exact, with the principal arctangent branches.
    pi = 16*a5 - 4*a239
    roots = {}
    for n in (1, 2, 3, 4):
        r = isqrt(n*SCALE*SCALE)
        roots[n] = Interval(F(r, SCALE), F(r+1, SCALE))
    assert 3 < pi.lo < pi.hi < F(22, 7)
    return pi, roots, {'arctan_terms': [n5, n239]}


def finite_model(pi, roots):
    occ = []
    for modes in itertools.combinations_with_replacement(range(8), 3):
        n = [0]*8
        for mode in modes:
            n[mode] += 1
        if sum(n[1::2]) == 1:
            occ.append(tuple(n))
    assert len(occ) == 40
    lookup = {n: i for i, n in enumerate(occ)}
    # A = -i*pi*H, where H has u=16 and phi=-2*pi/3.  Every pair
    # contains exact rational enclosures for real and imaginary parts.
    matrix = [[(Interval(0), Interval(0)) for _ in occ] for _ in occ]
    for col, n in enumerate(occ):
        ns = [n[2*j] + n[2*j+1] for j in range(4)]
        w = sum(v*(v-1)//2 for v in ns)
        matrix[col][col] = (Interval(0), -16*w*pi)
        for j in range(4):
            k = (j+1) % 4
            for spin in range(2):
                for src, dst, sign in ((2*k+spin, 2*j+spin, 1),
                                       (2*j+spin, 2*k+spin, -1)):
                    if n[src] == 0:
                        continue
                    new = list(n)
                    new[src] -= 1
                    new[dst] += 1
                    row = lookup[tuple(new)]
                    root = roots[n[src]*(n[dst]+1)]
                    real = sign*pi*root/2
                    imag = pi*roots[3]*root/2
                    oldr, oldi = matrix[row][col]
                    matrix[row][col] = oldr+real, oldi+imag
    sparse = []
    max_model_error_units = 0
    for row in matrix:
        terms = []
        row_norm_units = 0
        row_error_units = 0
        for col, (real, imag) in enumerate(row):
            if real.lo == real.hi == imag.lo == imag.hi == 0:
                continue
            rl, rh = fixed_interval(real)
            il, ih = fixed_interval(imag)
            re, im = (rl+rh)//2, (il+ih)//2
            terms.append((col, re, im))
            row_norm_units += abs(re)+abs(im)
            row_error_units += max(re-rl, rh-re) + max(im-il, ih-im)
        assert row_norm_units + row_error_units <= NORM_BOUND*SCALE
        max_model_error_units = max(max_model_error_units, row_error_units)
        sparse.append(terms)
    inputs = []
    for spins in ((0, 0, 1), (0, 1, 0), (1, 0, 0)):
        n = [0]*8
        for site, spin in enumerate(spins, 1):
            n[2*site+spin] = 1
        inputs.append(lookup[tuple(n)])
    return occ, sparse, inputs, max_model_error_units


def taylor_action(matrix, inputs, model_error_units):
    term = [[(0, 0) for _ in inputs] for _ in matrix]
    for j, i in enumerate(inputs):
        term[i][j] = (SCALE, 0)
    total = [row[:] for row in term]
    for k in range(1, ORDER+1):
        newterm = []
        denominator = SCALE*k
        for row, entries in enumerate(matrix):
            newrow = []
            for j in range(len(inputs)):
                real = imag = 0
                for col, ar, ai in entries:
                    xr, xi = term[col][j]
                    real += ar*xr - ai*xi
                    imag += ar*xi + ai*xr
                value = (real//denominator, imag//denominator)
                newrow.append(value)
                tr, ti = total[row][j]
                total[row][j] = (tr+value[0], ti+value[1])
            newterm.append(newrow)
        term = newterm
    # e^B < 3^B.  These rational bounds include the enclosure of A, all
    # complex arithmetic rounding, and the discarded exponential tail.
    exp_bound = 3**NORM_BOUND
    model_and_rounding = F(model_error_units + 2*ORDER, SCALE)
    tail = F(NORM_BOUND**(ORDER+1), factorial(ORDER+1))
    error = exp_bound*(model_and_rounding + tail)
    assert error < F(1, 10**170)
    return total, error


def effects(occ, action, action_error):
    result = []
    for p in PATTERNS:
        rows = [i for i, n in enumerate(occ)
                if tuple((n[2*j]+n[2*j+1]) % 2 for j in range(4)) == p]
        effect = []
        for a in range(3):
            erow = []
            for b in range(3):
                re = im = 0
                for i in rows:
                    ar, ai = action[i][a]
                    br, bi = action[i][b]
                    re += ar*br + ai*bi
                    im += ar*bi - ai*br
                erow.append((re//SCALE, im//SCALE))
            effect.append(erow)
        result.append(effect)
    # True propagator entries have modulus <=1.  Product perturbations are
    # bounded by 2*epsilon+epsilon^2; at most 40 terms occur in an effect.
    error = 40*(2*action_error + action_error*action_error) + F(2, SCALE)
    return result, error


def calibration_intervals(es, error):
    cycle = ((0, 0, 1), (1, 0, 0), (0, 1, 0))
    dmatrix = ((1, -1, 0), (-1, 0, 1), (0, 1, -1))
    cs, ds = [], []
    for e in es:
        c = F(2, 3*SCALE)*sum((cycle[i][j]-cycle[j][i])*e[j][i][1]
                              for i in range(3) for j in range(3))
        d = F(1, 6*SCALE)*sum(dmatrix[i][j]*e[j][i][0]
                              for i in range(3) for j in range(3))
        cs.append(Interval(c).centered_widen(4*error))
        ds.append(Interval(d).centered_widen(error))
    return cs, ds


def correction(cs, ds, b, target=(1, 0)):
    p, q = PIVOTS
    determinant = cs[p]*ds[q] - cs[q]*ds[p]
    assert not determinant.lo <= 0 <= determinant.hi
    residual_c = target[0] - sum((bp*cp for bp, cp in zip(b, cs)), Interval(0))
    residual_d = target[1] - sum((bp*dp for bp, dp in zip(b, ds)), Interval(0))
    alpha = (residual_c*ds[q] - cs[q]*residual_d)/determinant
    beta = (cs[p]*residual_d - residual_c*ds[p])/determinant
    square_change = max(2*abs(b[index])*z.abs_upper() + z.abs_upper()**2
                        for index, z in ((p, alpha), (q, beta)))
    return determinant, alpha, beta, square_change


def second_moment(es, error, b):
    def entry(i, j):
        # Reflection R exchanges indices 0 and 2 in this input sector.
        return sum(bp*bp*F(e[i][j][0]+e[2-i][2-j][0], 2*SCALE)
                   for bp, e in zip(b, es))
    a, bb, c = entry(0, 0), entry(0, 1), entry(0, 2)
    err = sum(bp*bp for bp in b)*error
    # Exact SU(2) invariance and reflection imply three fixed eigenvectors.
    # Their eigenvalues are linear functions of these three entries.
    eigenvalues = {
        'symmetric': Interval(a+bb+c).centered_widen(3*err),
        'reflection_odd': Interval(a-c).centered_widen(2*err),
        'reflection_even_standard': Interval(a+c-2*bb).centered_widen(4*err),
    }
    return eigenvalues


def sqrt_upper(x):
    """Strict rational square-root upper bound, by integer square root."""
    x = F(x)
    assert x >= 0
    root = isqrt((x.numerator*SCALE*SCALE)//x.denominator)
    return F(root+1, SCALE)


def strong_neighborhood(es, effect_error, cs, ds, b, alpha, beta, pi, roots):
    """All-outcome calibration and sectorwise quadratic variance enclosure."""
    scores = [Interval(v) for v in b]
    scores[PIVOTS[0]] += alpha
    scores[PIVOTS[1]] += beta
    g = [cs, [4*v for v in ds]]
    b0 = [[bp, F(second)] for bp, second in zip(b, RIGHT_INVERSE_SECOND)]
    calibration = [[sum((g[i][p]*b0[p][j] for p in range(8)), Interval(0))
                    for j in range(2)] for i in range(2)]
    det = calibration[0][0]*calibration[1][1]-calibration[0][1]*calibration[1][0]
    assert det.lo > F(999, 1000) and det.hi < F(1001, 1000)
    inverse = [[(row[0]*calibration[1][1]-row[1]*calibration[1][0])/det,
                (-row[0]*calibration[0][1]+row[1]*calibration[0][0])/det]
               for row in b0]
    row_squares = [sum(v.abs_upper()**2 for v in row) for row in inverse]
    sum_row_norms = sum(sqrt_upper(v) for v in row_squares)
    inverse_frobenius = sqrt_upper(sum(row_squares))
    score_l1 = sum(v.abs_upper() for v in scores)
    score_l2 = sqrt_upper(sum(v.abs_upper()**2 for v in scores))
    sectors = [[], [], []]
    for e in es:
        a = F(e[0][0][0]+e[2][2][0], 2*SCALE)
        bb = F(e[0][1][0]+e[2][1][0], 2*SCALE)
        c = F(e[0][2][0]+e[2][0][0], 2*SCALE)
        for dest, value, multiplier in zip(sectors, (a+bb+c, a-c, a+c-2*bb), (3, 2, 4)):
            enclosed = Interval(value).centered_widen(multiplier*effect_error)
            # These are eigenvalues of genuine effects, hence lie in [0,1].
            dest.append(Interval(max(F(0), enclosed.lo), min(F(1), enclosed.hi)))
    effect_radius = F(15, 10000)
    column_radius = 4*effect_radius/roots[3].lo
    contraction = column_radius*sum_row_norms
    assert contraction < 1
    correction_radius = column_radius*score_l1/(1-contraction)
    effect_moment_error = effect_radius*(score_l2+inverse_frobenius*correction_radius)**2
    sector_report = []
    sector_bounds = []
    for values in sectors:
        baseline = sum((v*a*a for v, a in zip(values, scores)), Interval(0))
        linear = [sum((v*a*row[j] for v, a, row in zip(values, scores, inverse)), Interval(0))
                  for j in range(2)]
        linear_norm = sqrt_upper(sum(v.abs_upper()**2 for v in linear))
        # The quadratic coefficient matrix is positive semidefinite, so its
        # trace bounds its largest eigenvalue without an eigenvalue algorithm.
        quadratic_trace = sum(v.hi*row_square for v, row_square in zip(values, row_squares))
        upper = (baseline.hi+2*correction_radius*linear_norm
                 +correction_radius**2*quadratic_trace+effect_moment_error)
        sector_bounds.append(upper)
        sector_report.append({
            'baseline': display(baseline),
            'linear_vector_norm_upper': decimal_bound(linear_norm, 12, True),
            'quadratic_matrix_trace_upper': decimal_bound(quadratic_trace, 12, True),
            'perturbed_second_moment_upper': decimal_bound(upper, 12, True),
        })
    second_upper = max(sector_bounds)
    assert second_upper < F(54, 100) < F(9, 16)
    parameter_radius = F(9, 100000)
    parameter_effect_bound = 6*parameter_radius+3*(pi.hi+parameter_radius)*parameter_radius
    assert parameter_effect_bound < effect_radius
    return {
        'effect_operator_norm_radius': '0.0015',
        'requires': 'complete SU(2)-invariant effects, flux reflection, exact recalibration',
        'right_inverse_definition': 'B=B0 (G B0)^(-1), G=(c;4d)',
        'right_inverse_reference_second_column': list(RIGHT_INVERSE_SECOND),
        'right_inverse_normalizer_determinant': display(det),
        'score_l1_upper': decimal_bound(score_l1, 12, True),
        'score_l2_upper': decimal_bound(score_l2, 12, True),
        'sum_right_inverse_row_norms_upper': decimal_bound(sum_row_norms, 12, True),
        'right_inverse_frobenius_norm_upper': decimal_bound(inverse_frobenius, 12, True),
        'neumann_contraction_upper': decimal_bound(contraction, 12, True),
        'two_component_correction_norm_upper': decimal_bound(correction_radius, 12, True),
        'sectors': dict(zip(['symmetric', 'reflection_odd', 'reflection_even_standard'], sector_report)),
        'second_moment_operator_norm_upper': decimal_bound(second_upper, 12, True),
        'certified_less_than': '0.54 < 9/16',
        'uniform_ring_parameter_box': {
            'strict_coordinate_radius': '0.00009',
            'coordinates': ['u', 'tau', 'phi'],
            'effect_operator_norm_upper': decimal_bound(parameter_effect_bound, 15, True),
            'requires': 'opposite flux signs and exact recalibration at actual parameters',
        },
        'scope': 'A sufficient stability theorem, not an apparatus tolerance estimate.',
    }


def positive_exponential(x, order=20):
    """Rational enclosure for exp(x), x>=0, by a positive series tail."""
    x = F(x)
    assert 0 <= x < order+2
    term = total = F(1)
    for k in range(1, order+1):
        term = term*x/k
        total += term
    first_omitted = term*x/(order+1)
    # Subsequent ratios are at most x/(order+2).
    remainder = first_omitted/(1-x/(order+2))
    return Interval(total, total+remainder)


def noise_certificate(b, alpha_correction, beta_correction, mixed, t_upper):
    scores = [Interval(v) for v in b]
    scores[PIVOTS[0]] += alpha_correction
    scores[PIVOTS[1]] += beta_correction
    a_upper = max(v.abs_upper() for v in scores)
    singletons = [v for p, v in zip(PATTERNS, scores) if sum(p) == 1]
    c1 = Interval((max(v.lo for v in singletons)-min(v.hi for v in singletons))/2,
                  (max(v.hi for v in singletons)-min(v.lo for v in singletons))/2)
    a_bound, c1_bound = F(2337, 1000), F(1499, 1000)
    t_bound, mixed_bound = F(390, 1000), F(285, 1000)
    assert a_upper < a_bound and c1.hi < c1_bound
    assert t_upper < t_bound and mixed.hi < mixed_bound
    exp_lambda = positive_exponential(F(1, 100))
    exp_3lambda = positive_exponential(F(3, 100))
    z = exp_lambda-1
    q = F(3, 1000)
    no_flip = (1-q)**4
    even_flip = 6*q*q*(1-q)**2 + q**4
    odd_flip = 4*q*(1-q)**3 + 4*q**3*(1-q)
    bias = (3*c1_bound*z*z
            + a_bound*even_flip/no_flip*(1+3*z*z)
            + 3*a_bound*odd_flip/no_flip*z)

    def variance_budget(q0):
        numerator = (no_flip*q0 + even_flip*a_bound*a_bound
                     + 3*z*z*(no_flip+even_flip)*a_bound*a_bound
                     + (3*z+z*z*z)*odd_flip*a_bound*a_bound)
        return exp_3lambda*numerator/(no_flip*no_flip)

    all_input = variance_budget(t_bound)
    mixed_input = variance_budget(mixed_bound)
    assert bias.hi < F(1430, 10**6)
    assert all_input.hi < F(411, 1000)
    assert mixed_input.hi < F(302, 1000)
    return {
        'score_absolute_upper': decimal_bound(a_upper, 15, True),
        'singleton_score_half_range': display(c1),
        'coarse_bounds': {'A': '2.337', 'C1': '1.499',
                          'T_operator_norm': '0.390', 'Q_maximally_mixed': '0.285'},
        'uniform_window': {'lambda': ['0', '0.01'], 'bit_error_q': ['0', '0.003']},
        'requires': ('equal spin-independent loss, independent symmetric bit errors, '
                     'known survival-and-no-flip correction, exactly calibrated nominal scores'),
        'exponential_series_order': 20,
        'exp_lambda_endpoint': display(exp_lambda, 18),
        'exp_3lambda_endpoint': display(exp_3lambda, 18),
        'bias_upper': decimal_bound(bias.hi, 15, True),
        'all_input_second_moment_upper': decimal_bound(all_input.hi, 15, True),
        'maximally_mixed_second_moment_upper': decimal_bound(mixed_input.hi, 15, True),
        'coarse_certified_bounds': {'bias': '0.001430', 'all_input_second_moment': '0.411',
                                    'maximally_mixed_second_moment': '0.302'},
        'scope': 'Calibration and noise-rate uncertainty require separate budgets.',
    }


def channel_recalibration_certificate(es, error, cs, ds, roots):
    """Exact auxiliary-D reconstruction and rate-free two-feature readout."""
    reference = tuple(4*F(v) for v in RIGHT_INVERSE_SECOND)
    det, alpha, beta, correction_norm = correction(cs, ds, reference, target=(0, 1))
    scores = [Interval(v) for v in reference]
    scores[PIVOTS[0]] += alpha
    scores[PIVOTS[1]] += beta
    auxiliary_eigenvalues = second_moment(es, error, reference)
    auxiliary_upper = max(v.hi for v in auxiliary_eigenvalues.values())+correction_norm
    auxiliary_a = max(v.abs_upper() for v in scores)
    singleton = [v for p, v in zip(PATTERNS, scores) if sum(p) == 1]
    auxiliary_c1 = (max(v.hi for v in singleton)-min(v.lo for v in singleton))/2
    assert auxiliary_a < 194 and auxiliary_c1 < 107 and auxiliary_upper < 646
    z = positive_exponential(F(1, 100))-1
    exp_3lambda = positive_exponential(F(3, 100))
    q = F(3, 1000)
    no_flip = (1-q)**4
    even_flip = 6*q*q*(1-q)**2+q**4
    odd_flip = 4*q*(1-q)**3+4*q**3*(1-q)

    def budgets(a, c1, t):
        bias = (3*c1*z*z+a*even_flip/no_flip*(1+3*z*z)
                +3*a*odd_flip/no_flip*z)
        second = exp_3lambda*(no_flip*t+even_flip*a*a
                             +3*z*z*(no_flip+even_flip)*a*a
                             +(3*z+z*z*z)*odd_flip*a*a)/(no_flip*no_flip)
        return bias.hi, second.hi

    error_a, variance_a = budgets(F(2337, 1000), F(1499, 1000), F(390, 1000))
    error_v, variance_v = budgets(F(194), F(107), F(646))
    # Normalized response M has rows (1+dc_a, dd_a), (dc_v,1+dd_v).
    # The response norm constraints follow from the anticommuting odd basis.
    dc_a, dd_a = 4*error_a/roots[3].lo, error_a/roots[3].lo
    dc_v, dd_v = 4*error_v/roots[3].lo, error_v/roots[3].lo
    assert dd_v < 1
    schur_lower = 1-dc_a-dd_a*dc_v/(1-dd_v)
    assert schur_lower > 0
    determinant_lower = (1-dd_v)*schur_lower
    assert determinant_lower > F(93, 100)
    # Use the Schur complement before bounding; bounding the inverse's
    # numerator and determinant independently loses their cancellation.
    wa = 1/schur_lower
    wv = wa*dd_a/(1-dd_v)
    corrected_second = (wa*sqrt_upper(variance_a)+wv*sqrt_upper(variance_v))**2
    assert corrected_second < F(446, 1000) < F(9, 16)
    raw_determinant_lower = (no_flip/exp_3lambda.hi)**2*determinant_lower
    absolute_score_upper = (wa*F(2337, 1000)+wv*194)*exp_3lambda.hi/no_flip
    assert absolute_score_upper < F(263, 100)
    return {
        'auxiliary_target': 'D=P12-P23',
        'rational_reference_scores': [str(v) for v in reference],
        'exact_two_outcome_corrections': {'alpha': display(alpha, 18), 'beta': display(beta, 18)},
        'auxiliary_score_absolute_upper': decimal_bound(auxiliary_a, 12, True),
        'auxiliary_singleton_half_range_upper': decimal_bound(auxiliary_c1, 12, True),
        'auxiliary_second_moment_norm_upper': decimal_bound(auxiliary_upper, 12, True),
        'auxiliary_coarse_bounds': {'A_D': '194', 'C1_D': '107', 'T_D_norm': '646'},
        'uniform_window': {'lambda': ['0', '0.01'], 'bit_error_q': ['0', '0.003']},
        'normalized_feature_budgets': {
            'chi_mean_error_upper': decimal_bound(error_a, 15, True),
            'D_mean_error_upper': decimal_bound(error_v, 15, True),
            'chi_second_moment_upper': decimal_bound(variance_a, 15, True),
            'D_second_moment_upper': decimal_bound(variance_v, 15, True),
        },
        'normalized_response_determinant_lower': decimal_bound(determinant_lower, 15),
        'raw_response_determinant_lower': decimal_bound(raw_determinant_lower, 15),
        'normalized_chi_weight_upper': decimal_bound(wa, 15, True),
        'normalized_D_weight_absolute_upper': decimal_bound(wv, 15, True),
        'channel_recalibrated_score_absolute_upper': decimal_bound(absolute_score_upper, 12, True),
        'channel_recalibrated_second_moment_upper': decimal_bound(corrected_second, 12, True),
        'certified_less_than': '0.446 < 9/16',
        'requires': ('exact full-channel response coefficients; equal spin-independent loss; '
                     'independent symmetric parity-bit errors; stable symmetry-preserving channel'),
        'operational_rates_required': False,
        'operational_calibration': 'Use raw two-feature response; common survival/no-flip factor cancels.',
        'finite_data_scope': 'Finite empirical calibration needs CHANNEL_CALIBRATION.md validation bounds.',
    }


def calibration_moment_certificate(es, error, b, alpha, beta):
    """Nominal product-probe variance profiles from exact score moments."""
    scores = [Interval(v) for v in b]
    scores[PIVOTS[0]] += alpha
    scores[PIVOTS[1]] += beta
    sectors = [[], [], []]
    for e in es:
        aa = F(e[0][0][0]+e[2][2][0], 2*SCALE)
        bb = F(e[0][1][0]+e[2][1][0], 2*SCALE)
        cc = F(e[0][2][0]+e[2][0][0], 2*SCALE)
        for dest, value, multiple in zip(sectors, (aa+bb+cc, aa-cc, aa+cc-2*bb), (3, 2, 4)):
            enclosed = Interval(value).centered_widen(multiple*error)
            dest.append(Interval(max(F(0), enclosed.lo), min(F(1), enclosed.hi)))
    second = [sum((w*v*v for w, v in zip(sector, scores)), Interval(0))
              for sector in sectors]
    fourth = [sum((w*v*v*v*v for w, v in zip(sector, scores)), Interval(0))
              for sector in sectors]
    specifications = {
        'xyz': (F(1, 8), (F(1, 2), F(1, 4), F(1, 4))),
        '001': (F(0), (F(1, 3), F(1, 2), F(1, 6))),
        '000': (F(0), (F(1), F(0), F(0))),
        '010': (F(0), (F(1, 3), F(0), F(2, 3))),
    }
    probes = {}
    for name, (mean, weights) in specifications.items():
        q = sum((w*v for w, v in zip(weights, second)), Interval(0))
        fourth_moment = sum((w*v for w, v in zip(weights, fourth)), Interval(0))
        variance = q-mean*mean
        squared_score_variance = fourth_moment-q*q
        assert variance.lo > 0 and squared_score_variance.lo > 0
        probes[name] = {
            'mean_exact': str(mean),
            'sector_weights_exact': [str(w) for w in weights],
            'second_moment': display(q),
            'fourth_moment': display(fourth_moment),
            'variance_X': display(variance),
            'variance_X_squared': display(squared_score_variance),
        }
    return {
        'scope': ('Exact nominal ideal-model planning profiles; not measured variances '
                  'of an implemented channel or a pilot-fitted score.'),
        'flux_allocation': 'independent equal-probability randomization',
        'score_definition': 'nominal rational reference scores plus exact two-outcome correction',
        'score_absolute_upper': decimal_bound(max(v.abs_upper() for v in scores), 15, True),
        'sector_order': ['symmetric', 'reflection_odd', 'reflection_even_standard'],
        'second_moment_eigenvalues': [display(v) for v in second],
        'fourth_moment_eigenvalues': [display(v) for v in fourth],
        'probes': probes,
    }


def run():
    if not __debug__:
        raise RuntimeError("Run the certificate without Python optimization.")
    pi, roots, constant_report = constants()
    occ, matrix, inputs, model_error_units = finite_model(pi, roots)
    action, action_error = taylor_action(matrix, inputs, model_error_units)
    es, effect_error = effects(occ, action, action_error)
    cs, ds = calibration_intervals(es, effect_error)
    b = tuple(F(value) for value in SCORE_DECIMALS)
    det, alpha, beta, correction_norm = correction(cs, ds, b)
    eig = second_moment(es, effect_error, b)
    calibrated_eig = {k: v.centered_widen(correction_norm) for k, v in eig.items()}
    upper = max(v.hi for v in calibrated_eig.values())
    assert upper < F(39, 100)  # A stronger nominal bound than 0.4.
    mixed = (2*calibrated_eig['symmetric'] + calibrated_eig['reflection_odd']
             + calibrated_eig['reflection_even_standard'])/4
    noise_report = noise_certificate(b, alpha, beta, mixed, upper)
    strong_report = strong_neighborhood(es, effect_error, cs, ds, b, alpha, beta, pi, roots)
    channel_report = channel_recalibration_certificate(es, effect_error, cs, ds, roots)
    moment_report = calibration_moment_certificate(es, effect_error, b, alpha, beta)
    # A finite open parameter box, all parameters varied while preserving the
    # uniform spin-independent ring and using recalibrated outcome scores.
    radius = F(1, 10**7)
    effect_variation = 6*radius + 3*(pi.hi+radius)*radius
    c_variation = 4*effect_variation/roots[3].lo
    d_variation = effect_variation/roots[3].lo
    box_c = [v.centered_widen(c_variation) for v in cs]
    box_d = [v.centered_widen(d_variation) for v in ds]
    box_det, box_alpha, box_beta, box_correction = correction(box_c, box_d, b)
    box_fixed_score_norm = sum(bp*bp for bp in b)*effect_variation
    box_second_upper = max(v.hi for v in eig.values()) + box_fixed_score_norm + box_correction
    assert box_second_upper < F(2, 5)
    # A larger, coordinate-free neighborhood of effects; the symmetry
    # assumptions guaranteeing the two-dimensional odd space remain required.
    effect_radius = F(1, 10**5)
    ball_c = [v.centered_widen(4*effect_radius/roots[3].lo) for v in cs]
    ball_d = [v.centered_widen(effect_radius/roots[3].lo) for v in ds]
    ball_det, ball_alpha, ball_beta, ball_correction = correction(ball_c, ball_d, b)
    ball_upper = (max(v.hi for v in eig.values())
                  + sum(bp*bp for bp in b)*effect_radius + ball_correction)
    assert ball_upper < F(9, 20)
    return {
        'certificate': 'PASS: exact integer and rational enclosures; no floating arithmetic',
        'parameters': {'u': 16, 'tau': 'pi', 'phi': '-2*pi/3'},
        'dimension': 40,
        'input_sector': ['001', '010', '100'],
        'fixed_point_bits': BITS,
        'taylor_order': ORDER,
        'matrix_infinity_norm_upper': NORM_BOUND,
        **constant_report,
        'matrix_rounding_error_units': model_error_units,
        'action_error_strict_upper': '1e-170',
        'patterns': [''.join(map(str, p)) for p in PATTERNS],
        'rational_reference_scores': list(SCORE_DECIMALS),
        'calibration_c': [display(v) for v in cs],
        'calibration_d': [display(v) for v in ds],
        'pivot_patterns': [''.join(map(str, PATTERNS[j])) for j in PIVOTS],
        'pivot_determinant': display(det),
        'exact_score_correction_intervals': {'alpha': display(alpha, 18), 'beta': display(beta, 18)},
        'calibrated_second_moment_eigenvalues': {k: display(v) for k, v in calibrated_eig.items()},
        'calibrated_second_moment_norm_upper': decimal_bound(upper, 12, True),
        'maximally_mixed_second_moment': display(mixed, 12),
        'noise_window': noise_report,
        'all_outcome_robustness': strong_report,
        'channel_recalibration': channel_report,
        'calibration_moments': moment_report,
        'effect_neighborhood': {
            'uniform_effect_operator_norm_radius': '0.00001',
            'requires': 'SU(2) invariance, flux reflection, exact recalibration, complete POVM',
            'pivot_determinant': display(ball_det),
            'score_correction_intervals': {'alpha': display(ball_alpha), 'beta': display(ball_beta)},
            'second_moment_norm_upper': decimal_bound(ball_upper, 12, True),
            'certified_less_than': '0.45',
        },
        'parameter_box': {
            'strict_coordinate_radius': '0.0000001',
            'coordinates': ['u', 'tau', 'phi'],
            'requires': 'uniform spin-independent ring, exact recalibration at actual parameters',
            'effect_variation_upper': decimal_bound(effect_variation, 15, True),
            'pivot_determinant': display(box_det),
            'score_correction_intervals': {'alpha': display(box_alpha), 'beta': display(box_beta)},
            'second_moment_norm_upper': decimal_bound(box_second_upper, 12, True),
            'certified_less_than': '0.4',
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional JSON report outside the repository.')
    args = parser.parse_args()
    report = run()
    text = json.dumps(report, indent=2) + '\n'
    if args.output:
        repo = Path(__file__).resolve().parents[1]
        target = args.output.resolve()
        assert not target.is_relative_to(repo), 'Generated evidence must remain outside the repository.'
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('x') as report_file:
            report_file.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
