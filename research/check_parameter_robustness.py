#!/usr/bin/env python3
"""Exact arithmetic for the primary-score robustness consequences.

The operator and statistical arguments are in PARAMETER_ROBUSTNESS.md.
This checker verifies their rational consequences, using the nominal exact
parity certificate as a hashed dependency; it is not an experimental record.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from math import isqrt
import json
from pathlib import Path


V = F('0.3287')
A = F('1.623')
VD = F(646)
AD = F(194)
ROOT_SCALE = 10**30


def sqrt_bounds(value):
    value = F(value)
    assert value >= 0
    integer = isqrt(value.numerator*ROOT_SCALE**2//value.denominator)
    lo, hi = F(integer, ROOT_SCALE), F(integer+1, ROOT_SCALE)
    assert lo*lo <= value < hi*hi
    return lo, hi


def decimal(value, places=24, upper=True):
    value = F(value)
    scale = 10**places
    scaled = value*scale
    integer = -((-scaled.numerator)//scaled.denominator) if upper else scaled.__floor__()
    sign = '-' if integer < 0 else ''
    integer = abs(integer)
    return f'{sign}{integer//scale}.{integer%scale:0{places}d}'


def bind_inputs(path):
    raw = path.read_bytes()
    report = json.loads(raw)
    assert report['certificate'] == (
        'PASS: exact integer and rational enclosures; no floating arithmetic'
    )
    assert report['parameters'] == {'u': 16, 'tau': 'pi', 'phi': '-2*pi/3'}
    primary = report['passive_separation']
    auxiliary = report['channel_recalibration']
    assert auxiliary['auxiliary_target'] == 'D=P12-P23'
    fields = {
        'primary_second_moment_upper': primary['calibrated_second_moment_norm_upper'],
        'primary_score_absolute_upper': primary['score_absolute_upper'],
        'auxiliary_second_moment_upper': auxiliary['auxiliary_second_moment_norm_upper'],
        'auxiliary_score_absolute_upper': auxiliary['auxiliary_score_absolute_upper'],
    }
    for value, bound in zip(fields.values(), (V, A, VD, AD)):
        assert 0 < F(value) < bound
    return {
        'certificate_sha256': hashlib.sha256(raw).hexdigest(),
        'certificate_file': path.name,
        'source_objects': ['passive_separation', 'channel_recalibration'],
        'verified_nominal_inputs': fields,
        'coarse_bounds_used': {'v': str(V), 'A': str(A), 'v_D': str(VD), 'A_D': str(AD)},
    }


def dilation_bounds(eta, moment, amplitude):
    root = sqrt_bounds(moment)[1]
    epsilon = 2*root*eta+amplitude*eta*eta
    qroot = root+amplitude*eta
    return epsilon, qroot


def recalibrated_profile(eta, claimed):
    epsilon, qroot = dilation_bounds(eta, V, A)
    epsilon_d, qroot_d = dilation_bounds(eta, VD, AD)
    h_lower = sqrt_bounds(F(3))[0]/4
    x, y = epsilon/h_lower, epsilon/(4*h_lower)
    z, w = epsilon_d/h_lower, epsilon_d/(4*h_lower)
    assert 0 <= w < 1
    kappa = 1-x-y*z/(1-w)
    assert kappa > 0
    second = ((qroot+y*qroot_d/(1-w))/kappa)**2
    assert second < claimed < F(13, 32)
    return {
        'dilation_distance_upper': str(eta),
        'primary_mean_error_upper': decimal(epsilon),
        'auxiliary_mean_error_upper': decimal(epsilon_d),
        'response_schur_lower': decimal(kappa, upper=False),
        'recalibrated_second_moment_upper': decimal(second),
        'second_moment_less_than': str(claimed),
    }


def recalibrated_bounds():
    profiles = [
        recalibrated_profile(F('0.0025'), F('0.394824826')),
        recalibrated_profile(F('0.0015'), F('0.366080282')),
    ]
    information, prior_upper = F(775, 297), F(4000, 3)
    conservative_coefficient = F('0.367')
    assert F(profiles[1]['second_moment_less_than']) < conservative_coefficient
    assert conservative_coefficient*information < 1
    threshold = conservative_coefficient*prior_upper/(1-conservative_coefficient*information)
    assert threshold < 12_000
    return {
        'group': 'exact_response_recalibrated_bounds',
        'profiles': profiles,
        'passive_unbiased_reference': '13/32',
        'small_region_all_estimator_comparison': {
            'passive_information_upper': str(information),
            'prior_information_strict_upper': str(prior_upper),
            'conservative_interacting_coefficient': str(conservative_coefficient),
            'strict_sample_threshold_upper': decimal(threshold),
            'sufficient_integer_sample_size': 12_000,
        },
    }


def physical_budget():
    pi_upper = F(22, 7)
    assert pi_upper*pi_upper < 10
    root10_upper = sqrt_bounds(F(10))[1]
    assert root10_upper < F('3.163')
    # |Phi-Phi_0| <= 0.001. The four-mode hopping spectrum gives
    # max(|cos(Phi/4)|,|sin(Phi/4)|) <= sqrt(3)/2 + 1/8000.
    trigonometric_upper = sqrt_bounds(F(3))[1]/2+F(1, 8000)
    assert F(1, 2)+F(1, 4000) < trigonometric_upper
    hopping_upper, flux_derivative_upper = F('5.197'), F('1.300')
    assert 6*trigonometric_upper < hopping_upper
    assert F(3, 2)*trigonometric_upper < flux_derivative_upper
    u_weight, time_weight, flux_weight = map(F, ('1.643', '3.163', '4.086'))
    u_coefficient = pi_upper*(time_weight+hopping_upper)/F('15.998')
    assert u_coefficient < u_weight
    assert pi_upper*flux_derivative_upper < flux_weight
    profiles = []
    for eta, radius in ((F('0.0025'), F('0.00028')),
                        (F('0.0015'), F('0.00016'))):
        assert eta/u_weight < F('0.002')
        assert eta/flux_weight < F('0.001')
        equal_box_budget = (u_weight+time_weight+flux_weight)*radius
        assert equal_box_budget < eta
        profiles.append({
            'weighted_budget_upper': str(eta),
            'sufficient_equal_coordinate_radius': str(radius),
            'equal_box_weighted_budget': decimal(equal_box_budget),
        })
    return {
        'group': 'physical_parameter_budget_arithmetic',
        'weights': {'abs_u_minus_16': str(u_weight),
                    'abs_tau_minus_pi': str(time_weight),
                    'abs_phi_plus_2pi_over_3': str(flux_weight)},
        'input_action_norm_squared_exact': 10,
        'three_particle_hopping_norm_upper': str(hopping_upper),
        'three_particle_flux_derivative_norm_upper': str(flux_derivative_upper),
        'repulsion_path_coefficient_upper': decimal(u_coefficient),
        'profiles': profiles,
    }


def fixed_score_risk():
    eta, n = F('0.0001'), 100_000
    beta_exact_upper, qroot = dilation_bounds(eta, V, A)
    beta, second = F('0.000114682'), F('0.328887')
    assert beta_exact_upper < beta and qroot*qroot < second
    independent_risk = second/n+beta*beta
    history_risk = (sqrt_bounds(second/n)[1]+beta)**2
    passive_lower = 1/(F(97600, 38907)*n+F(16000, 3))
    assert independent_risk < F('0.000003302022')
    assert history_risk < F('0.000003717980')
    assert passive_lower > F('0.000003903384') > history_risk
    information, prior_upper = F(775, 297), F(4000, 3)
    quadratic = (
        beta*beta*information,
        beta*beta*prior_upper+second*information-1,
        second*prior_upper,
    )
    assert quadratic[0] > 0
    endpoints = {}
    for endpoint, claimed in ((4000, F(-128)), (4_000_000, F(-17558))):
        value = quadratic[0]*endpoint**2+quadratic[1]*endpoint+quadratic[2]
        assert value < claimed < 0
        endpoints[str(endpoint)] = decimal(value)
    return {
        'group': 'fixed_score_finite_risk',
        'dilation_distance_upper': str(eta),
        'uniform_bias_upper': str(beta),
        'uniform_second_moment_upper': str(second),
        'sample_size': n,
        'independent_trial_mse_upper': decimal(independent_risk),
        'history_conditional_mse_upper': decimal(history_risk),
        'passive_all_estimator_mse_lower': decimal(passive_lower, upper=False),
        'independent_trial_interval': {
            'sample_size_inclusive': [4000, 4_000_000],
            'convex_comparison_polynomial_coefficients': list(map(str, quadratic)),
            'endpoint_polynomial_upper': endpoints,
        },
        'history_scope': ('Every conditional instrument must satisfy the same dilation bound. '
                          'The sample-size interval above is only for independent trials.'),
    }


def run(certificate):
    dependency = bind_inputs(certificate)
    groups = [recalibrated_bounds(), physical_budget(), fixed_score_risk()]
    return {
        'status': 'PASS',
        'scientific_groups': len(groups),
        'arithmetic': 'Exact fractions and enclosing integer square roots; no floating arithmetic.',
        'dependency': dependency,
        'groups': groups,
        'assumptions': [
            'Known three-particle input; ideal total-site parity detection.',
            'The calibrated parameter region preserves spin and paired-reflection covariance.',
            'The two flux settings use the same u and tau and opposite fluxes.',
            'Recalibrated bounds use exact true response coefficients, not finite empirical estimates.',
            'The characterized response remains fixed, or its exact coefficients are known per trial.',
            'The fixed-score bounds allow symmetry breaking under the stated uniform dilation bound.',
            'The region is separate from loss/bit-error windows and finite-reference calibration budgets.',
            'Arithmetic verifies consequences of the documented operator and statistical proofs.',
        ],
    }


def main():
    assert __debug__, 'Run without Python optimization.'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, required=True)
    parser.add_argument('--output', type=Path,
                        help='Optional new JSON report outside the repository.')
    args = parser.parse_args()
    report = run(args.certificate.resolve())
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        repo = Path(__file__).resolve().parents[1]
        target = args.output.resolve()
        assert not target.is_relative_to(repo), 'Generated evidence must remain outside the repository.'
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('x') as output:
            output.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
