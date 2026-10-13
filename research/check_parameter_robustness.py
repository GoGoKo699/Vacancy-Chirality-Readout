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
C1 = F('1.123')
VD = F(646)
AD = F(194)
C1D = F(107)
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
        'primary_singleton_half_range_upper': primary['singleton_half_range_upper'],
        'auxiliary_singleton_half_range_upper': auxiliary['auxiliary_singleton_half_range_upper'],
    }
    for value, bound in zip(fields.values(), (V, A, VD, AD, C1, C1D)):
        assert 0 < F(value) < bound
    return {
        'certificate_sha256': hashlib.sha256(raw).hexdigest(),
        'certificate_file': path.name,
        'source_objects': ['passive_separation', 'channel_recalibration'],
        'verified_nominal_inputs': fields,
        'coarse_bounds_used': {'v': str(V), 'A': str(A), 'C1': str(C1),
                               'v_D': str(VD), 'A_D': str(AD), 'C1_D': str(C1D)},
        'nominal_finite_calibration_inputs': primary['finite_calibration_inputs'],
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


def exp_upper(value, order=30):
    """Positive Taylor series plus a geometric bound on its omitted tail."""
    value = F(value)
    assert 0 <= value < order+2
    term = total = F(1)
    for k in range(1, order+1):
        term *= value/k
        total += term
    return total+term*value/(order+1)/(1-value/(order+2))


def noisy_feature_bounds(amplitude, half_range, clean_second, loss, bit_error):
    z = exp_upper(loss)-1
    alpha = (1-bit_error)**4
    odds = bit_error/(1-bit_error)
    even_ratio, odd_ratio = 6*odds**2+odds**4, 4*odds+4*odds**3
    error = (3*half_range*z*z+amplitude*even_ratio*(1+3*z*z)
             +3*amplitude*odd_ratio*z)
    second = exp_upper(3*loss)/alpha*(
        clean_second+even_ratio*amplitude**2
        +3*z*z*(1+even_ratio)*amplitude**2
        +(3*z+z**3)*odd_ratio*amplitude**2)
    return error, second


def joint_profile(eta, loss, bit_error, claimed):
    h_lower = sqrt_bounds(F(3))[0]/4
    features = []
    for moment, amplitude, half_range in ((V, A, C1), (VD, AD, C1D)):
        coherent_error, qroot = dilation_bounds(eta, moment, amplitude)
        noise_error, second = noisy_feature_bounds(
            amplitude, half_range, qroot*qroot, loss, bit_error)
        features.append((coherent_error+noise_error, second))
    (error_a, second_a), (error_d, second_d) = features
    x, y = error_a/h_lower, error_a/(4*h_lower)
    z, w = error_d/h_lower, error_d/(4*h_lower)
    assert 0 <= w < 1
    kappa = 1-x-y*z/(1-w)
    assert kappa > 0
    corrected = ((sqrt_bounds(second_a)[1]
                  +y*sqrt_bounds(second_d)[1]/(1-w))/kappa)**2
    assert corrected < claimed < F(13, 32)
    report = {
        'dilation_distance_upper': decimal(eta),
        'uniform_window': {'lambda': ['0', decimal(loss)],
                           'bit_error_q': ['0', decimal(bit_error)]},
        'normalized_chi_mean_error_upper': decimal(error_a),
        'normalized_D_mean_error_upper': decimal(error_d),
        'normalized_chi_second_moment_upper': decimal(second_a),
        'normalized_D_second_moment_upper': decimal(second_d),
        'response_schur_lower': decimal(kappa, upper=False),
        'channel_recalibrated_second_moment_upper': decimal(corrected),
        'second_moment_less_than': decimal(claimed),
        'operational_rates_required': False,
    }
    return report, (error_a, error_d, second_a, second_d, corrected)


def joint_bounds(nominal_inputs):
    eta = F('0.001')
    full, _ = joint_profile(eta, F('0.01'), F('0.003'), F('0.395873255'))
    loss, bit_error = F('0.005'), F('0.0015')
    small, data = joint_profile(eta, loss, bit_error, F('0.366684486'))
    error_a, error_d, second_a, second_d, corrected = data
    h_lower = sqrt_bounds(F(3))[0]/4
    row1, row2 = error_a/h_lower, error_d/(4*h_lower)
    error_norm = sqrt_bounds(row1*row1+row2*row2)[1]
    assert row1 < 1 and error_norm < 1
    rmin = (1-bit_error)**4/exp_upper(3*loss)
    b = 1-rmin+rmin*error_norm
    first_row = 1-rmin+rmin*row1
    mu, trace = 1+error_norm, second_a+second_d/16
    assert b < 1
    degree_bounds = [second_a, ((1+first_row)*sqrt_bounds(second_a)[1]
                               +row1*sqrt_bounds(second_d/16)[1])**2]
    for degree in (3, 4):
        tail = first_row*b**(degree-1)/(1-b)
        degree_bounds.append((sqrt_bounds(corrected)[1]+sqrt_bounds(trace)[1]*tail)**2)
    inputs = {
        'feature_definition': 'Raw Z=(X_primary,X_D/4), target x=(<chi>,<D>/4)',
        'raw_clean_factor_lower': decimal(rmin, upper=False),
        'normalized_response_row_error_upper': list(map(decimal, (row1, row2))),
        'normalized_response_error_norm_upper': decimal(error_norm),
        'response_residual_norm_upper': decimal(b),
        'first_response_residual_row_norm_upper': decimal(first_row),
        'response_norm_upper': decimal(mu),
        'feature_second_moment_trace_upper': decimal(trace),
        'mean_weight_second_moment_upper': list(map(decimal, degree_bounds)),
        'polynomial_degrees': [1, 2, 3, 4],
    }
    assert inputs['feature_definition'] == nominal_inputs['feature_definition']
    assert inputs['polynomial_degrees'] == nominal_inputs['polynomial_degrees']
    for key in ('response_residual_norm_upper', 'first_response_residual_row_norm_upper',
                'response_norm_upper', 'feature_second_moment_trace_upper'):
        assert 0 < F(inputs[key]) < F(nominal_inputs[key])
    assert all(0 < F(value) < F(nominal)
               for value, nominal in zip(inputs['mean_weight_second_moment_upper'],
                                          nominal_inputs['mean_weight_second_moment_upper']))
    equal_radius = F('0.00011')
    equal_budget = sum(map(F, ('1.643', '3.163', '4.086')))*equal_radius
    assert equal_budget < eta
    information, prior_upper, coefficient = F(775, 297), F(4000, 3), F('0.367')
    assert F(small['second_moment_less_than']) < coefficient
    threshold = coefficient*prior_upper/(1-coefficient*information)
    assert threshold < 12000
    return {
        'group': 'joint_parameter_and_noise_bounds',
        'profiles': [full, small],
        'sufficient_equal_coordinate_radius': str(equal_radius),
        'equal_coordinate_weighted_budget_upper': decimal(equal_budget),
        'small_window_all_estimator_sample_size_sufficient': 12000,
        'small_window_all_estimator_strict_sample_threshold_upper': decimal(threshold),
        'finite_calibration_source_window': {
            'dilation_distance_upper': '0.001',
            'lambda': ['0', '0.005'], 'bit_error_q': ['0', '0.0015'],
        },
        'finite_calibration_inputs': inputs,
        'finite_calibration_domination': (
            'All four norm bounds and four polynomial moment bounds are strictly below '
            'the dependency certificate finite_calibration_inputs. The existing finite '
            'risk profiles apply with their original reference, task, and drift assumptions.'),
        'requires': (
            'Spin and paired-reflection covariance; equal spin-independent Markov loss; '
            'independent identical parity-bit errors. Lambda uses the actual duration. '
            'Exact response coefficients are needed for the unbiased coefficient only; '
            'finite-reference polynomial calibration has its separate unconditional risk.'),
    }


def microscopic_and_loading_bounds():
    pi_upper = F(22, 7)
    pair_weight = pi_upper*(sqrt_bounds(F(10))[1]+3*sqrt_bounds(F(3))[1])/16
    one_body_weight = 3*pi_upper
    assert one_body_weight < F('9.429') and pair_weight < F('1.643')
    eta, loading_distance, n = F('0.0001'), F('0.0001'), 100_000
    beta, qroot = dilation_bounds(eta, V, A)
    assert beta < F('0.000114682') and qroot*qroot < F('0.328887')
    # Use the published outward bounds to check the displayed example.
    total_bias = F('0.000114682')+2*A*loading_distance
    total_second = F('0.328887')+A*A*loading_distance
    independent_risk = total_second/n+total_bias*total_bias
    history_risk = (sqrt_bounds(total_second/n)[1]+total_bias)**2
    passive_lower = 1/(F(97600, 38907)*n+F(16000, 3))
    assert independent_risk < F('0.000003484473')
    assert independent_risk < F('0.000003903384') < passive_lower
    return {
        'group': 'microscopic_control_and_loading_bounds',
        'one_body_integrated_weight_upper': decimal(one_body_weight),
        'onsite_pair_integrated_weight_upper': decimal(pair_weight),
        'sufficient_uniform_weights': {'one_body': '9.429', 'onsite_pair': '1.643'},
        'dilation_distance_upper': str(eta),
        'input_half_trace_distance_upper': str(loading_distance),
        'sample_size': n,
        'combined_bias_upper': decimal(total_bias),
        'combined_second_moment_upper': decimal(total_second),
        'independent_trial_mse_upper': decimal(independent_risk),
        'history_conditional_mse_upper': decimal(history_risk),
        'passive_all_estimator_mse_lower': decimal(passive_lower, upper=False),
        'scope': (
            'Microscopic dynamics stay in four orbitals and two internal levels; '
            'pair perturbations preserve occupation blocks. Loading distance is '
            'a uniform half-trace-distance assumption, not a leakage probability. '
            'All parity patterns are scored without postselection. The displayed '
            'passive separation is for independent trials only.'),
    }


def run(certificate):
    dependency = bind_inputs(certificate)
    groups = [recalibrated_bounds(), physical_budget(), fixed_score_risk(),
              joint_bounds(dependency['nominal_finite_calibration_inputs']),
              microscopic_and_loading_bounds()]
    return {
        'status': 'PASS',
        'scientific_groups': len(groups),
        'arithmetic': ('Exact fractions, enclosing integer square roots, and positive '
                       'Taylor exponential tails; no floating arithmetic.'),
        'dependency': dependency,
        'groups': groups,
        'assumptions': [
            'Known three-particle input except for the explicit loading-distance allowance; '
            'ideal detection except in the explicit joint-noise group.',
            'The calibrated parameter region preserves spin and paired-reflection covariance.',
            'The two flux settings use the same u and tau and opposite fluxes.',
            'Recalibrated bounds use exact true response coefficients, not finite empirical estimates.',
            'The characterized response remains fixed, or its exact coefficients are known per trial.',
            'The fixed-score bounds allow symmetry breaking under the stated uniform dilation bound.',
            'Only the explicit joint group combines parameter variation with loss and bit errors.',
            'The smaller joint window supplies finite-reference inputs dominated by the nominal bounds.',
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
