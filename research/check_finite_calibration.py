#!/usr/bin/env python3
"""Exact finite-calibration risk and total-preparation certificates.

Consumes the integer/rational parity certificate; no floating arithmetic,
empirical calibration record, or confidence-event substitution is used.
See FINITE_CALIBRATION.md for the matrix covariance and statistical proofs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import parity_certificate as pc


OUTWARD = {
    'response_residual_norm_upper': F('0.103947447477219'),
    'first_response_residual_row_norm_upper': F('0.043400571947183'),
    'response_norm_upper': F('1.065494885645340'),
    'feature_second_moment_trace_upper': F('44.219368969990164'),
}
V = tuple(map(F, ('0.344784300185690', '0.395119522625782',
                  '0.370837757092059', '0.366919338686148')))
STABLE_STAGES = (1_800_000, 18_000_000, 180_000_000)
ROBUST_STAGES = (2_000_000, 20_000_000, 200_000_000)
H2 = F(3, 16)
EPSILON_SQUARED = F(1, 1_000_000)


def upper(value: F, places: int = 24) -> str:
    return pc.decimal_bound(value, places, upper=True)


def bind_inputs(path: Path):
    raw = path.read_bytes()
    report = json.loads(raw)
    assert report['certificate'] == (
        'PASS: exact integer and rational enclosures; no floating arithmetic'
    ), 'The dependency is not a passing exact parity certificate.'
    assert report['parameters'] == {'u': 16, 'tau': 'pi', 'phi': '-2*pi/3'}
    source = report['passive_separation']
    assert source['full_channel']['uniform_window'] == {
        'lambda': ['0', '0.01'], 'bit_error_q': ['0', '0.003']
    }
    inputs = source['finite_calibration_inputs']
    assert inputs['polynomial_degrees'] == [1, 2, 3, 4]
    assert inputs['feature_definition'] == (
        'Raw Z=(X_primary,X_D/4), target x=(<chi>,<D>/4)'
    )
    for key, bound in OUTWARD.items():
        supplied = F(inputs[key])
        assert 0 <= supplied <= bound, (key, supplied, bound)
    supplied_v = inputs['mean_weight_second_moment_upper']
    assert len(supplied_v) == len(V)
    for supplied, bound in zip(supplied_v, V):
        assert 0 <= F(supplied) <= bound
    return {
        'certificate_sha256': hashlib.sha256(raw).hexdigest(),
        'certificate_file': path.name,
        'source_object': 'passive_separation.finite_calibration_inputs',
        'verified_supplied_inputs': inputs,
        'outward_bounds_used': {key: str(value) for key, value in OUTWARD.items()},
        'mean_weight_second_moment_bounds_used': list(map(str, V)),
    }


def passive_lower(pi):
    # Uniform pointwise FI bound over [-3h/40,3h/40].  A prior-averaged
    # per-copy bound cannot be substituted under adaptive task allocation.
    radius = F(3, 40)
    information = 16*(6+radius)/(3*(1-radius*radius)*(13+radius))
    assert information == F(2073600, 832093)
    assert 1/(H2*radius*radius) == F(25600, 27)
    prior_information = pi*pi/(H2*radius*radius)
    lower = (1/EPSILON_SQUARED-prior_information)/information
    assert lower.lo > 0
    return lower


def count_comparison(per_task, stages, n, tasks):
    total = sum(stages)+n*tasks
    bound = per_task*tasks
    low_ceil = pc.ceil_fraction(bound.lo)
    high_ceil = pc.ceil_fraction(bound.hi)
    assert low_ceil == high_ceil, 'The integer passive bound is unresolved.'
    assert F(total) < bound.lo, 'No preparation-count separation certified.'
    return {
        'tasks': tasks,
        'target_preparations_per_task': n,
        'reference_preparations': sum(stages),
        'interacting_total_preparations': total,
        'passive_expected_total_necessary_interval': pc.display(bound, 12),
        'passive_fixed_integer_budget_necessary': low_ceil,
        'scope': ('The ceiling applies to fixed or almost-sure integer budgets; '
                  'the unrounded interval applies to expected budgets.'),
    }


def stable_profile(per_task):
    b = OUTWARD['response_residual_norm_upper']
    r = OUTWARD['first_response_residual_row_norm_upper']
    mu = OUTWARD['response_norm_upper']
    h = OUTWARD['feature_second_moment_trace_upper']
    c = F(0)
    covariance = [c]
    allocations = []
    for nstage, v in zip(STABLE_STAGES, V):
        assert nstage % 5 == 0
        nx, nd = 4*nstage//5, nstage//5
        assert F(64, nx) == F(16, nd) == F(80, nstage)
        allocations.append({'x_probe': nx, 'd_probe': nd})
        c = b*b*c+F(80, nstage)*(v+h*c)
        covariance.append(c)
    bias_squared = H2*r*r*b**6
    mean_squared = H2*mu*mu
    n = 382_210
    risk = (V[3]+h*c)/n+mean_squared*c+bias_squared
    assert risk < F('0.0000009999796322') < EPSILON_SQUARED
    comparison = count_comparison(per_task, STABLE_STAGES, n, 13_100)
    assert comparison['interacting_total_preparations'] == 5_206_751_000
    assert comparison['passive_fixed_integer_budget_necessary'] == 5_207_568_409
    # The reference design attains the isotropic homoscedastic bound I/80.
    assert F(4, 5)*F(1, 8)**2 == F(1, 5)*F(1, 4)**2 == F(1, 80)
    return {
        'group': 'unconditional_finite_calibration_and_reuse',
        'reference_stage_totals': list(STABLE_STAGES),
        'reference_stage_allocations': allocations,
        'weight_covariance_operator_norm_bounds': list(map(upper, covariance)),
        'deterministic_bias_squared_upper': upper(bias_squared),
        'squared_target_mean_norm_bound': upper(mean_squared),
        'unconditional_mse_upper': upper(risk),
        'unconditional_mse_less_than': '0.0000009999796322',
        'comparison': comparison,
    }


def sqrt_bound(value):
    bound = pc.sqrt_upper(value)
    assert bound >= 0 and bound*bound >= value
    return bound


def robust_profile(per_task):
    b = OUTWARD['response_residual_norm_upper']
    r = OUTWARD['first_response_residual_row_norm_upper']
    mu = OUTWARD['response_norm_upper']
    h = OUTWARD['feature_second_moment_trace_upper']
    e, d, t = map(F, ('0.000081', '0.00002', '0.001'))
    q = b+e
    s0 = (F(1), 1+r, 1+r*(1+b), 1+r*(1+b+b*b))
    delta = F(0)
    c = F(0)
    deltas, covariance, moments = [delta], [c], []
    for j, nstage in enumerate(ROBUST_STAGES):
        vprime = (sqrt_bound(V[j])+sqrt_bound(h)*delta)**2
        moments.append(vprime)
        c = q*q*c+F(80, nstage)*(vprime+h*c)
        delta = q*delta+e*s0[j]
        covariance.append(c)
        deltas.append(delta)
    vprime4 = (sqrt_bound(V[3])+sqrt_bound(h)*delta)**2
    moments.append(vprime4)
    s4 = 1+(r+e)*(1+q+q*q)
    bias_squared = H2*((r+e)*q**3+s4*(e+d))**2
    n = 384_220
    sampling = (vprime4+h*c+t*(s4*s4+2*c))/n
    calibration = H2*(mu+d)**2*c
    risk = sampling+calibration+bias_squared
    assert risk < F('0.0000009999754581') < EPSILON_SQUARED
    angle = F('0.005')
    response_error = mu*(3*angle*angle+angle**3)
    assert response_error < F('0.000080045305') < e
    assert angle*angle/4 <= 3*angle*angle+angle**3
    assert 1-3*angle*angle-angle**3 > 0
    comparison = count_comparison(per_task, ROBUST_STAGES, n, 16_800)
    assert comparison['interacting_total_preparations'] == 6_676_896_000
    assert comparison['passive_fixed_integer_budget_necessary'] == 6_678_408_341
    return {
        'group': 'preparation_error_and_channel_drift',
        'reference_response_error_norm_upper': str(e),
        'target_response_drift_norm_upper': str(d),
        'feature_second_moment_matrix_drift_norm_upper': str(t),
        'reference_stage_totals': list(ROBUST_STAGES),
        'mean_weight_deviation_bounds': list(map(upper, deltas)),
        'mean_weight_second_moment_bounds': list(map(upper, moments)),
        'weight_covariance_operator_norm_bounds': list(map(upper, covariance)),
        'final_mean_weight_norm_upper': upper(s4),
        'deterministic_bias_squared_upper': upper(bias_squared),
        'target_sampling_mse_upper': upper(sampling),
        'random_calibration_mse_upper': upper(calibration),
        'unconditional_mse_upper': upper(risk),
        'unconditional_mse_less_than': '0.0000009999754581',
        'angular_preparation': {
            'pure_bloch_vector_angle_radians_upper': str(angle),
            'response_error_norm_upper': str(response_error),
            'requires': ('Cyclic permutation of the same actual pure triad; '
                         'second probe has identical first two Bloch vectors; '
                         'errors remain within the same orbital embedding.'),
        },
        'comparison': comparison,
    }


def run(certificate: Path):
    dependency = bind_inputs(certificate)
    pi, _, _ = pc.constants()
    per_task = passive_lower(pi)
    assert pc.ceil_fraction(per_task.lo) == pc.ceil_fraction(per_task.hi) == 397_525
    groups = [stable_profile(per_task), robust_profile(per_task)]
    return {
        'status': 'PASS',
        'scientific_groups': len(groups),
        'arithmetic': 'Exact fractions and outward rational intervals; no floating arithmetic.',
        'dependency': dependency,
        'target_rmse': '0.001',
        'passive_prior_radius_in_units_of_h': '3/40',
        'passive_uniform_fisher_information_bound': '2073600/832093',
        'passive_prior_information_pi_squared_factor': '25600/27',
        'passive_per_task_expected_count_necessary_interval': pc.display(per_task, 15),
        'groups': groups,
        'assumptions': [
            'Spin and paired-reflection covariance; fixed predetermined two-feature score.',
            'Independent reference stages and target trials; every preparation is counted.',
            'Targets are fixed independently of the shared calibration randomness.',
            'Separate arbitrary-state tasks; no promised common target model.',
            'Equal preparation costs for the total-count comparison.',
            'Noise, preparation, and drift windows are model assumptions, not apparatus data.',
        ],
    }


def main():
    assert __debug__, 'Run without Python optimization.'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, required=True,
                        help='Produced exact complete-parity certificate JSON.')
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
