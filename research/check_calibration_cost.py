#!/usr/bin/env python3
"""Three cross-check groups and one explicitly hypothetical calibration forecast."""
from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
import complete_parity as cp
from check_channel_calibration import PROBES
from calibration_cost import (ProbeSummary, bernstein_radius, certificate,
                              cycle_reuse_threshold, operator_allowances, stable_reuse_capacity,
                              leading_plan, preparation_budget, separable_required,
                              strict_reuse_threshold)


def midpoint(interval):
    return float((Fraction(interval[0])+Fraction(interval[1]))/2)


def run(exact_report):
    assert exact_report['certificate'] == 'PASS: exact integer and rational enclosures; no floating arithmetic'
    moments = exact_report['calibration_moments']
    cal = cp.calibration()
    fourth = np.einsum('p,pij->ij', cal['scores']**4, cal['even'])
    moment_errors = []
    for name, rho in zip(('xyz', '001', '000'), PROBES):
        p = moments['probes'][name]
        mean = np.trace(rho@cal['estimate']).real
        second = np.trace(rho@cal['second_moment']).real
        fourth_value = np.trace(rho@fourth).real
        for numerical, certified in ((mean, float(Fraction(p['mean_exact']))),
                                     (second, midpoint(p['second_moment'])),
                                     (fourth_value, midpoint(p['fourth_moment'])),
                                     (second-mean**2, midpoint(p['variance_X'])),
                                     (fourth_value-second**2, midpoint(p['variance_X_squared']))):
            moment_errors.append(abs(numerical-certified))
    assert max(moment_errors) < 2e-10

    # An imperfect mean operator and every corner of a simultaneous error box.
    # This checks interval propagation against independent operator norms;
    # it is not a simulated proof of statistical coverage.
    mean_op = 1.004*cp.fr.CHI-.002*cp.fr.D
    true_bias = float(np.linalg.norm(mean_op-cp.fr.CHI, 2))
    true_second = float(np.linalg.norm(cal['second_moment'], 2))
    count, bound, delta = 1000000, 2.337, .05
    base = {key: ProbeSummary(count, float(np.trace(rho@mean_op).real),
                             float(np.trace(rho@cal['second_moment']).real), .4, .5)
            for key, rho in zip(('x', 'd', '0'), PROBES)}
    initial = certificate(base, bound, delta)
    for signs in itertools.product((-1, 1), repeat=5):
        observed = dict(base)
        for i, key in enumerate(('x', 'd')):
            observed[key] = replace(observed[key], mean=base[key].mean
                                    + signs[i]*initial['mean_radii'][key])
        for i, key in enumerate(('x', 'd', '0')):
            observed[key] = replace(observed[key], second=base[key].second
                                    + signs[i+2]*initial['moment_radii'][key])
        validated = certificate(observed, bound, delta)
        assert validated['beta']+1e-13 >= true_bias
        assert validated['second_moment_upper']+1e-13 >= true_second
    # The Bernstein range term survives even a zero observed sample variance.
    assert bernstein_radius(0., 2*bound, 100, delta) > 0
    group1 = dict(nominal_moment_max_crosscheck_error=max(moment_errors),
                  operator_error_box_corners=32, true_imperfect_bias=true_bias,
                  scope='Numerical checks of exact moment and confidence-box algebra; not empirical coverage evidence.')

    # The separable endpoint distribution saturates the information bound;
    # independent quadrature checks prior normalization and information.
    b = 4/3
    theta_grid = np.linspace(-.21, .21, 31)
    chord_errors = []
    for theta in theta_grid:
        for x in np.linspace(-b, b, 51):
            f = x*x/(1+theta*x)
            chord = ((b+x)*b*b/(1+theta*b)+(b-x)*b*b/(1-theta*b))/(2*b)
            assert f <= chord+1e-14
        endpoint = .5*b*b/(1+theta*b)+.5*b*b/(1-theta*b)
        chord_errors.append(abs(endpoint-b*b/(1-b*b*theta*theta)))
    epsilon = .01
    a = math.sqrt(math.pi*epsilon/b)
    nodes, weights = np.polynomial.legendre.leggauss(64)
    normalization = float(weights@np.cos(math.pi*nodes/2)**2)
    prior_information = float((math.pi/a)**2*(weights@np.sin(math.pi*nodes/2)**2))
    assert abs(normalization-1) < 1e-14
    assert abs(prior_information-math.pi**2/a**2) < 1e-10
    def necessary_at(radius):
        return (epsilon**-2-math.pi**2/radius**2)*(1-b*b*radius**2)/(b*b)
    optimum = (1/(b*epsilon)-math.pi)**2
    assert abs(necessary_at(a)-optimum) < 1e-10
    assert necessary_at(.99*a) < optimum and necessary_at(1.01*a) < optimum
    assert separable_required(epsilon) == 5164
    group2 = dict(endpoint_information_max_error=max(chord_errors),
                  prior_integral=normalization, prior_information=prior_information,
                  epsilon=epsilon, optimized_necessary_count=optimum,
                  required_integer_copies=separable_required(epsilon),
                  scope='Analytical all-biased bound is proved in SEPARABLE_BENCHMARK.md; grid and quadrature are cross-checks.')

    # ONE planning profile. Population moments are substituted for hypothetical
    # observed summaries. This is not a generated validation record, and does
    # not predict the probability of a real batch passing the certificate.
    counts = {'x': 5500000, 'd': 3000000, '0': 100000}
    profile = {}
    for key, label in (('x', 'xyz'), ('d', '001'), ('0', '000')):
        p = moments['probes'][label]
        # If hypothetical empirical first/second/fourth moments equal this
        # profile, the unbiased sample variances include N/(N-1).
        variance_factor = counts[key]/(counts[key]-1)
        profile[key] = ProbeSummary(counts[key], float(Fraction(p['mean_exact'])),
                                    midpoint(p['second_moment']),
                                    variance_factor*midpoint(p['variance_X']),
                                    variance_factor*midpoint(p['variance_X_squared']))
    validation = certificate(profile, bound, delta)
    budget = preparation_budget(validation, epsilon, pilot_count=10000)
    assert budget['certified_conditional_mse'] <= epsilon**2
    assert budget['ideal_pauli_copies'] == 5625
    assert budget['total_reference_copies'] == 8610000
    for field, comparator in (('strict_reuse_vs_pauli', budget['ideal_pauli_copies']),
                              ('strict_reuse_vs_biased_lower_bound', budget['biased_separable_necessary_copies'])):
        k = budget[field]
        assert k is not None
        charge = budget['total_reference_copies']
        saving = comparator-budget['unknown_copies']
        assert charge < k*saving
        assert k == 1 or charge >= (k-1)*saving
    assert strict_reuse_threshold(100, 50, 100) == 3  # strict, not non-strict
    assert strict_reuse_threshold(100, 100, 100) is None
    assert strict_reuse_threshold(100, 101, 100) is None
    assert strict_reuse_threshold(0, 99, 100) == 1
    assert strict_reuse_threshold(100, 50, 100, .1) == 1
    try:
        preparation_budget(validation, validation['beta'], 0)
    except ValueError:
        pass
    else:
        raise AssertionError('A calibration floor must forbid unattainable target error.')
    vx = midpoint(moments['probes']['xyz']['variance_X'])
    vd = midpoint(moments['probes']['001']['variance_X'])
    nominal_plan = leading_plan(vx, vd, .389156263913, delta)
    # The worst-state ellipse, rather than a sum of independent maxima,
    # yields the plug-in calibration allocation.
    fraction_x = 4*vx/(4*vx+vd)
    assert abs(12*vx/fraction_x-3*vd/(1-fraction_x)) < 1e-13
    assert abs(12*vx/fraction_x-nominal_plan['plugin_mse_coefficient']) < 1e-13
    noisy_plan = leading_plan(.445, .445, .445, delta)
    m, n = budget['total_reference_copies'], budget['unknown_copies']
    # Full-cycle accounting reduces to copy accounting only at equal costs.
    cycle_rows = []
    for comparator in (budget['ideal_pauli_copies'], budget['biased_separable_necessary_copies']):
        assert cycle_reuse_threshold(m, n, comparator, 1., 1., 1.) == strict_reuse_threshold(m, n, comparator)
        cycle_rows.append(dict(comparator_copies=comparator,
                               native_to_comparator_cost_must_be_less_than=comparator/n,
                               reference_cost_ratio_must_be_less_than={
                                   str(k): k*(comparator-n)/m for k in (1, 10, 100, 1000)}))
    assert cycle_reuse_threshold(100, 50, 100, 1., 2., 1.) is None
    assert cycle_reuse_threshold(100, 50, 100, 0., 3., 1.) is None
    # Check the capacity boundary for a complete sequential stable epoch.
    required_k = budget['strict_reuse_vs_biased_lower_bound']
    required_copies = m+required_k*n
    assert required_copies == 61497869
    assert stable_reuse_capacity(m, n, required_copies, 1., 1.) == required_k
    assert stable_reuse_capacity(m, n, required_copies-1, 1., 1.) == required_k-1
    assert stable_reuse_capacity(m, n, m-1, 1., 1.) == 0
    assert stable_reuse_capacity(10, 5, 30., 1., 1., 5., 5.) == 1
    drifted = preparation_budget(operator_allowances(validation, 1e-4), epsilon, 10000)
    assert drifted['unknown_copies'] == 4476
    assert drifted['strict_reuse_vs_biased_lower_bound'] == 12515
    slack = epsilon**2-validation['second_moment_upper']/n-validation['beta']**2
    assert slack > 0
    mean_drift_limit = math.sqrt(validation['beta']**2+slack)-validation['beta']

    # Construct two valid, covariant three-outcome input POVMs with different
    # response strengths. Matching their mixture calibrates every input;
    # correlating input states with the channel labels invalidates that result.
    p, q = np.array([.3, .7]), np.array([.6, .4])
    strengths = np.array([.5, 1.])
    mix_score = 1/float(p@strengths)
    channels = [np.array([a*(cp.fr.I8+cp.fr.CHI)/2,
                          a*(cp.fr.I8-cp.fr.CHI)/2, (1-a)*cp.fr.I8]) for a in strengths]
    scores = np.array([mix_score, -mix_score, 0.])
    for channel in channels:
        assert np.linalg.norm(channel.sum(axis=0)-cp.fr.I8, 2) < 1e-14
        assert min(np.linalg.eigvalsh(effect)[0] for effect in channel) >= -1e-14
    means = np.array([np.einsum('p,pij->ij', scores, channel) for channel in channels])
    seconds = np.array([np.einsum('p,pij->ij', scores**2, channel) for channel in channels])
    matched_mean = np.einsum('i,ijk->jk', p, means)
    assert np.linalg.norm(matched_mean-cp.fr.CHI, 2) < 1e-14
    tv = float(np.sum(abs(p-q))/2)
    mismatch_mean = float(np.linalg.norm(np.einsum('i,ijk->jk', q-p, means), 2))
    mismatch_second = float(np.linalg.norm(np.einsum('i,ijk->jk', q-p, seconds), 2))
    assert 0 < mismatch_mean <= 2*mix_score*tv
    assert 0 < mismatch_second <= mix_score**2*tv
    values, vectors = np.linalg.eigh(cp.fr.CHI)
    rho_plus = np.outer(vectors[:, -1], vectors[:, -1].conj())
    rho_minus = np.outer(vectors[:, 0], vectors[:, 0].conj())
    correlated_output = sum(pi*np.trace(ai@rho).real
                            for pi, ai, rho in zip(p, means, (rho_plus, rho_minus)))
    average_target = np.trace(cp.fr.CHI@(p[0]*rho_plus+p[1]*rho_minus)).real
    assert abs(correlated_output-average_target) > .05
    # A common random setting/noise component is enough to invalidate the
    # independent-pixel variance formula, even for independent local noise.
    pixels, frames, correlation = 100, 20, .02
    covariance = (1-correlation)*np.eye(pixels)+correlation*np.ones((pixels, pixels))
    frame_variance = float(np.ones(pixels)@covariance@np.ones(pixels)/pixels**2)
    factor = 1+(pixels-1)*correlation
    assert abs(frame_variance/frames-factor/(pixels*frames)) < 1e-15
    use_case = dict(full_cycle_requirements=cycle_rows,
                    independent_samples_to_equal_cost_biased_break_even=required_copies,
                    baseline_mse_slack=slack, separate_mean_drift_limit=mean_drift_limit,
                    separate_second_drift_limit=n*slack,
                    illustrative_mean_drift_1e_minus4=drifted,
                    mixture_crosscheck=dict(total_variation=tv,
                        mean_mismatch=mismatch_mean, second_mismatch=mismatch_second,
                        state_channel_correlation_error=abs(correlated_output-average_target)),
                    shared_frame_variance_inflation=factor,
                    scope='Analytical requirements and constructed cross-checks, not measured apparatus '
                          'costs, rates, spatial noise, or stability. The class-wide time comparison '
                          'additionally needs a valid comparator cost lower bound.')
    group3 = dict(planning_profile={k: vars(v) for k, v in profile.items()},
                  score_bound=bound, confidence_failure=delta, validation=validation,
                  budget=budget, nominal_leading_plan=nominal_plan,
                  conservative_noisy_leading_plan=noisy_plan,
                  repeated_state_scan=use_case,
                  scope='HYPOTHETICAL nominal summaries, not experimental data, guaranteed-pass counts, '
                        'or an unconditional finite-risk advantage. Pilot count is an explicit illustrative '
                        'charge, not a claim that 10000 preparations suffice to select this score. '
                        'No preparation error or channel drift is assumed in this profile.')
    return dict(status='PASS', scientific_groups=3,
                groups=dict(moment_and_confidence_algebra=group1,
                            biased_separable_bound=group2,
                            preparation_cost_accounting=group3))


def main():
    if not __debug__:
        raise RuntimeError('Run with Python assertions enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--certificate', type=Path)
    args = parser.parse_args()
    target = args.output.resolve()
    if target.is_relative_to(Path(__file__).resolve().parents[1]):
        parser.error('Use a new report outside the repository.')
    certificate_path = args.certificate or target.parent/'complete-parity-certificate.json'
    report = run(json.loads(certificate_path.read_text()))
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS: three calibration-cost cross-check groups.')


if __name__ == '__main__':
    main()
