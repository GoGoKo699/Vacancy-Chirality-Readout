"""Finite calibration accounting; floating calculations, not interval proofs.

Confidence claims require fixed counts and scores frozen independently of the
validation data. Model moments used for planning are not empirical data.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ProbeSummary:
    count: int
    mean: float
    second: float
    variance: float
    squared_variance: float


def bernstein_radius(variance, width, count, delta, statistics=5):
    """Two-sided Maurer--Pontil radius with a union over statistics.

variance is the UNBIASED sample variance, not a fitted population variance.
Correlated statistics from one batch are allowed by the union bound.
"""
    if (not isinstance(count, int) or count < 2 or statistics < 1
            or not 0 < delta < 1 or variance < 0 or width < 0
            or not all(map(math.isfinite, (variance, width, delta)))):
        raise ValueError("Invalid variance, range, count, or confidence budget.")
    log_factor = math.log(4*statistics/delta)
    return (math.sqrt(2*variance*log_factor/count)
            + 7*width*log_factor/(3*(count-1)))


def certificate(probes, score_bound, delta=.05):
    """Certificate from the two means and three second moments of x,d,0.

Returns confidence-event bounds, not an unconditional risk or a pass guarantee.
No preparation or drift allowance is implicitly included.
"""
    if not math.isfinite(score_bound) or score_bound <= 0:
        raise ValueError("A positive finite bound on the frozen score is required.")
    if set(probes) != {"x", "d", "0"}:
        raise ValueError("Provide exactly the three reference summaries x,d,0.")
    mean_radii, moment_radii = {}, {}
    for name, p in probes.items():
        if not all(map(math.isfinite, (p.mean, p.second, p.variance, p.squared_variance))):
            raise ValueError("Non-finite probe summary.")
        if abs(p.mean) > score_bound or not 0 <= p.second <= score_bound**2:
            raise ValueError("Probe summary is outside the frozen score range.")
        # Reject impossible negative variances; callers must provide empirical
        # unbiased variances, or explicitly label their output as a forecast.
        moment_radii[name] = bernstein_radius(
            p.squared_variance, score_bound**2, p.count, delta)
        if name != "0":
            mean_radii[name] = bernstein_radius(p.variance, 2*score_bound, p.count, delta)
    bc = abs(8*probes["x"].mean-1)+8*mean_radii["x"]
    bd = abs(probes["d"].mean)+mean_radii["d"]
    beta = math.sqrt(3*bc*bc/16+3*bd*bd)
    qx, qd, q0 = (probes[k].second for k in ("x", "d", "0"))
    nx, nd, n0 = (moment_radii[k] for k in ("x", "d", "0"))
    upper_blocks = [q0+n0, 3*qd-2*qx+3*nd+2*nx,
                    6*qx-3*qd-2*q0+6*nx+3*nd+2*n0]
    upper = min(score_bound**2, max(0., *upper_blocks))
    return dict(beta=beta, second_moment_upper=upper,
                mean_radii=mean_radii, moment_radii=moment_radii,
                block_upper_s_o_e=upper_blocks,
                reference_validation=sum(p.count for p in probes.values()),
                confidence_failure=delta)


def separable_required(epsilon):
    """Necessary unknown-copy count for uniform MSE, including biased estimates."""
    if not 0 < epsilon < 1/(16*math.pi):
        raise ValueError("The optimized prior requires 0 < epsilon < 1/(16*pi).")
    return math.ceil((3/(4*epsilon)-math.pi)**2)


def strict_reuse_threshold(reference_count, target_count, comparator_count, weight=1.):
    """Least integer K with target_count + weight*reference_count/K < comparator.

None means there is no finite strict saving for these certified counts.
"""
    if (not isinstance(reference_count, int) or reference_count < 0
            or not isinstance(target_count, int) or target_count < 1
            or not isinstance(comparator_count, int) or comparator_count < 1
            or not math.isfinite(weight) or weight < 0):
        raise ValueError("Counts must be nonnegative integers and costs nonnegative.")
    saving = comparator_count-target_count
    if saving <= 0:
        return None
    return math.floor(weight*reference_count/saving)+1


def preparation_budget(validation, epsilon, pilot_count, reference_weight=1.):
    """Count pilot, all reference trials, and unknowns; retain calibration bias."""
    if not isinstance(pilot_count, int) or pilot_count < 0:
        raise ValueError("Count every pilot preparation, using a nonnegative integer.")
    beta, upper = validation["beta"], validation["second_moment_upper"]
    if not math.isfinite(epsilon) or not 0 <= beta < epsilon or upper < 0:
        raise ValueError("A positive error budget beyond the calibration bias is required.")
    count = max(1, math.ceil(upper/(epsilon**2-beta**2)))
    references = pilot_count+validation["reference_validation"]
    pauli = math.ceil((9/16)/epsilon**2)
    separable = separable_required(epsilon)
    return dict(epsilon=epsilon, unknown_copies=count, pilot_copies=pilot_count,
                validation_copies=validation["reference_validation"],
                total_reference_copies=references, reference_weight=reference_weight,
                certified_conditional_mse=upper/count+beta**2,
                ideal_pauli_copies=pauli, biased_separable_necessary_copies=separable,
                strict_reuse_vs_pauli=strict_reuse_threshold(references, count, pauli, reference_weight),
                strict_reuse_vs_biased_lower_bound=strict_reuse_threshold(
                    references, count, separable, reference_weight))


def leading_plan(vx, vd, worst_second, delta=.05):
    """Separate confidence-radius and regular plug-in asymptotic forecasts."""
    if min(vx, vd) <= 0 or not 0 < worst_second < 9/16 or not 0 < delta < 1:
        raise ValueError("Positive variances, confidence, and a separable margin are required.")
    confidence_c = 6*math.log(20/delta)*(2*math.sqrt(vx)+math.sqrt(vd))**2
    plugin_c = 12*vx+3*vd
    margin = (.75-math.sqrt(worst_second))**2
    return dict(confidence_radius_coefficient=confidence_c,
                confidence_allocation_x_over_d=2*math.sqrt(vx/vd),
                confidence_leading_reuse_threshold=confidence_c/margin,
                plugin_mse_coefficient=plugin_c,
                plugin_allocation_x_over_d=4*vx/vd,
                plugin_leading_reuse_threshold=plugin_c/margin,
                scope="Sufficient leading planning only; excludes pilot/moment-only overhead, "
                      "finite remainders and acceptance margins; not a necessary reuse count.")
