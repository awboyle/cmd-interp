"""Age grid and posterior quantiles."""

from __future__ import annotations

import numpy as np

AGE_FLOOR_MYR = 0.0
AGE_CEILING_MYR = 4000.0


def age_grid(reference_ages_myr, age_step_myr=1.0, floor_to=AGE_FLOOR_MYR, ceiling_to=AGE_CEILING_MYR, coarse_step_myr=5.0):
    """Segmented age grid (Myr).

    ``age_step_myr`` nodes between consecutive reference ages, the same step below the youngest
    reference down to ``floor_to``, and ``coarse_step_myr`` nodes above the oldest reference up to
    ``ceiling_to``.  Outside the reference ages nothing constrains the likelihood, so a star's
    posterior is held flat there at its value at the nearest fitted node.
    """
    reference_ages_myr = np.asarray(reference_ages_myr, dtype=float)
    age_segments = []
    floor_nodes = np.arange(reference_ages_myr[0] - age_step_myr, floor_to - 1e-9, -age_step_myr)[::-1]
    if len(floor_nodes):
        age_segments.append(floor_nodes.astype(float))
    for segment_index, (segment_min, segment_max) in enumerate(zip(reference_ages_myr[:-1], reference_ages_myr[1:])):
        segment = np.arange(segment_min, segment_max, age_step_myr, dtype=float)
        if not len(segment) or not np.isclose(segment[-1], segment_max):
            segment = np.append(segment, segment_max)
        else:
            segment[-1] = segment_max
        if segment_index:
            segment = segment[1:]
        age_segments.append(segment)
    ceiling_nodes = np.arange(reference_ages_myr[-1] + coarse_step_myr, ceiling_to + 1e-9, coarse_step_myr)
    if len(ceiling_nodes):
        age_segments.append(ceiling_nodes.astype(float))
    return np.concatenate(age_segments)


def posterior_quantile(age_grid_myr, density, probability):
    """Trapezoid-rule quantile of one posterior per column of ``density`` (n_age, n_star)."""
    cdf = np.zeros_like(density)
    cdf[1:] = np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(age_grid_myr)[:, None], axis=0)
    cdf /= cdf[-1]
    upper = np.argmax(cdf >= probability, axis=0)
    upper = np.maximum(upper, 1)
    lower = upper - 1
    column = np.arange(density.shape[1])
    cdf_lower = cdf[lower, column]
    cdf_upper = cdf[upper, column]
    fraction = np.divide(
        probability - cdf_lower,
        cdf_upper - cdf_lower,
        out=np.zeros_like(cdf_lower),
        where=cdf_upper > cdf_lower,
    )
    return age_grid_myr[lower] + fraction * (age_grid_myr[upper] - age_grid_myr[lower])


def node_widths(age_grid_myr):
    """Trapezoid weight (Myr) of every grid node."""
    age_grid_myr = np.asarray(age_grid_myr, dtype=float)
    widths = np.zeros_like(age_grid_myr)
    step = np.diff(age_grid_myr)
    widths[:-1] += 0.5 * step
    widths[1:] += 0.5 * step
    return widths


def posterior_density(posterior, age_grid_myr):
    """Convert unit-sum node weights (rows = stars) into probability densities per Myr.

    The grid mixes 1 Myr and 5 Myr steps, so node weights must not be plotted or averaged
    as densities directly.  Each returned row integrates to one with the trapezoid rule.
    """
    posterior = np.atleast_2d(np.asarray(posterior, dtype=float))
    age_grid_myr = np.asarray(age_grid_myr, dtype=float)
    norm = np.trapezoid(posterior, age_grid_myr, axis=1) if hasattr(np, "trapezoid") else np.trapz(posterior, age_grid_myr, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return posterior / norm[:, None]
