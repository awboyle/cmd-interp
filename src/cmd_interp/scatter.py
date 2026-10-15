"""Local CMD scatter of the reference sequences."""

from __future__ import annotations

import numpy as np
from scipy.stats import median_abs_deviation

# A sliding-window scatter needs at least this many members, otherwise it is filled along colour.
MIN_MEMBERS_SCATTER = 8


def exact_local_scatter(query_color, reference):
    """Normal-scaled MAD of the member residuals within half_color_width of each query colour."""
    member_color = reference["member_color"]
    residual = reference["member_residual"]
    half_width = reference["half_color_width"]
    query_color = np.asarray(query_color, dtype=float)
    left = np.searchsorted(member_color, query_color - half_width, side="left")
    right = np.searchsorted(member_color, query_color + half_width, side="right")
    pair_code = left.astype(np.int64) * (len(member_color) + 1) + right
    unique_code, inverse = np.unique(pair_code, return_inverse=True)
    unique_scatter = np.full(len(unique_code), np.nan, dtype=float)
    for index, code in enumerate(unique_code):
        first = int(code // (len(member_color) + 1))
        last = int(code % (len(member_color) + 1))
        if last - first < 2:
            continue
        scatter = median_abs_deviation(residual[first:last], scale="normal", nan_policy="omit")
        if np.isfinite(scatter) and scatter > 0:
            unique_scatter[index] = scatter
    return unique_scatter[inverse.reshape(-1)]


def window_counts(query_color, reference):
    member_color = reference["member_color"]
    half_width = reference["half_color_width"]
    return (np.searchsorted(member_color, query_color + half_width, side="right")
            - np.searchsorted(member_color, query_color - half_width, side="left"))


def fill_scatter(query_color, reference, scatter, min_members=None):
    """Replace NaN scatters, and windows holding fewer than ``min_members`` members, by the
    reference's scatter at the nearest colour where it is defined from at least that many."""
    min_members = MIN_MEMBERS_SCATTER if min_members is None else min_members
    query_color = np.asarray(query_color, dtype=float)
    scatter = np.asarray(scatter, dtype=float).copy()
    scatter[window_counts(query_color, reference) < min_members] = np.nan
    bad = ~np.isfinite(scatter)
    if not bad.any():
        return scatter
    grid = np.arange(reference["color"][0], reference["color"][-1] + 0.025, 0.05)
    table = exact_local_scatter(grid, reference)
    table[window_counts(grid, reference) < min_members] = np.nan
    good = np.isfinite(table)
    if not good.any():
        return scatter
    scatter[bad] = np.interp(query_color[bad], grid[good], table[good])
    return scatter


def nearest_two_scatter(abs_g, reference_mag, reference_scatter):
    """Per-star scatter: mean of the two reference sequences nearest in magnitude at the star's
    colour.  Above the youngest sequence only the youngest's scatter is used; below the oldest
    only the oldest's.  reference_mag / reference_scatter: (n_ref, n_star), youngest first."""
    n_ref, n_star = reference_mag.shape
    out = np.full(n_star, np.nan)
    for j in range(n_star):
        mags, scat = reference_mag[:, j], reference_scatter[:, j]
        ok = np.isfinite(scat)
        if not ok.any():
            continue
        if abs_g[j] <= mags[0]:
            cand = [0] if ok[0] else [int(np.flatnonzero(ok)[0])]
        elif abs_g[j] >= mags[-1]:
            cand = [n_ref - 1] if ok[-1] else [int(np.flatnonzero(ok)[-1])]
        else:
            order = np.argsort(np.abs(mags - abs_g[j]))
            cand = [int(k) for k in order if ok[k]][:2]
        out[j] = float(np.mean(scat[cand]))
    return out
