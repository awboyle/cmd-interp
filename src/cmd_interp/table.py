"""Ages for a table of stars."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .core import DEFAULT_SCATTER_MODE, DEFAULT_SYSTEMATIC_MAG, infer_ages
from .references import load_references

DEFAULT_ABS_G_ERR = 0.03   # mag; used where the table has no (valid) per-star absolute-G error

# output columns
AGE_COLUMN = "cmd_age_myr"                  # adopted age: MAP, or the reference age at a limit
LOWER_COLUMN = "cmd_age_lower_myr"          # 16th percentile (the youngest reference age for 'old' limits)
UPPER_COLUMN = "cmd_age_upper_myr"          # 84th percentile (the youngest reference age for 'young' limits)
MEDIAN_COLUMN = "cmd_age_median_myr"
P025_COLUMN = "cmd_age_p2p5_myr"
P975_COLUMN = "cmd_age_p97p5_myr"
LIMIT_COLUMN = "cmd_age_limit"              # 'young' (< youngest reference), 'old' (> oldest reference) or ''
NREF_COLUMN = "cmd_age_n_references"
SIGMA_COLUMN = "cmd_age_scatter_mag"
MH_COLUMN = "cmd_age_mh"
CLIPPED_COLUMN = "cmd_age_metallicity_clipped"
FITTED_COLUMN = "cmd_age_fitted"
COLOR_RANGE_COLUMN = "cmd_age_color_out_of_range"   # (BP-RP)_0 outside the fitted range (default: redder than 3.0)


def read_table(path):
    """Read a .csv, .feather, .parquet or .fits table into a DataFrame."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".feather":
        return pd.read_feather(path)
    if suffix in (".parquet", ".pq"):
        return pd.read_parquet(path)
    if suffix in (".csv", ".txt"):
        return pd.read_csv(path)
    if suffix in (".fits", ".fit"):
        from astropy.table import Table  # optional dependency
        return Table.read(path).to_pandas()
    raise ValueError(f"Unsupported table format: {path}")


def write_table(frame, path):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".feather":
        frame.reset_index(drop=True).to_feather(path)
    elif suffix in (".parquet", ".pq"):
        frame.to_parquet(path, index=False)
    elif suffix in (".csv", ".txt"):
        frame.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported table format: {path}")


def ages_for_table(
    data,
    *,
    color_col="BpmRp0",
    abs_g_col="phot_g_mean_mag_0",
    color_err_col="BpmRp0_err",
    abs_g_err_col="phot_g_mean_mag_0_err",
    abs_g_err=DEFAULT_ABS_G_ERR,
    mh=None,
    mh_col=None,
    fit_mask=None,
    systematic_mag=DEFAULT_SYSTEMATIC_MAG,
    scatter_mode=DEFAULT_SCATTER_MODE,
    references=None,
    reference_ages=None,
    exclude_references=None,
    return_posterior=True,
    posterior_out=None,
    **kwargs,
):
    """Append CMD ages to a table of stars.

    data: DataFrame or path (.csv/.feather/.parquet/.fits) holding extinction-corrected
        (BP-RP)_0 and absolute G, their errors, and optionally [M/H].
    abs_g_err: fallback where ``abs_g_err_col`` is missing or not finite.
    mh / mh_col: a single [M/H] for every star, or the name of a per-star [M/H] column
        (NaN rows are treated as solar).  Default: solar.
    fit_mask: boolean array / column name selecting the stars to fit (e.g. a pre-main-sequence
        selection); other rows get NaN ages.  Default: every row with finite photometry.
        Main-sequence stars get broad, uninformative posteriors from this method.
    reference_ages / exclude_references: new ages (Myr) for any reference cluster, e.g.
        ``{"beta_Pic": 21.0}``, and clusters to leave out; see infer_ages().
    Other keyword arguments are passed to infer_ages() (e.g. max_color, binary_fraction).

    Returns (table, result): a copy of the table with the cmd_age_* columns appended, and the
    infer_ages() output for the fitted rows (``result['row']`` gives their positions in the
    table; ``result['posterior']`` is (n_fitted, n_age) if return_posterior).
    """
    frame = read_table(data) if isinstance(data, (str, Path)) else data.copy()
    n_row = len(frame)
    missing = [col for col in (color_col, abs_g_col, color_err_col) if col not in frame]
    if missing:
        raise KeyError(f"Missing column(s): {', '.join(missing)}")

    color = frame[color_col].to_numpy(dtype=float)
    mag = frame[abs_g_col].to_numpy(dtype=float)
    color_err = frame[color_err_col].to_numpy(dtype=float)
    if abs_g_err_col in frame:
        g_err = frame[abs_g_err_col].to_numpy(dtype=float)
        g_err = np.where(np.isfinite(g_err) & (g_err >= 0), g_err, abs_g_err)
    else:
        g_err = np.full(n_row, float(abs_g_err))
    if mh_col is not None:
        if mh is not None:
            raise ValueError("Pass either mh or mh_col, not both")
        mh_values = frame[mh_col].to_numpy(dtype=float)
    else:
        mh_values = np.full(n_row, 0.0 if mh is None else float(mh))

    if fit_mask is None:
        selected = np.ones(n_row, dtype=bool)
    elif isinstance(fit_mask, str):
        selected = frame[fit_mask].fillna(False).to_numpy(dtype=bool)
    else:
        selected = np.asarray(fit_mask, dtype=bool)
        if selected.shape != (n_row,):
            raise ValueError("fit_mask must have one entry per row")
    finite = np.isfinite(color) & np.isfinite(mag) & np.isfinite(color_err) & (color_err >= 0)
    rows = np.flatnonzero(selected & finite)

    if references is None:
        references = load_references(reference_ages=reference_ages, exclude=exclude_references)
    elif reference_ages is not None or exclude_references is not None:
        raise ValueError("Pass reference_ages / exclude_references, or a preloaded references list, not both")
    result = infer_ages(
        color[rows], mag[rows], color_err[rows], g_err[rows],
        mh=mh_values[rows] if np.any(np.nan_to_num(mh_values[rows]) != 0) else None,
        systematic_mag=systematic_mag, scatter_mode=scatter_mode, references=references,
        return_posterior=return_posterior, posterior_out=posterior_out, **kwargs,
    )
    result["row"] = rows

    age = np.full(n_row, np.nan)
    lower = np.full(n_row, np.nan)
    upper = np.full(n_row, np.nan)
    median = np.full(n_row, np.nan)
    p025 = np.full(n_row, np.nan)
    p975 = np.full(n_row, np.nan)
    limit = np.full(n_row, "", dtype=object)
    n_ref = np.zeros(n_row, dtype=int)
    sigma = np.full(n_row, np.nan)
    clipped = np.zeros(n_row, dtype=bool)
    fitted = np.zeros(n_row, dtype=bool)
    out_of_range = np.zeros(n_row, dtype=bool)

    young = result["limit"] == "young"
    old = result["limit"] == "old"
    success = result["has_posterior"] & (result["limit"] == "")
    age[rows[young]] = result["youngest_age"][young]
    upper[rows[young]] = result["youngest_age"][young]
    age[rows[old]] = result["oldest_age"][old]
    lower[rows[old]] = result["oldest_age"][old]
    age[rows[success]] = result["age_map_myr"][success]
    lower[rows[success]] = result["age_myr_16"][success]
    upper[rows[success]] = result["age_myr_84"][success]
    median[rows[success]] = result["age_myr_50"][success]
    p025[rows[success]] = result["age_myr_025"][success]
    p975[rows[success]] = result["age_myr_975"][success]
    limit[rows] = result["limit"]
    n_ref[rows] = result["n_references"]
    sigma[rows] = result["scatter_sigma_mag"]
    clipped[rows] = result["metallicity_clipped"]
    fitted[rows] = result["has_posterior"]
    out_of_range[rows] = result["color_out_of_range"]

    frame[AGE_COLUMN] = age
    frame[LOWER_COLUMN] = lower
    frame[UPPER_COLUMN] = upper
    frame[MEDIAN_COLUMN] = median
    frame[P025_COLUMN] = p025
    frame[P975_COLUMN] = p975
    frame[LIMIT_COLUMN] = limit
    frame[NREF_COLUMN] = n_ref
    frame[SIGMA_COLUMN] = sigma
    frame[MH_COLUMN] = np.where(np.isfinite(mh_values), mh_values, 0.0)
    frame[CLIPPED_COLUMN] = clipped
    frame[FITTED_COLUMN] = fitted
    frame[COLOR_RANGE_COLUMN] = out_of_range
    return frame, result
