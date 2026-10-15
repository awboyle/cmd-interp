"""Empirical reference sequences (de-extincted cluster CMD sequences of known age)."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd

# Adopted ages (Myr) of the bundled reference clusters, youngest first.
REFERENCE_AGES = {
    "UBC_17a": 4.6,
    "beta_Pic": 24.3,      # beta Pic moving group, LDB age (Binks & Jeffries 2014)
    "IC_2602": 52.5,
    "Melotte_20": 79.0,    # alpha Persei
    "Melotte_22": 127.4,   # Pleiades
    "Melotte_111": 675.0,  # Coma Berenices: terminal sequence, reaches ~3300 K
}

COLOR_BIN_WIDTH = 0.5        # mag; colour window of the local member scatter
SCATTER_RESIDUAL_CUT = 0.3   # mag; members farther than this from their sequence are not used for the scatter


def bundled_reference_dir():
    """Directory of the reference library shipped with the package."""
    return Path(str(resources.files("cmd_interp") / "data" / "references"))


def _sequence_points(frame, color_col="BpmRp0", mag_col="GMag0"):
    return (
        frame[[color_col, mag_col]]
        .dropna()
        .groupby(color_col, as_index=False)[mag_col]
        .mean()
        .sort_values(color_col)
    )


def load_references(
    reference_dir=None,
    reference_ages=None,
    *,
    color_bin_width=COLOR_BIN_WIDTH,
    scatter_residual_cut=SCATTER_RESIDUAL_CUT,
    include_member_ids=False,
):
    """Load the reference sequences and their member residuals.

    reference_dir: directory holding ``<name>_sequence.feather`` (columns BpmRp0, GMag0)
        and ``<name>_members.feather`` (columns BpmRp0, phot_g_mean_mag_0, and GaiaDR3 if
        member IDs are wanted).  Defaults to the bundled library.
    reference_ages: mapping of cluster name -> adopted age in Myr.  Defaults to
        REFERENCE_AGES.  Pass a subset to leave clusters out (at least two are needed).

    Returns a list of dicts sorted by age (youngest first), the input of infer_ages().
    """
    reference_dir = bundled_reference_dir() if reference_dir is None else Path(reference_dir).expanduser()
    reference_ages = REFERENCE_AGES if reference_ages is None else dict(reference_ages)
    if len(reference_ages) < 2:
        raise ValueError("At least two reference clusters are needed.")

    references = []
    for cluster, age_myr in sorted(reference_ages.items(), key=lambda item: item[1]):
        sequence = pd.read_feather(reference_dir / f"{cluster}_sequence.feather", columns=["BpmRp0", "GMag0"])
        points = _sequence_points(sequence)
        members = pd.read_feather(
            reference_dir / f"{cluster}_members.feather",
            columns=None if include_member_ids else ["BpmRp0", "phot_g_mean_mag_0"],
        ).dropna(subset=["BpmRp0", "phot_g_mean_mag_0"])
        if include_member_ids and ("GaiaDR3" not in members or not pd.api.types.is_integer_dtype(members["GaiaDR3"])):
            raise ValueError(f"{cluster} needs integer GaiaDR3 member IDs.")
        members = members.loc[members.BpmRp0.between(points.BpmRp0.min(), points.BpmRp0.max())].copy()
        fitted_mag = np.interp(members.BpmRp0, points.BpmRp0, points.GMag0)
        members["residual_g"] = members.phot_g_mean_mag_0 - fitted_mag
        if scatter_residual_cut is not None:
            members = members.loc[members.residual_g.abs().le(scatter_residual_cut)]
        members = members.sort_values("BpmRp0")
        reference = {
            "cluster": cluster,
            "age_myr": float(age_myr),
            "color": points.BpmRp0.to_numpy(dtype=float),
            "mag": points.GMag0.to_numpy(dtype=float),
            "member_color": members.BpmRp0.to_numpy(dtype=float),
            "member_mag": members.phot_g_mean_mag_0.to_numpy(dtype=float),
            "member_residual": members.residual_g.to_numpy(dtype=float),
            "half_color_width": color_bin_width / 2.0,
        }
        if include_member_ids:
            reference["member_source_ids"] = members["GaiaDR3"].to_numpy()
        references.append(reference)
    return references
