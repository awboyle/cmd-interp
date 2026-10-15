"""Empirical reference sequences (de-extincted cluster CMD sequences of known age)."""

from __future__ import annotations

import re
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

# Other accepted names for the bundled references (matched ignoring case, spaces, '_', '-' and '.').
REFERENCE_ALIASES = {
    "ubc17a": "UBC_17a",
    "betapic": "beta_Pic", "betapictoris": "beta_Pic", "bpmg": "beta_Pic", "βpic": "beta_Pic",
    "ic2602": "IC_2602",
    "melotte20": "Melotte_20", "alphaper": "Melotte_20", "alphapersei": "Melotte_20", "αper": "Melotte_20",
    "melotte22": "Melotte_22", "pleiades": "Melotte_22", "m45": "Melotte_22",
    "melotte111": "Melotte_111", "comaber": "Melotte_111", "comaberenices": "Melotte_111",
}

COLOR_BIN_WIDTH = 0.5        # mag; colour window of the local member scatter
SCATTER_RESIDUAL_CUT = 0.3   # mag; members farther than this from their sequence are not used for the scatter


def bundled_reference_dir():
    """Directory of the reference library shipped with the package."""
    return Path(str(resources.files("cmd_interp") / "data" / "references"))


def reference_name(name):
    """Canonical name of a bundled reference ('Pleiades', 'alpha Per', 'beta_Pic', ... -> file name);
    names that are not bundled references are returned unchanged."""
    key = re.sub(r"[\s_\-.]", "", str(name)).lower()
    if key in REFERENCE_ALIASES:
        return REFERENCE_ALIASES[key]
    for canonical in REFERENCE_AGES:
        if re.sub(r"[\s_\-.]", "", canonical).lower() == key:
            return canonical
    return str(name)


def resolve_reference_ages(reference_ages=None, exclude=None, use_default_ages=True):
    """The reference ages to use (Myr), youngest first.

    reference_ages: mapping of cluster -> age in Myr that replaces the adopted age of those
        clusters (e.g. ``{"beta_Pic": 21.0}`` or ``{"Pleiades": 112}``); all other clusters
        keep their default ages.  With use_default_ages=False it is the complete list instead
        (e.g. for a custom reference library).
    exclude: cluster names to leave out.
    """
    ages = dict(REFERENCE_AGES) if use_default_ages else {}
    for name, age in (reference_ages or {}).items():
        age = float(age)
        if not np.isfinite(age) or age <= 0:
            raise ValueError(f"Reference age of {name} must be a positive number of Myr, got {age}")
        ages[reference_name(name)] = age
    for name in exclude or ():
        canonical = reference_name(name)
        if canonical not in ages:
            raise ValueError(f"Cannot exclude unknown reference {name!r}; references: {', '.join(ages)}")
        ages.pop(canonical)
    if len(ages) < 2:
        raise ValueError("At least two reference clusters are needed.")
    values = sorted(ages.values())
    if np.any(np.diff(values) <= 0):
        raise ValueError(f"Reference ages must all differ: {ages}")
    return dict(sorted(ages.items(), key=lambda item: item[1]))


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
    exclude=None,
    use_default_ages=True,
    color_bin_width=COLOR_BIN_WIDTH,
    scatter_residual_cut=SCATTER_RESIDUAL_CUT,
    include_member_ids=False,
):
    """Load the reference sequences and their member residuals.

    reference_dir: directory holding ``<name>_sequence.feather`` (columns BpmRp0, GMag0)
        and ``<name>_members.feather`` (columns BpmRp0, phot_g_mean_mag_0, and GaiaDR3 if
        member IDs are wanted).  Defaults to the bundled library.
    reference_ages: ages (Myr) that replace the defaults of the named clusters, e.g.
        ``{"beta_Pic": 21.0}``; see resolve_reference_ages().  Aliases such as "Pleiades"
        or "alpha Per" are accepted.
    exclude: clusters to leave out (at least two must remain).
    use_default_ages: False to take ``reference_ages`` as the complete list (custom library).

    Returns a list of dicts sorted by age (youngest first), the input of infer_ages().
    """
    reference_dir = bundled_reference_dir() if reference_dir is None else Path(reference_dir).expanduser()
    ages = resolve_reference_ages(reference_ages, exclude, use_default_ages)

    references = []
    for cluster, age_myr in ages.items():
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
