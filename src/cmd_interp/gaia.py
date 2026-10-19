"""Ages straight from Gaia DR3 source IDs: catalogue query -> de-extinction and uncertainties -> CMD ages.

    table, result = cmd_interp.ages_from_gaia([6794047652729201024, ...], a0=0.0)
    star = cmd_interp.star_age_from_gaia(6794047652729201024, a0=0.0)   # AU Mic

The photometry and astrometry come from gaiadr3.gaia_source through the ESA Gaia archive
(anonymous queries; requires astroquery: pip install "cmd-interp[gaia]").  A_0 must be
given explicitly: a value (or one per star), or "dustmap" for the EXPLORE 3D dust maps.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .photometry import deredden_gaia
from .table import ages_for_table

GAIA_QUERY_COLUMNS = [
    "source_id", "ra", "dec", "pmra", "pmdec", "parallax", "parallax_error", "ruwe",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
    "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_bp_mean_flux", "phot_bp_mean_flux_error",
    "phot_rp_mean_flux", "phot_rp_mean_flux_error", "azero_gspphot",
]
CHUNK = 2000   # synchronous archive queries return at most 2000 rows


def query_gaia_dr3(source_ids, *, chunk_size=CHUNK, credentials_file=None):
    """gaia_source rows for the given DR3 source IDs, as a DataFrame with a dr3_source_id column.

    Anonymous by default.  credentials_file: optional astroquery credentials file (user name and
    password on two lines) to log in, e.g. for higher archive quotas.  Never put credentials in code.
    """
    try:
        from astroquery.gaia import Gaia
    except ImportError as error:  # pragma: no cover
        raise ImportError('Querying Gaia needs astroquery: pip install "cmd-interp[gaia]"') from error
    ids = np.unique(np.asarray(source_ids, dtype=np.int64))
    if credentials_file is not None:
        Gaia.login(credentials_file=str(credentials_file))
    frames = []
    try:
        for start in range(0, len(ids), chunk_size):
            block = ", ".join(str(int(i)) for i in ids[start:start + chunk_size])
            query = f"SELECT {', '.join(GAIA_QUERY_COLUMNS)} FROM gaiadr3.gaia_source WHERE source_id IN ({block})"
            frames.append(Gaia.launch_job(query, verbose=False).get_results().to_pandas())
    finally:
        if credentials_file is not None:
            Gaia.logout()
    rows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=GAIA_QUERY_COLUMNS)
    rows.columns = [c.lower() for c in rows.columns]
    return rows.rename(columns={"source_id": "dr3_source_id"}).astype({"dr3_source_id": np.int64})


def _per_star(value, n, name):
    """Scalar, string (column / 'dustmap') or one value per requested ID."""
    if value is None or isinstance(value, str) or np.ndim(value) == 0:
        return value
    value = np.asarray(value, dtype=float)
    if value.shape != (n,):
        raise ValueError(f"{name} must be a single value or one value per source ID ({n})")
    return value


def ages_from_gaia(
    source_ids,
    *,
    a0,
    a0_err=0.0,
    dust_map=None,
    mh=None,
    inflate_parallax_error=True,
    gaia_table=None,
    credentials_file=None,
    **age_kwargs,
):
    """Query Gaia DR3, de-extinct the photometry, propagate its uncertainties and infer CMD ages.

    source_ids: Gaia DR3 source ID(s).
    a0: required.  A_0 (mag at 550 nm) for every star, one value per source ID, or "dustmap"
        to integrate the EXPLORE 3D dust maps (dust_map=ExploreDustMap(...) or CMD_INTERP_DUST_DIR).
    a0_err: uncertainty of a supplied A_0 (single value or one per source ID).
    mh: stellar [M/H] (single value or one per source ID); default solar.
    gaia_table: optional DataFrame of gaia_source rows (with dr3_source_id or source_id) to use
        instead of querying the archive.
    Other keyword arguments go to ages_for_table / infer_ages (systematic_mag, reference_ages,
    exclude_references, max_color, return_posterior, ...).

    Returns (table, result): one row per requested ID in the order given (found_in_gaia flags
    IDs missing from DR3), holding the gaia_source columns, the de-extincted photometry and
    uncertainties (deredden_gaia) and the cmd_age_* columns; result is the infer_ages output
    for the fitted rows (result['row'] indexes the table).
    """
    ids = np.atleast_1d(np.asarray(source_ids, dtype=np.int64))
    n = len(ids)
    a0 = _per_star(a0, n, "a0")
    a0_err = _per_star(a0_err, n, "a0_err")
    mh = _per_star(mh, n, "mh")

    if gaia_table is None:
        rows = query_gaia_dr3(ids, credentials_file=credentials_file)
    else:
        rows = gaia_table.rename(columns={"source_id": "dr3_source_id", "SOURCE_ID": "dr3_source_id"})
    rows = rows.drop_duplicates("dr3_source_id")
    table = pd.DataFrame({"dr3_source_id": ids}).merge(rows, on="dr3_source_id", how="left")
    table["found_in_gaia"] = table.dr3_source_id.isin(rows.dr3_source_id).to_numpy()

    extra = {}
    if isinstance(a0, np.ndarray):
        table["_a0"] = a0; a0 = "_a0"
    if isinstance(a0_err, np.ndarray):
        table["_a0_err"] = a0_err; a0_err = "_a0_err"
    if isinstance(mh, np.ndarray):
        table["_mh"] = mh; extra["mh_col"] = "_mh"
    elif mh is not None:
        extra["mh"] = float(mh)

    table = deredden_gaia(table, a0=a0, a0_err=a0_err, dust_map=dust_map, inflate_parallax_error=inflate_parallax_error)
    age_kwargs.setdefault("_warn_stacklevel", 4)
    table, result = ages_for_table(table, **extra, **age_kwargs)
    return table.drop(columns=[c for c in ("_a0", "_a0_err", "_mh") if c in table]), result


def star_record(table, result, index=0):
    """One row of ages_from_gaia output as a star_age()-style dict (scalars plus 'posterior' and 'age_grid_myr')."""
    row = table.iloc[index]
    out = {key: (value.item() if hasattr(value, "item") else value) for key, value in row.items()}
    hit = np.flatnonzero(np.asarray(result["row"]) == index)
    for key, value in result.items():
        if key in ("age_grid_myr", "row"):
            continue
        if key == "posterior":
            out[key] = (np.asarray(value[hit[0]]) if len(hit) and value is not None
                        else np.zeros(len(result["age_grid_myr"]), np.float32))
        else:
            item = value[hit[0]] if len(hit) else np.nan
            out[key] = item.item() if hasattr(item, "item") else item
    out["age_grid_myr"] = result["age_grid_myr"]
    out.update(color=out.get("BpmRp0"), abs_g=out.get("phot_g_mean_mag_0"),
               color_err=out.get("BpmRp0_err"), abs_g_err=out.get("phot_g_mean_mag_0_err"))
    return out


def star_age_from_gaia(source_id, *, a0, a0_err=0.0, dust_map=None, mh=None, **kwargs):
    """ages_from_gaia for one star, returned like star_age(): a dict of scalars plus 'posterior'
    (n_age), 'age_grid_myr' and the star's Gaia / de-extincted photometry columns."""
    kwargs.setdefault("_warn_stacklevel", 5)
    table, result = ages_from_gaia([source_id], a0=a0, a0_err=a0_err, dust_map=dust_map, mh=mh, **kwargs)
    return star_record(table, result)
