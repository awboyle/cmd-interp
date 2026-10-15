#!/usr/bin/env python
"""Rebuild src/cmd_interp/data/parsec_mg_offsets.csv from PARSEC v1.2S isochrones (provenance only).

For each [M/H], age and (BP-RP)_0 on the grid below, M_G is read where the PMS/MS branch
(label <= 1) of the PARSEC isochrone nearest in log age crosses that colour, walking from the
lowest mass upward so that the first crossing is taken.  dMG = M_G([M/H]) - M_G(solar)
(positive = fainter).

Isochrones (CMD 3.7, Gaia EDR3 / YBC photometry, logAge 6.60-10.10 in 0.05 dex steps) are read
from --cache (parquet files named parsec_gaiaEDR3_MH<+0p00>.parquet with columns logAge, Mini,
label, Gmag, G_BPmag, G_RPmag) or, with --download, fetched with ezpadova
(pip install git+https://github.com/mfouesneau/ezpadova) and written to the cache.

Run:  python scripts/build_parsec_offsets.py --cache parsec_cache [--download] [--output table.csv]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

METALLICITIES = [-0.5, -0.3, -0.16, -0.1, -0.08, -0.05, -0.03, 0.0, 0.1, 0.15, 0.2, 0.3]
AGES_MYR = [5, 10, 20, 30, 50, 100, 200, 500, 1000]
COLOURS = [round(0.8 + 0.1 * k, 1) for k in range(22)]   # 0.8 ... 2.9
LOGAGE_MIN, LOGAGE_MAX, LOGAGE_STEP = 6.60, 10.10, 0.05
PHOTSYS = "YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat"
MAX_LABEL = 1   # PMS (0) and MS (1)
OUTPUT = Path(__file__).resolve().parents[1] / "src" / "cmd_interp" / "data" / "parsec_mg_offsets.csv"


def _cache_file(cache, mh):
    return Path(cache) / f"parsec_gaiaEDR3_MH{f'{mh:+.2f}'.replace('.', 'p')}.parquet"


def load_grid(cache, mh, download=False):
    path = _cache_file(cache, mh)
    if not path.exists():
        if not download:
            raise FileNotFoundError(f"{path} missing; pass --download to fetch it")
        from ezpadova import get_isochrones
        df = get_isochrones(logage=(LOGAGE_MIN, LOGAGE_MAX + 1e-4, LOGAGE_STEP), MH=(mh, mh, 0.0), photsys_file=PHOTSYS)
        df = df[["logAge", "MH", "Mini", "label", "Gmag", "G_BPmag", "G_RPmag"]]
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(path)
    df = pd.read_parquet(path)
    df = df[df.label <= MAX_LABEL].copy()
    df["colour"] = df.G_BPmag - df.G_RPmag
    return df


def nearest_logage(df, age_myr):
    grid = np.sort(df.logAge.unique())
    return grid[np.argmin(np.abs(grid - np.log10(age_myr * 1e6)))]


def mag_at_colour(iso, colour):
    """M_G where the PMS/MS branch first crosses ``colour`` (lowest mass upward); NaN if never."""
    iso = iso.sort_values("Mini")
    c = iso.colour.to_numpy()
    m = iso.Gmag.to_numpy()
    cross = np.flatnonzero((c[:-1] - colour) * (c[1:] - colour) <= 0)
    if len(cross) == 0:
        return np.nan
    i = cross[0]
    if c[i + 1] == c[i]:
        return m[i]
    return m[i] + (m[i + 1] - m[i]) * (colour - c[i]) / (c[i + 1] - c[i])


def build(cache, download=False):
    rows = []
    for mh in METALLICITIES:
        df = load_grid(cache, mh, download)
        for age in AGES_MYR:
            la = nearest_logage(df, age)
            iso = df[df.logAge == la]
            for colour in COLOURS:
                rows.append(dict(MH=mh, age_myr=age, logAge=la, colour=colour, MG=mag_at_colour(iso, colour)))
    table = pd.DataFrame(rows)
    solar = table[table.MH == 0.0].set_index(["age_myr", "colour"]).MG
    table["dMG"] = table.MG.to_numpy() - solar.reindex(list(zip(table.age_myr, table.colour))).to_numpy()
    return table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    table = build(args.cache, args.download)
    if not np.isfinite(table.dMG).all():
        raise SystemExit("Some (age, colour) cells have no PMS/MS crossing; the offset grid must be complete.")
    table.to_csv(args.output, index=False)
    print(f"{len(table)} rows -> {args.output}")


if __name__ == "__main__":
    main()
