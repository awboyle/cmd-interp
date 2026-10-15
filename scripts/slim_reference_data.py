#!/usr/bin/env python
"""Build the bundled reference library from the research archive (provenance only).

Reads (never modifies) the de-extincted empirical sequences and member lists of the six
production reference clusters and writes slim copies, holding only the columns the age
engine uses, to src/cmd_interp/data/references/:

    <cluster>_sequence.feather  BpmRp0, GMag0                       (empirical sequence)
    <cluster>_members.feather   GaiaDR3, BpmRp0, phot_g_mean_mag_0  (members, for the local scatter)

and copies the PARSEC metallicity-offset table to src/cmd_interp/data/parsec_mg_offsets.csv.

Run:  python scripts/slim_reference_data.py --source ~/python/isochrone_fitter
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CLUSTERS = ["UBC_17a", "beta_Pic", "IC_2602", "Melotte_20", "Melotte_22", "Melotte_111"]
PACKAGE_DATA = Path(__file__).resolve().parents[1] / "src" / "cmd_interp" / "data"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True,
                        help="isochrone_fitter directory holding deextincted_isochrones/ and metallicity_sweep/")
    args = parser.parse_args(argv)
    source = args.source.expanduser().resolve()
    out_dir = PACKAGE_DATA / "references"
    out_dir.mkdir(parents=True, exist_ok=True)

    for cluster in CLUSTERS:
        sequence = pd.read_feather(source / "deextincted_isochrones" / f"{cluster}_corrected_isochrone.feather",
                                   columns=["BpmRp0", "GMag0"])
        members = pd.read_feather(source / "deextincted_isochrones" / f"{cluster}_members.feather",
                                  columns=["GaiaDR3", "BpmRp0", "phot_g_mean_mag_0"])
        sequence.to_feather(out_dir / f"{cluster}_sequence.feather")
        members.reset_index(drop=True).to_feather(out_dir / f"{cluster}_members.feather")
        print(f"{cluster}: {len(sequence)} sequence points, {len(members)} members")

    offsets = pd.read_csv(source / "metallicity_sweep" / "metallicity_offsets_fine.csv")
    offsets[["MH", "age_myr", "logAge", "colour", "MG", "dMG"]].to_csv(PACKAGE_DATA / "parsec_mg_offsets.csv", index=False)
    print(f"PARSEC offsets: {len(offsets)} rows")


if __name__ == "__main__":
    main()
