#!/usr/bin/env python
"""Build the bundled reference library from de-extincted cluster sequences and member lists.

For each reference cluster the source directory must hold

    <cluster>_corrected_isochrone.feather   the empirical sequence, with columns BpmRp0, GMag0
    <cluster>_members.feather               the members, with columns GaiaDR3, BpmRp0, phot_g_mean_mag_0

(de-extincted (BP-RP)_0 and absolute M_G,0).  Slim copies holding only those columns are written to
src/cmd_interp/data/references/ as <cluster>_sequence.feather and <cluster>_members.feather.
The PARSEC metallicity offsets are built separately by scripts/build_parsec_offsets.py.

Run:  python scripts/slim_reference_data.py --source path/to/sequences
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CLUSTERS = ["UBC_17a", "beta_Pic", "IC_2602", "Melotte_20", "Melotte_22", "Melotte_111"]
OUTPUT = Path(__file__).resolve().parents[1] / "src" / "cmd_interp" / "data" / "references"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="Directory holding the sequence and member files.")
    parser.add_argument("--clusters", nargs="+", default=CLUSTERS, help="Reference clusters to include.")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    source = args.source.expanduser().resolve()
    args.output.mkdir(parents=True, exist_ok=True)

    for cluster in args.clusters:
        sequence = pd.read_feather(source / f"{cluster}_corrected_isochrone.feather", columns=["BpmRp0", "GMag0"])
        members = pd.read_feather(source / f"{cluster}_members.feather", columns=["GaiaDR3", "BpmRp0", "phot_g_mean_mag_0"])
        sequence.to_feather(args.output / f"{cluster}_sequence.feather")
        members.reset_index(drop=True).to_feather(args.output / f"{cluster}_members.feather")
        print(f"{cluster}: {len(sequence)} sequence points, {len(members)} members")


if __name__ == "__main__":
    main()
