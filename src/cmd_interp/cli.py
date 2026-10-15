"""Command-line interface: ``cmd-interp star ...`` and ``cmd-interp table ...``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .core import DEFAULT_SCATTER_MODE, DEFAULT_SYSTEMATIC_MAG, MAX_COLOR, SCATTER_MODES, star_age
from .grid import age_grid
from .references import REFERENCE_AGES, load_references
from .table import DEFAULT_ABS_G_ERR, ages_for_table, read_table, write_table


def _add_model_arguments(parser):
    parser.add_argument("--systematic-mag", type=float, default=DEFAULT_SYSTEMATIC_MAG,
                        help=f"Extra magnitude scatter added in quadrature (default {DEFAULT_SYSTEMATIC_MAG}).")
    parser.add_argument("--scatter-mode", choices=SCATTER_MODES, default=DEFAULT_SCATTER_MODE,
                        help=f"Empirical-scatter prescription (default {DEFAULT_SCATTER_MODE}).")
    parser.add_argument("--reference-age", action="append", default=[], metavar="NAME=MYR",
                        help="Use a different age for a reference cluster, e.g. --reference-age beta_Pic=21 "
                             "(repeatable; aliases such as Pleiades or alpha_Per work). Defaults: "
                             + ", ".join(f"{k}={v:g}" for k, v in REFERENCE_AGES.items()) + ".")
    parser.add_argument("--exclude-reference", action="append", default=[], metavar="NAME",
                        help=f"Leave out a reference cluster (repeatable). Choices: {', '.join(REFERENCE_AGES)}.")
    parser.add_argument("--max-color", type=float, default=MAX_COLOR,
                        help=f"Stars redder than this (BP-RP)_0 get no age (default {MAX_COLOR}).")
    parser.add_argument("--reference-dir", type=Path, default=None,
                        help="Use a custom reference library instead of the bundled one.")
    parser.add_argument("--binary-fraction", type=float, default=0.0,
                        help="Unresolved-binary mixture fraction (default 0 = off).")
    parser.add_argument("--binary-dm", type=float, nargs="+", default=[],
                        help="Brightening nodes (mag) of the binary component, e.g. 0.3 0.5 0.75.")


def _references(args):
    ages = {}
    for item in args.reference_age:
        name, sep, value = item.partition("=")
        if not sep:
            raise SystemExit(f"--reference-age needs NAME=MYR, got {item!r}")
        ages[name.strip()] = float(value)
    try:
        return load_references(args.reference_dir, ages or None, exclude=args.exclude_reference or None)
    except ValueError as error:
        raise SystemExit(str(error)) from error


def _model_kwargs(args):
    return dict(systematic_mag=args.systematic_mag, scatter_mode=args.scatter_mode, references=_references(args),
                binary_fraction=args.binary_fraction, binary_dm=args.binary_dm, max_color=args.max_color)


def build_parser():
    parser = argparse.ArgumentParser(prog="cmd-interp", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    star = sub.add_parser("star", help="Age of one star.")
    star.add_argument("--color", type=float, required=True, help="Extinction-corrected (BP-RP)_0.")
    star.add_argument("--abs-g", type=float, required=True, help="Extinction-corrected absolute G.")
    star.add_argument("--color-err", type=float, required=True)
    star.add_argument("--abs-g-err", type=float, required=True)
    star.add_argument("--mh", type=float, default=None, help="Stellar [M/H] in dex (default: solar).")
    star.add_argument("--plot", type=Path, default=None, help="Save a CMD + posterior figure here.")
    star.add_argument("--posterior", type=Path, default=None, help="Save age grid and posterior density as CSV.")
    star.add_argument("--json", action="store_true", help="Machine-readable output.")
    _add_model_arguments(star)

    table = sub.add_parser("table", help="Ages for every star in a table (.csv/.feather/.parquet/.fits).")
    table.add_argument("input", type=Path)
    table.add_argument("output", type=Path, help="Output table (.csv/.feather/.parquet).")
    table.add_argument("--color-col", default="BpmRp0")
    table.add_argument("--abs-g-col", default="phot_g_mean_mag_0")
    table.add_argument("--color-err-col", default="BpmRp0_err")
    table.add_argument("--abs-g-err-col", default="phot_g_mean_mag_0_err")
    table.add_argument("--abs-g-err", type=float, default=DEFAULT_ABS_G_ERR,
                       help=f"Fallback absolute-G error (default {DEFAULT_ABS_G_ERR}).")
    metal = table.add_mutually_exclusive_group()
    metal.add_argument("--mh", type=float, default=None, help="One [M/H] for every star (default: solar).")
    metal.add_argument("--mh-col", default=None, help="Per-star [M/H] column (NaN = solar).")
    table.add_argument("--fit-col", default=None, help="Boolean column selecting the stars to fit (default: all).")
    table.add_argument("--posterior-output", type=Path, default=None,
                       help="Write the posteriors (float32 .npy, unit-sum node weights, one row per fitted star) plus "
                            "<stem>_age_grid.npy and <stem>_index.csv (row of each posterior in the output table).")
    table.add_argument("--overwrite", action="store_true")
    table.add_argument("--chunk-size", type=int, default=2000)
    _add_model_arguments(table)

    dered = sub.add_parser("deredden", help="De-extincted absolute G, BP, RP, (BP-RP)_0 and uncertainties from gaia_source columns.")
    dered.add_argument("input", type=Path, help="Table with parallax, parallax_error, ruwe, phot_*_mean_mag, phot_*_mean_flux(_error).")
    dered.add_argument("output", type=Path)
    a0 = dered.add_mutually_exclusive_group()
    a0.add_argument("--a0", type=float, default=0.0, help="A_0 (mag at 550 nm) for every star (default 0).")
    a0.add_argument("--a0-col", default=None, help="Per-star A_0 column.")
    a0.add_argument("--dustmap", action="store_true", help="Integrate the EXPLORE dust maps to each star (needs ra, dec).")
    dered.add_argument("--a0-err", type=float, default=0.0, help="Uncertainty of a supplied A_0 (default 0).")
    dered.add_argument("--a0-err-col", default=None, help="Per-star A_0 uncertainty column.")
    dered.add_argument("--dust-values", type=Path, default=None, help="EXPLORE density cube (.h5).")
    dered.add_argument("--dust-errors", type=Path, default=None, help="EXPLORE density-error cube (.h5).")
    dered.add_argument("--no-parallax-inflation", action="store_true", help="Use the catalogue parallax errors as they are.")
    dered.add_argument("--overwrite", action="store_true")

    gaia = sub.add_parser("gaia", help="Query Gaia DR3 by source ID, de-extinct, and infer ages (needs astroquery).")
    gaia.add_argument("source_ids", nargs="*", type=int, help="Gaia DR3 source IDs.")
    gaia.add_argument("--ids-file", type=Path, default=None, help="Text file with one source ID per line.")
    gaia_a0 = gaia.add_mutually_exclusive_group(required=True)
    gaia_a0.add_argument("--a0", type=float, help="A_0 (mag at 550 nm) for every star.")
    gaia_a0.add_argument("--dustmap", action="store_true", help="A_0 from the EXPLORE 3D dust maps.")
    gaia.add_argument("--a0-err", type=float, default=0.0, help="Uncertainty of --a0 (default 0).")
    gaia.add_argument("--dust-values", type=Path, default=None, help="EXPLORE density cube (.h5).")
    gaia.add_argument("--dust-errors", type=Path, default=None, help="EXPLORE density-error cube (.h5).")
    gaia.add_argument("--mh", type=float, default=None, help="[M/H] for every star (default solar).")
    gaia.add_argument("--output", type=Path, default=None, help="Write the full table (.csv/.feather/.parquet).")
    gaia.add_argument("--plot", type=Path, default=None, help="CMD + posterior figure (single star only).")
    gaia.add_argument("--credentials-file", type=Path, default=None, help="Optional astroquery Gaia credentials file.")
    _add_model_arguments(gaia)
    return parser


def _run_gaia(args):
    from .gaia import ages_from_gaia, star_record
    ids = list(args.source_ids)
    if args.ids_file is not None:
        ids += [int(line.split()[0]) for line in args.ids_file.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not ids:
        raise SystemExit("Give at least one Gaia DR3 source ID (or --ids-file).")
    dust_map = None
    if args.dustmap:
        from .dustmap import ExploreDustMap
        dust_map = ExploreDustMap(args.dust_values, args.dust_errors)
    a0 = "dustmap" if args.dustmap else args.a0
    kwargs = _model_kwargs(args)
    if args.plot is not None and len(ids) != 1:
        raise SystemExit("--plot works for a single source ID")
    table, result = ages_from_gaia(ids, a0=a0, a0_err=args.a0_err, dust_map=dust_map, mh=args.mh,
                                   credentials_file=args.credentials_file, return_posterior=args.plot is not None,
                                   _warn_stacklevel=4, **kwargs)
    if args.plot is not None:
        from .plotting import plot_star
        fig, _ = plot_star(star_record(table, result), references=kwargs["references"])
        fig.savefig(args.plot, dpi=150, bbox_inches="tight")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_table(table, args.output)
    show = ["dr3_source_id", "BpmRp0", "BpmRp0_err", "phot_g_mean_mag_0", "phot_g_mean_mag_0_err", "extinction_a0",
            "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr", "cmd_age_limit"]
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(table[show].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    missing = table.loc[~table.found_in_gaia, "dr3_source_id"].tolist()
    if missing:
        print(f"Not found in Gaia DR3: {missing}")
    return 0


def _run_deredden(args):
    from .photometry import deredden_gaia
    if args.output.exists() and not args.overwrite:
        raise SystemExit(f"Output exists: {args.output}; pass --overwrite to replace it.")
    data = read_table(args.input)
    dust_map = None
    if args.dustmap:
        from .dustmap import ExploreDustMap
        dust_map = ExploreDustMap(args.dust_values, args.dust_errors)
        a0 = "dustmap"
    else:
        a0 = args.a0_col if args.a0_col else args.a0
    out = deredden_gaia(data, a0=a0, a0_err=args.a0_err_col if args.a0_err_col else args.a0_err, dust_map=dust_map,
                        inflate_parallax_error=not args.no_parallax_inflation)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_table(out, args.output)
    print(f"{len(out):,} rows; {int(np.isfinite(out.BpmRp0).sum()):,} with de-extincted photometry -> {args.output}")
    return 0


def _run_star(args):
    kwargs = _model_kwargs(args)
    star = star_age(args.color, args.abs_g, args.color_err, args.abs_g_err, mh=args.mh, **kwargs)
    summary = {key: star[key] for key in ("age_map_myr", "age_myr_16", "age_myr_50", "age_myr_84", "age_myr_025",
                                          "age_myr_975", "limit", "n_references", "mh", "metallicity_clipped",
                                          "color_out_of_range", "blue_caution")}
    summary = {key: (None if isinstance(value, float) and not np.isfinite(value) else value) for key, value in summary.items()}
    if args.posterior is not None:
        from .grid import posterior_density
        density = posterior_density(star["posterior"], star["age_grid_myr"])[0]
        np.savetxt(args.posterior, np.column_stack([star["age_grid_myr"], density]), delimiter=",",
                   header="age_myr,posterior_density_per_myr", comments="")
    if args.plot is not None:
        from .plotting import plot_star
        fig, _ = plot_star(star, references=kwargs["references"])
        fig.savefig(args.plot, dpi=150, bbox_inches="tight")
    if args.json:
        print(json.dumps(summary, indent=2))
    elif star["color_out_of_range"]:
        print(f"No age: (BP-RP)_0 = {args.color:g} is outside the fitted colour range (max {args.max_color:g}).")
    elif not star["has_posterior"]:
        print("No age: the star lies outside the colour range covered by at least two reference sequences.")
    else:
        if star["limit"] == "young":
            print(f"Age < {star['youngest_age']:g} Myr: brighter than the youngest covering sequence.")
        elif star["limit"] == "old":
            print(f"Age > {star['oldest_age']:g} Myr: fainter than the oldest covering sequence.")
        prefix = "Posterior (flat beyond the references): " if star["limit"] else ""
        print(f"{prefix}MAP age {star['age_map_myr']:.1f} Myr; median {star['age_myr_50']:.1f} Myr; "
              f"68% interval [{star['age_myr_16']:.1f}, {star['age_myr_84']:.1f}] Myr")
    return 0


def _run_table(args):
    if args.output.exists() and not args.overwrite:
        raise SystemExit(f"Output exists: {args.output}; pass --overwrite to replace it.")
    kwargs = _model_kwargs(args)
    data = read_table(args.input)
    posterior_out = None
    if args.posterior_output is not None:
        # one row per fitted star: the same selection as ages_for_table (fit column and finite photometry)
        color_err = data[args.color_err_col].to_numpy(float)
        fitted = (np.isfinite(data[args.color_col].to_numpy(float)) & np.isfinite(data[args.abs_g_col].to_numpy(float))
                  & np.isfinite(color_err) & (color_err >= 0))
        if args.fit_col is not None:
            fitted &= data[args.fit_col].fillna(False).to_numpy(bool)
        n_fit = int(fitted.sum())
        grid = age_grid([r["age_myr"] for r in kwargs["references"]])
        posterior_out = np.lib.format.open_memmap(args.posterior_output, mode="w+", dtype=np.float32, shape=(n_fit, len(grid)))
        posterior_out[:] = 0.0
    table, result = ages_for_table(
        data, color_col=args.color_col, abs_g_col=args.abs_g_col, color_err_col=args.color_err_col,
        abs_g_err_col=args.abs_g_err_col, abs_g_err=args.abs_g_err, mh=args.mh, mh_col=args.mh_col,
        fit_mask=args.fit_col, return_posterior=False, posterior_out=posterior_out, chunk_size=args.chunk_size,
        **kwargs,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_table(table, args.output)
    if posterior_out is not None:
        posterior_out.flush()
        stem = args.posterior_output.with_suffix("")
        np.save(f"{stem}_age_grid.npy", result["age_grid_myr"])
        pd.DataFrame({"row": result["row"], "has_posterior": result["has_posterior"], "limit": result["limit"]}).to_csv(
            f"{stem}_index.csv", index=False)
    n_fit = int(result["has_posterior"].sum())
    print(f"{len(table):,} rows; {len(result['row']):,} fitted; {n_fit:,} with a posterior "
          f"({int((result['limit'] == 'young').sum()):,} young limits, {int((result['limit'] == 'old').sum()):,} old limits) "
          f"-> {args.output}")
    return 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "star":
        return _run_star(args)
    if args.command == "deredden":
        return _run_deredden(args)
    if args.command == "gaia":
        return _run_gaia(args)
    return _run_table(args)


if __name__ == "__main__":
    raise SystemExit(main())
