# Changelog

## 0.1.0 (2026-10-15)

First release.

- `infer_ages`, `star_age` and `ages_for_table`: age posteriors from six empirical reference sequences
  (UBC 17a, β Pic, IC 2602, α Per, Pleiades, Coma Ber).
- Likelihood width: the empirical scatter of the reference sequences, interpolated in log age, plus a 0.05 mag
  systematic. Stars outside the envelope of the reference sequences keep a single width at every age; stars
  brighter than the youngest covering sequence use that sequence's width, so they get a flat young posterior.
- Optional per-star [M/H] correction from PARSEC v1.2S offsets at fixed color.
- Ages only for (BP-RP)_0 <= 3.0 by default (`max_color`); redder stars get NaN and `color_out_of_range`.
- `BlueColorWarning` and a `blue_caution` flag for stars bluer than (BP-RP)_0 = 0.8, where the reference
  sequences start to cross. Ages are still returned.
- `reference_ages=` / `exclude_references=` (and `--reference-age NAME=MYR`): change the age of any reference
  cluster, or leave clusters out. Aliases such as "Pleiades" and "alpha Per" are accepted.
- `deredden_gaia` (and `cmd-interp deredden`): de-extincted absolute G, BP, RP and (BP-RP)_0 with uncertainties
  from gaia_source columns. A_0 is user-supplied (default) or integrated from the EXPLORE 3D dust maps
  (`ExploreDustMap`; the FITS cubes from CDS J/A+A/664/A174 or the G-Tomo HDF5 cubes).
- `ages_from_gaia` / `star_age_from_gaia` (and `cmd-interp gaia`): Gaia DR3 source IDs -> anonymous archive query
  -> de-extinction and uncertainties -> CMD ages, with A_0 given explicitly or from the dust maps.
- `cmd-interp` command-line tool and diagnostic plots.
