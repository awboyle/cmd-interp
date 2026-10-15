# Changelog

## 0.1.0 (2026-10-15)

First release.

- `infer_ages`, `star_age` and `ages_for_table`: age posteriors from six empirical reference sequences
  (UBC 17a, β Pic, IC 2602, α Per, Pleiades, Coma Ber).
- Defaults: hybrid scatter and a 0.05 mag systematic. These reproduce the original research pipeline exactly,
  except for stars brighter than the youngest covering sequence.
- Hybrid scatter: stars brighter than the youngest covering sequence use that sequence's filled width at every age.
  This gives a flat young posterior instead of a spurious peak near 10 Myr. The original pipeline used the
  age-dependent width for these stars.
- Optional per-star [M/H] correction from PARSEC v1.2S offsets at fixed colour.
- Optional unresolved-binary likelihood mixture.
- Ages only for (BP-RP)_0 <= 3.0 by default (`max_color`); redder stars get NaN and `color_out_of_range`.
- `reference_ages=` / `exclude_references=` (and `--reference-age NAME=MYR`): change the age of any reference
  cluster, or leave clusters out. Aliases such as "Pleiades" and "alpha Per" are accepted.
- `deredden_gaia` (and `cmd-interp deredden`): de-extincted absolute G, BP, RP and (BP-RP)_0 with uncertainties
  from gaia_source columns. A_0 is user-supplied (default) or integrated from the EXPLORE 3D dust maps
  (`ExploreDustMap`). This reproduces the original calculate_dr3_temps.py pipeline to < 1e-12 mag.
- `cmd-interp` command-line tool and diagnostic plots.
