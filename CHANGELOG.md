# Changelog

## 0.1.0 (2026-10-15)

First release.

- `infer_ages`, `star_age` and `ages_for_table`: age posteriors from six empirical reference sequences
  (UBC 17a, β Pic, IC 2602, α Per, Pleiades, Coma Ber).
- Defaults: hybrid scatter and a 0.05 mag systematic. These reproduce the original research pipeline exactly.
- Optional per-star [M/H] correction from PARSEC v1.2S offsets at fixed colour.
- Optional unresolved-binary likelihood mixture.
- `cmd-interp` command-line tool and diagnostic plots.
