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
- `cmd-interp` command-line tool and diagnostic plots.
