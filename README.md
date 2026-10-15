# cmd-interp

Stellar ages from a star's position in the Gaia colour–magnitude diagram, interpolated
between **empirical sequences of benchmark clusters** rather than model isochrones.
The method is aimed at pre-main-sequence stars (≈5–700 Myr), where the
contraction of K and M dwarfs toward the main sequence is age-sensitive. Each star gets a
full age posterior.

```python
import cmd_interp

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
star["age_map_myr"], star["age_myr_16"], star["age_myr_84"]
# (74.5, 62.7, 83.9)
```

## Installation

```bash
pip install git+https://github.com/awboyle/cmd-interp
pip install "cmd-interp[plot] @ git+https://github.com/awboyle/cmd-interp"       # + matplotlib for the figures
pip install "cmd-interp[plot,dust] @ git+https://github.com/awboyle/cmd-interp"  # + h5py/astropy for the dust maps
```

Requires Python ≥ 3.10, numpy, scipy, pandas and pyarrow. The reference library (≈300 kB) ships with the package.

## Inputs

All photometry must be **extinction corrected**, matching the reference sequences:

| quantity | meaning |
|---|---|
| `color` | Gaia DR3 (G_BP − G_RP)₀ |
| `abs_g` | absolute, extinction-corrected Gaia G, M_G,0 |
| `color_err`, `abs_g_err` | 1σ uncertainties (mag); include the parallax error in `abs_g_err` |
| `mh` (optional) | stellar [M/H] in dex; default solar |

Ages are only computed for **(BP−RP)₀ ≤ 3.0**. Redder stars get NaN and `color_out_of_range = True`.
The limit can be changed with `max_color=`, and an optional blue limit set with `min_color=`.

**Caution bluer than (BP−RP)₀ = 0.8.** Here the reference sequences start to cross: they are no longer
ordered by age, especially near the Coma Ber turn-off. Ages are still returned. But every call that
includes such stars issues a `cmd_interp.BlueColorWarning`, and the stars are flagged `blue_caution`
(`cmd_age_blue_caution` in tables). To silence the warning:
`warnings.filterwarnings("ignore", category=cmd_interp.BlueColorWarning)`.

If you are starting from Gaia DR3 catalogue columns, `cmd_interp.deredden_gaia` computes these
inputs and their uncertainties (see [De-extincted photometry](#de-extincted-photometry)).

## Usage

### Many stars

```python
import numpy as np
import cmd_interp

result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err)        # arrays, one entry per star
result["age_map_myr"], result["age_myr_50"], result["limit"]
result["posterior"]      # (n_star, n_age) float32, unit-sum weights on result["age_grid_myr"]
```

From a table (DataFrame, or a .csv/.feather/.parquet/.fits path):

```python
table, result = cmd_interp.ages_for_table(
    "stars.parquet",
    color_col="BpmRp0", abs_g_col="phot_g_mean_mag_0",
    color_err_col="BpmRp0_err", abs_g_err_col="phot_g_mean_mag_0_err",
    mh_col="mh_xgboost",          # optional per-star [M/H]; NaN rows are treated as solar
    fit_mask="is_pms",            # optional: only fit these rows
)
```

This appends `cmd_age_myr` (MAP, or the reference age at a limit), `cmd_age_lower_myr` /
`cmd_age_upper_myr` (16th/84th percentiles), `cmd_age_median_myr`, `cmd_age_p2p5_myr` /
`cmd_age_p97p5_myr`, `cmd_age_limit`, `cmd_age_n_references`, `cmd_age_mh`,
`cmd_age_metallicity_clipped`, `cmd_age_color_out_of_range` and `cmd_age_fitted`. `result["row"]` maps
each posterior row to its table row.

### Changing the reference ages

The adopted ages (table under [Method](#method)) are used by default. To use a different age for any
reference, pass only the ones you want to change; the others keep their defaults:

```python
cmd_interp.star_age(2.0, 8.0, 0.01, 0.03, reference_ages={"beta_Pic": 21.0})
cmd_interp.ages_for_table(table, reference_ages={"Pleiades": 112, "alpha Per": 85})
cmd_interp.infer_ages(..., exclude_references=["Coma Ber"])                 # leave a reference out
cmd_interp.resolve_reference_ages({"IC 2602": 46})                         # inspect the resulting set
```

Names are matched loosely: `beta_Pic`, `beta pic`, `bPMG`, `Pleiades`, `Melotte_22`, `alpha Per`, `Coma Ber`,
`IC 2602` and similar all work. The age grid and the interpolation follow the new ages. For many
calls, load the library once with `refs = cmd_interp.load_references(reference_ages=...)` and pass
`references=refs`.

### De-extincted photometry

`deredden_gaia` takes gaia_source columns (`parallax`, `parallax_error`, `ruwe`, `phot_{g,bp,rp}_mean_mag`,
`phot_{g,bp,rp}_mean_flux`, `phot_{g,bp,rp}_mean_flux_error`, plus `ra`/`dec` for the dust map). It returns
absolute, extinction-corrected `phot_g_mean_mag_0`, `BpmRp0` (and BP, RP) with uncertainties: exactly the
columns `ages_for_table` reads by default.

```python
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12, a0_err=0.03)   # one A_0 for every star (default 0)
phot = cmd_interp.deredden_gaia(gaia_table, a0="a0_column")          # per-star A_0 column
dust = cmd_interp.ExploreDustMap("explore_cube_density_values_010pc_v2.h5",
                                 "explore_cube_density_errors_050pc_v2.h5")
phot = cmd_interp.deredden_gaia(gaia_table, a0="dustmap", dust_map=dust)   # A_0 from the EXPLORE 3D dust maps
table, result = cmd_interp.ages_for_table(phot)
```

- **Extinction law.** A_0 (extinction at 550 nm) is converted with the Fitzpatrick (2019) Gaia EDR3
  coefficients k_X((BP−RP)₀, A_0). The upper-HRD table is used for G − 5 log₁₀(d/10 pc) < 5, otherwise the
  main-sequence table. The intrinsic colour is iterated to 0.01 mag.
- **Dust map.** `a0="dustmap"` integrates the EXPLORE density cube (Vergely et al. 2022) to d = 1000/parallax in
  5 pc steps. It also returns the line-of-sight integral of the density-error cube (perfectly correlated
  voxels, so an upper bound) and the local density dA_0/dd. The cubes (~230 MB) are not bundled: download
  them from EXPLORE and pass their paths, or set `CMD_INTERP_DUST_DIR`.
- **Uncertainties.** These combine:
  - flux errors and the Vega zero points (Riello et al. 2021);
  - the parallax error, inflated with RUWE and G and propagated through both the distance modulus
    and A_0(d);
  - the A_0 error (dust-map, or your `a0_err`).

  In (BP−RP)₀ the distance modulus cancels and only E(BP−RP) carries the A_0 error.
- **Not included.** The full equations are in the `cmd_interp.photometry` docstring. Not included: the parallax
  zero point, a photometric error floor, and extinction-law uncertainty.
- **Bad parallaxes.** Stars with parallax ≤ 0 get NaN.

### Command line

```bash
cmd-interp star --color 2.0 --abs-g 8.0 --color-err 0.01 --abs-g-err 0.03 --mh 0.1 --plot star.png
cmd-interp table stars.parquet ages.parquet --mh-col mh_xgboost --posterior-output posteriors.npy
cmd-interp table stars.csv ages.csv --reference-age beta_Pic=21 --exclude-reference Coma_Ber
cmd-interp deredden gaia.csv phot.csv --a0 0.12 --a0-err 0.03          # or --a0-col, or --dustmap --dust-values ... --dust-errors ...
cmd-interp table --help
```

With `--posterior-output`, the posteriors are written as a float32 memmap (`posteriors.npy`) along with
`posteriors_age_grid.npy` and `posteriors_index.csv`, which gives the table row of each posterior. This
keeps memory bounded for large catalogues (1 M stars ≈ 5 GB).

### Plots

```python
from cmd_interp.plotting import plot_star
fig, axes = plot_star(cmd_interp.star_age(2.0, 8.0, 0.01, 0.03))
```

## Settings

| argument | default | meaning |
|---|---|---|
| `systematic_mag` | **0.05** | extra magnitude scatter added in quadrature at every age; it covers the field metallicity spread and the zero point of the reference sequences |
| `mh` | solar | stellar [M/H]; see below |
| `scatter_mode` | `"hybrid"` | how the empirical width of the reference sequences enters σ (`hybrid`, `interpolated_fill`, `interpolated`, `nearest_two`) |
| `binary_fraction`, `binary_dm` | off | optional unresolved-binary mixture: with probability f_b the star is brighter than the sequence by `binary_dm` |
| `reference_ages` | adopted ages | new ages (Myr) for any reference, e.g. `{"beta_Pic": 21}`; the rest keep their defaults |
| `exclude_references` | none | references to leave out |
| `max_color`, `min_color` | 3.0, none | (BP−RP)₀ range that gets ages; outside it the outputs are NaN |
| `references` | bundled | a preloaded `load_references(...)` list, e.g. a custom library (`use_default_ages=False`) |

## Method

**Reference library.** De-extincted empirical sequences of six clusters and moving groups:

| reference | adopted age (Myr) |
|---|---|
| UBC 17a | 4.6 |
| β Pictoris moving group | 24.3 (lithium depletion boundary) |
| IC 2602 | 52.5 |
| α Persei (Melotte 20) | 79.0 |
| Pleiades (Melotte 22) | 127.4 |
| Coma Berenices (Melotte 111) | 675.0 |

For each star, the references whose sequences cover its colour are used (at least two are
needed). Their magnitudes at that colour are interpolated in log age with a shape-preserving
cubic (PCHIP), giving a model magnitude m(t) on a 1 Myr grid between reference ages.

**Likelihood.**

L(t) = Σᵢ wᵢ N(M_G; m_t(cᵢ), σ_t)

The sum is over Gauss–Hermite nodes cᵢ of the colour uncertainty, renormalised over the colour
range the references share.

σ_t² = σ_G² + s_t² + σ_sys²

Here s_t is the local width of the reference sequences: the normal-scaled MAD of member
residuals in a 0.5 mag colour window, interpolated in log age. σ_G is each star's own
absolute-G error.

In `hybrid` mode (the default), a star outside the reference envelope keeps one width at every
age. Otherwise an age-dependent s_t would make the likelihood favour the broadest reference
rather than the nearest one.
- **Fainter than every covering sequence:** the mean width of its two nearest sequences.
- **Brighter than the youngest covering sequence:** that sequence's width at the star's colour. That is
  UBC 17a at every colour that gets an age. Where UBC 17a has fewer than 8 members in the
  window, the width is taken from the nearest colour that has enough. It equals s_t at the youngest
  age, so σ is continuous across the sequence. Such stars get a flat posterior below the youngest
  age with a falling tail above it, instead of a spurious peak at older ages.

**Prior and age grid.** The prior is uniform in linear age. The posterior is evaluated between the
youngest and oldest covering reference. It is held flat at its edge values down to 0 Myr and up to
4000 Myr (5 Myr steps above 675 Myr), so a star outside the library gets a plateau rather than a
spike. MAP ages come from the fitted range only.

**Limits.** A star brighter than the youngest covering sequence gets `limit = "young"`
(age < that reference), and one fainter than the oldest gets `limit = "old"` (age > that reference).
Main-sequence stars are mostly "old" limits with broad posteriors: the method has no
leverage once a star has reached the main sequence.

**Posterior conventions.** `posterior` holds unit-sum **node weights** on `age_grid_myr`. The
grid mixes 1 Myr and 5 Myr steps, so use `cmd_interp.posterior_density(posterior, age_grid)`
to get a density per Myr before plotting or multiplying with other age likelihoods.

### Metallicity correction

The reference clusters are treated as solar. For a star of known [M/H], every model sequence is
shifted by the PARSEC v1.2S offset at fixed colour,

dM_G([M/H], t, c) = M_G^PARSEC([M/H], t, c) − M_G^PARSEC(0, t, c),

before the comparison. This "de-metallicities" the star: it is moved to the CMD position
it would have at solar metallicity. The shift depends on age, so it is applied at every trial
age rather than once. The offsets are tabulated for [M/H] −0.5 to +0.3, 5–1000 Myr and
(BP−RP)₀ 0.8–2.9 (`scripts/build_parsec_offsets.py`), and are interpolated linearly in
([M/H], log t, colour). Values outside the table are taken at its edge, flagged in
`metallicity_clipped`, and trigger a warning. Ages between 4.6 and 5 Myr use the 5 Myr offset.

At (BP−RP)₀ ≈ 2, a 0.1 dex metallicity error changes the age by roughly 5–25%. Using
a measured [M/H] removes that bias; it does not account for the uncertainty in [M/H]
itself, which the 0.05 mag systematic term partly absorbs.

## Validation

**Agreement with the original pipeline.** The default settings reproduce the original
research implementation exactly for every star not brighter than the youngest covering sequence
(`tests/test_regression.py`). Stars above that sequence differ by design, because of the bright-side
width described under Method.

On 40,000 TARS pre-main-sequence stars, 7.6% are young limits. For these, the median fraction of
posterior at or below the youngest age rises from 0.83 to 0.93. For stars more than 1 mag above
UBC 17a it rises from 0.08 to 1.00.

**Leave-one-out tests.** Each reference cluster was removed from the library and its photometric
PMS members (BP−RP ≥ 1.5) were re-dated:

| held-out cluster | median MAP age (Myr) | bias (dex) | scatter (dex) | 68% / 95% coverage |
|---|---|---|---|---|
| β Pic (24.3) | 21 | −0.06 | 0.20 | 32 / 76% |
| IC 2602 (52.5) | 56 | +0.03 | 0.12 | 64 / 86% |
| α Per (79.0) | 75 | −0.02 | 0.10 | 57 / 84% |
| Pleiades (127.4) | 282 | +0.35 | 0.56 | 19 / 56% |

The Pleiades fold fails because removing it leaves a 79–675 Myr gap, so ages between the
Pleiades and Coma Ber are poorly constrained by construction. UBC 17a and Coma Ber are the
ends of the library and cannot be held out.

## Limitations

- Only stars still above the main sequence carry real age information; the sequences
  converge for F/G stars older than ~50 Myr.
- Unresolved binaries look younger (brighter). Remove them first, or enable the binary mixture.
- Ages between 127 and 675 Myr rest on interpolation between the Pleiades and Coma Ber.
- Stars redder than (BP−RP)₀ = 3.0 get no age by default. Bluer than 0.8 the reference sequences start
  to cross (flagged and warned, see [Inputs](#inputs)). By BP−RP ≲ 0.5 they are not ordered by age at
  all, so ages and young/old limits there are not meaningful.
- The reference sequences, and so the ages, inherit the adopted cluster ages and the
  extinction corrections of the members.

## Data credits

- Cluster membership and Gaia DR3 photometry: Hunt & Reffert (2023, A&A 673, A114).
- Empirical cluster sequences: Rottensteiner & Meingast (2024).
- β Pic members: the Montreal Open Clusters and Associations database (MOCAdb; Gagné 2024). β Pic age: Binks & Jeffries (2014).
- Metallicity offsets: PARSEC v1.2S (Bressan et al. 2012) through the CMD 3.7 web interface.
- Extinction law: Fitzpatrick et al. (2019), Gaia EDR3 coefficients (bundled in `data/extinction/`).
- 3D dust maps: EXPLORE / Vergely et al. (2022) (not bundled).
- Photometric zero points: Riello et al. (2021). Parallax-error inflation for RUWE ≥ 1.25: Eq. 3 and Fig. 3 of arXiv:2206.05439.

`scripts/` documents how the bundled data were made: `slim_reference_data.py` and `build_parsec_offsets.py`.

## Citation

If you use this package, please cite it (see `CITATION.cff`) and the data sources above.

## License

MIT (see `LICENSE`).
