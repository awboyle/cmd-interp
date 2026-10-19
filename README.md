# cmd-interp

**Stellar ages from a star's position in the Gaia colour–magnitude diagram (CMD), interpolated
between the empirical sequences of benchmark clusters rather than theoretical isochrones.**

The method is designed for pre-main-sequence stars, roughly 5–700 Myr old. K and M dwarfs are still
contracting toward the main sequence at these ages, so their CMD position is age-sensitive.
For every star the package returns a full age posterior. It also returns point estimates, credible
intervals, and flags for stars that lie beyond the youngest or oldest reference sequence.

```python
import cmd_interp

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
star["age_map_myr"], star["age_myr_16"], star["age_myr_84"]      # (74.5, 62.7, 83.9)

star = cmd_interp.star_age_from_gaia(6794047652729201024, a0="dustmap")   # AU Mic, from its Gaia DR3 ID
```

## Installation

```bash
pip install git+https://github.com/awboyle/cmd-interp
```

Optional extras:

```bash
pip install "cmd-interp[gaia] @ git+https://github.com/awboyle/cmd-interp"   # astroquery, to start from Gaia DR3 IDs
pip install "cmd-interp[dust] @ git+https://github.com/awboyle/cmd-interp"   # h5py + astropy, for the 3D dust maps
pip install "cmd-interp[plot] @ git+https://github.com/awboyle/cmd-interp"   # matplotlib, for the diagnostic figures
```

Requires Python ≥ 3.10, numpy, scipy, pandas and pyarrow. The reference sequences (≈300 kB) are included.

## Quick start

The fit itself needs four numbers per star: the extinction-corrected colour (G_BP − G_RP)₀, the
extinction-corrected absolute magnitude M_G,0, and their uncertainties. The examples below start from
progressively less: first those four numbers, then Gaia photometry plus an extinction, then a Gaia DR3
source ID alone, and finally the same with a known metallicity.

### 1. You have de-reddened photometry and uncertainties

```python
import cmd_interp

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
print(star["age_map_myr"], star["age_myr_16"], star["age_myr_84"], star["limit"])
# 74.5 62.7 83.9 ''
```

`star` also holds the full posterior (`star["posterior"]` on `star["age_grid_myr"]`), the median and
the 2.5/97.5 percentiles. `limit` is `""` for a star inside the reference library, `"young"` for one
brighter than the youngest reference sequence (so only an upper limit on the age), and `"old"` for one
fainter than the oldest.

For many stars pass arrays to `infer_ages`, or a table to `ages_for_table`:

```python
result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err)       # arrays, one entry per star
result["age_map_myr"], result["age_myr_50"], result["limit"]
result["posterior"]          # (n_star, n_age) unit-sum weights on result["age_grid_myr"]

table, result = cmd_interp.ages_for_table("stars.parquet")               # or a DataFrame, .csv, .feather, .fits
```

`ages_for_table` reads the columns `BpmRp0`, `phot_g_mean_mag_0`, `BpmRp0_err` and
`phot_g_mean_mag_0_err` by default (`color_col=` etc. to rename), and returns the table with
`cmd_age_*` columns added (see [Tables](#tables)).

### 2. You have Gaia photometry and an extinction, but no uncertainties

Starting from the Gaia DR3 `gaia_source` columns of your stars (magnitudes, fluxes and their errors,
parallax and its error, RUWE), `deredden_gaia` computes the four fit inputs and *derives the
uncertainties for you*. It propagates the flux errors, the photometric zero points, and the parallax
error (inflated as a function of RUWE and magnitude) into both the absolute magnitude and the
reddening. The only thing you supply is the extinction A_0 (at 550 nm):

```python
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12)                 # one A_0 for every star
phot = cmd_interp.deredden_gaia(gaia_table, a0="a0")                 # a per-star A_0 column
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12, a0_err=0.03)    # with an A_0 uncertainty, if you have one
table, result = cmd_interp.ages_for_table(phot)
```

`phot` is a copy of your table with the columns `BpmRp0`, `phot_g_mean_mag_0` and their `_err`
columns added (plus `distance`, `extinction_a0` and the inflated parallax error), ready for
`ages_for_table`. If you only have source IDs, `ages_from_gaia` fetches the Gaia rows first:

```python
table, result = cmd_interp.ages_from_gaia([64895139073954944, 6794047652729201024], a0=0.12)
star = cmd_interp.star_age_from_gaia(64895139073954944, a0=0.12)     # a Pleiades M dwarf
print(star["BpmRp0"], star["phot_g_mean_mag_0"], star["phot_g_mean_mag_0_err"], star["age_map_myr"])
# 2.463 9.371 0.010 117.0
```

### 3. You have only a Gaia DR3 source ID

With the EXPLORE 3D dust maps, one source ID is enough: the package queries the Gaia archive
(anonymously), integrates the dust map to the star's distance to get A_0 and its uncertainty, corrects
the photometry, propagates every uncertainty, and fits the age. This needs the `gaia` and `dust`
extras and a one-time download of the two dust cubes (~230 MB) from the
[EXPLORE G-Tomo service](https://explore-platform.eu): `explore_cube_density_values_010pc_v2.h5` and
`explore_cube_density_errors_050pc_v2.h5`. Put them in a folder and point the package at it:

```bash
export CMD_INTERP_DUST_DIR=/path/to/explore_cubes
```

```python
star = cmd_interp.star_age_from_gaia(64895139073954944, a0="dustmap")
print(star["extinction_a0"], star["extinction_a0_err"], star["age_map_myr"], star["age_myr_16"])
# 0.165 0.003 119.0 130.1
```

Or pass the cube paths explicitly:

```python
dust = cmd_interp.ExploreDustMap("explore_cube_density_values_010pc_v2.h5",
                                 "explore_cube_density_errors_050pc_v2.h5")
table, result = cmd_interp.ages_from_gaia(ids, a0="dustmap", dust_map=dust)
```

The same works from the command line:

```bash
cmd-interp gaia 64895139073954944 --dustmap --plot pleiades_star.png
```

The cubes cover about 3 kpc around the Sun; beyond their edge A_0 is held at the edge value. Without
the dust maps you must give A_0 yourself (example 2); A_0 = 0 is a fine approximation for stars
within a few tens of parsecs.

### 4. You know the metallicity

The reference sequences are solar-metallicity clusters. A metal-poor star sits below the solar sequence
at fixed colour, so it would be read as older than it is, and a metal-rich star as younger. If you have
[M/H] for the star, pass it and the sequences are shifted to that metallicity before the comparison
(see [Metallicity](#metallicity) for how):

```python
star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03, mh=-0.3)
print(star["age_map_myr"])                       # 34.3, against 74.5 at solar metallicity

star = cmd_interp.star_age_from_gaia(64895139073954944, a0="dustmap", mh=-0.1)    # 93 Myr; +0.1 gives 201 Myr
table, result = cmd_interp.ages_from_gaia(ids, a0="dustmap", mh=mh_values)          # one [M/H] per star
table, result = cmd_interp.ages_for_table(phot, mh_col="mh")                        # from a table column
```

Missing or NaN values are treated as solar. The offsets are tabulated for [M/H] from −0.5 to +0.3 and
(BP−RP)₀ from 0.8 to 2.9; outside that range the edge value is used and the star is flagged
`metallicity_clipped`. At (BP−RP)₀ ≈ 2 a 0.1 dex change in [M/H] moves the age by roughly 20%, so for
young K and M dwarfs the metallicity matters whenever it is known to better than that.

## Inputs

| quantity | meaning |
|---|---|
| `color` | (G_BP − G_RP)₀ |
| `abs_g` | absolute, extinction-corrected G magnitude, M_G,0 |
| `color_err`, `abs_g_err` | 1σ uncertainties in mag; `abs_g_err` should include the parallax error |
| `mh` (optional) | stellar [M/H] in dex; solar by default |

**Colour range.**
- **Redder than (BP−RP)₀ = 3.0:** no age (NaN), flagged `color_out_of_range`. The limit can be changed
  with `max_color=`; `min_color=` sets an optional blue limit.
- **Bluer than (BP−RP)₀ = 0.8:** the reference sequences start to cross and are no longer ordered by age.
  Ages are still returned, but such stars trigger a `cmd_interp.BlueColorWarning` and are flagged
  `blue_caution`. To silence the warning:
  `warnings.filterwarnings("ignore", category=cmd_interp.BlueColorWarning)`.

## Usage in more detail

### Tables

`ages_for_table` takes a DataFrame, or a path to a .csv, .feather, .parquet or .fits file:

```python
table, result = cmd_interp.ages_for_table(
    "stars.parquet",
    color_col="BpmRp0", abs_g_col="phot_g_mean_mag_0",          # these are the defaults
    color_err_col="BpmRp0_err", abs_g_err_col="phot_g_mean_mag_0_err",
    mh_col="mh",                  # optional per-star [M/H]; missing values are treated as solar
    fit_mask="young_candidate",   # optional boolean column: only fit these rows
)
```

The returned table gains these columns:
- `cmd_age_myr`: MAP age, or the reference age for a limit;
- `cmd_age_median_myr`;
- `cmd_age_lower_myr` / `cmd_age_upper_myr`: 16th/84th percentiles;
- `cmd_age_p2p5_myr` / `cmd_age_p97p5_myr`;
- `cmd_age_limit`;
- `cmd_age_n_references`;
- `cmd_age_mh` and `cmd_age_metallicity_clipped`;
- `cmd_age_color_out_of_range` and `cmd_age_blue_caution`;
- `cmd_age_fitted`.

`result["row"]` gives the table row of each posterior.

### Extinction correction

`deredden_gaia` converts Gaia DR3 `gaia_source` columns into absolute, extinction-corrected photometry
with uncertainties. Its output has exactly the column names `ages_for_table` reads by default.
- **Required columns:** `parallax`, `parallax_error`, `ruwe`, `phot_{g,bp,rp}_mean_mag`,
  `phot_{g,bp,rp}_mean_flux` and `phot_{g,bp,rp}_mean_flux_error`.
- **For the dust map:** also `ra` and `dec`.

```python
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12, a0_err=0.03)   # one A_0 for every star (default 0)
phot = cmd_interp.deredden_gaia(gaia_table, a0="a0_column")          # per-star A_0 from a column
phot = cmd_interp.deredden_gaia(gaia_table, a0="dustmap", dust_map=dust)   # A_0 from the EXPLORE 3D dust maps
```

- **Extinction law.** A_0 (the extinction at 550 nm) is converted to A_G, A_BP and A_RP with the
  Fitzpatrick et al. (2019) Gaia EDR3 coefficients k_X((BP−RP)₀, A_0). The intrinsic colour is iterated
  to convergence. Stars brighter than G − 5 log₁₀(d / 10 pc) = 5 use the upper-HRD coefficients.
- **3D dust maps.** `a0="dustmap"` integrates the EXPLORE extinction-density cube (Vergely et al. 2022)
  along the line of sight to d = 1000/parallax.
  - It also integrates the density-error cube. This treats the voxel errors as fully correlated, so it
    is an upper bound on the dust-map uncertainty.
  - It records the local extinction density, so the distance error can be propagated into A_0.
  - The cubes are not included: download them from the [EXPLORE](https://explore-platform.eu) G-Tomo
    service. Pass their paths to `ExploreDustMap`, set `CMD_INTERP_DUST_DIR` to the folder that holds
    them, or set `CMD_INTERP_EXPLORE_VALUES` and `CMD_INTERP_EXPLORE_ERRORS` to the two files.
- **Uncertainties.** These combine:
  - the flux errors;
  - the Gaia photometric zero points (Riello et al. 2021);
  - the parallax error, inflated as a function of RUWE and G and propagated through both the
    distance modulus and the distance dependence of A_0;
  - the A_0 uncertainty.

  In (BP−RP)₀ the distance modulus cancels, so only the reddening E(BP−RP) carries the A_0 error. The
  full equations are in the `cmd_interp.photometry` docstring.
- **Not included:** the Gaia parallax zero point, a photometric error floor, and uncertainty in the
  extinction law.
- **Bad parallaxes:** stars with parallax ≤ 0 get NaN.

### From Gaia DR3 source IDs

`ages_from_gaia` runs the whole chain:
1. Queries `gaiadr3.gaia_source` in the ESA Gaia archive (anonymously).
2. Corrects for extinction and propagates the uncertainties.
3. Fits the ages.

```python
table, result = cmd_interp.ages_from_gaia(ids, a0=0.1)                       # one A_0 for all stars
table, result = cmd_interp.ages_from_gaia(ids, a0=a0_values, mh=mh_values)    # one A_0 and [M/H] per star
table, result = cmd_interp.ages_from_gaia(ids, a0="dustmap", dust_map=dust)   # A_0 from the 3D dust maps
star = cmd_interp.star_age_from_gaia(6794047652729201024, a0=0.0)             # a single star
```

- **A_0 is required:** a value, one value per ID, or `"dustmap"`.
- **Row order:** rows follow the order of the IDs. IDs not found in DR3 keep their row, with
  `found_in_gaia = False` and NaN ages.
- **Other options** (`systematic_mag`, `reference_ages`, `max_color`, …) are passed to the age fit.
- **Your own Gaia rows:** pass them with `gaia_table=` to skip the archive query.
- **Logging in:** `credentials_file=` logs in with an astroquery credentials file.

### Command line

```bash
cmd-interp star --color 2.0 --abs-g 8.0 --color-err 0.01 --abs-g-err 0.03 --mh 0.1 --plot star.png
cmd-interp table stars.parquet ages.parquet --mh-col mh --posterior-output posteriors.npy
cmd-interp table stars.csv ages.csv --reference-age beta_Pic=21 --exclude-reference Coma_Ber
cmd-interp deredden gaia.csv phot.csv --a0 0.12 --a0-err 0.03      # or --a0-col, or --dustmap
cmd-interp gaia 6794047652729201024 --a0 0 --plot aumic.png         # also --ids-file, --dustmap, --output
cmd-interp table --help
```

`--posterior-output` writes three files:
- the posteriors, as a float32 array (`posteriors.npy`);
- the age grid (`posteriors_age_grid.npy`);
- an index giving the table row of each posterior (`posteriors_index.csv`).

The posteriors are written to disk as they are computed, so memory stays bounded for large catalogues
(about 5 GB per million stars).

### Plots

```python
from cmd_interp.plotting import plot_star
fig, axes = plot_star(cmd_interp.star_age(2.0, 8.0, 0.01, 0.03))   # CMD position + posterior
```

`examples/quickstart.py` runs the four quick-start examples.

### Changing the reference ages

The default ages are listed under [Method](#method). To use a different age for any reference cluster,
pass only the ones you want to change; the others keep their defaults:

```python
cmd_interp.star_age(2.0, 8.0, 0.01, 0.03, reference_ages={"beta_Pic": 21.0})
cmd_interp.ages_for_table(table, reference_ages={"Pleiades": 112, "alpha Per": 85})
cmd_interp.infer_ages(..., exclude_references=["Coma Ber"])     # leave a reference out
cmd_interp.resolve_reference_ages({"IC 2602": 46})             # inspect the resulting set
```

Names are matched loosely, so `beta_Pic`, `beta pic`, `BPMG`, `Pleiades`, `Melotte_22`, `alpha Per`,
`Coma Ber` and `IC 2602` all work. For many calls, load the library once with
`refs = cmd_interp.load_references(reference_ages=...)` and pass `references=refs`.

## Settings

| argument | default | meaning |
|---|---|---|
| `systematic_mag` | 0.05 | extra magnitude scatter added in quadrature at every age (field metallicity spread, zero point of the reference sequences) |
| `mh` | solar | stellar [M/H]; see [Metallicity](#metallicity) |
| `scatter_mode` | `"hybrid"` | how the intrinsic width of the reference sequences enters the likelihood (`hybrid`, `interpolated_fill`, `interpolated`, `nearest_two`) |
| `binary_fraction`, `binary_dm` | off | optional unresolved-binary mixture: with probability f_b the star is brighter than the sequence by `binary_dm` |
| `reference_ages` | see below | new ages (Myr) for any reference cluster |
| `exclude_references` | none | reference clusters to leave out |
| `max_color`, `min_color` | 3.0, none | (BP−RP)₀ range that receives ages |
| `references` | bundled | a preloaded `load_references(...)` list, or a custom library (`use_default_ages=False`) |

## Method

### Reference sequences

The package uses de-extincted empirical CMD sequences of six clusters and moving groups:

| reference | default age (Myr) |
|---|---|
| UBC 17a | 4.6 |
| β Pictoris moving group | 24.3 (lithium depletion boundary) |
| IC 2602 | 52.5 |
| α Persei (Melotte 20) | 79.0 |
| Pleiades (Melotte 22) | 127.4 |
| Coma Berenices (Melotte 111) | 675.0 |

For each star, the reference sequences that cover its colour are used (at least two are needed). Their
magnitudes at that colour are interpolated in log age with a shape-preserving cubic (PCHIP). This gives
a model magnitude m(t) on a 1 Myr grid between the reference ages.

### Likelihood

L(t) = Σᵢ wᵢ N(M_G; m_t(cᵢ), σ_t)

The sum runs over Gauss–Hermite nodes cᵢ of the colour uncertainty, renormalised over the colour range
the sequences share. The width is

σ_t² = σ_G² + s_t² + σ_sys²

- σ_G is the star's absolute-magnitude uncertainty.
- σ_sys is `systematic_mag`.
- s_t is the intrinsic width of the reference sequences: a robust (normal-scaled MAD) scatter of the
  cluster members in a 0.5 mag colour window, interpolated in log age.

In the default `hybrid` mode, a star outside the envelope of the reference sequences keeps a single
width at every age. Otherwise the likelihood would favour whichever reference sequence is broadest,
rather than the nearest one.
- **Fainter than every covering sequence:** the mean width of the two nearest sequences.
- **Brighter than the youngest covering sequence:** the width of that sequence. This matches s_t at the
  youngest age, so σ is continuous across the sequence. Such stars get a flat posterior below the
  youngest age.

### Prior, age grid and limits

- **Prior.** The prior is uniform in linear age.
- **Fitted range.** The likelihood is evaluated between the youngest and oldest reference sequence covering
  the star.
- **Outside that range.** The posterior is held flat at its edge values down to 0 Myr and up to 4000 Myr
  (5 Myr steps above 675 Myr). A star beyond the library therefore gets a plateau rather than a spike.
  MAP ages are taken from the fitted range only.
- **Limits.** A star brighter than the youngest covering sequence is flagged `limit = "young"` (younger than
  that reference). One fainter than the oldest is flagged `limit = "old"` (older than that reference).
- **Main-sequence stars** are mostly old limits with broad posteriors. The method has little leverage once
  a star has reached the main sequence.

### Posterior conventions

`posterior` holds unit-sum **node weights** on `age_grid_myr`. The grid mixes 1 Myr and 5 Myr steps, so
use `cmd_interp.posterior_density(posterior, age_grid)` to get a density per Myr before plotting or
combining with other age likelihoods.

### Metallicity

The reference clusters are treated as solar-metallicity. For a star of known [M/H], every model sequence
is first shifted by the offset predicted by PARSEC v1.2S isochrones at fixed colour:

dM_G([M/H], t, c) = M_G^PARSEC([M/H], t, c) − M_G^PARSEC(0, t, c)

The offset is positive (fainter) for metal-poor stars and negative for metal-rich ones; at
(BP−RP)₀ = 2 and 100 Myr it is +0.30 mag at [M/H] = −0.2 and −0.21 mag at +0.2. Shifting the
sequences is equivalent to moving the star to the CMD position it would have at solar metallicity.
Because the shift depends on age, it is applied at every trial age rather than once, and the same
shifted sequences decide whether the star counts as above or below the library.
- **Tabulated range:** [M/H] −0.5 to +0.3, 5–1000 Myr and (BP−RP)₀ 0.8–2.9, interpolated linearly in
  [M/H], log age and colour.
- **Outside the table:** values are taken at the table edge, flagged `metallicity_clipped`, and trigger a
  warning.

At (BP−RP)₀ ≈ 2, a 0.1 dex metallicity error changes the inferred age by roughly 20%. Supplying [M/H]
removes this bias. The uncertainty of the [M/H] value itself is not propagated, but the systematic term
partly absorbs it.

## Limitations

- **Main-sequence stars:** only stars still above the main sequence carry real age information. The
  sequences converge for F and G stars older than about 50 Myr, and for K dwarfs beyond the Pleiades,
  so an M dwarf's posterior at a few hundred Myr can have a long tail toward the old plateau.
- **Binaries:** unresolved binaries are brighter and therefore look younger. Remove them beforehand, or
  enable the binary mixture.
- **127–675 Myr:** ages here rest on interpolation between the Pleiades and Coma Ber.
- **Colour range:** stars redder than (BP−RP)₀ = 3.0 get no age by default. Bluer than 0.8 the sequences
  start to cross, and bluer than about 0.5 they are not ordered by age at all.
- **Inherited systematics:** the ages carry the adopted reference-cluster ages and the extinction
  corrections of their members.

## Data and credits

- **Cluster membership and Gaia DR3 photometry:** Hunt & Reffert (2023, A&A 673, A114).
- **Empirical cluster sequences:** Rottensteiner & Meingast (2024).
- **β Pic members:** the Montreal Open Clusters and Associations database (MOCAdb; Gagné 2024).
- **β Pic age:** Binks & Jeffries (2014).
- **Metallicity offsets:** PARSEC v1.2S (Bressan et al. 2012), via the CMD 3.7 web interface.
- **Extinction law:** Fitzpatrick et al. (2019) Gaia EDR3 coefficients (included).
- **3D dust maps:** EXPLORE, Vergely et al. (2022) (not included).
- **Photometric zero points:** Riello et al. (2021).
- **Parallax-uncertainty inflation:** arXiv:2206.05439.

## Development

```bash
git clone https://github.com/awboyle/cmd-interp
cd cmd-interp
pip install -e ".[dev]"
pytest
```

## Citation

If you use this package, please cite it (see `CITATION.cff`) together with the data sources above.

## License

MIT (see `LICENSE`).
