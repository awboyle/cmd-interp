# cmd-interp

**Stellar ages from a star's position in the Gaia colour–magnitude diagram (CMD), interpolated
between the empirical sequences of benchmark clusters rather than theoretical isochrones.**

The method is designed for pre-main-sequence stars, roughly 5–700 Myr old. K and M dwarfs are still
contracting toward the main sequence at these ages, so their CMD position is age-sensitive. For every
star the package returns a full age posterior, point estimates, credible intervals, and a flag for
stars that lie beyond the youngest or oldest reference sequence. The method, its calibration and its
limitations are described in the accompanying paper (Boyle et al., in preparation); this page covers
how to run the code.

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

The fit needs four numbers per star: the extinction-corrected colour (G_BP − G_RP)₀, the
extinction-corrected absolute magnitude M_G,0, and their uncertainties. The examples below start from
progressively less: those four numbers, then Gaia photometry plus an extinction, then a Gaia DR3
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
brighter than the youngest reference sequence (only an upper limit on the age), and `"old"` for one
fainter than the oldest.

For many stars pass arrays to `infer_ages`, or a table to `ages_for_table`:

```python
result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err)       # arrays, one entry per star
result["age_map_myr"], result["age_myr_50"], result["limit"]
result["posterior"]          # (n_star, n_age) unit-sum weights on result["age_grid_myr"]

table, result = cmd_interp.ages_for_table("stars.parquet")               # or a DataFrame, .csv, .feather, .fits
```

`ages_for_table` reads the columns `BpmRp0`, `phot_g_mean_mag_0`, `BpmRp0_err` and
`phot_g_mean_mag_0_err` by default (`color_col=` etc. to rename) and returns the table with
`cmd_age_*` columns added (MAP, median, 16th/84th and 2.5th/97.5th percentiles, limit, flags).

### 2. You have Gaia photometry and an extinction, but no uncertainties

Starting from the Gaia DR3 `gaia_source` columns of your stars (magnitudes, fluxes and their errors,
parallax and its error, RUWE), `deredden_gaia` computes the four fit inputs and *derives the
uncertainties for you* from the flux errors, the photometric zero points and the parallax error. The
only thing you supply is the extinction A_0 (at 550 nm):

```python
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12)                 # one A_0 for every star
phot = cmd_interp.deredden_gaia(gaia_table, a0="a0")                 # a per-star A_0 column
phot = cmd_interp.deredden_gaia(gaia_table, a0=0.12, a0_err=0.03)    # with an A_0 uncertainty, if you have one
table, result = cmd_interp.ages_for_table(phot)
```

`phot` is a copy of your table with `BpmRp0`, `phot_g_mean_mag_0` and their `_err` columns added,
ready for `ages_for_table`. If you only have source IDs, `ages_from_gaia` fetches the Gaia rows first:

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

Without the dust maps you must give A_0 yourself (example 2); A_0 = 0 is a fine approximation for
stars within a few tens of parsecs.

### 4. You know the metallicity

The reference sequences are solar-metallicity clusters. A metal-poor star sits below the solar sequence
at fixed colour, so it would be read as older than it is, and a metal-rich star as younger. If you have
[M/H] for the star, pass it and the sequences are shifted to that metallicity (using PARSEC isochrone
offsets) before the comparison:

```python
star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03, mh=-0.3)
print(star["age_map_myr"])                       # 34.3, against 74.5 at solar metallicity

star = cmd_interp.star_age_from_gaia(64895139073954944, a0="dustmap", mh=-0.1)    # 93 Myr; +0.1 gives 201 Myr
table, result = cmd_interp.ages_from_gaia(ids, a0="dustmap", mh=mh_values)          # one [M/H] per star
table, result = cmd_interp.ages_for_table(phot, mh_col="mh")                        # from a table column
```

Missing or NaN values are treated as solar. The offsets are tabulated for [M/H] from −0.5 to +0.3 and
(BP−RP)₀ from 0.8 to 2.9; outside that range the edge value is used and the star is flagged
`metallicity_clipped`. At (BP−RP)₀ ≈ 2 a 0.1 dex change in [M/H] moves the age by roughly 20%.

## Things to know

- **Colour range.** Stars redder than (BP−RP)₀ = 3.0 get no age (NaN, flagged `color_out_of_range`;
  change with `max_color=`). Bluer than 0.8 the reference sequences start to cross and ages there should
  be used with care: such stars are flagged `blue_caution` and raise a `cmd_interp.BlueColorWarning`
  (silence it with `warnings.filterwarnings("ignore", category=cmd_interp.BlueColorWarning)`).
- **Main-sequence stars.** Only stars still above the main sequence carry real age information. Once a
  star has settled on the main sequence it comes back as an old limit with a broad posterior.
- **Binaries.** Unresolved binaries are brighter and therefore look younger. Remove them beforehand, or
  enable the binary mixture (`binary_fraction=`, `binary_dm=`).
- **Posteriors** are unit-sum node weights on an age grid that mixes 1 Myr and 5 Myr steps; use
  `cmd_interp.posterior_density(posterior, age_grid)` to get a density per Myr before plotting or
  combining with other age likelihoods.

## Command line

```bash
cmd-interp star --color 2.0 --abs-g 8.0 --color-err 0.01 --abs-g-err 0.03 --mh 0.1 --plot star.png
cmd-interp table stars.parquet ages.parquet --mh-col mh --posterior-output posteriors.npy
cmd-interp deredden gaia.csv phot.csv --a0 0.12 --a0-err 0.03      # or --a0-col, or --dustmap
cmd-interp gaia 6794047652729201024 --a0 0 --plot aumic.png         # also --ids-file, --dustmap, --output
cmd-interp table --help
```

`--posterior-output` writes the posteriors as a float32 array, the age grid, and an index giving the
table row of each posterior. They are written as they are computed, so memory stays bounded for large
catalogues (about 5 GB per million stars).

## Plots

```python
from cmd_interp.plotting import plot_star
fig, axes = plot_star(cmd_interp.star_age(2.0, 8.0, 0.01, 0.03))   # CMD position + posterior
```

`examples/quickstart.py` runs the four quick-start examples.

## Reference clusters and options

The ages are interpolated between the empirical sequences of six clusters and moving groups:

| reference | default age (Myr) |
|---|---|
| UBC 17a | 4.6 |
| β Pictoris moving group | 24.3 |
| IC 2602 | 52.5 |
| α Persei (Melotte 20) | 79.0 |
| Pleiades (Melotte 22) | 127.4 |
| Coma Berenices (Melotte 111) | 675.0 |

To use a different age for any of them, pass only the ones you want to change; names are matched
loosely (`beta_Pic`, `BPMG`, `Pleiades`, `Melotte_22`, `alpha Per`, `Coma Ber` all work):

```python
cmd_interp.star_age(2.0, 8.0, 0.01, 0.03, reference_ages={"beta_Pic": 21.0})
cmd_interp.ages_for_table(table, reference_ages={"Pleiades": 112, "alpha Per": 85})
cmd_interp.infer_ages(..., exclude_references=["Coma Ber"])     # leave a reference out
```

Other options accepted by every entry point:

| argument | default | meaning |
|---|---|---|
| `systematic_mag` | 0.05 | extra magnitude scatter added in quadrature at every age |
| `mh` | solar | stellar [M/H] |
| `binary_fraction`, `binary_dm` | off | unresolved-binary mixture: with probability f_b the star is brighter than the sequence by `binary_dm` |
| `reference_ages`, `exclude_references` | defaults | change or drop reference clusters |
| `max_color`, `min_color` | 3.0, none | (BP−RP)₀ range that receives ages |
| `scatter_mode` | `"hybrid"` | how the width of the reference sequences enters the likelihood |

How the likelihood, the sequence widths, the prior and the flat extension outside the library are
defined, and how the ages compare with independent cluster ages, is in the paper.

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

If you use this package, please cite the paper and the package (see `CITATION.cff`) together with the
data sources above.

## License

MIT (see `LICENSE`).
