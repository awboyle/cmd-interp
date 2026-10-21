# cmd-interp

**Stellar ages from a star's position in the Gaia color–magnitude diagram (CMD), interpolated
between the empirical sequences of benchmark clusters rather than theoretical isochrones.**

The method is designed for pre-main-sequence stars, roughly 5–700 Myr old. K and M dwarfs are still
contracting toward the main sequence at these ages, so their CMD position is age-sensitive. For every
star the package returns a full age posterior, point estimates, credible intervals, and a flag for
stars that lie beyond the youngest or oldest reference sequence. The method, its calibration and its
limitations are described in the accompanying paper (Boyle et al. in prep., "The Fraction of Bound Stars as a Function of Age within 500 pc"); this page covers
how to run the code.


```python
import cmd_interp

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
print(star["age_map_myr"], star["age_myr_16"], star["age_myr_84"])
# 74.5 62.7 83.9
```

## Installation

```bash
pip install git+https://github.com/awboyle/cmd-interp
```

Optional extras:

```bash
pip install "cmd-interp[gaia] @ git+https://github.com/awboyle/cmd-interp"   # astroquery, to start from Gaia DR3 IDs
pip install "cmd-interp[dust] @ git+https://github.com/awboyle/cmd-interp"   # astropy, for the 3D dust maps
pip install "cmd-interp[plot] @ git+https://github.com/awboyle/cmd-interp"   # matplotlib, for the diagnostic figures
```

Requires Python ≥ 3.10, numpy, scipy, pandas and pyarrow. The reference sequences (≈300 kB) are included.

## Quick start

The fit needs four numbers per star: the extinction-corrected color (G_BP − G_RP)₀, the
extinction-corrected absolute magnitude M_G,0, and their uncertainties. The examples below start from
progressively less: those four numbers, then Gaia photometry plus an extinction, then a Gaia DR3
source ID alone, and finally the same with a known metallicity. Every code block is self-contained
and can be copied and run as it is. Examples that query the Gaia archive need the `gaia` extra and a
network connection; the ones marked *dust maps* also need the cubes described in example 3.

### 1. You have de-reddened photometry and uncertainties

One star:

```python
import cmd_interp

color = 2.0          # (G_BP - G_RP)_0, extinction-corrected
abs_g = 8.0          # M_G,0, absolute and extinction-corrected
color_err = 0.01     # 1-sigma uncertainty of the color (mag)
abs_g_err = 0.03     # 1-sigma uncertainty of the absolute magnitude (mag)

star = cmd_interp.star_age(color=color, abs_g=abs_g, color_err=color_err, abs_g_err=abs_g_err)

print(star["age_map_myr"])                            # 74.5    most probable age (Myr)
print(star["age_myr_16"], star["age_myr_84"])         # 62.7 83.9    68% credible interval
print(star["age_myr_50"])                             # 73.7    median
print(star["age_myr_025"], star["age_myr_975"])       # 52.1 92.8    95% credible interval
print(star["limit"])                                  # ''      inside the reference library
print(star["posterior"].shape, star["age_grid_myr"][:3])   # (1343,) [0.6 1.6 2.6]
```

`star` is a dict with these keys:

| key | meaning |
|---|---|
| `age_map_myr` | most probable age: the peak of the posterior (Myr) |
| `age_myr_16`, `age_myr_84` | 16th and 84th percentiles of the posterior, the 68% credible interval (Myr) |
| `age_myr_50` | median of the posterior (Myr) |
| `age_myr_025`, `age_myr_975` | 2.5th and 97.5th percentiles, the 95% credible interval (Myr) |
| `limit` | `""` inside the reference library; `"young"` if the star is brighter than the youngest reference sequence (the age is only an upper limit); `"old"` if fainter than the oldest |
| `posterior` | the full posterior as unit-sum weights on the age grid, shape `(1343,)` |
| `age_grid_myr` | the age grid (Myr), shape `(1343,)`, from 0.6 to 4000 Myr in mixed 1 Myr and 5 Myr steps |
| `has_posterior` | `True` when an age was computed; `False` if fewer than two reference sequences reach the star's color, or the color is out of range |
| `n_references` | number of reference sequences that cover the star's color (at most 6) |
| `youngest_age`, `oldest_age` | ages (Myr) of the youngest and oldest reference sequences covering the star |
| `youngest_mag`, `oldest_mag` | M_G of those two sequences at the star's color, i.e. the bright and faint edges of the library |
| `scatter_sigma_mag` | likelihood width (mag) used for a star outside the library; NaN inside it, where the width varies with age |
| `mh` | the [M/H] used in the fit (0 = solar) |
| `metallicity_clipped` | `True` if [M/H] or the color fell outside the tabulated PARSEC offsets and the edge value was used |
| `color_out_of_range` | `True` if the color is outside `min_color`–`max_color` (default 3.0); no age is returned |
| `blue_caution` | `True` for (BP−RP)₀ < 0.8, where the reference sequences start to cross and ages are less reliable |
| `color`, `abs_g`, `color_err`, `abs_g_err` | the inputs, echoed back |

Many stars, as arrays:

```python
import numpy as np
import cmd_interp

color = np.array([1.2, 2.0, 2.6])            # (G_BP - G_RP)_0 of each star
abs_g = np.array([5.5, 8.0, 10.4])           # M_G,0 of each star
color_err = np.array([0.01, 0.01, 0.02])     # color uncertainties (mag)
abs_g_err = np.array([0.03, 0.03, 0.05])     # magnitude uncertainties (mag)

result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err)

print(result["age_map_myr"])      # [ 21.6  74.5 675. ]
print(result["age_myr_16"])       # [ 16.9  62.7 926.3]
print(result["age_myr_84"])       # [ 29.3  83.9 3414.5]
print(result["limit"])            # ['' '' 'old']    the third star is fainter than the oldest sequence
print(result["posterior"].shape)  # (3, 1343)        one row per star on result["age_grid_myr"]
```

`result` has the same keys as `star_age` (except the echoed inputs), each now an array:

| key | meaning |
|---|---|
| `age_map_myr`, `age_myr_16`, `age_myr_84`, `age_myr_50`, `age_myr_025`, `age_myr_975` | the point estimates and percentiles of the single-star case, one per star, shape `(n_stars,)`; NaN where no age could be computed |
| `limit` | `""`, `"young"` or `"old"` per star |
| `posterior` | unit-sum weights, shape `(n_stars, 1343)`, one row per star |
| `age_grid_myr` | the common age grid (Myr), shape `(1343,)` |
| `has_posterior`, `n_references`, `youngest_age`, `oldest_age`, `youngest_mag`, `oldest_mag`, `scatter_sigma_mag`, `mh`, `metallicity_clipped`, `color_out_of_range`, `blue_caution` | as for `star_age`, one entry per star |

Many stars, as a table. `ages_for_table` reads the columns `BpmRp0`, `phot_g_mean_mag_0`,
`BpmRp0_err` and `phot_g_mean_mag_0_err` by default (`color_col=`, `abs_g_col=`, `color_err_col=`,
`abs_g_err_col=` to rename them) and returns a copy of the table with `cmd_age_*` columns added:

```python
import pandas as pd
import cmd_interp

stars = pd.DataFrame({
    "name": ["star A", "star B", "star C"],
    "BpmRp0": [1.2, 2.0, 2.6],                   # (G_BP - G_RP)_0
    "phot_g_mean_mag_0": [5.5, 8.0, 10.4],       # M_G,0
    "BpmRp0_err": [0.01, 0.01, 0.02],            # color uncertainty (mag)
    "phot_g_mean_mag_0_err": [0.03, 0.03, 0.05], # magnitude uncertainty (mag)
})

table, result = cmd_interp.ages_for_table(stars)

print(table[["name", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr", "cmd_age_limit"]].round(1))
#      name  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr cmd_age_limit
# 0  star A         21.6               16.9               29.3
# 1  star B         74.5               62.7               83.9
# 2  star C        675.0              675.0                NaN           old
```

`table` is a copy of the input with these columns added:

| column | meaning |
|---|---|
| `cmd_age_myr` | the adopted age (Myr): the most probable age, or at a limit the age of the edge reference sequence (4.6 Myr for `young`, 675 Myr for `old`) |
| `cmd_age_lower_myr` | 16th percentile (Myr); for an `old` limit the oldest reference age, NaN for a `young` limit |
| `cmd_age_upper_myr` | 84th percentile (Myr); for a `young` limit the youngest reference age, NaN for an `old` limit |
| `cmd_age_median_myr` | median (Myr); NaN at a limit |
| `cmd_age_p2p5_myr`, `cmd_age_p97p5_myr` | 2.5th and 97.5th percentiles (Myr); NaN at a limit |
| `cmd_age_limit` | `""`, `"young"` or `"old"`, as `limit` above |
| `cmd_age_n_references` | number of reference sequences covering the star's color |
| `cmd_age_scatter_mag` | likelihood width (mag) for a star outside the library; NaN inside it |
| `cmd_age_mh` | the [M/H] used (0 = solar) |
| `cmd_age_metallicity_clipped` | `True` if [M/H] or the color was outside the tabulated offsets |
| `cmd_age_fitted` | `True` if a posterior was computed for the row; `False` for rows left out by `fit_mask`, with non-finite photometry, out of the color range, or with fewer than two covering sequences |
| `cmd_age_color_out_of_range` | `True` if the color is outside `min_color`–`max_color` |
| `cmd_age_blue_caution` | `True` for (BP−RP)₀ < 0.8 |

The second return value `result` is the `infer_ages` output for the fitted rows only (the arrays of the
previous table, including `result["posterior"]` and `result["age_grid_myr"]`), plus `result["row"]`,
the row index in `table` of each fitted star. Instead of a DataFrame you can pass a file path (`.csv`,
`.feather`, `.parquet` or `.fits`):

```python
import cmd_interp

table, result = cmd_interp.ages_for_table("stars.parquet")
table.to_parquet("stars_with_ages.parquet")
```

### 2. You have Gaia photometry and an extinction, but no uncertainties

Starting from the Gaia DR3 `gaia_source` columns of your stars (magnitudes, fluxes and their errors,
parallax and its error, RUWE), `deredden_gaia` computes the four fit inputs and *derives the
uncertainties for you* from the flux errors, the photometric zero points and the parallax error. The
only thing you supply is the extinction A_0 (at 550 nm). The table below holds the real `gaia_source`
rows of AU Mic and a Pleiades M dwarf, so this runs without a network connection:

```python
import pandas as pd
import cmd_interp

gaia_table = pd.DataFrame({
    "dr3_source_id":           [6794047652729201024, 64895139073954944],   # AU Mic, a Pleiades M dwarf
    "parallax":                [102.9432, 7.309302],                       # mas
    "parallax_error":          [0.02311219, 0.03100616],
    "ruwe":                    [0.9273568, 1.100511],
    "phot_g_mean_mag":         [7.843440, 15.13400],
    "phot_bp_mean_mag":        [8.909956, 16.50113],
    "phot_rp_mean_mag":        [6.798632, 13.98515],
    "phot_g_mean_flux":        [1.372685e7, 1.664744e4],
    "phot_g_mean_flux_error":  [1.100688e4, 22.51610],
    "phot_bp_mean_flux":       [3.727644e6, 3427.412],
    "phot_bp_mean_flux_error": [8432.241, 29.77980],
    "phot_rp_mean_flux":       [1.512535e7, 20188.24],
    "phot_rp_mean_flux_error": [1.871785e4, 83.50732],
})
a0 = 0.12            # extinction A_0 at 550 nm (mag), here the same value for every star

phot = cmd_interp.deredden_gaia(gaia_table, a0=a0)
table, result = cmd_interp.ages_for_table(phot)

print(table[["dr3_source_id", "BpmRp0", "BpmRp0_err", "phot_g_mean_mag_0", "phot_g_mean_mag_0_err", "cmd_age_myr"]])
#          dr3_source_id    BpmRp0  BpmRp0_err  phot_g_mean_mag_0  phot_g_mean_mag_0_err  cmd_age_myr
# 0  6794047652729201024  2.063841    0.005469           7.819516               0.002930         40.3
# 1    64895139073954944  2.463277    0.011456           9.370757               0.010248        117.0
```

`phot` is a copy of your table with these columns added, ready for `ages_for_table`:

| column | meaning |
|---|---|
| `distance` | 1000 / parallax (pc); NaN for a non-positive parallax |
| `extinction_a0` | the A_0 used for the star (mag) |
| `extinction_a0_err_cube` | the A_0 uncertainty you supplied (`a0_err`, 0 here), or the dust-map error in example 3 |
| `extinction_a0_gradient` | dA_0/dd at the star (mag/pc); 0 when A_0 is supplied, non-zero with the dust maps |
| `extinction_a0_err` | total A_0 uncertainty: `extinction_a0_err_cube` plus the distance error times the gradient, in quadrature |
| `parallax_error_inflated` | the parallax error after the RUWE-dependent inflation (equal to `parallax_error` if `inflate_parallax_error=False`) |
| `phot_g_mean_mag_0`, `phot_bp_mean_mag_0`, `phot_rp_mean_mag_0` | absolute, de-extincted G, BP and RP (mag) |
| `BpmRp0` | the de-extincted color (G_BP − G_RP)₀ (mag) |
| `phot_g_mean_mag_0_err`, `phot_bp_mean_mag_0_err`, `phot_rp_mean_mag_0_err` | magnitude uncertainties (mag): flux error, zero point, parallax and A_0 terms in quadrature |
| `BpmRp0_err` | color uncertainty (mag): BP and RP flux errors, zero points and the A_0 term in quadrature |

`table` then carries the `cmd_age_*` columns of example 1 as well. A per-star extinction can be given
as a column name, and an extinction uncertainty is propagated if you have one:

```python
import pandas as pd
import cmd_interp

gaia_table = pd.DataFrame({
    "dr3_source_id":           [6794047652729201024, 64895139073954944],
    "parallax":                [102.9432, 7.309302],
    "parallax_error":          [0.02311219, 0.03100616],
    "ruwe":                    [0.9273568, 1.100511],
    "phot_g_mean_mag":         [7.843440, 15.13400],
    "phot_bp_mean_mag":        [8.909956, 16.50113],
    "phot_rp_mean_mag":        [6.798632, 13.98515],
    "phot_g_mean_flux":        [1.372685e7, 1.664744e4],
    "phot_g_mean_flux_error":  [1.100688e4, 22.51610],
    "phot_bp_mean_flux":       [3.727644e6, 3427.412],
    "phot_bp_mean_flux_error": [8432.241, 29.77980],
    "phot_rp_mean_flux":       [1.512535e7, 20188.24],
    "phot_rp_mean_flux_error": [1.871785e4, 83.50732],
    "a0":                      [0.01, 0.17],       # per-star A_0 (mag)
    "a0_err":                  [0.01, 0.03],       # per-star A_0 uncertainty (mag)
})

phot = cmd_interp.deredden_gaia(gaia_table, a0="a0", a0_err="a0_err")   # column names
table, result = cmd_interp.ages_for_table(phot)

print(table[["dr3_source_id", "extinction_a0", "phot_g_mean_mag_0_err", "cmd_age_myr"]])
#          dr3_source_id  extinction_a0  phot_g_mean_mag_0_err  cmd_age_myr
# 0  6794047652729201024           0.01               0.007792         40.3
# 1    64895139073954944           0.17               0.023085        120.0
```

| column | meaning |
|---|---|
| `extinction_a0` | the per-star A_0 read from the `a0` column (mag) |
| `extinction_a0_err` | the per-star uncertainty read from the `a0_err` column (mag); it is also stored in `extinction_a0_err_cube` |
| `phot_g_mean_mag_0_err` | the M_G,0 uncertainty, now larger because the A_0 uncertainty is included |
| `cmd_age_myr` | the most probable age (Myr); the Pleiades star moves from 117.0 to 120.0 Myr with the higher A_0 |

If you only have source IDs, `ages_from_gaia` fetches the `gaia_source` rows first (needs the `gaia`
extra and a network connection; the archive is queried anonymously):

```python
import cmd_interp

source_ids = [64895139073954944, 6794047652729201024]   # Gaia DR3 source IDs
a0 = 0.12                                                # A_0 (mag) for every star

table, result = cmd_interp.ages_from_gaia(source_ids, a0=a0)

print(table[["dr3_source_id", "found_in_gaia", "BpmRp0", "phot_g_mean_mag_0", "cmd_age_myr", "cmd_age_limit"]])
#          dr3_source_id  found_in_gaia    BpmRp0  phot_g_mean_mag_0  cmd_age_myr cmd_age_limit
# 0    64895139073954944           True  2.463275           9.370754        117.0
# 1  6794047652729201024           True  2.063842           7.819516         40.3
```

`table` has one row per requested ID, in the order given, with:

| column | meaning |
|---|---|
| `dr3_source_id` | the Gaia DR3 source ID you asked for |
| `found_in_gaia` | `True` if the ID exists in DR3; a missing ID keeps its row with NaN everywhere else |
| `ra`, `dec`, `pmra`, `pmdec`, `parallax`, `parallax_error`, `ruwe`, `phot_{g,bp,rp}_mean_mag`, `phot_{g,bp,rp}_mean_flux`, `phot_{g,bp,rp}_mean_flux_error`, `azero_gspphot` | the `gaia_source` columns fetched from the archive (`azero_gspphot` is Gaia's own GSP-Phot A_0, for comparison; it is not used) |
| `distance`, `extinction_a0`, `extinction_a0_err`, …, `BpmRp0`, `BpmRp0_err` | the de-reddened photometry and uncertainties of `deredden_gaia`, as in the table above |
| `cmd_age_myr`, `cmd_age_lower_myr`, …, `cmd_age_blue_caution` | the age columns of example 1 |

`result` is the `infer_ages` output for the fitted rows with `result["row"]`, as in example 1.

For a single star, `star_age_from_gaia` returns a dict like `star_age` with the Gaia and de-reddened
photometry columns included:

```python
import cmd_interp

source_id = 64895139073954944   # a Pleiades M dwarf
a0 = 0.12                       # A_0 (mag)

star = cmd_interp.star_age_from_gaia(source_id, a0=a0)

print(star["BpmRp0"], star["phot_g_mean_mag_0"], star["phot_g_mean_mag_0_err"])   # 2.463 9.371 0.010
print(star["age_map_myr"])                                                        # 117.0
```

`star` holds every key of the `star_age` dict (example 1) plus every column of the `ages_from_gaia`
table above as a scalar:

| key | meaning |
|---|---|
| `BpmRp0`, `phot_g_mean_mag_0`, `BpmRp0_err`, `phot_g_mean_mag_0_err` | the four fit inputs derived from the Gaia photometry; also available as `color`, `abs_g`, `color_err`, `abs_g_err` |
| `age_map_myr`, `age_myr_16`, `age_myr_84`, `limit`, `posterior`, `age_grid_myr`, … | the age results of `star_age` |
| `cmd_age_myr`, `cmd_age_lower_myr`, … | the same results in their table-column form |
| `ra`, `dec`, `parallax`, `distance`, `extinction_a0`, … | the Gaia and de-reddened photometry columns |

### 3. You have only a Gaia DR3 source ID

With the EXPLORE 3D dust maps, one source ID is enough: the package queries the Gaia archive
(anonymously), integrates the dust map to the star's distance to get A_0 and its uncertainty, corrects
the photometry, propagates every uncertainty, and fits the age. This needs the `gaia` and `dust`
extras and a one-time download of two FITS cubes (~260 MB) from the CDS record of
[Vergely et al. (2022), J/A+A/664/A174](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174#/browse):
`explore_cube_density_values_010pc_v2.fits` (the extinction density, 10 pc resolution) and
`explore_cube_density_errors_050pc_v2.fits` (its error, 50 pc resolution), both in the `fits/`
subdirectory. Put them in a folder and point the package at it, either in the shell before starting
Python:

```bash
export CMD_INTERP_DUST_DIR=/path/to/explore_cubes
```

or from Python:

```python
import os
import cmd_interp

os.environ["CMD_INTERP_DUST_DIR"] = "/path/to/explore_cubes"   # folder holding the two .fits cubes
source_id = 64895139073954944                                  # a Pleiades M dwarf

star = cmd_interp.star_age_from_gaia(source_id, a0="dustmap")

print(star["extinction_a0"], star["extinction_a0_err"])   # 0.165 0.003    A_0 from the dust map
print(star["age_map_myr"])                                # 119.0
```

The dust-map keys in `star` (the rest is as in example 2):

| key | meaning |
|---|---|
| `distance` | 1000 / parallax (pc), the distance the dust map is integrated to |
| `extinction_a0` | A_0 (mag) integrated along the line of sight from the Sun to the star |
| `extinction_a0_err_cube` | the dust-map uncertainty of A_0 (mag): the error cube integrated along the same line of sight |
| `extinction_a0_gradient` | dA_0/dd at the star (mag/pc), the local dust density |
| `extinction_a0_err` | total A_0 uncertainty (mag): `extinction_a0_err_cube` and the parallax error times `extinction_a0_gradient`, in quadrature; this is the value propagated into the photometry |
| `age_map_myr` | the most probable age (Myr) with this A_0 and its uncertainty included |

Or pass the cube paths explicitly:

```python
import cmd_interp

dust = cmd_interp.ExploreDustMap(
    "/path/to/explore_cubes/explore_cube_density_values_010pc_v2.fits",
    "/path/to/explore_cubes/explore_cube_density_errors_050pc_v2.fits",
)
source_ids = [64895139073954944, 6794047652729201024]   # a Pleiades M dwarf, AU Mic

table, result = cmd_interp.ages_from_gaia(source_ids, a0="dustmap", dust_map=dust)

print(table[["dr3_source_id", "extinction_a0", "extinction_a0_err", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr"]])
#          dr3_source_id  extinction_a0  extinction_a0_err  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr
# 0    64895139073954944        0.16478           0.003259        119.0         130.060208        2976.067306
# 1  6794047652729201024        0.00746           0.000100         40.3          29.316139          45.512950
```

| column | meaning |
|---|---|
| `extinction_a0`, `extinction_a0_err` | A_0 (mag) from the dust map and its total uncertainty, as above |
| `cmd_age_myr` | the most probable age (Myr) |
| `cmd_age_lower_myr`, `cmd_age_upper_myr` | 16th and 84th percentiles (Myr); for the Pleiades star the posterior is broad and skewed to old ages, so its most probable age lies below the 16th percentile |

The same works from the command line:

```bash
cmd-interp gaia 64895139073954944 --dustmap --plot pleiades_star.png
```

The HDF5 versions of the same cubes from the [EXPLORE G-Tomo service](https://explore-platform.eu)
(`.h5`, needs `h5py`) are accepted too. Without the dust maps you must give A_0 yourself (example 2);
A_0 = 0 is a fine approximation for stars within a few tens of parsecs.

### 4. You know the metallicity

The reference sequences are solar-metallicity clusters. A metal-poor star sits below the solar sequence
at fixed color, so it would be read as older than it is, and a metal-rich star as younger. If you have
[M/H] for the star, pass it as `mh` and the sequences are shifted to that metallicity (using PARSEC
isochrone offsets) before the comparison:

```python
import cmd_interp

color = 2.0          # (G_BP - G_RP)_0
abs_g = 8.0          # M_G,0
color_err = 0.01     # mag
abs_g_err = 0.03     # mag

solar = cmd_interp.star_age(color, abs_g, color_err, abs_g_err)
metal_poor = cmd_interp.star_age(color, abs_g, color_err, abs_g_err, mh=-0.3)   # [M/H] = -0.3 dex
metal_rich = cmd_interp.star_age(color, abs_g, color_err, abs_g_err, mh=+0.1)   # [M/H] = +0.1 dex

print(solar["age_map_myr"], metal_poor["age_map_myr"], metal_rich["age_map_myr"])   # 74.5 34.3 88.0
```

| key | meaning |
|---|---|
| `age_map_myr` | the most probable age (Myr) with the sequences shifted to the given [M/H] |
| `mh` | the [M/H] used: 0.0, −0.3 and 0.1 for the three fits |
| `metallicity_clipped` | `True` if [M/H] was outside −0.5 to +0.3 or the color outside 0.8 to 2.9, so the edge offset was used; `False` for all three here |

In a table, give one [M/H] per star in a column (`mh_col=`) or one value for every star (`mh=`).
Missing or NaN values are treated as solar:

```python
import numpy as np
import pandas as pd
import cmd_interp

stars = pd.DataFrame({
    "name": ["star A", "star B", "star C"],
    "BpmRp0": [1.2, 2.0, 2.6],                   # (G_BP - G_RP)_0
    "phot_g_mean_mag_0": [5.5, 8.0, 10.4],       # M_G,0
    "BpmRp0_err": [0.01, 0.01, 0.02],            # mag
    "phot_g_mean_mag_0_err": [0.03, 0.03, 0.05], # mag
    "mh": [np.nan, -0.1, 0.1],                   # [M/H] (dex); NaN = solar
})

table, result = cmd_interp.ages_for_table(stars, mh_col="mh")

print(table[["name", "mh", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr", "cmd_age_limit"]].round(1))
#      name   mh  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr cmd_age_limit
# 0  star A  NaN         21.6               16.9               29.3
# 1  star B -0.1         57.5               48.2               68.6
# 2  star C  0.1        675.0              675.0                NaN           old
```

| column | meaning |
|---|---|
| `mh` | your input column, left as it is (NaN for star A) |
| `cmd_age_mh` | the [M/H] actually used: 0.0 for star A (NaN means solar), −0.1 and 0.1 for the others |
| `cmd_age_metallicity_clipped` | `True` where the edge offset was used; `False` for all three |
| `cmd_age_myr`, `cmd_age_lower_myr`, `cmd_age_upper_myr`, `cmd_age_limit` | as in example 1; star B is now 57.5 Myr instead of 74.5 Myr because a metal-poor star of its brightness is younger than a solar one |

From source IDs, pass one value per star (needs the `gaia` extra and the dust maps, or `a0=` as in
example 2):

```python
import cmd_interp

source_ids = [64895139073954944, 6794047652729201024]   # a Pleiades M dwarf, AU Mic
mh_values = [-0.1, 0.0]                                 # [M/H] (dex) of each star

table, result = cmd_interp.ages_from_gaia(source_ids, a0="dustmap", mh=mh_values)

print(table[["dr3_source_id", "cmd_age_mh", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr"]])
#          dr3_source_id  cmd_age_mh  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr
# 0    64895139073954944        -0.1         93.0          82.991828          110.94843
# 1  6794047652729201024         0.0         40.3          29.316139           45.51295

star = cmd_interp.star_age_from_gaia(64895139073954944, a0="dustmap", mh=-0.1)
print(star["age_map_myr"])    # 93.0, against 119.0 at solar metallicity and 201.4 at [M/H] = +0.1
```

| column / key | meaning |
|---|---|
| `cmd_age_mh` | the [M/H] used for each star, taken from `mh_values` in the order of `source_ids` |
| `cmd_age_myr`, `cmd_age_lower_myr`, `cmd_age_upper_myr` | the most probable age and 68% interval (Myr) at that metallicity |
| `age_map_myr` | the same most probable age for the single-star call |

The offsets are tabulated for [M/H] from −0.5 to +0.3 and (BP−RP)₀ from 0.8 to 2.9; outside that
range the edge value is used and the star is flagged `metallicity_clipped`. At (BP−RP)₀ ≈ 2 a 0.1 dex
change in [M/H] moves the age by roughly 20%, and metal-poor offsets are two to three times larger
than metal-rich ones.

## Things to know

- **Color range.** Stars redder than (BP−RP)₀ = 3.0 get no age (NaN, flagged `color_out_of_range`;
  change with `max_color=`). Bluer than 0.8 the reference sequences start to cross and ages there should
  be used with care: such stars are flagged `blue_caution` and raise a `cmd_interp.BlueColorWarning`
  (silence it with `warnings.filterwarnings("ignore", category=cmd_interp.BlueColorWarning)`).
- **Main-sequence stars.** Only stars still above the main sequence carry real age information. Once a
  star has settled on the main sequence it comes back as an old limit with a broad posterior.
- **Binaries.** Unresolved binaries are brighter and therefore look younger. Remove them beforehand
  (for example with RUWE or the Gaia non-single-star tables); the fit treats every star as single.
- **Posteriors** are unit-sum node weights on an age grid that mixes 1 Myr and 5 Myr steps; use
  `posterior_density` to get a density per Myr before plotting or combining with other age likelihoods:

```python
import numpy as np
import cmd_interp

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
age_grid = star["age_grid_myr"]                                  # (1343,) ages in Myr
weights = star["posterior"]                                      # (1343,) unit-sum node weights

density = cmd_interp.posterior_density(weights, age_grid)        # (1, 1343) probability per Myr

print(np.sum(weights), np.trapz(density[0], age_grid))           # 0.99999994 1.0  (float32 weights)
```

| variable | meaning |
|---|---|
| `age_grid` | the age grid (Myr), shape `(1343,)` |
| `weights` | the posterior as node weights that sum to 1, shape `(1343,)`; a weight is the probability of the grid node, not a density, so wider 5 Myr nodes carry more weight than 1 Myr nodes at the same density |
| `density` | probability per Myr on the same grid, shape `(1, 1343)` (one row per posterior; `posterior_density` accepts a 2-D array of many posteriors); its integral over the grid is 1 |

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
catalogs (about 5 GB per million stars).

## Plots

`plot_star` draws the star on the CMD with the reference sequences next to its age posterior (needs
the `plot` extra):

```python
import cmd_interp
from cmd_interp.plotting import plot_star

star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03, mh=-0.3)

fig, axes = plot_star(star)
fig.savefig("star.png", dpi=150, bbox_inches="tight")
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
loosely (`beta_Pic`, `BPMG`, `Pleiades`, `Melotte_22`, `alpha Per`, `Coma Ber` all work). A reference
can also be left out entirely:

```python
import cmd_interp

color, abs_g, color_err, abs_g_err = 2.0, 8.0, 0.01, 0.03

default = cmd_interp.star_age(color, abs_g, color_err, abs_g_err)
younger_bpic = cmd_interp.star_age(color, abs_g, color_err, abs_g_err, reference_ages={"beta_Pic": 21.0})
older_pair = cmd_interp.star_age(color, abs_g, color_err, abs_g_err, reference_ages={"Pleiades": 112, "alpha Per": 85})
no_coma = cmd_interp.star_age(color, abs_g, color_err, abs_g_err, exclude_references=["Coma Ber"])

print(default["age_map_myr"], younger_bpic["age_map_myr"], older_pair["age_map_myr"], no_coma["age_map_myr"])
# 74.5 74.5 80.5 74.5
```

| variable | meaning |
|---|---|
| `default["age_map_myr"]` | the most probable age (Myr) with the default reference ages |
| `younger_bpic["age_map_myr"]` | the same with β Pic at 21 Myr: unchanged, because this star is bracketed by IC 2602 and α Per |
| `older_pair["age_map_myr"]` | with α Per at 85 Myr and the Pleiades at 112 Myr: the age moves with the bracketing sequence |
| `no_coma["age_map_myr"]` | with Coma Ber left out (`n_references` drops to 5): unchanged for this star |

Other options accepted by every entry point:

| argument | default | meaning |
|---|---|---|
| `systematic_mag` | 0.05 | extra magnitude scatter added in quadrature at every age |
| `mh` | solar | stellar [M/H] |
| `reference_ages`, `exclude_references` | defaults | change or drop reference clusters |
| `max_color`, `min_color` | 3.0, none | (BP−RP)₀ range that receives ages |

For example, a fit with a larger systematic term:

```python
import cmd_interp

star = cmd_interp.star_age(
    color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03,
    systematic_mag=0.08,      # extra scatter (mag) added in quadrature at every age; default 0.05
)
print(star["age_map_myr"], star["age_myr_16"], star["age_myr_84"])   # 74.5 60.0 86.0
```

| key | meaning |
|---|---|
| `age_map_myr` | the most probable age (Myr): unchanged, since the extra scatter widens the likelihood without moving its peak |
| `age_myr_16`, `age_myr_84` | the 68% credible interval (Myr): wider than the default 62.7–83.9 Myr |

How the likelihood, the sequence widths, the prior and the flat extension outside the library are
defined, and how the ages compare with independent cluster ages, is in the paper.

## Data and credits

- **Cluster membership and Gaia DR3 photometry:** Hunt & Reffert (2023, A&A 673, A114).
- **Empirical cluster sequences:** Rottensteiner & Meingast (2024).
- **β Pic members:** the Montreal Open Clusters and Associations database (MOCAdb; Gagné 2024).
- **β Pic age:** Binks & Jeffries (2014).
- **Metallicity offsets:** PARSEC v1.2S (Bressan et al. 2012), via the CMD 3.7 web interface.
- **Extinction law:** Fitzpatrick et al. (2019) Gaia EDR3 coefficients (included).
- **3D dust maps:** EXPLORE, Vergely et al. (2022, A&A 664, A174), [CDS J/A+A/664/A174](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174) (not included).
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
