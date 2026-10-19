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
pip install "cmd-interp[dust] @ git+https://github.com/awboyle/cmd-interp"   # h5py + astropy, for the 3D dust maps
pip install "cmd-interp[plot] @ git+https://github.com/awboyle/cmd-interp"   # matplotlib, for the diagnostic figures
```

Requires Python ≥ 3.10, numpy, scipy, pandas and pyarrow. The reference sequences (≈300 kB) are included.

## Quick start

The fit needs four numbers per star: the extinction-corrected colour (G_BP − G_RP)₀, the
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
color_err = 0.01     # 1-sigma uncertainty of the colour (mag)
abs_g_err = 0.03     # 1-sigma uncertainty of the absolute magnitude (mag)

star = cmd_interp.star_age(color=color, abs_g=abs_g, color_err=color_err, abs_g_err=abs_g_err)

print(star["age_map_myr"])                            # 74.5    most probable age (Myr)
print(star["age_myr_16"], star["age_myr_84"])         # 62.7 83.9    68% credible interval
print(star["age_myr_50"])                             # 73.7    median
print(star["age_myr_025"], star["age_myr_975"])       # 52.1 92.8    95% credible interval
print(star["limit"])                                  # ''      inside the reference library
print(star["posterior"].shape, star["age_grid_myr"][:3])   # (1343,) [0.6 1.6 2.6]
```

`star["posterior"]` is the full posterior as unit-sum weights on the age grid `star["age_grid_myr"]`.
`limit` is `""` for a star inside the reference library, `"young"` for one brighter than the youngest
reference sequence (only an upper limit on the age), and `"old"` for one fainter than the oldest.

Many stars, as arrays:

```python
import numpy as np
import cmd_interp

color = np.array([1.2, 2.0, 2.6])            # (G_BP - G_RP)_0 of each star
abs_g = np.array([5.5, 8.0, 10.4])           # M_G,0 of each star
color_err = np.array([0.01, 0.01, 0.02])     # colour uncertainties (mag)
abs_g_err = np.array([0.03, 0.03, 0.05])     # magnitude uncertainties (mag)

result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err)

print(result["age_map_myr"])      # [ 21.6  74.5 675. ]
print(result["age_myr_16"])       # [ 16.9  62.7 926.3]
print(result["age_myr_84"])       # [ 29.3  83.9 3414.5]
print(result["limit"])            # ['' '' 'old']    the third star is fainter than the oldest sequence
print(result["posterior"].shape)  # (3, 1343)        one row per star on result["age_grid_myr"]
```

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
    "BpmRp0_err": [0.01, 0.01, 0.02],            # colour uncertainty (mag)
    "phot_g_mean_mag_0_err": [0.03, 0.03, 0.05], # magnitude uncertainty (mag)
})

table, result = cmd_interp.ages_for_table(stars)

print(table[["name", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr", "cmd_age_limit"]].round(1))
#      name  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr cmd_age_limit
# 0  star A         21.6               16.9               29.3
# 1  star B         74.5               62.7               83.9
# 2  star C        675.0              675.0                NaN           old
```

The columns added are `cmd_age_myr` (most probable age, or the edge reference age at a limit),
`cmd_age_lower_myr` and `cmd_age_upper_myr` (16th and 84th percentiles), `cmd_age_median_myr`,
`cmd_age_p2p5_myr` and `cmd_age_p97p5_myr`, `cmd_age_limit`, `cmd_age_n_references`,
`cmd_age_scatter_mag`, `cmd_age_mh`, `cmd_age_metallicity_clipped`, `cmd_age_fitted`,
`cmd_age_color_out_of_range` and `cmd_age_blue_caution`. The second return value `result` is the
`infer_ages` output for the fitted rows (`result["posterior"]`, `result["age_grid_myr"]`, and
`result["row"]` giving each fitted star's row in the table). Instead of a DataFrame you can pass a
file path (`.csv`, `.feather`, `.parquet` or `.fits`):

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

`phot` is a copy of your table with `BpmRp0`, `phot_g_mean_mag_0` and their `_err` columns added,
ready for `ages_for_table`. A per-star extinction can be given as a column name, and an extinction
uncertainty is propagated if you have one:

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

### 3. You have only a Gaia DR3 source ID

With the EXPLORE 3D dust maps, one source ID is enough: the package queries the Gaia archive
(anonymously), integrates the dust map to the star's distance to get A_0 and its uncertainty, corrects
the photometry, propagates every uncertainty, and fits the age. This needs the `gaia` and `dust`
extras and a one-time download of the two dust cubes (~230 MB) from the
[EXPLORE G-Tomo service](https://explore-platform.eu): `explore_cube_density_values_010pc_v2.h5` and
`explore_cube_density_errors_050pc_v2.h5`. Put them in a folder and point the package at it, either
in the shell before starting Python:

```bash
export CMD_INTERP_DUST_DIR=/path/to/explore_cubes
```

or from Python:

```python
import os
import cmd_interp

os.environ["CMD_INTERP_DUST_DIR"] = "/path/to/explore_cubes"   # folder holding the two .h5 cubes
source_id = 64895139073954944                                  # a Pleiades M dwarf

star = cmd_interp.star_age_from_gaia(source_id, a0="dustmap")

print(star["extinction_a0"], star["extinction_a0_err"])   # 0.165 0.003    A_0 from the dust map
print(star["age_map_myr"])                                # 119.0
```

Or pass the cube paths explicitly:

```python
import cmd_interp

dust = cmd_interp.ExploreDustMap(
    "/path/to/explore_cubes/explore_cube_density_values_010pc_v2.h5",
    "/path/to/explore_cubes/explore_cube_density_errors_050pc_v2.h5",
)
source_ids = [64895139073954944, 6794047652729201024]   # a Pleiades M dwarf, AU Mic

table, result = cmd_interp.ages_from_gaia(source_ids, a0="dustmap", dust_map=dust)

print(table[["dr3_source_id", "extinction_a0", "extinction_a0_err", "cmd_age_myr", "cmd_age_lower_myr", "cmd_age_upper_myr"]])
#          dr3_source_id  extinction_a0  extinction_a0_err  cmd_age_myr  cmd_age_lower_myr  cmd_age_upper_myr
# 0    64895139073954944        0.16478           0.003259        119.0         130.060208        2976.067306
# 1  6794047652729201024        0.00746           0.000100         40.3          29.316139          45.512950
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

The offsets are tabulated for [M/H] from −0.5 to +0.3 and (BP−RP)₀ from 0.8 to 2.9; outside that
range the edge value is used and the star is flagged `metallicity_clipped`. At (BP−RP)₀ ≈ 2 a 0.1 dex
change in [M/H] moves the age by roughly 20%, and metal-poor offsets are two to three times larger
than metal-rich ones.

## Things to know

- **Colour range.** Stars redder than (BP−RP)₀ = 3.0 get no age (NaN, flagged `color_out_of_range`;
  change with `max_color=`). Bluer than 0.8 the reference sequences start to cross and ages there should
  be used with care: such stars are flagged `blue_caution` and raise a `cmd_interp.BlueColorWarning`
  (silence it with `warnings.filterwarnings("ignore", category=cmd_interp.BlueColorWarning)`).
- **Main-sequence stars.** Only stars still above the main sequence carry real age information. Once a
  star has settled on the main sequence it comes back as an old limit with a broad posterior.
- **Binaries.** Unresolved binaries are brighter and therefore look younger. Remove them beforehand, or
  enable the binary mixture (`binary_fraction=`, `binary_dm=`; see the options table below).
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

Other options accepted by every entry point:

| argument | default | meaning |
|---|---|---|
| `systematic_mag` | 0.05 | extra magnitude scatter added in quadrature at every age |
| `mh` | solar | stellar [M/H] |
| `binary_fraction`, `binary_dm` | off | unresolved-binary mixture: with probability f_b the star is brighter than the sequence by `binary_dm` |
| `reference_ages`, `exclude_references` | defaults | change or drop reference clusters |
| `max_color`, `min_color` | 3.0, none | (BP−RP)₀ range that receives ages |
| `scatter_mode` | `"hybrid"` | how the width of the reference sequences enters the likelihood |

For example, a fit that allows 30% of stars to be equal-mass binaries (0.75 mag brighter) with a
larger systematic term:

```python
import cmd_interp

star = cmd_interp.star_age(
    color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03,
    binary_fraction=0.3,      # prior probability that the star is an unresolved binary
    binary_dm=(0.75,),        # brightening (mag) of the binary component(s) of the mixture
    systematic_mag=0.08,      # extra scatter (mag) added in quadrature at every age
)
print(star["age_map_myr"])
```

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
