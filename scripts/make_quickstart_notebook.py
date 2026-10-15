"""Build examples/cmd_interp_quickstart.ipynb (run from the repository root)."""
import nbformat as nbf

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip()))  # noqa: E731
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip()))  # noqa: E731

md("""
# cmd-interp quick start

Stellar ages from empirical cluster sequences in the Gaia CMD.

1. One star from its de-extincted photometry
2. One star straight from its Gaia DR3 source ID (queries the ESA archive; needs internet)
3. Many Gaia IDs at once, with A₀ from the EXPLORE dust maps, and their posteriors
4. Settings: metallicity, systematic scatter, reference ages
5. Your own table → ages, and saving the results
6. The colour limits (redder than 3.0 → NaN; bluer than 0.8 → warning)

Defaults: hybrid scatter, 0.05 mag systematic, solar [M/H], the six bundled references (UBC 17a → Coma Ber).
""")
code("""
%matplotlib inline
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cmd_interp as ci
from cmd_interp.plotting import plot_star

print("cmd-interp", ci.__version__)
print("reference ages (Myr):", ci.REFERENCE_AGES)


def summary(star):
    \"\"\"One-line summary of a star_age / star_age_from_gaia result.\"\"\"
    if not star["has_posterior"]:
        return "no age (outside the reference colour range)"
    text = (f"MAP {star['age_map_myr']:.1f} Myr, median {star['age_myr_50']:.1f} Myr, "
            f"68% [{star['age_myr_16']:.1f}, {star['age_myr_84']:.1f}] Myr")
    if star["limit"] == "young":
        text += f"  -> limit: younger than {star['youngest_age']:g} Myr"
    elif star["limit"] == "old":
        text += f"  -> limit: older than {star['oldest_age']:g} Myr"
    return text
""")

md("""
## 1. One star from de-extincted photometry

Inputs: (BP−RP)₀, absolute M_G,0 and their 1σ errors (mag). Optionally `mh=` for [M/H].
""")
code("""
star = ci.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
print(summary(star))
fig, axes = plot_star(star)
plt.show()
""")

md("""
## 2. One star from its Gaia DR3 source ID

This queries `gaiadr3.gaia_source`, de-extincts the photometry, propagates the uncertainties and fits
the age. `a0` must be given explicitly: a number, or `"dustmap"` (section 3).
Example: AU Mic, an M1 member of the β Pic moving group at 9.7 pc, so A₀ ≈ 0.
""")
code("""
aumic = ci.star_age_from_gaia(6794047652729201024, a0=0.0)
print(f"(BP-RP)_0 = {aumic['BpmRp0']:.3f} ± {aumic['BpmRp0_err']:.3f},  "
      f"M_G,0 = {aumic['phot_g_mean_mag_0']:.3f} ± {aumic['phot_g_mean_mag_0_err']:.3f}")
print(summary(aumic))
fig, axes = plot_star(aumic)
plt.show()
""")

md("""
## 3. Many Gaia IDs, A₀ from the EXPLORE dust maps

Point `DUST_DIR` at the folder holding the two EXPLORE cubes. If they are not found, this cell falls back
to a single A₀ for every star.
""")
code("""
DUST_DIR = Path("~/Downloads").expanduser()
VALUES = DUST_DIR / "explore_cube_density_values_010pc_v2.h5"
ERRORS = DUST_DIR / "explore_cube_density_errors_050pc_v2.h5"

ids = [6794047652729201024,   # AU Mic
       6665284743653466112, 5973004833797555584, 1118791671544467840, 6723345211080957696,
       5411128214187804672, 357574689027564928]

if VALUES.exists() and ERRORS.exists():
    dust = ci.ExploreDustMap(VALUES, ERRORS)
    table, result = ci.ages_from_gaia(ids, a0="dustmap", dust_map=dust)
else:
    print("EXPLORE cubes not found -> using A0 = 0.1 for every star")
    table, result = ci.ages_from_gaia(ids, a0=0.1)

cols = ["dr3_source_id", "BpmRp0", "BpmRp0_err", "phot_g_mean_mag_0", "phot_g_mean_mag_0_err", "extinction_a0",
        "cmd_age_myr", "cmd_age_median_myr", "cmd_age_lower_myr", "cmd_age_upper_myr", "cmd_age_limit"]
table[cols].round(3)
""")
md("""
### The posteriors

`result["posterior"]` holds unit-sum **node weights** on `result["age_grid_myr"]`. The grid mixes 1 Myr and
5 Myr steps, so convert with `ci.posterior_density` before plotting. `result["row"]` maps each posterior to its table row.
""")
code("""
density = ci.posterior_density(result["posterior"], result["age_grid_myr"])
fig, ax = plt.subplots(figsize=(8, 4.5))
for k, row in enumerate(result["row"]):
    if not result["has_posterior"][k]:
        continue
    y = density[k] / np.nanmax(density[k])
    ax.plot(result["age_grid_myr"][1:], y[1:], lw=1.8,
            label=f"{table.dr3_source_id[row]}  (BP-RP)$_0$={table.BpmRp0[row]:.2f}")
ax.set(xscale="log", xlim=(0.6, 4000), ylim=(0, 1.05), xlabel="Age [Myr]", ylabel="Relative posterior density")
ax.legend(fontsize=8, frameon=False)
plt.show()
""")

md("""
## 4. Settings

Every function takes the same options:
- `mh` — stellar [M/H]. Each sequence is shifted by the PARSEC offset at that metallicity before the comparison.
- `systematic_mag` — extra scatter added in quadrature; default 0.05 mag.
- `reference_ages` — new ages for any reference cluster; the others keep their defaults.
- `exclude_references` — clusters to leave out.
""")
code("""
variants = {
    "default": {},
    "[M/H] = -0.2": dict(mh=-0.2),
    "[M/H] = +0.2": dict(mh=0.2),
    "systematic 0.10 mag": dict(systematic_mag=0.10),
    "no systematic": dict(systematic_mag=0.0),
    "beta Pic at 21 Myr": dict(reference_ages={"beta_Pic": 21.0}),
    "Pleiades 112, alpha Per 85": dict(reference_ages={"Pleiades": 112, "alpha Per": 85}),
    "without Coma Ber": dict(exclude_references=["Coma Ber"]),
}
rows = []
for name, kwargs in variants.items():
    s = ci.star_age(2.0, 8.0, 0.01, 0.03, **kwargs)
    rows.append(dict(setting=name, MAP=s["age_map_myr"], median=s["age_myr_50"], p16=s["age_myr_16"], p84=s["age_myr_84"]))
pd.DataFrame(rows).round(1)
""")
md("""
**Reading the numbers.** The posterior is flat beyond the reference ages: below 4.6 Myr, and from 675 up to 4000 Myr.
For a star close to the main sequence (here, after shifting to [M/H] = +0.2), most of the probability can sit on
that old plateau. The median and 68% interval then lie above the MAP, because the MAP is taken from the fitted
range only. Check `limit` and the posterior plot, not just one number.
""")
code("""
# the resulting reference set for a given override
ci.resolve_reference_ages({"IC 2602": 46})
""")

md("""
## 5. Your own table

`ages_for_table` takes a DataFrame or a .csv/.feather/.parquet/.fits path. Default columns are `BpmRp0`,
`phot_g_mean_mag_0`, `BpmRp0_err` and `phot_g_mean_mag_0_err`; rename them with `color_col=` etc.

If you start from raw gaia_source columns, run `ci.deredden_gaia(table, a0=...)` first. Its output
already has those column names.
""")
code("""
# example: reuse the de-extincted photometry from section 3 (any DataFrame with these columns works)
mine = table[["dr3_source_id", "BpmRp0", "phot_g_mean_mag_0", "BpmRp0_err", "phot_g_mean_mag_0_err"]].copy()
mine["mh"] = [0.0, 0.1, -0.1, 0.0, np.nan, 0.0, 0.0]          # optional per-star [M/H]; NaN = solar

ages, res = ci.ages_for_table(mine, mh_col="mh", systematic_mag=0.05)
ages.filter(regex="dr3_source_id|BpmRp0$|cmd_age_(myr|lower|upper|limit|mh)").round(2)
""")
code("""
# save the table and the posteriors
out_dir = Path("cmd_interp_output"); out_dir.mkdir(exist_ok=True)
ages.to_csv(out_dir / "ages.csv", index=False)
np.save(out_dir / "posteriors.npy", res["posterior"])          # rows = res["row"], unit-sum node weights
np.save(out_dir / "age_grid_myr.npy", res["age_grid_myr"])
sorted(p.name for p in out_dir.iterdir())
""")

md("""
## 6. Colour limits

- **(BP−RP)₀ > 3.0:** no age (NaN) and `color_out_of_range = True`. Change it with `max_color=`.
- **(BP−RP)₀ < 0.8:** the reference sequences start to cross. Ages are still returned, but a
  `BlueColorWarning` is printed and the star is flagged `blue_caution`.
""")
code("""
check = ci.infer_ages(color=[0.6, 2.0, 3.3], abs_g=[3.4, 8.0, 11.0], color_err=0.01, abs_g_err=0.03)
pd.DataFrame({"BP-RP": [0.6, 2.0, 3.3], "age_map_myr": check["age_map_myr"],
              "blue_caution": check["blue_caution"], "color_out_of_range": check["color_out_of_range"]})
""")

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3 (ipykernel)", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
})
nbf.write(nb, "examples/cmd_interp_quickstart.ipynb")
print("written")
