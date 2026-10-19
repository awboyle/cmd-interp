"""Quick start: the four README examples.

1. de-reddened photometry and uncertainties in hand;
2. Gaia photometry and an extinction, uncertainties derived for you (needs network access);
3. a Gaia DR3 source ID alone, extinction from the EXPLORE 3D dust maps (needs network access and the
   dust cubes, see README);
4. the same with a known metallicity.

Examples 2-4 are skipped when the archive or the dust maps are not available.
"""

import warnings

import cmd_interp

PLEIADES_M_DWARF = 64895139073954944   # Gaia DR3 source ID

# ---- 1. you have (BP-RP)_0, M_G,0 and their uncertainties
star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03)
print(f"1. MAP {star['age_map_myr']:.1f} Myr, 68% [{star['age_myr_16']:.1f}, {star['age_myr_84']:.1f}] Myr, limit '{star['limit']}'")

# ---- 2. Gaia DR3 source ID plus an extinction A_0; the uncertainties come from the Gaia errors
try:
    star = cmd_interp.star_age_from_gaia(PLEIADES_M_DWARF, a0=0.12)
except Exception as error:   # no network, or astroquery not installed
    print(f"2. skipped ({error})")
else:
    print(f"2. (BP-RP)_0 = {star['BpmRp0']:.3f}, M_G,0 = {star['phot_g_mean_mag_0']:.3f} +/- {star['phot_g_mean_mag_0_err']:.3f}, "
          f"MAP {star['age_map_myr']:.1f} Myr")

# ---- 3. source ID alone: A_0 and its uncertainty from the EXPLORE dust maps
try:
    star = cmd_interp.star_age_from_gaia(PLEIADES_M_DWARF, a0="dustmap")
except Exception as error:   # dust cubes not found (set CMD_INTERP_DUST_DIR), no network, or extras missing
    print(f"3. skipped ({error})")
else:
    print(f"3. A_0 = {star['extinction_a0']:.3f} +/- {star['extinction_a0_err']:.3f}, MAP {star['age_map_myr']:.1f} Myr")

# ---- 4. with a known metallicity
star = cmd_interp.star_age(color=2.0, abs_g=8.0, color_err=0.01, abs_g_err=0.03, mh=-0.3)
print(f"4. [M/H] = -0.3: MAP {star['age_map_myr']:.1f} Myr (74.5 Myr at solar metallicity)")

try:
    from cmd_interp.plotting import plot_star
except ImportError:
    pass
else:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig, _ = plot_star(cmd_interp.star_age(2.0, 8.0, 0.01, 0.03, mh=-0.3))
    fig.savefig("quickstart_star.png", dpi=150, bbox_inches="tight")
    print("wrote quickstart_star.png")
