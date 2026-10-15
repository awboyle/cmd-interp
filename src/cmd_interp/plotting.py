"""Diagnostic plots for a single star (requires matplotlib: pip install "cmd-interp[plot]")."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator

from .grid import posterior_density
from .metallicity import parsec_offset
from .references import load_references


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:  # pragma: no cover
        raise ImportError('Plotting needs matplotlib: pip install "cmd-interp[plot]"') from error
    return plt


def plot_cmd(star, references=None, ax=None, color_half_width=0.6, mag_half_width=2.0, n_grid_ages=25, show_members=True):
    """The star on the CMD with the reference sequences and the interpolated sequences between them.

    star: output of star_age().  If the star has a non-solar [M/H], the reference sequences are
    drawn as they appear at that metallicity (PARSEC-shifted, dashed) next to the solar ones.
    """
    plt = _pyplot()
    references = load_references() if references is None else references
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    color, abs_g = star["color"], star["abs_g"]
    color_grid = np.linspace(color - color_half_width, color + color_half_width, 400)
    covering = [r for r in references if r["color"][0] <= color <= r["color"][-1]]
    cmap = plt.colormaps["viridis"]
    shades = np.linspace(0.08, 0.92, len(references))
    mh = star.get("mh", 0.0) or 0.0
    offset = parsec_offset() if mh != 0 else None

    if len(covering) >= 2:
        lo = max(r["color"][0] for r in covering)
        hi = min(r["color"][-1] for r in covering)
        cgrid = color_grid[(color_grid >= lo) & (color_grid <= hi)]
        mags = np.vstack([np.interp(cgrid, r["color"], r["mag"]) for r in covering])
        log_ages = np.log10([r["age_myr"] for r in covering])
        between = np.linspace(log_ages[0], log_ages[-1], n_grid_ages)
        for k, curve in enumerate(PchipInterpolator(log_ages, mags, axis=0)(between)):
            ax.plot(cgrid, curve, color="0.75", lw=0.6, zorder=1,
                    label="interpolated sequences" if k == 0 else None)

    for reference, shade in zip(references, shades):
        line_color = cmap(shade)
        if show_members:
            ax.scatter(reference["member_color"], reference["member_mag"], s=5, alpha=0.15,
                       color=line_color, rasterized=True, zorder=0)
        in_range = (color_grid >= reference["color"][0]) & (color_grid <= reference["color"][-1])
        seq = np.interp(color_grid[in_range], reference["color"], reference["mag"])
        ax.plot(color_grid[in_range], seq, color=line_color, lw=2,
                label=f"{reference['cluster']} ({reference['age_myr']:g} Myr)", zorder=3)
        if offset is not None and in_range.any():
            shift = offset(mh, np.log10(reference["age_myr"]), color_grid[in_range])
            ax.plot(color_grid[in_range], seq + shift, color=line_color, lw=1.2, ls="--", zorder=3)

    ax.errorbar(color, abs_g, xerr=star["color_err"], yerr=star["abs_g_err"], marker="*", markersize=13,
                color="k", ecolor="k", elinewidth=1.2, capsize=3, ls="none", zorder=10,
                label="star" + (f" ([M/H] = {mh:+.2f}; dashed = shifted sequences)" if mh else ""))
    ax.set(xlabel=r"$(G_{BP}-G_{RP})_0$ [mag]", ylabel=r"$M_{G,0}$ [mag]",
           xlim=(color - color_half_width, color + color_half_width),
           ylim=(abs_g + mag_half_width, abs_g - mag_half_width))
    ax.legend(frameon=False, fontsize="small")
    return ax


def plot_posterior(star, ax=None, log_age=True):
    """The star's posterior density per Myr with its MAP, median and 68% interval."""
    plt = _pyplot()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    age = star["age_grid_myr"]
    density = posterior_density(star["posterior"], age)[0]
    if not np.isfinite(density).any():
        ax.text(0.5, 0.5, "no posterior (star outside the reference colour range)",
                ha="center", va="center", transform=ax.transAxes)
        return ax
    density = density / np.nanmax(density)
    ax.plot(age, density, color="k", lw=1.5)
    in_68 = (age >= star["age_myr_16"]) & (age <= star["age_myr_84"])
    ax.fill_between(age[in_68], density[in_68], color="C0", alpha=0.25, label="68% interval")
    ax.axvline(star["age_myr_50"], color="C0", ls="--", lw=1.2, label="median")
    ax.axvline(star["age_map_myr"], color="C3", ls=":", lw=1.2, label="MAP")
    if log_age:
        ax.set_xscale("log")
        ax.set_xlim(max(age[1], 0.5), age[-1])
    ax.set(xlabel="Age [Myr]", ylabel="Relative posterior density", ylim=(0, 1.05))
    title = (f"MAP {star['age_map_myr']:.1f} Myr; median {star['age_myr_50']:.1f} Myr; "
             f"68% [{star['age_myr_16']:.1f}, {star['age_myr_84']:.1f}] Myr")
    if star.get("limit"):
        title += f"\nlimit: {star['limit']}"
    ax.set_title(title, fontsize="medium")
    ax.legend(frameon=False)
    return ax


def plot_star(star, references=None, figsize=(12, 5)):
    """Two-panel summary: CMD position and posterior."""
    plt = _pyplot()
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    plot_cmd(star, references=references, ax=axes[0])
    plot_posterior(star, ax=axes[1])
    fig.tight_layout()
    return fig, axes
