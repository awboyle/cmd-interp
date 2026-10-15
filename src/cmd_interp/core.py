"""Age posteriors from the empirical reference sequences."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator

from .grid import age_grid as _age_grid
from .grid import posterior_quantile
from .metallicity import parsec_offset
from .references import load_references
from .scatter import exact_local_scatter, fill_scatter, nearest_two_scatter

SCATTER_MODES = ("hybrid", "nearest_two", "interpolated", "interpolated_fill")
DEFAULT_SCATTER_MODE = "hybrid"
DEFAULT_SYSTEMATIC_MAG = 0.05


def infer_ages(
    color,
    abs_g,
    color_err,
    abs_g_err,
    *,
    mh=None,
    systematic_mag=DEFAULT_SYSTEMATIC_MAG,
    scatter_mode=DEFAULT_SCATTER_MODE,
    references=None,
    binary_fraction=0.0,
    binary_dm=(),
    binary_dm_weights=None,
    age_step_myr=1.0,
    n_color_quadrature=15,
    chunk_size=2000,
    return_posterior=True,
    posterior_out=None,
    metallicity_offset=None,
    verbose=False,
):
    """Age posterior of every star from its position in the (BP-RP)_0 vs M_G CMD.

    Parameters
    ----------
    color, abs_g : array_like
        Extinction-corrected Gaia (BP-RP)_0 colour and absolute G magnitude M_G,0.
    color_err, abs_g_err : array_like or float
        1-sigma uncertainties (mag).  The likelihood is marginalised over the true colour
        with Gauss-Hermite quadrature (``n_color_quadrature`` nodes).
    mh : None, float or array_like
        Stellar [M/H] (dex).  None or 0 = solar (no correction).  Otherwise every model
        sequence is shifted by the PARSEC offset dM_G([M/H], age, colour) before the star is
        compared with it, i.e. the star is moved to its solar-metallicity CMD position at
        each trial age.  NaN entries are treated as solar.
    systematic_mag : float
        Extra magnitude scatter added in quadrature at every age (field metallicity spread
        and zero point of the reference sequences).  Default 0.05.
    scatter_mode : {'hybrid', 'nearest_two', 'interpolated', 'interpolated_fill'}
        How the empirical scatter of the reference sequences enters sigma.  'interpolated'
        interpolates the references' local scatter in log age; 'interpolated_fill' also
        fills colours where a reference has too few members; 'nearest_two' uses one sigma
        per star (mean of the two bracketing references); 'hybrid' (default) uses
        'interpolated_fill' inside the reference envelope, 'nearest_two' for stars
        fainter than every covering reference, and the youngest covering reference's
        (filled) width for stars brighter than it.  The star's own abs_g_err and the
        systematic always enter in quadrature.
    references : list, optional
        Output of load_references(); defaults to the bundled six-cluster library.
    binary_fraction, binary_dm, binary_dm_weights :
        Optional unresolved-binary mixture: at every age the magnitude likelihood is
        (1 - f_b) N(M_G; m_t, s) + f_b sum_k w_k N(M_G; m_t - dm_k, s).  Off by default.
    age_step_myr : float
        Grid step between reference ages (Myr).
    return_posterior : bool
        Return the (n_star, n_age) posterior array (float32).  Set False for large
        catalogues, or pass ``posterior_out`` (e.g. a np.memmap) to receive it instead.
    metallicity_offset : callable, optional
        dM_G(mh, log10_age_myr, colour); defaults to the bundled PARSEC table.

    Returns
    -------
    dict of arrays, one entry per star unless noted:
        age_map_myr, age_myr_025, age_myr_16, age_myr_50, age_myr_84, age_myr_975,
        limit ('young' brighter than the youngest covering sequence, 'old' fainter than the
        oldest, '' otherwise), has_posterior, n_references, scatter_sigma_mag,
        youngest_age, oldest_age, youngest_mag, oldest_mag (metallicity-shifted),
        mh, metallicity_clipped, posterior (n_star, n_age; unit-sum node weights, zero for
        stars without a fit) and age_grid_myr (n_age).

    Notes
    -----
    Each star is fitted with the references whose sequences cover its colour (at least
    two).  The likelihood is evaluated between the youngest and oldest covering reference;
    below and above, the posterior is held flat at its value at the nearest fitted node
    (down to 0 Myr and up to 4000 Myr).  The prior is uniform in linear age.  MAP ages come
    from the fitted range only.
    """
    color = np.asarray(color, dtype=float).reshape(-1)
    n_star = len(color)
    abs_g = np.asarray(abs_g, dtype=float).reshape(-1)
    color_err = np.broadcast_to(np.asarray(color_err, dtype=float), color.shape).copy()
    abs_g_err = np.broadcast_to(np.asarray(abs_g_err, dtype=float), color.shape).copy()
    if abs_g.shape != color.shape:
        raise ValueError("color and abs_g must have the same length")
    if scatter_mode not in SCATTER_MODES:
        raise ValueError(f"scatter_mode must be one of {SCATTER_MODES}")
    if not np.isfinite(systematic_mag) or systematic_mag < 0:
        raise ValueError("systematic_mag must be finite and non-negative")
    if int(n_color_quadrature) != n_color_quadrature or n_color_quadrature < 1:
        raise ValueError("n_color_quadrature must be a positive integer")
    if references is None:
        references = load_references()
    if len(references) < 2:
        raise ValueError("At least two references are needed")

    binary_dm = np.atleast_1d(np.asarray(binary_dm, dtype=float))
    if binary_fraction > 0:
        if not (0 < binary_fraction < 1) or not len(binary_dm):
            raise ValueError("binary_fraction must be in (0, 1) with at least one binary_dm node")
        bw = np.ones(len(binary_dm)) if binary_dm_weights is None else np.asarray(binary_dm_weights, dtype=float)
        log_binary_w = np.log(binary_fraction * bw / bw.sum())

    # metallicity: per-star [M/H], NaN -> solar; all-solar skips the shift entirely
    mh_star = np.zeros(n_star) if mh is None else np.broadcast_to(np.asarray(mh, dtype=float), color.shape).copy()
    mh_star[~np.isfinite(mh_star)] = 0.0
    use_mh = bool(np.any(mh_star != 0))
    offset = None
    metallicity_clipped = np.zeros(n_star, dtype=bool)
    if use_mh:
        offset = parsec_offset() if metallicity_offset is None else metallicity_offset
        if hasattr(offset, "clipped"):
            metallicity_clipped = np.asarray(offset.clipped(mh_star, color), dtype=bool) & (mh_star != 0)
            if metallicity_clipped.any():
                warnings.warn(
                    f"{int(metallicity_clipped.sum())} star(s) have [M/H] or (BP-RP)_0 outside the PARSEC offset "
                    "table; their offsets were taken at the table edge (see 'metallicity_clipped').",
                    stacklevel=2,
                )

    ref_age = np.array([reference["age_myr"] for reference in references])
    full_age_grid = _age_grid(ref_age, age_step_myr=age_step_myr)
    usable = (np.isfinite(color) & np.isfinite(abs_g) & np.isfinite(color_err) & (color_err >= 0)
              & np.isfinite(abs_g_err) & (abs_g_err >= 0))

    cover = np.vstack([(color >= r["color"][0]) & (color <= r["color"][-1]) & usable for r in references])
    hermite_nodes, hermite_weights = np.polynomial.hermite.hermgauss(int(n_color_quadrature))
    hermite_weights = hermite_weights / np.sqrt(np.pi)

    map_age = np.full(n_star, np.nan)
    age_025 = np.full(n_star, np.nan)
    age_16 = np.full(n_star, np.nan)
    age_50 = np.full(n_star, np.nan)
    age_84 = np.full(n_star, np.nan)
    age_975 = np.full(n_star, np.nan)
    sigma_used = np.full(n_star, np.nan)
    n_refs_used = cover.sum(axis=0)
    youngest_mag = np.full(n_star, np.nan)
    oldest_mag = np.full(n_star, np.nan)
    youngest_age = np.full(n_star, np.nan)
    oldest_age = np.full(n_star, np.nan)
    if posterior_out is None and return_posterior:
        posterior_out = np.zeros((n_star, len(full_age_grid)), dtype=np.float32)
    elif posterior_out is not None and posterior_out.shape != (n_star, len(full_age_grid)):
        raise ValueError(f"posterior_out must have shape {(n_star, len(full_age_grid))}")

    patterns = {}
    for j in range(n_star):
        patterns.setdefault(tuple(np.flatnonzero(cover[:, j])), []).append(j)
    n_groups_done = 0
    for pattern, members in patterns.items():
        members = np.asarray(members)
        refs = [references[k] for k in pattern]
        if len(refs) == 0:
            continue
        sub_ages = np.array([r["age_myr"] for r in refs])
        log_sub_ages = np.log10(sub_ages)
        # sequence magnitudes / scatter of the covering references at the stars' colours
        ref_mag = np.vstack([np.interp(color[members], r["color"], r["mag"]) for r in refs])
        if use_mh:   # the sequences as they would appear at each star's metallicity
            ref_mag = ref_mag + offset(mh_star[members][None, :], log_sub_ages[:, None], color[members][None, :])
        ref_scatter = np.vstack([exact_local_scatter(color[members], r) for r in refs])
        youngest_mag[members], oldest_mag[members] = ref_mag[0], ref_mag[-1]
        youngest_age[members], oldest_age[members] = refs[0]["age_myr"], refs[-1]["age_myr"]
        if len(refs) < 2:
            continue
        if scatter_mode == "nearest_two":
            star_sigma = nearest_two_scatter(abs_g[members], ref_mag, ref_scatter)
        elif scatter_mode == "interpolated":
            star_sigma = np.where(np.all(np.isfinite(ref_scatter), axis=0), 1.0, np.nan)
        elif scatter_mode == "interpolated_fill":
            # sigma(t) depends on the age only, never on which side of a sequence the star falls
            ref_scatter = np.vstack([fill_scatter(color[members], r, ref_scatter[k]) for k, r in enumerate(refs)])
            star_sigma = np.where(np.all(np.isfinite(ref_scatter), axis=0), 1.0, np.nan)
        else:  # hybrid
            # interpolated_fill inside the reference envelope; one age-independent width per star outside it, because
            # there an age-dependent sigma makes the likelihood prefer the broadest reference instead of the nearest one.
            # Below every covering sequence: nearest_two.  Above the youngest covering sequence: that sequence's filled
            # width, i.e. the value sigma(t) has at the youngest age, so sigma is continuous across the sequence.
            n2_sigma = nearest_two_scatter(abs_g[members], ref_mag, ref_scatter)
            ref_scatter = np.vstack([fill_scatter(color[members], r, ref_scatter[k]) for k, r in enumerate(refs)])
            fainter = abs_g[members] > ref_mag.max(axis=0)
            brighter = abs_g[members] < ref_mag[0]
            bright_sigma = np.where(np.isfinite(ref_scatter[0]), ref_scatter[0], n2_sigma)
            outside = fainter | brighter
            star_sigma = np.where(fainter, n2_sigma, np.where(
                brighter, bright_sigma, np.where(np.all(np.isfinite(ref_scatter), axis=0), 1.0, np.nan)))
        if scatter_mode != "hybrid":
            outside = np.zeros(len(members), dtype=bool)
        valid = np.isfinite(star_sigma)
        grid_mask = (full_age_grid >= sub_ages[0]) & (full_age_grid <= sub_ages[-1])
        age_grid = full_age_grid[grid_mask]
        grid_index = np.flatnonzero(grid_mask)
        floor_index = np.flatnonzero(full_age_grid < sub_ages[0])
        old_index = np.flatnonzero(full_age_grid > sub_ages[-1])
        grid_ext = np.concatenate([full_age_grid[floor_index], age_grid, full_age_grid[old_index]])
        ext_index = np.concatenate([floor_index, grid_index, old_index])
        log_age_grid = np.log10(age_grid)
        cmin = max(r["color"][0] for r in refs)
        cmax = min(r["color"][-1] for r in refs)
        valid_members = members[valid]
        for start in range(0, len(valid_members), chunk_size):
            index = valid_members[start:start + chunk_size]
            local = np.flatnonzero(np.isin(members, index))
            sample_color = color[index][:, None] + np.sqrt(2.0) * color_err[index][:, None] * hermite_nodes[None, :]
            supported = (sample_color >= cmin) & (sample_color <= cmax)
            weights = supported * hermite_weights[None, :]
            wsum = weights.sum(axis=1, keepdims=True)
            good = wsum[:, 0] > 0
            if not good.all():
                index, local, sample_color, weights, wsum = index[good], local[good], sample_color[good], weights[good], wsum[good]
            if not len(index):
                continue
            weights = weights / wsum
            reference_sample_mag = np.empty((len(refs), len(index), int(n_color_quadrature)))
            for k, r in enumerate(refs):
                reference_sample_mag[k] = np.interp(sample_color.ravel(), r["color"], r["mag"]).reshape(sample_color.shape)
            model_sample_mag = PchipInterpolator(log_sub_ages, reference_sample_mag, axis=0, extrapolate=False)(log_age_grid)
            if use_mh:   # de-metallicity: shift every model sequence by the star's PARSEC offset at that age and colour
                model_sample_mag = model_sample_mag + offset(
                    mh_star[index][None, :, None], log_age_grid[:, None, None], sample_color[None, :, :]
                )
            if scatter_mode == "nearest_two":
                sigma = np.sqrt(abs_g_err[index] ** 2 + star_sigma[local] ** 2 + systematic_mag ** 2)[None, :, None]
            else:
                scatter_grid = PchipInterpolator(log_sub_ages, ref_scatter[:, local], axis=0, extrapolate=False)(log_age_grid)
                sigma = np.sqrt(abs_g_err[index][None, :] ** 2 + scatter_grid ** 2 + systematic_mag ** 2)[:, :, None]
                if outside[local].any():
                    sigma_n2 = np.sqrt(abs_g_err[index] ** 2 + star_sigma[local] ** 2 + systematic_mag ** 2)[None, :, None]
                    sigma = np.where(outside[local][None, :, None], sigma_n2, sigma)
            residual = abs_g[index][None, :, None] - model_sample_mag
            sample_ll = -0.5 * (residual / sigma) ** 2 - np.log(sigma)
            if binary_fraction > 0:
                sample_ll = sample_ll + np.log1p(-binary_fraction)
                for dm_k, lw_k in zip(binary_dm, log_binary_w):
                    sample_ll = np.logaddexp(sample_ll, lw_k - 0.5 * ((residual + dm_k) / sigma) ** 2 - np.log(sigma))
            sample_ll = np.where(weights[None, :, :] > 0, sample_ll, -np.inf)
            sample_max = np.max(sample_ll, axis=2)
            log_likelihood = sample_max + np.log(np.sum(weights[None, :, :] * np.exp(sample_ll - sample_max[:, :, None]), axis=2))
            finite_ll = np.all(np.isfinite(log_likelihood), axis=0)
            if not finite_ll.all():
                index, local, log_likelihood = index[finite_ll], local[finite_ll], log_likelihood[:, finite_ll]
            if not len(index):
                continue
            density = np.exp(log_likelihood - np.max(log_likelihood, axis=0))
            map_age[index] = age_grid[np.argmax(log_likelihood, axis=0)]
            # flat below the youngest / above the oldest covering reference (unconstrained there)
            dens_ext = np.vstack([np.repeat(density[:1], len(floor_index), axis=0), density,
                                  np.repeat(density[-1:], len(old_index), axis=0)])
            age_025[index] = posterior_quantile(grid_ext, dens_ext, 0.025)
            age_16[index] = posterior_quantile(grid_ext, dens_ext, 0.16)
            age_50[index] = posterior_quantile(grid_ext, dens_ext, 0.50)
            age_84[index] = posterior_quantile(grid_ext, dens_ext, 0.84)
            age_975[index] = posterior_quantile(grid_ext, dens_ext, 0.975)
            sigma_used[index] = star_sigma[local] if scatter_mode == "nearest_two" else np.where(outside[local], star_sigma[local], np.nan)
            if posterior_out is not None:
                posterior_out[np.ix_(index, ext_index)] = (dens_ext / dens_ext.sum(axis=0)).T.astype(np.float32)
        n_groups_done += 1
        if verbose:
            print(f"coverage group {n_groups_done}/{len(patterns)} ({len(refs)} references, {len(members):,} stars) done", flush=True)

    covered = n_refs_used >= 2
    brighter = covered & (abs_g < youngest_mag)
    fainter = covered & (abs_g > oldest_mag)
    limit = np.where(fainter, "old", np.where(brighter, "young", "")).astype(object)

    result = {
        "age_map_myr": map_age,
        "age_myr_025": age_025,
        "age_myr_16": age_16,
        "age_myr_50": age_50,
        "age_myr_84": age_84,
        "age_myr_975": age_975,
        "limit": limit,
        "has_posterior": np.isfinite(map_age),
        "n_references": n_refs_used,
        "scatter_sigma_mag": sigma_used,
        "youngest_age": youngest_age,
        "oldest_age": oldest_age,
        "youngest_mag": youngest_mag,
        "oldest_mag": oldest_mag,
        "mh": mh_star,
        "metallicity_clipped": metallicity_clipped,
        "age_grid_myr": full_age_grid,
    }
    if return_posterior or posterior_out is not None:
        result["posterior"] = posterior_out
    return result


def star_age(color, abs_g, color_err, abs_g_err, mh=None, **kwargs):
    """Age of a single star: infer_ages() for one star, returned as scalars plus its posterior.

    Returns a dict with age_map_myr, age_myr_16/50/84 (and 2.5/97.5), limit, n_references,
    mh, metallicity_clipped, posterior (n_age,) and age_grid_myr.  Keyword arguments are
    passed to infer_ages().
    """
    result = infer_ages([color], [abs_g], [color_err], [abs_g_err], mh=None if mh is None else [mh], **kwargs)
    out = {}
    for key, value in result.items():
        if key == "age_grid_myr":
            out[key] = value
        elif key == "posterior":
            out[key] = np.asarray(value[0])
        else:
            item = value[0]
            out[key] = item.item() if hasattr(item, "item") else item
    out.update(color=float(color), abs_g=float(abs_g), color_err=float(color_err), abs_g_err=float(abs_g_err))
    return out
