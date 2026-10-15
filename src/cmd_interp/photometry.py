"""De-extincted absolute Gaia DR3 photometry and its uncertainties.

``deredden_gaia`` turns gaia_source columns into the inputs of the age fit: absolute,
extinction-corrected M_G,0 and (BP-RP)_0 with uncertainties (columns phot_g_mean_mag_0,
BpmRp0, phot_g_mean_mag_0_err, BpmRp0_err), plus the BP and RP magnitudes.

Extinction
==========
A_0 (extinction at 550 nm) is either supplied (a single value, one per star, or a column) or
taken from the EXPLORE 3D dust maps (``a0="dustmap"``, see cmd_interp.dustmap).  A_0 is
converted to A_G, A_BP and A_RP with the Fitzpatrick (2019) / Gaia EDR3 extinction-law
polynomials k_X((BP-RP)_0, A_0).  The upper-HRD table is used for stars with
G - 5 log10(d / 10 pc) < 5 and the main-sequence table otherwise.  The intrinsic colour is
iterated (up to 5 times, until it changes by < 0.01 mag).

Uncertainty budget
==================
X is one of G, BP, RP; F_X and sigma_F,X are the mean flux and its error; plx and sigma_plx
are the parallax and its error (mas); d = 1000/plx (pc); k_X = A_X / A_0; rho = dA_0/dd is
the extinction density (mag/pc) at the star; sigma_A0,map is the dust-map (or user-supplied)
A_0 error; sigma_zp,X is the Vega zero-point uncertainty (Riello et al. 2021).

    sigma_plx  = f(RUWE, G, plx) * sigma_plx,catalogue          (inflate_parallax_uncertainty)
    sigma_A0^2 = sigma_A0,map^2 + (rho d sigma_plx / plx)^2

    sigma^2(m_X,0) = (2.5 / ln 10)^2 (sigma_F,X / F_X)^2              [flux]
                   + sigma_zp,X^2                                       [zero point]
                   + [(5 / ln 10 + k_X rho d) sigma_plx / plx]^2        [parallax: DM + A_0(d)]
                   + k_X^2 sigma_A0,map^2                               [dust map]

    sigma^2((BP-RP)_0) = (2.5 / ln 10)^2 [(sigma_F,BP / F_BP)^2 + (sigma_F,RP / F_RP)^2]
                       + sigma_zp,BP^2 + sigma_zp,RP^2 + (k_BP - k_RP)^2 sigma_A0^2

The parallax enters each magnitude twice with the same sign: through the distance modulus,
and through A_0, which grows with distance at the local rate rho.  So the two derivatives
are summed before squaring.  In the colour the distance modulus cancels.  With a
user-supplied A_0 there is no dust gradient (rho = 0).

Not included: the Gaia parallax zero point (Lindegren et al. 2021) and its ~10 uas floor; a
photometric error floor; uncertainty of the extinction-law coefficients; dust beyond the
edge of the dust map; the covariance between M_G,0 and (BP-RP)_0 induced by sigma_A0.
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd

# Vega zero-point uncertainties, Riello et al. (2021)
ZP_ERR = {"G": 0.0027553202, "BP": 0.0027901700, "RP": 0.0037793818}
COEFFICIENT_NAMES = ["Intercept", "X", "X2", "X3", "A", "A2", "A3", "XA", "AX2", "XA2"]
UPPER_HRD_ABS_G = 5.0          # stars with G - 5 log10(d/10) brighter than this use the upper-HRD extinction law
N_ITERATIONS, TOLERANCE = 5, 0.01

GAIA_COLUMNS = ["parallax", "parallax_error", "ruwe",
                "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
                "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_bp_mean_flux", "phot_bp_mean_flux_error",
                "phot_rp_mean_flux", "phot_rp_mean_flux_error"]


@lru_cache(maxsize=2)
def _extinction_law(kind):
    """{band: 10 coefficients} of the Fitzpatrick (2019) EDR3 k_X(BP-RP, A_0) polynomial."""
    name = {"top": "Fitz19_EDR3_HRDTop.csv", "main": "Fitz19_EDR3_MainSequence.csv"}[kind]
    with resources.as_file(resources.files("cmd_interp") / "data" / "extinction" / name) as path:
        table = pd.read_csv(path)
    rows = table[table.Xname == "BPRP"].set_index("Kname")
    return {band: rows.loc[f"k{band}", COEFFICIENT_NAMES].to_numpy(dtype=float) for band in ("G", "BP", "RP")}


def extinction_coefficient(band, bprp0, a0, upper_hrd):
    """k_X = A_X / A_0 for band 'G', 'BP' or 'RP' (Fitzpatrick 2019, Gaia EDR3); vectorised.

    upper_hrd: boolean array selecting the upper-HRD coefficients (else main sequence)."""
    x, a, upper_hrd = np.broadcast_arrays(np.asarray(bprp0, dtype=float), np.asarray(a0, dtype=float), np.asarray(upper_hrd, dtype=bool))
    terms = np.stack([np.ones_like(x), x, x ** 2, x ** 3, a, a ** 2, a ** 3, a * x, a * x ** 2, x * a ** 2])
    k_main = np.tensordot(_extinction_law("main")[band], terms, axes=1)
    k_top = np.tensordot(_extinction_law("top")[band], terms, axes=1)
    return np.where(upper_hrd, k_top, k_main)


def _al_uncertainty_per_ccd(g_mag):
    """Along-scan uncertainty per CCD (mas), 'EDR3 adjusted' curve of Fig. 3 of arXiv:2206.05439."""
    g_vals = [4, 5, 6, 7, 8.2, 8.4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    sigma_eta = [0.4, 0.35, 0.15, 0.17, 0.23, 0.13, 0.13, 0.135, 0.125, 0.13, 0.15, 0.23, 0.36, 0.63, 1.05, 2.05, 4.1]
    return np.interp(g_mag, g_vals, sigma_eta)


def _f_high_ruwe(ruwe, parallax, g_mag, f0=3.73, alpha=2.77, beta=0.065, gamma=-0.056):
    """Parallax-error inflation factor for RUWE >= 1.25 (Eq. 3 of arXiv:2206.05439)."""
    sigma_eta = _al_uncertainty_per_ccd(g_mag)
    fmax = f0 * (parallax / 10) ** beta * (sigma_eta / 0.1) ** gamma
    return 1 + (fmax - 1) * (1 - np.exp(-alpha * (ruwe - 1)))


def inflate_parallax_uncertainty(ruwe, parallax, parallax_error, g_mag):
    """Catalogue parallax errors inflated as a function of RUWE, G and parallax (vectorised).

    RUWE < 1: unchanged.  1 < RUWE < 1.25: the G-dependent factor
    A exp(-(G - G0)^2 / b^2) + p0 + p1 G + p2 G^2.  Otherwise: _f_high_ruwe."""
    ruwe, parallax, parallax_error, g_mag = np.broadcast_arrays(*(np.asarray(v, dtype=float)
                                                                  for v in (ruwe, parallax, parallax_error, g_mag)))
    with np.errstate(invalid="ignore", divide="ignore"):
        f_low = 0.21 * np.exp(-((g_mag - 12.65) ** 2) / 0.90 ** 2) + 1.141 + 0.0040 * g_mag - 0.00062 * g_mag ** 2
        f_high = _f_high_ruwe(ruwe, parallax, g_mag)
    f = np.where(ruwe < 1, 1.0, np.where((ruwe > 1) & (ruwe < 1.25), f_low, f_high))
    return f * parallax_error


def _column_or_value(frame, value, name):
    if isinstance(value, str):
        if value not in frame:
            raise KeyError(f"{name} column {value!r} not in the table")
        return frame[value].to_numpy(dtype=float)
    return np.broadcast_to(np.asarray(value, dtype=float), (len(frame),)).astype(float)


def deredden_gaia(data, a0=0.0, *, a0_err=0.0, dust_map=None, inflate_parallax_error=True):
    """De-extincted absolute Gaia photometry and uncertainties for a table of stars.

    data: DataFrame with the gaia_source columns parallax, parallax_error, ruwe,
        phot_{g,bp,rp}_mean_mag, phot_{g,bp,rp}_mean_flux and phot_{g,bp,rp}_mean_flux_error
        (plus ra, dec for the dust map).
    a0: extinction A_0 (mag at 550 nm): a number for every star, an array, a column name, or
        "dustmap" to integrate the EXPLORE dust maps to each star's distance.  Default 0.
    a0_err: uncertainty of a user-supplied A_0 (number, array or column name); default 0.
        With "dustmap" it is replaced by the dust-map error.
    dust_map: an ExploreDustMap (or None to locate the cubes from the environment).
    inflate_parallax_error: apply inflate_parallax_uncertainty (needs ruwe).

    Returns a copy of the table with distance, extinction_a0, extinction_a0_err_cube,
    extinction_a0_gradient, extinction_a0_err, parallax_error_inflated,
    phot_{g,bp,rp}_mean_mag_0 (absolute, de-extincted), BpmRp0 and their _err columns.
    Stars with parallax <= 0 get NaN.
    """
    frame = data.copy()
    need = [c for c in GAIA_COLUMNS if c not in frame and not (c == "ruwe" and not inflate_parallax_error)]
    if need:
        raise KeyError(f"Missing gaia_source column(s): {', '.join(need)}")
    n = len(frame)
    parallax = frame.parallax.to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        distance = np.where(parallax > 0, 1000.0 / parallax, np.nan)

    if isinstance(a0, str) and a0 == "dustmap":
        if dust_map is None:
            from .dustmap import ExploreDustMap
            dust_map = ExploreDustMap()
        if "ra" not in frame or "dec" not in frame:
            raise KeyError("The dust map needs ra and dec columns")
        a0_values, a0_err_map, gradient = dust_map.extinction(frame.ra.to_numpy(float), frame.dec.to_numpy(float), distance)
    else:
        a0_values = _column_or_value(frame, a0, "a0")
        a0_err_map = _column_or_value(frame, a0_err, "a0_err")
        gradient = np.zeros(n)

    g, bp, rp = (frame[f"phot_{b}_mean_mag"].to_numpy(dtype=float) for b in ("g", "bp", "rp"))
    with np.errstate(invalid="ignore", divide="ignore"):
        dm = 5 * (np.log10(distance) - 1)
        upper_hrd = (g - 5 * np.log10(distance / 10)) < UPPER_HRD_ABS_G

        # intrinsic colour, iterated with the BP and RP coefficients
        bprp = bp - rp
        bprp0 = np.full(n, np.nan)
        active = np.ones(n, dtype=bool)
        for _ in range(N_ITERATIONS):
            new = (bp - extinction_coefficient("BP", bprp, a0_values, upper_hrd) * a0_values) - (
                rp - extinction_coefficient("RP", bprp, a0_values, upper_hrd) * a0_values)
            converged = active & (np.abs(new - bprp) < TOLERANCE)
            bprp0[converged] = new[converged]
            active &= ~converged
            bprp = np.where(active, new, bprp)
        bprp0[active] = bprp[active]

        if inflate_parallax_error:
            plx_err = inflate_parallax_uncertainty(frame.ruwe.to_numpy(float), parallax, frame.parallax_error.to_numpy(float), g)
        else:
            plx_err = frame.parallax_error.to_numpy(dtype=float)
        frac_plx_err = plx_err / parallax
        a0_err_total = np.sqrt(a0_err_map ** 2 + (gradient * distance * frac_plx_err) ** 2)

        k, mag0, mag_err, flux_term = {}, {}, {}, {}
        for band, mag in (("G", g), ("BP", bp), ("RP", rp)):
            k[band] = extinction_coefficient(band, bprp0, a0_values, upper_hrd)
            flux = frame[f"phot_{band.lower()}_mean_flux"].to_numpy(dtype=float)
            flux_err = frame[f"phot_{band.lower()}_mean_flux_error"].to_numpy(dtype=float)
            flux_term[band] = (2.5 / (np.log(10) * flux)) ** 2 * flux_err ** 2
            plx_term = ((5 / np.log(10) + k[band] * gradient * distance) * frac_plx_err) ** 2
            mag_err[band] = np.sqrt(flux_term[band] + ZP_ERR[band] ** 2 + plx_term + k[band] ** 2 * a0_err_map ** 2)
            mag0[band] = mag - dm - k[band] * a0_values
        color_err = np.sqrt(flux_term["BP"] + flux_term["RP"] + ZP_ERR["BP"] ** 2 + ZP_ERR["RP"] ** 2
                            + ((k["BP"] - k["RP"]) * a0_err_total) ** 2)

    frame["distance"] = distance
    frame["extinction_a0"] = a0_values
    frame["extinction_a0_err_cube"] = a0_err_map
    frame["extinction_a0_gradient"] = gradient
    frame["extinction_a0_err"] = a0_err_total
    frame["parallax_error_inflated"] = plx_err
    frame["phot_g_mean_mag_0"] = mag0["G"]
    frame["phot_bp_mean_mag_0"] = mag0["BP"]
    frame["phot_rp_mean_mag_0"] = mag0["RP"]
    frame["BpmRp0"] = mag0["BP"] - mag0["RP"]
    frame["phot_g_mean_mag_0_err"] = mag_err["G"]
    frame["phot_bp_mean_mag_0_err"] = mag_err["BP"]
    frame["phot_rp_mean_mag_0_err"] = mag_err["RP"]
    frame["BpmRp0_err"] = color_err
    return frame
