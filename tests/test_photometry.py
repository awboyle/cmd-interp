import numpy as np
import pandas as pd
import pytest

import cmd_interp as ci
from cmd_interp.photometry import ZP_ERR


def _stars():
    return pd.DataFrame({
        "ra": [56.75, 130.1, 200.0], "dec": [24.1, 19.6, -30.0],
        "parallax": [7.4, 5.5, 2.0], "parallax_error": [0.02, 0.05, 0.1], "ruwe": [0.95, 1.1, 2.0],
        "phot_g_mean_mag": [12.0, 4.0, 15.0], "phot_bp_mean_mag": [12.6, 4.2, 15.9], "phot_rp_mean_mag": [11.3, 3.7, 14.1],
        "phot_g_mean_flux": [3e5, 5e8, 2e4], "phot_g_mean_flux_error": [60, 2e5, 15],
        "phot_bp_mean_flux": [1e5, 2e8, 5e3], "phot_bp_mean_flux_error": [80, 1e5, 20],
        "phot_rp_mean_flux": [2e5, 3e8, 2e4], "phot_rp_mean_flux_error": [90, 1.5e5, 30],
    })


def test_zero_extinction_is_the_distance_modulus():
    s = _stars()
    out = ci.deredden_gaia(s, a0=0.0)
    dm = 5 * (np.log10(1000 / s.parallax) - 1)
    np.testing.assert_allclose(out.phot_g_mean_mag_0, s.phot_g_mean_mag - dm)
    np.testing.assert_allclose(out.BpmRp0, s.phot_bp_mean_mag - s.phot_rp_mean_mag)
    # error budget: flux + zero point + parallax (no dust terms)
    frac = out.parallax_error_inflated / s.parallax
    expected = np.sqrt((2.5 / np.log(10)) ** 2 * (s.phot_g_mean_flux_error / s.phot_g_mean_flux) ** 2 + ZP_ERR["G"] ** 2
                       + (5 / np.log(10) * frac) ** 2)
    np.testing.assert_allclose(out.phot_g_mean_mag_0_err, expected)
    colour = np.sqrt((2.5 / np.log(10)) ** 2 * ((s.phot_bp_mean_flux_error / s.phot_bp_mean_flux) ** 2
                                                  + (s.phot_rp_mean_flux_error / s.phot_rp_mean_flux) ** 2)
                     + ZP_ERR["BP"] ** 2 + ZP_ERR["RP"] ** 2)
    np.testing.assert_allclose(out.BpmRp0_err, colour)


def test_extinction_reddens_and_dims_and_a0_error_propagates():
    s = _stars()
    clean, dusty = ci.deredden_gaia(s, a0=0.0), ci.deredden_gaia(s, a0=0.5, a0_err=0.1)
    assert np.all(dusty.phot_g_mean_mag_0 < clean.phot_g_mean_mag_0)
    assert np.all(dusty.BpmRp0 < clean.BpmRp0)
    assert np.all(dusty.BpmRp0_err > clean.BpmRp0_err) and np.all(dusty.phot_g_mean_mag_0_err > clean.phot_g_mean_mag_0_err)
    k_g = ci.extinction_coefficient("G", dusty.BpmRp0, 0.5, np.array([False, True, False]))
    # k_G is evaluated at the iterated colour, which differs slightly from the final output BpmRp0
    np.testing.assert_allclose(clean.phot_g_mean_mag_0 - dusty.phot_g_mean_mag_0, k_g * 0.5, rtol=1e-3)
    assert 0.6 < k_g[0] < 1.0                                           # A_G / A_0 for a K dwarf
    per_star = ci.deredden_gaia(s.assign(a0=[0.0, 0.5, 0.0]), a0="a0")
    np.testing.assert_allclose(per_star.phot_g_mean_mag_0[[0, 2]], clean.phot_g_mean_mag_0[[0, 2]])


def test_parallax_inflation_regimes():
    err = ci.inflate_parallax_uncertainty([0.9, 1.1, 1.0, 2.0], [5.0] * 4, [0.1] * 4, [12.0] * 4)
    assert err[0] == 0.1                                               # RUWE < 1: unchanged
    assert 0.1 < err[1] < 0.2                                          # 1 < RUWE < 1.25: G-dependent factor
    assert np.isclose(err[2], 0.1)                                     # RUWE = 1: high-RUWE formula gives 1
    assert err[3] > err[1]                                             # RUWE 2: strongly inflated
    out = ci.deredden_gaia(_stars(), a0=0.0, inflate_parallax_error=False)
    np.testing.assert_array_equal(out.parallax_error_inflated, _stars().parallax_error)


def test_negative_parallax_gives_nan_and_missing_columns_raise():
    s = _stars(); s.loc[0, "parallax"] = -1.0
    out = ci.deredden_gaia(s, a0=0.1)
    assert np.isnan(out.loc[0, "phot_g_mean_mag_0"]) and np.isfinite(out.loc[1, "phot_g_mean_mag_0"])
    with pytest.raises(KeyError):
        ci.deredden_gaia(s.drop(columns="phot_bp_mean_flux"), a0=0.1)


def test_output_feeds_the_age_fit():
    out = ci.deredden_gaia(_stars(), a0=0.05)
    table, _ = ci.ages_for_table(out)
    assert "cmd_age_myr" in table


def test_dust_map_on_a_uniform_cube(tmp_path):
    h5py = pytest.importorskip("h5py")
    pytest.importorskip("astropy")
    density, err_density = 1e-3, 2e-4                                  # mag/pc
    for name, value in (("values.h5", density), ("errors.h5", err_density)):
        with h5py.File(tmp_path / name, "w") as hf:
            d = hf.create_dataset("explore/cube_datas", data=np.full((101, 101, 101), value, dtype=np.float32))
            d.attrs["sun_position"] = np.array([50.0, 50.0, 50.0])     # Sun at the centre, cube spans +-1000 pc
            d.attrs["gridstep_values"] = np.array([20, 20, 20])
    dust = ci.ExploreDustMap(tmp_path / "values.h5", tmp_path / "errors.h5")
    s = _stars().iloc[[0, 2]].reset_index(drop=True)                  # 135 pc and 500 pc
    out = ci.deredden_gaia(s, a0="dustmap", dust_map=dust)
    # the line-of-sight sum includes the sample at d = 0, so a uniform cube gives rho (d + step)
    np.testing.assert_allclose(out.extinction_a0, density * (out.distance + 5), rtol=1e-3)
    np.testing.assert_allclose(out.extinction_a0_err_cube, err_density * (out.distance + 5), rtol=1e-3)
    np.testing.assert_allclose(out.extinction_a0_gradient, density, rtol=1e-5)
    assert np.all(out.extinction_a0_err > out.extinction_a0_err_cube)  # distance term added
