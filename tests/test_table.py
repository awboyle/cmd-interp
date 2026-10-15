import numpy as np
import pandas as pd
import pytest

import cmd_interp as ci


@pytest.fixture
def table():
    return pd.DataFrame({
        "BpmRp0": [2.0, 2.0, 2.0, 1.8, np.nan, 5.0],
        "phot_g_mean_mag_0": [8.0, 4.0, 11.0, 7.4, 8.0, 14.0],
        "BpmRp0_err": [0.01, 0.01, 0.01, 0.02, 0.01, 0.01],
        "phot_g_mean_mag_0_err": [0.03, np.nan, 0.03, 0.02, 0.03, 0.03],
        "feh": [0.0, np.nan, 0.1, -0.2, 0.0, 0.0],
        "is_pms": [True, True, True, False, True, True],
    })


def test_columns_and_limits(table):
    out, result = ci.ages_for_table(table)
    assert len(out) == len(table)
    assert list(result["row"]) == [0, 1, 2, 3, 5]                        # NaN colour row is not fitted
    assert out.loc[0, "cmd_age_limit"] == "" and np.isfinite(out.loc[0, "cmd_age_median_myr"])
    assert out.loc[1, "cmd_age_limit"] == "young" and out.loc[1, "cmd_age_myr"] == 4.6
    assert out.loc[2, "cmd_age_limit"] == "old" and out.loc[2, "cmd_age_lower_myr"] == 675.0
    assert np.isnan(out.loc[4, "cmd_age_myr"]) and np.isnan(out.loc[5, "cmd_age_myr"])
    assert out.loc[5, "cmd_age_n_references"] == 0
    assert result["posterior"].shape == (5, len(result["age_grid_myr"]))
    star = ci.star_age(2.0, 8.0, 0.01, 0.03)
    assert out.loc[0, "cmd_age_myr"] == star["age_map_myr"]


def test_g_error_fallback(table):
    out, _ = ci.ages_for_table(table, abs_g_err=0.03)
    alt = table.copy()
    alt.loc[1, "phot_g_mean_mag_0_err"] = 0.03
    out2, _ = ci.ages_for_table(alt)
    assert out.loc[1, "cmd_age_scatter_mag"] == out2.loc[1, "cmd_age_scatter_mag"] or np.isnan(out.loc[1, "cmd_age_scatter_mag"])
    np.testing.assert_array_equal(out.cmd_age_myr, out2.cmd_age_myr)


def test_fit_mask_and_metallicity_column(table):
    out, result = ci.ages_for_table(table, fit_mask="is_pms", mh_col="feh")
    assert list(result["row"]) == [0, 1, 2, 5]
    assert np.isnan(out.loc[3, "cmd_age_myr"])
    assert list(out.cmd_age_mh) == [0.0, 0.0, 0.1, -0.2, 0.0, 0.0]
    with pytest.raises(ValueError):
        ci.ages_for_table(table, mh=0.1, mh_col="feh")


def test_systematic_is_user_settable(table):
    narrow, _ = ci.ages_for_table(table, systematic_mag=0.0)
    wide, _ = ci.ages_for_table(table, systematic_mag=0.2)
    w_narrow = narrow.cmd_age_upper_myr[0] - narrow.cmd_age_lower_myr[0]
    w_wide = wide.cmd_age_upper_myr[0] - wide.cmd_age_lower_myr[0]
    assert w_wide > w_narrow
