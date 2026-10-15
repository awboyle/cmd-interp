import numpy as np
import pandas as pd

import cmd_interp as ci


def test_red_stars_get_nan():
    result = ci.infer_ages([2.9, 3.0, 3.01, 3.5], [9.0, 9.2, 9.3, 11.0], 0.01, 0.03, mh=0.1)
    assert list(result["color_out_of_range"]) == [False, False, True, True]
    assert np.isfinite(result["age_map_myr"][:2]).all() and np.isnan(result["age_map_myr"][2:]).all()
    assert list(result["limit"][2:]) == ["", ""] and not result["metallicity_clipped"][2:].any()
    assert np.all(result["posterior"][2:] == 0)


def test_limit_is_adjustable():
    star = ci.star_age(3.2, 11.0, 0.01, 0.03, max_color=3.5)
    assert not star["color_out_of_range"] and np.isfinite(star["age_map_myr"])
    blue = ci.star_age(0.5, 3.0, 0.01, 0.03, min_color=0.8)
    assert blue["color_out_of_range"] and np.isnan(blue["age_map_myr"])


def test_table_flags_red_stars():
    data = pd.DataFrame({"BpmRp0": [2.0, 3.4], "phot_g_mean_mag_0": [8.0, 11.0], "BpmRp0_err": 0.01, "phot_g_mean_mag_0_err": 0.03})
    out, _ = ci.ages_for_table(data)
    assert list(out.cmd_age_color_out_of_range) == [False, True]
    assert np.isfinite(out.cmd_age_myr[0]) and np.isnan(out.cmd_age_myr[1]) and out.cmd_age_limit[1] == ""
