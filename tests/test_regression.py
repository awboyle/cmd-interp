"""Regression test of the default model (hybrid scatter, 0.05 mag systematic).

tests/data/default_model_regression.npz holds frozen reference outputs for 200 synthetic stars. Stars brighter
than the youngest covering sequence are excluded here and tested in test_bright_side.py."""

from pathlib import Path

import numpy as np

import cmd_interp as ci

FIXTURE = Path(__file__).parent / "data" / "default_model_regression.npz"


def test_defaults_reproduce_frozen_outputs():
    ref = np.load(FIXTURE)
    result = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"])
    red = ref["color"] > ci.MAX_COLOR                                  # no ages redward of (BP-RP)_0 = 3
    assert red.sum() > 5 and np.all(np.isnan(result["age_map_myr"][red])) and result["color_out_of_range"][red].all()
    above = (result["n_references"] >= 2) & (ref["abs_g"] < result["youngest_mag"])
    keep = ~above & ~red
    assert keep.sum() > 100 and above.sum() > 10
    for key in ("age_map_myr", "age_myr_025", "age_myr_16", "age_myr_50", "age_myr_84", "age_myr_975"):
        np.testing.assert_allclose(result[key][keep], ref[key][keep], rtol=1e-9, equal_nan=True, err_msg=key)
    np.testing.assert_allclose(result["posterior"][keep], ref["posterior"][keep], rtol=1e-5, atol=1e-8)


def test_explicit_settings_match_defaults():
    ref = np.load(FIXTURE)
    a = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"], return_posterior=False)
    b = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"], return_posterior=False,
                      scatter_mode="hybrid", systematic_mag=0.05, mh=0.0)
    np.testing.assert_array_equal(a["age_myr_50"], b["age_myr_50"])
