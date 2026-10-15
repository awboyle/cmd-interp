"""The package reproduces the original research pipeline (hybrid scatter, 0.05 mag systematic).

tests/data/hyb05_regression.npz was produced by the original batch_empirical_isochrone_ages.infer_batch
on 200 synthetic stars."""

from pathlib import Path

import numpy as np

import cmd_interp as ci

FIXTURE = Path(__file__).parent / "data" / "hyb05_regression.npz"


def test_defaults_reproduce_reference_pipeline():
    ref = np.load(FIXTURE)
    result = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"])
    for key in ("age_map_myr", "age_myr_025", "age_myr_16", "age_myr_50", "age_myr_84", "age_myr_975"):
        np.testing.assert_allclose(result[key], ref[key], rtol=1e-9, equal_nan=True, err_msg=key)
    np.testing.assert_allclose(result["posterior"], ref["posterior"], rtol=1e-5, atol=1e-8)
    assert np.isfinite(ref["age_map_myr"]).sum() > 100


def test_explicit_settings_match_defaults():
    ref = np.load(FIXTURE)
    a = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"], return_posterior=False)
    b = ci.infer_ages(ref["color"], ref["abs_g"], ref["color_err"], ref["abs_g_err"], return_posterior=False,
                      scatter_mode="hybrid", systematic_mag=0.05, mh=0.0)
    np.testing.assert_array_equal(a["age_myr_50"], b["age_myr_50"])
