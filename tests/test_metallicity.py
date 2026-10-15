import warnings

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

import cmd_interp as ci

STARS = dict(color=np.array([1.2, 1.6, 2.0, 2.4, 2.8, 3.0]), abs_g=np.array([5.5, 6.5, 7.6, 8.6, 9.6, 10.0]),
             color_err=np.full(6, 0.01), abs_g_err=np.full(6, 0.03))


def _infer(**kwargs):
    return ci.infer_ages(STARS["color"], STARS["abs_g"], STARS["color_err"], STARS["abs_g_err"], **kwargs)


def test_solar_offset_is_zero():
    offset = ci.parsec_offset()
    assert np.all(offset(0.0, np.linspace(0.7, 3.0, 20), np.linspace(0.8, 2.9, 20)) == 0.0)


def test_zero_metallicity_is_identical_to_none():
    a, b = _infer(), _infer(mh=np.zeros(6))
    np.testing.assert_array_equal(a["posterior"], b["posterior"])


def test_per_star_metallicity_matches_separate_runs():
    mh = np.array([-0.2, 0.1, -0.2, 0.1, -0.2, 0.1])
    mixed = _infer(mh=mh)
    for z in (-0.2, 0.1):
        sel = mh == z
        single = ci.infer_ages(STARS["color"][sel], STARS["abs_g"][sel], STARS["color_err"][sel],
                               STARS["abs_g_err"][sel], mh=z)
        np.testing.assert_array_equal(mixed["age_myr_50"][sel], single["age_myr_50"])


@pytest.mark.parametrize("age,colour,mh", [(10, 2.0, -0.3), (20, 2.5, 0.2), (40, 2.0, 0.2), (100, 2.0, -0.3), (20, 1.5, -0.3)])
def test_shifted_star_recovers_its_age(age, colour, mh):
    """A star on the solar sequence of age t, moved by the PARSEC offset for its [M/H], is dated back to t."""
    refs = ci.load_references()
    mags = [np.interp(colour, r["color"], r["mag"]) for r in refs]
    m_solar = PchipInterpolator(np.log10([r["age_myr"] for r in refs]), mags)(np.log10(age))
    m_obs = m_solar + ci.parsec_offset()(mh, np.log10(age), colour)
    corrected = ci.star_age(colour, m_obs, 0.01, 0.02, mh=mh)
    assert corrected["age_myr_16"] <= age <= corrected["age_myr_84"]
    naive = ci.star_age(colour, m_obs, 0.01, 0.02)
    assert abs(np.log10(corrected["age_map_myr"] / age)) < abs(np.log10(naive["age_map_myr"] / age)) + 0.1


def test_clipping_is_flagged():
    with pytest.warns(UserWarning, match="outside the PARSEC offset table"):
        result = _infer(mh=np.array([0.1, 0.1, 0.1, 0.1, 0.1, -0.8]))
    assert list(result["metallicity_clipped"]) == [False, False, False, False, False, True]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _infer(mh=np.array([0.1, 0.1, 0.1, 0.1, 0.1, np.nan]))  # NaN = solar, and the 3.0 star is solar -> no warning
