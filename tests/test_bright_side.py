"""Stars brighter than the youngest covering sequence get that sequence's (filled) width at every age."""

import numpy as np
import pytest

import cmd_interp as ci
from cmd_interp.scatter import exact_local_scatter, fill_scatter


def _mass_at_or_below_youngest(star):
    return star["posterior"][star["age_grid_myr"] <= star["youngest_age"]].sum()


def _youngest_seq(color):
    youngest = ci.load_references()[0]
    return float(np.interp(color, youngest["color"], youngest["mag"]))


@pytest.mark.parametrize("color,height", [(2.0, 2.9), (1.5, 1.5), (2.5, 1.2), (2.8, 3.0)])
def test_far_above_youngest_is_flat_young(color, height):
    star = ci.star_age(color, _youngest_seq(color) - height, 0.01, 0.03)
    assert star["limit"] == "young"
    assert star["age_map_myr"] == star["youngest_age"]
    assert _mass_at_or_below_youngest(star) >= 0.95
    age, post = star["age_grid_myr"], star["posterior"]
    floor = post[age < star["youngest_age"]]
    np.testing.assert_allclose(floor, floor[0], rtol=1e-6)            # flat plateau below the youngest age
    tail = post[age >= star["youngest_age"]]
    assert np.all(np.diff(tail) <= 1e-12)                            # and only a falling tail above it


def test_just_above_youngest_keeps_older_tail():
    star = ci.star_age(2.0, _youngest_seq(2.0) - 0.05, 0.01, 0.03)
    assert star["limit"] == "young"
    assert 0.02 < 1.0 - _mass_at_or_below_youngest(star) < 0.98
    assert star["age_myr_84"] > star["youngest_age"]


def test_bright_side_uses_filled_youngest_width_and_own_error():
    color = 1.0                                                       # UBC 17a window has < 8 members here
    youngest = ci.load_references()[0]
    filled = fill_scatter([color], youngest, exact_local_scatter([color], youngest))[0]
    star = ci.star_age(color, _youngest_seq(color) - 1.0, 0.01, 0.03)
    assert np.isclose(star["scatter_sigma_mag"], filled)
    a = ci.star_age(color, _youngest_seq(color) - 0.1, 0.01, 0.01)
    b = ci.star_age(color, _youngest_seq(color) - 0.1, 0.01, 0.10)  # larger own G error -> broader posterior
    assert (b["age_myr_84"] - b["age_myr_16"]) > (a["age_myr_84"] - a["age_myr_16"])
