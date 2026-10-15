import numpy as np
import pytest

import cmd_interp as ci


def test_bundled_library():
    refs = ci.load_references()
    assert [r["cluster"] for r in refs] == list(ci.REFERENCE_AGES)
    ages = [r["age_myr"] for r in refs]
    assert ages == sorted(ages)
    for r in refs:
        assert np.all(np.diff(r["color"]) > 0)
        assert len(r["member_color"]) > 10
        assert np.all(np.abs(r["member_residual"]) <= 0.3)


def test_subset_and_ids():
    refs = ci.load_references(reference_ages={"IC_2602": 52.5, "Melotte_22": 127.4}, use_default_ages=False,
                              include_member_ids=True)
    assert [r["cluster"] for r in refs] == ["IC_2602", "Melotte_22"]
    assert refs[0]["member_source_ids"].dtype.kind == "i"
    with pytest.raises(ValueError):
        ci.load_references(reference_ages={"IC_2602": 52.5}, use_default_ages=False)


def test_excluding_a_reference_changes_the_grid():
    star = ci.star_age(2.0, 8.0, 0.01, 0.03, exclude_references=["Coma Ber"])
    assert star["oldest_age"] == 127.4
    with pytest.raises(ValueError):
        ci.load_references(exclude=["NGC 9999"])


def test_override_one_age_keeps_the_others():
    refs = ci.load_references(reference_ages={"beta_Pic": 21.0})
    assert {r["cluster"]: r["age_myr"] for r in refs} == {**ci.REFERENCE_AGES, "beta_Pic": 21.0}
    assert ci.REFERENCE_AGES["beta_Pic"] == 24.3                       # defaults untouched


@pytest.mark.parametrize("alias,canonical", [("Pleiades", "Melotte_22"), ("alpha Per", "Melotte_20"),
                                             ("beta pic", "beta_Pic"), ("Coma_Ber", "Melotte_111"), ("ic 2602", "IC_2602")])
def test_aliases(alias, canonical):
    ages = ci.resolve_reference_ages({alias: 99.0})
    assert ages[canonical] == 99.0 and len(ages) == len(ci.REFERENCE_AGES)


def test_reordered_ages_are_sorted_and_validated():
    ages = ci.resolve_reference_ages({"IC_2602": 90.0})               # now older than alpha Per (79)
    assert list(ages)[2:4] == ["Melotte_20", "IC_2602"]
    for bad in ({"IC_2602": -1}, {"IC_2602": 79.0}):                   # non-positive / duplicate age
        with pytest.raises(ValueError):
            ci.resolve_reference_ages(bad)


def test_new_ages_change_the_inferred_age():
    base = ci.star_age(2.0, 8.0, 0.01, 0.03)
    shifted = ci.star_age(2.0, 8.0, 0.01, 0.03, reference_ages={"Pleiades": 112, "alpha Per": 85})
    assert shifted["age_map_myr"] > base["age_map_myr"]
    with pytest.raises(ValueError):
        ci.infer_ages([2.0], [8.0], [0.01], [0.03], references=ci.load_references(), reference_ages={"Pleiades": 112})
