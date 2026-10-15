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
    refs = ci.load_references(reference_ages={"IC_2602": 52.5, "Melotte_22": 127.4}, include_member_ids=True)
    assert [r["cluster"] for r in refs] == ["IC_2602", "Melotte_22"]
    assert refs[0]["member_source_ids"].dtype.kind == "i"
    with pytest.raises(ValueError):
        ci.load_references(reference_ages={"IC_2602": 52.5})


def test_excluding_a_reference_changes_the_grid():
    refs = ci.load_references(reference_ages={k: v for k, v in ci.REFERENCE_AGES.items() if k != "Melotte_111"})
    star = ci.star_age(2.0, 8.0, 0.01, 0.03, references=refs)
    assert star["oldest_age"] == 127.4
