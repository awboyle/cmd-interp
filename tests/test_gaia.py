"""ages_from_gaia without the network (gaia_table / patched query); one optional live test."""
import os
import warnings

import numpy as np
import pandas as pd
import pytest

import cmd_interp as ci
import cmd_interp.gaia as gaia_module


def _gaia_rows():
    return pd.DataFrame({
        "SOURCE_ID": np.array([111, 222, 333], dtype=np.int64), "ra": [56.75, 130.1, 200.0], "dec": [24.1, 19.6, -30.0],
        "parallax": [7.4, 5.5, 2.0], "parallax_error": [0.02, 0.05, 0.1], "ruwe": [0.95, 1.1, 1.05],
        "phot_g_mean_mag": [13.9, 12.5, 15.0], "phot_bp_mean_mag": [15.1, 13.2, 15.9], "phot_rp_mean_mag": [12.8, 11.7, 14.1],
        "phot_g_mean_flux": [5e4, 2e5, 2e4], "phot_g_mean_flux_error": [20, 40, 15],
        "phot_bp_mean_flux": [1e4, 6e4, 5e3], "phot_bp_mean_flux_error": [30, 50, 20],
        "phot_rp_mean_flux": [4e4, 2e5, 2e4], "phot_rp_mean_flux_error": [40, 60, 30],
    })


def test_pipeline_matches_the_individual_steps():
    rows = _gaia_rows()
    table, result = ci.ages_from_gaia([333, 111, 999], a0=[0.1, 0.05, 0.2], mh=[0.0, 0.1, 0.0], gaia_table=rows)
    assert list(table.dr3_source_id) == [333, 111, 999]                     # input order, missing ID kept
    assert list(table.found_in_gaia) == [True, True, False]
    assert np.isnan(table.cmd_age_myr[2]) and list(result["row"]) == [0, 1]
    phot = ci.deredden_gaia(rows.rename(columns={"SOURCE_ID": "dr3_source_id"}).set_index("dr3_source_id").loc[[333, 111]]
                            .reset_index().assign(a0=[0.1, 0.05]), a0="a0")
    np.testing.assert_allclose(table.BpmRp0[:2], phot.BpmRp0)
    direct = ci.infer_ages(phot.BpmRp0, phot.phot_g_mean_mag_0, phot.BpmRp0_err, phot.phot_g_mean_mag_0_err, mh=[0.0, 0.1])
    np.testing.assert_array_equal(result["age_myr_50"], direct["age_myr_50"])
    np.testing.assert_array_equal(table.cmd_age_mh[:2], [0.0, 0.1])
    assert "_a0" not in table and "_mh" not in table


def test_a0_is_required_and_lengths_are_checked():
    with pytest.raises(TypeError):
        ci.ages_from_gaia([111], gaia_table=_gaia_rows())
    with pytest.raises(ValueError):
        ci.ages_from_gaia([111, 222], a0=[0.1, 0.2, 0.3], gaia_table=_gaia_rows())


def test_single_star_and_patched_query(monkeypatch):
    calls = []
    monkeypatch.setattr(gaia_module, "query_gaia_dr3",
                        lambda ids, **kw: calls.append(list(ids)) or _gaia_rows().rename(columns={"SOURCE_ID": "dr3_source_id"}))
    star = ci.star_age_from_gaia(111, a0=0.05, systematic_mag=0.08)
    assert calls == [[111]]
    assert star["dr3_source_id"] == 111 and star["found_in_gaia"]
    assert np.isclose(star["posterior"].sum(), 1.0, atol=1e-5) and np.isfinite(star["age_map_myr"])
    assert star["color"] == star["BpmRp0"]
    missing = ci.star_age_from_gaia(999, a0=0.05)
    assert not missing["found_in_gaia"] and np.isnan(missing["age_map_myr"]) and missing["posterior"].sum() == 0


def test_cli_gaia(monkeypatch, capsys, tmp_path):
    from cmd_interp.cli import main
    monkeypatch.setattr(gaia_module, "query_gaia_dr3", lambda ids, **kw: _gaia_rows().rename(columns={"SOURCE_ID": "dr3_source_id"}))
    (tmp_path / "ids.txt").write_text("# ids\n222\n")
    assert main(["gaia", "111", "--ids-file", str(tmp_path / "ids.txt"), "--a0", "0.05", "--output", str(tmp_path / "out.csv")]) == 0
    out = pd.read_csv(tmp_path / "out.csv")
    assert list(out.dr3_source_id) == [111, 222] and "cmd_age_myr" in out
    assert "111" in capsys.readouterr().out


@pytest.mark.skipif(os.environ.get("CMD_INTERP_NETWORK_TESTS") != "1", reason="set CMD_INTERP_NETWORK_TESTS=1 to query the Gaia archive")
def test_live_archive_query():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table, _ = ci.ages_from_gaia([6794047652729201024], a0=0.0)    # AU Mic (beta Pic MG)
    assert table.found_in_gaia[0] and np.isfinite(table.BpmRp0[0])
