import json

import numpy as np
import pandas as pd

from cmd_interp.cli import main


def test_star_json(capsys, tmp_path):
    assert main(["star", "--color", "2.0", "--abs-g", "8.0", "--color-err", "0.01", "--abs-g-err", "0.03",
                 "--json", "--posterior", str(tmp_path / "p.csv")]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert 50 < summary["age_map_myr"] < 100 and summary["limit"] == ""
    posterior = np.loadtxt(tmp_path / "p.csv", delimiter=",", skiprows=1)
    assert posterior.shape[1] == 2


def test_table_with_posteriors(tmp_path):
    data = pd.DataFrame({"BpmRp0": [2.0, 1.5, np.nan], "phot_g_mean_mag_0": [8.0, 6.0, 7.0],
                         "BpmRp0_err": [0.01, 0.01, 0.01], "phot_g_mean_mag_0_err": [0.03, 0.03, 0.03]})
    data.to_csv(tmp_path / "in.csv", index=False)
    out = tmp_path / "out.parquet"
    post = tmp_path / "post.npy"
    assert main(["table", str(tmp_path / "in.csv"), str(out), "--posterior-output", str(post), "--mh", "0.1",
                 "--systematic-mag", "0.08"]) == 0
    table = pd.read_parquet(out)
    assert "cmd_age_myr" in table and table.cmd_age_mh.eq(0.1).all()
    posterior = np.load(post)
    grid = np.load(tmp_path / "post_age_grid.npy")
    index = pd.read_csv(tmp_path / "post_index.csv")
    assert posterior.shape == (2, len(grid)) and list(index.row) == [0, 1]
    np.testing.assert_allclose(posterior.sum(axis=1), 1.0, atol=1e-5)


def test_reference_age_and_color_limit_flags(capsys):
    assert main(["star", "--color", "2.0", "--abs-g", "8.0", "--color-err", "0.01", "--abs-g-err", "0.03", "--json",
                 "--reference-age", "Pleiades=112", "--reference-age", "alpha_Per=85"]) == 0
    shifted = json.loads(capsys.readouterr().out)
    assert shifted["age_map_myr"] > 74.5
    assert main(["star", "--color", "3.3", "--abs-g", "11.0", "--color-err", "0.01", "--abs-g-err", "0.03", "--json"]) == 0
    red = json.loads(capsys.readouterr().out)
    assert red["color_out_of_range"] and red["age_map_myr"] is None


def test_deredden_then_table(tmp_path):
    gaia = pd.DataFrame({"parallax": [7.4, 5.0], "parallax_error": [0.02, 0.05], "ruwe": [0.95, 1.1],
                         "phot_g_mean_mag": [12.0, 13.5], "phot_bp_mean_mag": [12.6, 14.5], "phot_rp_mean_mag": [11.3, 12.5],
                         "phot_g_mean_flux": [3e5, 7e4], "phot_g_mean_flux_error": [60, 30],
                         "phot_bp_mean_flux": [1e5, 2e4], "phot_bp_mean_flux_error": [80, 40],
                         "phot_rp_mean_flux": [2e5, 6e4], "phot_rp_mean_flux_error": [90, 50]})
    gaia.to_csv(tmp_path / "gaia.csv", index=False)
    assert main(["deredden", str(tmp_path / "gaia.csv"), str(tmp_path / "phot.csv"), "--a0", "0.1", "--a0-err", "0.02"]) == 0
    assert main(["table", str(tmp_path / "phot.csv"), str(tmp_path / "ages.csv")]) == 0
    ages = pd.read_csv(tmp_path / "ages.csv")
    assert {"BpmRp0", "BpmRp0_err", "cmd_age_myr"} <= set(ages.columns) and np.isfinite(ages.cmd_age_myr).any()
