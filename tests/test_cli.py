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
