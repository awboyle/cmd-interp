import matplotlib

matplotlib.use("Agg")

import cmd_interp as ci  # noqa: E402
from cmd_interp.plotting import plot_star  # noqa: E402


def test_plot_star(tmp_path):
    fig, axes = plot_star(ci.star_age(2.0, 8.0, 0.01, 0.03, mh=0.1))
    fig.savefig(tmp_path / "star.png")
    assert len(axes) == 2
