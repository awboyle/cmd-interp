import numpy as np

import cmd_interp as ci
from cmd_interp.grid import posterior_quantile


def test_age_grid_segments():
    grid = ci.age_grid([4.6, 24.3, 675.0])
    assert grid[0] == np.float64(4.6 - 4)                      # floor nodes at 1 Myr down towards 0
    assert {4.6, 24.3, 675.0} <= set(np.round(grid, 6))       # every reference age is a node
    assert np.all(np.diff(grid) > 0)
    assert grid[-1] == 4000.0 and np.isclose(grid[-1] - grid[-2], 5.0)


def test_posterior_quantile_uniform():
    grid = np.linspace(0.0, 100.0, 1001)
    density = np.ones((len(grid), 2))
    np.testing.assert_allclose(posterior_quantile(grid, density, 0.16), [16.0, 16.0])
    np.testing.assert_allclose(posterior_quantile(grid, density, 0.5), [50.0, 50.0])


def test_posterior_density_normalised():
    star = ci.star_age(2.0, 8.0, 0.01, 0.03)
    density = ci.posterior_density(star["posterior"], star["age_grid_myr"])[0]
    integral = np.sum(0.5 * (density[1:] + density[:-1]) * np.diff(star["age_grid_myr"]))
    assert np.isclose(integral, 1.0)
    assert np.isclose(star["posterior"].sum(), 1.0, atol=1e-5)
    assert np.isclose(ci.node_widths(star["age_grid_myr"]).sum(), star["age_grid_myr"][-1] - star["age_grid_myr"][0])
