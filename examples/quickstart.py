"""Quick start: ages for a few synthetic stars, one with a measured metallicity, and a diagnostic figure."""

import numpy as np

import cmd_interp

color = np.array([1.6, 2.0, 2.4, 2.0])
abs_g = np.array([6.2, 8.0, 8.9, 10.5])
color_err = np.full(4, 0.01)
abs_g_err = np.full(4, 0.03)
mh = np.array([0.0, 0.0, -0.2, 0.0])

result = cmd_interp.infer_ages(color, abs_g, color_err, abs_g_err, mh=mh)
for k in range(len(color)):
    print(f"star {k}: MAP {result['age_map_myr'][k]:7.1f} Myr, 68% [{result['age_myr_16'][k]:7.1f}, "
          f"{result['age_myr_84'][k]:7.1f}] Myr, limit '{result['limit'][k]}'")

try:
    from cmd_interp.plotting import plot_star
except ImportError:
    pass
else:
    fig, _ = plot_star(cmd_interp.star_age(2.4, 8.9, 0.01, 0.03, mh=-0.2))
    fig.savefig("quickstart_star.png", dpi=150, bbox_inches="tight")
    print("wrote quickstart_star.png")
