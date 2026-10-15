"""A_0 from the EXPLORE 3D dust maps (Vergely et al. 2022).

The extinction density cube (A_0 at 550 nm per pc) is integrated along the line of sight
from the Sun to the star in 5 pc steps.  The same is done with the density *error* cube,
which is integrated linearly, i.e. voxel errors are taken as perfectly correlated along the
sightline.  That makes it an upper bound on the dust-map uncertainty.  The local density
dA_0/dd at the star is also returned, so the distance error can be propagated into A_0.
Beyond the edge of a cube, A_0 is frozen at its edge value and the edge density is used as
the gradient.

The cubes are not shipped with the package (~200 MB and ~30 MB).  Download them from the
EXPLORE G-Tomo service (https://explore-platform.eu), then point to them with:

    ExploreDustMap(values_path, errors_path)

or set the environment variables CMD_INTERP_EXPLORE_VALUES and CMD_INTERP_EXPLORE_ERRORS,
or CMD_INTERP_DUST_DIR for a directory that holds the default file names
(explore_cube_density_values_010pc_v2.h5, explore_cube_density_errors_050pc_v2.h5).

Requires h5py and astropy: pip install "cmd-interp[dust]".
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import scipy.interpolate as spi

DEFAULT_VALUES_FILE = "explore_cube_density_values_010pc_v2.h5"
DEFAULT_ERRORS_FILE = "explore_cube_density_errors_050pc_v2.h5"
STEP_PC = 5


def _require(module):
    try:
        return __import__(module)
    except ImportError as error:  # pragma: no cover
        raise ImportError(f'The dust-map option needs {module}: pip install "cmd-interp[dust]"') from error


def load_cube(hdf5file):
    """Load an EXPLORE cube and the coordinates of its axes, with the Sun at (0, 0, 0) (pc).

    Returns (headers, cube, axes, min_axes, max_axes).  Original reader by N. Leclerc,
    G. Plum and S. Ferron.
    """
    h5py = _require("h5py")
    with h5py.File(hdf5file, "r") as hf:
        cube = hf["explore/cube_datas"][:]
        headers = {k: v for k, v in hf["explore/cube_datas"].attrs.items()}
    sun_position = headers["sun_position"]
    gridstep_values = headers["gridstep_values"]
    new_sun_position = np.append(sun_position[1:], sun_position[0])
    headers["new_sun_position"] = new_sun_position
    headers["values_unit"] = "A0(550nm)/parsec"
    min_axes = -1 * new_sun_position * gridstep_values
    max_axes = np.abs(min_axes)
    axes = tuple(np.linspace(min_axes[k], max_axes[k], cube.shape[k]) for k in range(3))
    return headers, cube, axes, min_axes, max_axes


def _galactic_unit_vectors(ra_deg, dec_deg):
    """Galactic cartesian unit vectors (x, y, z) of ICRS directions; shape (n, 3)."""
    _require("astropy")
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    sc = SkyCoord(ra=np.atleast_1d(ra_deg) * u.deg, dec=np.atleast_1d(dec_deg) * u.deg,
                  distance=np.ones(np.size(ra_deg)) * u.pc, frame="icrs")
    return sc.transform_to("galactic").represent_as("cartesian").get_xyz().value.T


def _line_of_sight(unit_xyz, target_distance, cube, axes, max_axes, step_pc=STEP_PC):
    """Integrated cube value and local density at ``target_distance`` (pc) along one sightline."""
    nonzero = unit_xyz != 0
    max_pc = np.amin(np.abs(np.asarray(max_axes)[nonzero] / unit_xyz[nonzero]))
    distances = np.arange(0, max_pc, step_pc)
    density = spi.interpn(axes, cube, distances[:, None] * unit_xyz[None, :], method="linear")
    xvalues = np.arange(0, len(density) * step_pc, step_pc)
    yvalues = np.nancumsum(density) * step_pc
    if target_distance >= np.nanmax(xvalues):
        extinction, gradient = yvalues[-1], density[-1]
    else:
        extinction = spi.interp1d(xvalues, yvalues)(target_distance)
        gradient = spi.interp1d(xvalues, density)(target_distance)
    return float(np.around(extinction, decimals=5)), float(gradient)


class ExploreDustMap:
    """The EXPLORE density and density-error cubes, loaded once and queried star by star."""

    def __init__(self, values_path=None, errors_path=None, step_pc=STEP_PC):
        dust_dir = os.environ.get("CMD_INTERP_DUST_DIR")
        if values_path is None:
            values_path = os.environ.get("CMD_INTERP_EXPLORE_VALUES") or (dust_dir and Path(dust_dir) / DEFAULT_VALUES_FILE)
        if errors_path is None:
            errors_path = os.environ.get("CMD_INTERP_EXPLORE_ERRORS") or (dust_dir and Path(dust_dir) / DEFAULT_ERRORS_FILE)
        if not values_path or not errors_path:
            raise FileNotFoundError(
                "EXPLORE cubes not found: pass ExploreDustMap(values_path, errors_path) or set CMD_INTERP_DUST_DIR "
                f"(holding {DEFAULT_VALUES_FILE} and {DEFAULT_ERRORS_FILE})"
            )
        self.values_path, self.errors_path, self.step_pc = Path(values_path), Path(errors_path), step_pc
        _, self._values, self._values_axes, _, self._values_max = load_cube(self.values_path)
        _, self._errors, self._errors_axes, _, self._errors_max = load_cube(self.errors_path)

    def extinction(self, ra_deg, dec_deg, distance_pc):
        """A_0, its dust-map error (line-of-sight integral of the error cube) and dA_0/dd (mag/pc).

        Returns three arrays.  Stars with a non-finite or non-positive distance get NaN.
        """
        ra_deg, dec_deg, distance_pc = (np.atleast_1d(np.asarray(x, dtype=float)) for x in (ra_deg, dec_deg, distance_pc))
        a0 = np.full(len(ra_deg), np.nan)
        a0_err = np.full(len(ra_deg), np.nan)
        gradient = np.full(len(ra_deg), np.nan)
        ok = np.isfinite(ra_deg) & np.isfinite(dec_deg) & np.isfinite(distance_pc) & (distance_pc > 0)
        if not ok.any():
            return a0, a0_err, gradient
        unit = _galactic_unit_vectors(ra_deg[ok], dec_deg[ok])
        for k, i in enumerate(np.flatnonzero(ok)):
            a0[i], gradient[i] = _line_of_sight(unit[k], distance_pc[i], self._values, self._values_axes, self._values_max, self.step_pc)
            a0_err[i], _ = _line_of_sight(unit[k], distance_pc[i], self._errors, self._errors_axes, self._errors_max, self.step_pc)
        return a0, a0_err, gradient
