"""A_0 from the EXPLORE 3D dust maps (Vergely et al. 2022).

The extinction density cube (A_0 at 550 nm per pc) is integrated along the line of sight
from the Sun to the star in 5 pc steps.  The same is done with the density *error* cube,
which is integrated linearly, i.e. voxel errors are taken as perfectly correlated along the
sightline.  That makes it an upper bound on the dust-map uncertainty.  The local density
dA_0/dd at the star is also returned, so the distance error can be propagated into A_0.
Beyond the edge of a cube, A_0 is frozen at its edge value and the edge density is used as
the gradient.

The cubes are not shipped with the package (~220 MB and ~40 MB).  Download the FITS cubes
explore_cube_density_values_010pc_v2.fits and explore_cube_density_errors_050pc_v2.fits from
the CDS record of Vergely et al. (2022), J/A+A/664/A174
(https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174), then point to them with:

    ExploreDustMap(values_path, errors_path)

or set the environment variables CMD_INTERP_EXPLORE_VALUES and CMD_INTERP_EXPLORE_ERRORS,
or CMD_INTERP_DUST_DIR for a directory that holds the default file names.  The HDF5 versions
of the same cubes distributed by the EXPLORE G-Tomo service (.h5) are also accepted and need
h5py.

Requires astropy: pip install "cmd-interp[dust]".
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import numpy as np
import scipy.interpolate as spi

DEFAULT_VALUES_STEM = "explore_cube_density_values_010pc_v2"
DEFAULT_ERRORS_STEM = "explore_cube_density_errors_050pc_v2"
CUBE_SUFFIXES = (".fits", ".fits.gz", ".h5")
DEFAULT_VALUES_FILE = DEFAULT_VALUES_STEM + CUBE_SUFFIXES[0]
DEFAULT_ERRORS_FILE = DEFAULT_ERRORS_STEM + CUBE_SUFFIXES[0]
STEP_PC = 5


def _require(module, hint='pip install "cmd-interp[dust]"'):
    try:
        return importlib.import_module(module)
    except ImportError as error:  # pragma: no cover
        raise ImportError(f"The dust-map option needs {module}: {hint}") from error


def _read_fits(path):
    """The CDS FITS cube: NAXIS1/2/3 = X/Y/Z, so the array is transposed to (X, Y, Z)."""
    fits = _require("astropy.io.fits")
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header
        cube = np.ascontiguousarray(hdul[0].data.T, dtype=np.float32)
    headers = {k: header[k] for k in header}
    headers["new_sun_position"] = np.array([header["SUN_POSX"], header["SUN_POSY"], header["SUN_POSZ"]], dtype=float)
    headers["gridstep_values"] = np.full(3, float(header["STEP"]))
    return headers, cube


def _read_hdf5(path):
    """The EXPLORE G-Tomo HDF5 cube (explore/cube_datas, already (X, Y, Z))."""
    h5py = _require("h5py", hint="pip install h5py")
    with h5py.File(path, "r") as hf:
        cube = hf["explore/cube_datas"][:]
        headers = {k: v for k, v in hf["explore/cube_datas"].attrs.items()}
    sun_position = headers["sun_position"]
    headers["new_sun_position"] = np.append(sun_position[1:], sun_position[0])
    headers["gridstep_values"] = np.asarray(headers["gridstep_values"], dtype=float)
    return headers, cube


def load_cube(path):
    """Load an EXPLORE cube (FITS from CDS, or HDF5 from the G-Tomo service) and the coordinates
    of its axes, with the Sun at (0, 0, 0) (pc).

    Returns (headers, cube, axes, min_axes, max_axes).  Original HDF5 reader by N. Leclerc,
    G. Plum and S. Ferron.
    """
    path = Path(path)
    headers, cube = _read_hdf5(path) if path.suffix.lower() in (".h5", ".hdf5") else _read_fits(path)
    new_sun_position = headers["new_sun_position"]
    gridstep_values = headers["gridstep_values"]
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


def _default_cube(dust_dir, stem):
    """<dust_dir>/<stem> with the first suffix of CUBE_SUFFIXES that exists (the .fits name if none does)."""
    if not dust_dir:
        return None
    for suffix in CUBE_SUFFIXES:
        candidate = Path(dust_dir) / (stem + suffix)
        if candidate.exists():
            return candidate
    return Path(dust_dir) / (stem + CUBE_SUFFIXES[0])


class ExploreDustMap:
    """The EXPLORE density and density-error cubes, loaded once and queried star by star."""

    def __init__(self, values_path=None, errors_path=None, step_pc=STEP_PC):
        dust_dir = os.environ.get("CMD_INTERP_DUST_DIR")
        if values_path is None:
            values_path = os.environ.get("CMD_INTERP_EXPLORE_VALUES") or _default_cube(dust_dir, DEFAULT_VALUES_STEM)
        if errors_path is None:
            errors_path = os.environ.get("CMD_INTERP_EXPLORE_ERRORS") or _default_cube(dust_dir, DEFAULT_ERRORS_STEM)
        if not values_path or not errors_path:
            raise FileNotFoundError(
                "EXPLORE cubes not found: pass ExploreDustMap(values_path, errors_path) or set CMD_INTERP_DUST_DIR "
                f"to a folder holding {DEFAULT_VALUES_FILE} and {DEFAULT_ERRORS_FILE} (download them from "
                "https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/664/A174)"
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
