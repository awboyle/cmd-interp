"""PARSEC metallicity offsets of the CMD sequence at fixed colour.

The empirical reference sequences are treated as solar.  For a star of known [M/H], the
PARSEC v1.2S shift of the pre-main-sequence/main-sequence locus at fixed (BP-RP)_0,

    dM_G(Z, t, c) = M_G(PARSEC, [M/H]=Z, age t, colour c) - M_G(PARSEC, [M/H]=0, t, c)

(positive = fainter), is added to every model sequence before the star is compared with it.
That is the same as moving the star to its solar-metallicity position, but the shift
depends on age, so it is applied at every trial age rather than once.

The bundled table (data/parsec_mg_offsets.csv; built by scripts/build_parsec_offsets.py)
spans [M/H] -0.5 to +0.3, 5-1000 Myr and (BP-RP)_0 0.8-2.9 and is interpolated linearly
in ([M/H], log10 age, colour).  Inputs outside these ranges are clipped to the edge (no
extrapolation); infer_ages() flags stars whose [M/H] or colour had to be clipped.
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator


def _table_path():
    return resources.files("cmd_interp") / "data" / "parsec_mg_offsets.csv"


class MetallicityOffset:
    """Callable dM_G(mh, log10_age_myr, colour); broadcasts its three arguments."""

    def __init__(self, table):
        table = table.copy()
        table["colour"] = table.colour.round(3)
        mh = np.sort(table.MH.unique())
        log_age = np.log10(np.sort(table.age_myr.unique()))
        colour = np.sort(table.colour.unique())
        cube = (
            table.assign(logAge=np.log10(table.age_myr))
            .pivot_table(index=["MH", "logAge", "colour"], values="dMG")
            .dMG.to_numpy()
            .reshape(len(mh), len(log_age), len(colour))
        )
        if not (np.isfinite(cube).all() and cube.shape == (len(mh), len(log_age), len(colour))):
            raise ValueError("The metallicity offset table must be a complete, finite ([M/H], age, colour) grid.")
        self._rgi = RegularGridInterpolator((mh, log_age, colour), cube, method="linear", bounds_error=True)
        self._lo = np.array([mh[0], log_age[0], colour[0]])
        self._hi = np.array([mh[-1], log_age[-1], colour[-1]])
        self.mh_range = (float(mh[0]), float(mh[-1]))
        self.age_range_myr = (float(10 ** log_age[0]), float(10 ** log_age[-1]))
        self.colour_range = (float(colour[0]), float(colour[-1]))

    def __call__(self, mh, log10_age_myr, colour):
        mh, log10_age_myr, colour = np.broadcast_arrays(
            np.asarray(mh, float), np.asarray(log10_age_myr, float), np.asarray(colour, float)
        )
        xi = np.stack([mh, log10_age_myr, colour], axis=-1).reshape(-1, 3)
        xi = np.clip(xi, self._lo, self._hi)
        return self._rgi(xi).reshape(mh.shape)

    def clipped(self, mh, colour):
        """True where [M/H] or the (central) colour lies outside the table and is clipped."""
        mh = np.asarray(mh, float)
        colour = np.asarray(colour, float)
        return ((mh < self.mh_range[0]) | (mh > self.mh_range[1])
                | (colour < self.colour_range[0]) | (colour > self.colour_range[1]))


@lru_cache(maxsize=1)
def parsec_offset():
    """The bundled PARSEC offset interpolator (cached)."""
    with resources.as_file(_table_path()) as path:
        return MetallicityOffset(pd.read_csv(path))
