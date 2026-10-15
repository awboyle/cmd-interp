"""cmd-interp: stellar ages from empirical cluster sequences in the Gaia colour-magnitude diagram."""

from .core import (BLUE_CAUTION_COLOR, DEFAULT_SCATTER_MODE, DEFAULT_SYSTEMATIC_MAG, MAX_COLOR, BlueColorWarning,
                   infer_ages, star_age)
from .dustmap import ExploreDustMap
from .grid import age_grid, node_widths, posterior_density
from .metallicity import parsec_offset
from .photometry import deredden_gaia, extinction_coefficient, inflate_parallax_uncertainty
from .references import REFERENCE_AGES, bundled_reference_dir, load_references, resolve_reference_ages
from .table import ages_for_table

__version__ = "0.1.0"

__all__ = [
    "BLUE_CAUTION_COLOR",
    "BlueColorWarning",
    "DEFAULT_SCATTER_MODE",
    "DEFAULT_SYSTEMATIC_MAG",
    "ExploreDustMap",
    "MAX_COLOR",
    "REFERENCE_AGES",
    "age_grid",
    "ages_for_table",
    "bundled_reference_dir",
    "deredden_gaia",
    "extinction_coefficient",
    "infer_ages",
    "inflate_parallax_uncertainty",
    "load_references",
    "node_widths",
    "parsec_offset",
    "posterior_density",
    "resolve_reference_ages",
    "star_age",
]
