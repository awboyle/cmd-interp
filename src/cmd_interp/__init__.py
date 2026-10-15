"""cmd-interp: stellar ages from empirical cluster sequences in the Gaia colour-magnitude diagram."""

from .core import DEFAULT_SCATTER_MODE, DEFAULT_SYSTEMATIC_MAG, infer_ages, star_age
from .grid import age_grid, node_widths, posterior_density
from .metallicity import parsec_offset
from .references import REFERENCE_AGES, bundled_reference_dir, load_references
from .table import ages_for_table

__version__ = "0.1.0"

__all__ = [
    "DEFAULT_SCATTER_MODE",
    "DEFAULT_SYSTEMATIC_MAG",
    "REFERENCE_AGES",
    "age_grid",
    "ages_for_table",
    "bundled_reference_dir",
    "infer_ages",
    "load_references",
    "node_widths",
    "parsec_offset",
    "posterior_density",
    "star_age",
]
