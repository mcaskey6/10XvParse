"""Plotting and analysis helpers for the 10x-vs-Parse comparison notebooks.

The read-to-matrix preprocessing that used to live under ``XvP_utils.preprocessing``
has been replaced by the Snakemake workflow in ``workflow/``; this package now holds
only the downstream plotting/analysis code used by the notebooks in ``Notebooks/``.

Order matters below: ``processing`` is imported first so the names it exports (e.g.
``init_processing``) are available to the modules that import them.
"""
from .processing import *        # noqa: F401,F403
from .basic_plots import *       # noqa: F401,F403
from .parse_plots import *       # noqa: F401,F403
from .combo_plots import *       # noqa: F401,F403
from .cross_comparison import *  # noqa: F401,F403

# Backwards compatibility for notebooks written against the old plotting subpackage
# (`from XvP_utils import plotting` and `from XvP_utils.plotting import ...`): alias
# this package as `plotting`, both as an attribute (for the former) and as a
# sys.modules entry (so the latter's submodule lookup resolves here). Remove once
# those notebook imports are updated to `import XvP_utils as plotting`.
import sys as _sys
plotting = _sys.modules[__name__]
_sys.modules[__name__ + ".plotting"] = plotting
