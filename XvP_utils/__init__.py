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
