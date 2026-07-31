'''Get Parse and 10X H5AD files for comparison for Analysis 4.'''

from __future__ import annotations
from XvP_utils.preprocessing import load_10x, load_parse, subsample_parse, subsample_10x, get_subsample_num, setup_logger, get_genebody_plot
from XvP_utils.preprocessing import RunSettings
from pathlib import Path
from logging import Logger
import os


def load_all(settings: RunSettings, config_file: str, logger: Logger) -> None:
    '''Run the full loading pipeline for both 10X and Parse.'''
    load_10x(settings, config_file, "10x", logger)
    load_parse(settings, config_file, "parse", logger)


def subsample_all(settings: RunSettings, config_file: str, subsample_num: int, logger: Logger) -> None:
    '''Subsample both datasets for cross-technology comparison.'''
    subsample_10x(settings, config_file, "10x", subsample_num, logger)
    subsample_parse(settings, config_file, "parse", subsample_num, logger)

def genebody_plots(settings: RunSettings, config_file: str, logger: Logger) -> None:
    get_genebody_plot(
        settings=settings,
        config_file=config_file,
        tenx_assay = "10x",
        parse_assay = "parse",
        logger = logger, 
    )


if __name__ == "__main__":
    settings = RunSettings(
        root_dir=Path(__file__).parent.parent,
        config_name="analysis4.yaml",
        overwrite=True,
        run_kb=True,
        threads=16,
        max_workers=4,
    )

    config_file = settings.root_dir / "Configs" / settings.config_name
    os.makedirs(settings.root_dir / "Logs", exist_ok=True)
    logger = setup_logger(settings.root_dir / "Logs" / "analysis4.txt")

    #load_all(settings, config_file, logger)
    subsample_num = get_subsample_num(settings, config_file, ten_x_assays=["10x"], parse_assays=["parse"], logger=logger)
    subsample_all(settings, config_file, subsample_num, logger)
    genebody_plots(settings, config_file, logger)

