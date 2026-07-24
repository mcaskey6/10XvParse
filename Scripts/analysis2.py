'''Get Parse and 10X H5AD files for comparison for Analysis 2.'''

from __future__ import annotations
from XvP_utils.preprocessing import load_10x, load_parse, load_10x_hashtags, subsample_parse, subsample_10x, get_subsample_num, setup_logger, get_genebody_plot
from XvP_utils.preprocessing import RunSettings
from pathlib import Path
from logging import Logger
import os


def load_all(settings: RunSettings, config_file: str, logger: Logger) -> None:
    '''Run the pipeline to get H5AD files for both 10X and Parse data from the same analysis'''
    load_10x(settings, config_file, "10x", logger)
    load_10x_hashtags(settings, config_file, "10x_hashtags", logger)
    load_parse(settings, config_file, "parse", logger)
    load_parse(settings, config_file, "parse_mini", logger)

def subsample_standard(settings: RunSettings, config_file: str, logger: Logger) -> None:
    '''Run the pipeline to subsample 10X and Parse data to the same number of reads and
    get H5AD files for the 'standard' set of datasets'''

    subsample_num_standard = get_subsample_num(
        settings, config_file,
        ten_x_assays=["10x"],
        parse_assays=["parse"],
        logger=logger,
    )
    subsample_10x(settings, config_file, "10x", subsample_num_standard, logger, tag="standard")
    subsample_parse(settings, config_file, "parse", subsample_num_standard, logger)

def subsample_mini(settings: RunSettings, config_file: str, logger: Logger) -> None:
    '''Run the pipeline to subsample 10X and Parse data to the same number of reads and
    get H5AD files for the 'mini' set of datasets'''

    subsample_num_mini = get_subsample_num(
        settings, config_file,
        ten_x_assays=["10x"],
        parse_assays=["parse_mini"],
        logger=logger,
    )
    subsample_10x(settings, config_file, "10x", subsample_num_mini, logger, tag="mini")
    subsample_parse(settings, config_file, "parse_mini", subsample_num_mini, logger)


def subsample_all(settings: RunSettings, config_file: str, logger: Logger) -> None:
    '''Run the pipeline to subsample 10X and Parse data to the same number of reads and
    get H5AD files'''

    subsample_standard(settings, config_file, logger)
    subsample_mini(settings, config_file, logger)

def genebody_plots(settings: RunSettings, config_file: str, logger: Logger) -> None:
    get_genebody_plot(
        settings=settings,
        config_file=config_file,
        tenx_assay = "10x",
        parse_assay = "parse",
        logger = logger, 
        tag = "standard"
    )

    get_genebody_plot(
        settings=settings,
        config_file=config_file,
        tenx_assay = "10x",
        parse_assay = "parse_mini",
        logger = logger, 
        tag = "mini"
    )


if __name__ == "__main__":
    settings = RunSettings(
        root_dir=Path(__file__).parent.parent,
        config_name="analysis2.yaml",
        overwrite=False,
        run_kb=True,
        threads=16,
        max_workers=4
    )

    config_file = settings.root_dir / "Configs" / settings.config_name
    os.makedirs(settings.root_dir / "Logs", exist_ok=True)
    logger = setup_logger(settings.root_dir / "Logs" / "analysis2.txt")

    load_all(settings, config_file, logger)
    subsample_all(settings, config_file, logger)
    genebody_plots(settings, config_file, logger)
