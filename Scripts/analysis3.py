'''Get Parse and 10X H5AD files for comparison for Analysis 3.'''

from __future__ import annotations
from XvP_utils.preprocessing import load_10x, load_parse, subsample_parse, subsample_10x, get_subsample_num, setup_logger
from XvP_utils.preprocessing import RunSettings
from pathlib import Path
import os


def load_all(settings: RunSettings, config_file: str, logger) -> None:
    '''Run the pipeline to get H5AD files for both 10X and Parse data from the same analysis'''
    load_10x(settings, config_file, "10x_H1", logger)
    load_10x(settings, config_file, "10x_H2", logger)
    load_parse(settings, config_file, "parse_H1", logger)
    load_parse(settings, config_file, "parse_H2", logger)


def subsample_all(settings: RunSettings, config_file: str, subsample_num: int, logger) -> None:
    '''Run the pipeline to get H5AD files for subsampled 10X and Parse data from the same analysis'''
    subsample_10x(settings, config_file, "10x_H1", subsample_num, logger)
    subsample_10x(settings, config_file, "10x_H2", subsample_num, logger)
    subsample_parse(settings, config_file, "parse_H1", subsample_num, logger)
    subsample_parse(settings, config_file, "parse_H2", subsample_num, logger)


if __name__ == "__main__":
    settings = RunSettings(
        root_dir=Path(__file__).parent.parent,
        config_name="analysis3.yaml",
        overwrite=False,
        threads=16,
        max_workers=4
    )

    config_file = settings.root_dir / "Configs" / settings.config_name
    os.makedirs(settings.root_dir / "Logs", exist_ok=True)
    logger = setup_logger(settings.root_dir / "Logs" / "analysis3.txt")

    load_all(settings, config_file, logger)
    subsample_num = get_subsample_num(
        settings, config_file,
        ten_x_assays=["10x_H1", "10x_H2"],
        parse_assays=["parse_H1", "parse_H2"],
        logger=logger,
    )
    subsample_all(settings, config_file, subsample_num, logger)
