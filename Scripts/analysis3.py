'''Get Parse and 10X H5AD files for comparison for Analysis 2.'''

from __future__ import annotations
from XvP_utils.fasta_utils import load10X, loadParse, subsampleParse, subsample10X
from XvP_utils.fasta_utils import RunSettings
from pathlib import Path

def loadAll(settings: RunSettings, config_file: str) -> None:
    '''Run the pipeline to get H5AD files for both 10X and Parse data from the same analysis'''
    load10X(settings, config_file, "10x_H1")
    load10X(settings, config_file, "10x_H2")
    loadParse(settings, config_file, "parse_H1")
    loadParse(settings, config_file, "parse_H2")

def subsampleAll(settings: RunSettings, config_file: str) -> None:
    '''Run the pipeline to get H5AD files for subsampled 10X and Parse data from the same analysis'''
    subsample10X(settings, config_file, "10x_H1")
    subsample10X(settings, config_file, "10x_H2")
    subsampleParse(settings, config_file, "parse_H1")
    subsampleParse(settings, config_file, "parse_H2")

if __name__ == "__main__":
    # Specify the files to be run and set pipeline parameteres
    settings = settings = RunSettings(
        root_dir=Path("/home/mcaskey/10XvParse"),
        config_name="analysis3.yaml",
        overwrite=False,
        subsample_num=89116454,
        threads=16,
        max_workers=4
    )

    # Get the config_file path
    config_file = settings.root_dir / "Configs" / settings.config_name

    # Run the pipeline
    loadAll(settings, config_file)
    subsampleAll(settings, config_file)