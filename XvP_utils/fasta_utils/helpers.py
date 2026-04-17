from __future__ import annotations

from pathlib import Path
import logging
import subprocess
import requests
from .classes import LibraryFiles, PipelinePaths, AnalysisConfig

def download_file(url: str, destination: Path, logger: logging.Logger) -> None:
    '''Download a file to a destination given a url'''

    logger.debug("Downloading %s -> %s", url, destination)
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()

    with open(destination, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)


def setup_logger(log_file: Path) -> logging.Logger:
    '''Set up logger'''

    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    handler = logging.FileHandler(log_file, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    logger.addHandler(handler)

    return logger


def build_libraries(config: AnalysisConfig, paths: PipelinePaths) -> list[LibraryFiles]:
    '''Generate names for FASTA files to be dumped with fasterq-dump'''

    libraries: list[LibraryFiles] = []
    for i in range(len(config.sra)):
        lib_name = f"Lib{i}"
        libraries.append(
            LibraryFiles(
                name=lib_name,
                read1_fasta=paths.dumped_dir / f"{lib_name}_{config.r1_num}.fasta",
                read2_fasta=paths.dumped_dir / f"{lib_name}_{config.r2_num}.fasta",
            )
        )
    return libraries
def run_command(cmd: list[str], logger: logging.Logger) -> None:
    '''A wrapper for subproccess.run to set up logging'''

    logger.debug("Running command: %s", " ".join(cmd))
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    if result.stdout:
        logger.debug(result.stdout)
    if result.stderr:
        logger.error(result.stderr)