from __future__ import annotations

from pathlib import Path
import logging
import subprocess
import urllib.request
import requests
from .classes import LibraryFiles, BasePaths, AnalysisConfig

def download_file(url: str, destination: Path, logger: logging.Logger) -> None:
    '''Download a file over HTTP/HTTPS to a destination path.'''
    logger.debug("Downloading %s -> %s", url, destination)
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()

    with open(destination, "wb") as f:
        for chunk in response.iter_content(chunk_size=65536):
            if chunk:
                f.write(chunk)


def download_ftp_file(url: str, destination: Path, logger: logging.Logger) -> None:
    '''Download a file over FTP to a destination path.'''
    logger.debug("Downloading %s -> %s", url, destination)
    with urllib.request.urlopen(url) as response, open(destination, "wb") as out_f:
        while chunk := response.read(65536):
            out_f.write(chunk)


def setup_logger(log_file: Path) -> logging.Logger:
    '''Create a named logger that writes to log_file and stdout.'''
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(log_file.stem)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    return logger


def build_era_libraries(config: AnalysisConfig, paths: BasePaths) -> list[LibraryFiles]:
    '''Generate paths for FASTQ files downloaded from ENA. Files arrive pre-gzipped.'''
    libraries: list[LibraryFiles] = []
    for i in range(len(config.era)):
        lib_name = f"Lib{i}"
        libraries.append(
            LibraryFiles(
                name=lib_name,
                read1_fasta=paths.dumped_dir / f"{lib_name}_{config.r1_num}.fastq.gz",
                read2_fasta=paths.dumped_dir / f"{lib_name}_{config.r2_num}.fastq.gz",
                is_gzipped=True,
            )
        )
    return libraries


def build_libraries(config: AnalysisConfig, paths: BasePaths) -> list[LibraryFiles]:
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


def build_local_libraries(config: AnalysisConfig, paths: BasePaths) -> list[LibraryFiles]:
    '''Discover pre-existing gzipped FASTQ files in dumped_dir. Expects Lib{i}_{read_num}.fastq.gz naming.'''
    libraries: list[LibraryFiles] = []
    i = 0
    while True:
        lib_name = f"Lib{i}"
        r1 = paths.dumped_dir / f"{lib_name}_{config.r1_num}.fastq.gz"
        if not r1.exists():
            break
        r2 = paths.dumped_dir / f"{lib_name}_{config.r2_num}.fastq.gz"
        libraries.append(LibraryFiles(name=lib_name, read1_fasta=r1, read2_fasta=r2, is_gzipped=True))
        i += 1
    return libraries


def run_command(cmd: list[str], logger: logging.Logger) -> None:
    '''A wrapper for subproccess.run to set up logging'''

    logger.debug("Running command: %s", " ".join(cmd))
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    if result.stdout:
        logger.debug(result.stdout)
    if result.stderr:
        logger.debug(result.stderr)