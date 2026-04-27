from __future__ import annotations

import logging

from .classes import RunSettings, AnalysisConfig, TenXPaths, HashtagsPaths, ParsePaths
from . import utils


def _download_reference(settings: RunSettings, paths: TenXPaths, config: AnalysisConfig, logger: logging.Logger) -> None:
    '''Download reference genome. Builds a barnyard (dual-species) reference when config.is_barnyard.'''
    ref_exists = paths.genome_file.is_file() and paths.gtf_file.is_file()
    if not ref_exists or settings.overwrite:
        if config.is_barnyard:
            utils.make_barnyard_reference(paths, config, logger)
        else:
            utils.get_reference(paths.genome_file, paths.gtf_file, config.genome_url, config.gtf_url, logger)
    else:
        logger.info("Reference files already exist. Skipping reference download.")


def _run_core_pipeline(settings: RunSettings, paths, config: AnalysisConfig, assay: str, logger: logging.Logger) -> None:
    '''Route to ERA or SRA download pipeline based on what accessions the config provides.'''
    if config.era:
        utils.era_core_pipeline(settings, paths, config, assay, logger)
    else:
        utils.core_pipeline(settings, paths, config, assay, logger)


def load_10x(settings: RunSettings, config_file: str, assay: str, logger: logging.Logger) -> None:
    '''Get H5AD files for 10X data. Supports SRA and ERA sources, single-species and barnyard references.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = TenXPaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    _download_reference(settings, paths, config, logger)
    _run_core_pipeline(settings, paths, config, assay, logger)

    utils.pseudoalign_10x(
        paths=paths,
        fastq_files=paths.multiplexed_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )


def load_10x_hashtags(settings: RunSettings, config_file: str, hashtags: str, logger: logging.Logger) -> None:
    '''Get H5AD files for 10X Hashtags for Demultiplexing'''
    config = AnalysisConfig.from_yaml(config_file, hashtags)
    paths = HashtagsPaths.build(settings, config, hashtags)
    paths.ensure_dirs(logger)

    utils.core_pipeline(settings, paths, config, hashtags, logger)

    utils.pseudoalign_10x_hashtags(
        paths=paths,
        fastq_files=paths.multiplexed_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )


def load_parse(settings: RunSettings, config_file: str, assay: str, logger: logging.Logger) -> None:
    '''Get H5AD files for Parse data. Supports SRA and ERA sources, single-species and barnyard references.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = ParsePaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    _download_reference(settings, paths, config, logger)
    _run_core_pipeline(settings, paths, config, assay, logger)

    filter_parse_fastqs_exist = all(p.is_file() for p in paths.filtered_files)
    if not filter_parse_fastqs_exist or settings.overwrite:
        utils.filter_parse_fastqs(
            paths=paths,
            threads=settings.threads,
            logger=logger
        )

    # Split parse FASTQ into polyT and randO files
    randOpolyT_exist = all(p.is_file() for p in paths.randO_files + paths.polyT_files)
    if not randOpolyT_exist or settings.overwrite:
        utils.extract_rando_polyt(
            paths=paths,
            threads=settings.threads,
            logger=logger
        )

    # Pseudoalign reads separated by barcode type
    utils.pseudoalign_parse(
        paths=paths,
        fastq_files=paths.filtered_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger)


def subsample_parse(settings: RunSettings, config_file: str, assay: str, subsample_num: int, logger: logging.Logger) -> None:
    '''Subsample the parse, polyT and randO FASTQ files for comparison.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = ParsePaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    sampled_parse_exist = all(p.is_file() for p in paths.sampled_files)
    if not sampled_parse_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.filtered_files,
            output_files=paths.sampled_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_polyT_exist = all(p.is_file() for p in paths.sampled_polyT_files)
    if not sampled_polyT_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.polyT_files,
            output_files=paths.sampled_polyT_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_randO_exist = all(p.is_file() for p in paths.sampled_randO_files)
    if not sampled_randO_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.randO_files,
            output_files=paths.sampled_randO_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    utils.pseudoalign_parse(
        paths=paths,
        fastq_files=paths.sampled_files,
        kb_out_dir=paths.kb_sub_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )

    utils.pseudoalign_parse(
        paths=paths,
        fastq_files=paths.sampled_polyT_files,
        kb_out_dir=paths.kb_sub_polyT_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger
    )

    utils.pseudoalign_parse(
        paths=paths,
        fastq_files=paths.sampled_randO_files,
        kb_out_dir=paths.kb_sub_randO_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger
    )


def subsample_10x(settings: RunSettings, config_file: str, assay: str, subsample_num: int, logger: logging.Logger) -> None:
    '''Subsample the 10X FASTQ files for comparison.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = TenXPaths.build(settings, config, assay)
    paths.ensure_dirs(logger)

    sampled_10x_exist = all(p.is_file() for p in paths.sampled_files)
    if not sampled_10x_exist or settings.overwrite:
        utils.subsample_fastqs(
            fastq_files=paths.multiplexed_files,
            output_files=paths.sampled_files,
            num_reads=subsample_num,
            threads=settings.threads,
            logger=logger
        )

    utils.pseudoalign_10x(
        paths=paths,
        fastq_files=paths.sampled_files,
        kb_out_dir=paths.kb_sub_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )


def get_subsample_num(
    settings: RunSettings,
    config_file: str,
    ten_x_assays: list[str],
    parse_assays: list[str],
    logger: logging.Logger,
) -> int:
    '''Count reads in processed FASTQ files and return the minimum across all assays.'''
    counts: dict[str, int] = {}

    for assay in ten_x_assays:
        config = AnalysisConfig.from_yaml(config_file, assay)
        paths = TenXPaths.build(settings, config, assay)
        f = paths.multiplexed_files[0]
        n = utils._count_reads(f)
        counts[str(f)] = n
        logger.info("%s: %d reads", f, n)

    for assay in parse_assays:
        config = AnalysisConfig.from_yaml(config_file, assay)
        paths = ParsePaths.build(settings, config, assay)
        for f in [paths.filtered_files[0], paths.polyT_files[0], paths.randO_files[0]]:
            n = utils._count_reads(f)
            counts[str(f)] = n
            logger.info("%s: %d reads", f, n)

    minimum = min(counts.values())
    logger.info("Minimum reads across all assays: %d", minimum)
    return minimum
