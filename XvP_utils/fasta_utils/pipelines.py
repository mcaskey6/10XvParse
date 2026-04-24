from __future__ import annotations

from .classes import RunSettings, AnalysisConfig, PipelinePaths
from . import helpers
from . import utils

def load10X(settings: RunSettings, config_file:str, assay:str) -> None:
    '''Get H5AD files for 10X data.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = PipelinePaths.build_10x(settings, config, assay)
    logger = helpers.setup_logger(paths.log_file)
    paths.ensure_dirs(logger)
    utils.core_pipeline(settings, paths, config, assay, logger)
    
    # Pseudalign multiplexed files
    utils.pseudoalign10X(
        paths=paths,
        fastq_files=paths.multiplexed_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )

def load10xHashtags(settings: RunSettings, config_file:str, hashtags:str) -> None:
    '''Get H5AD files for 10X Hashtags for Demultiplexing'''
    config = AnalysisConfig.from_yaml(config_file, hashtags)
    paths = PipelinePaths.build_10xhashtags(settings, config, hashtags)
    logger = helpers.setup_logger(paths.log_file)
    paths.ensure_dirs(logger)
    utils.core_pipeline(settings, paths, config, hashtags, logger)
    
    utils.pseudoalign10XHashtags(
        paths=paths,
        fastq_files=paths.multiplexed_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )


def loadParse(settings: RunSettings, config_file:str, assay:str) -> None:
    '''Get H5AD files for Parse data. Includes files separated by barcode type (polyT vs. randO).'''

    # Run core pipeline shared with 10X
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = PipelinePaths.build_parse(settings, config, assay)
    logger = helpers.setup_logger(paths.log_file)
    paths.ensure_dirs(logger)
    utils.core_pipeline(settings, paths, config, assay, logger)
    
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
        utils.extractRandOPolyT(
            paths=paths, 
            threads=settings.threads, 
            logger=logger
        )

    # Pseudoalign reads separated by barcode type
    utils.pseudoalignParse(
        paths=paths,
        fastq_files=paths.filtered_files,
        kb_out_dir=paths.kb_all_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger)
    
def subsampleParse(settings: RunSettings, config_file:str, assay:str) -> None:
    '''Subsample the parse, polyT and randO FASTQ files for comparison.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = PipelinePaths.build_parse(settings, config, assay)
    logger = helpers.setup_logger(paths.log_file)
    paths.ensure_dirs(logger)

    sampled_parse_exist = all(p.is_file() for p in paths.sampled_files)
    if not sampled_parse_exist or settings.overwrite:
        # subsample parse FASTQs 
        utils.subsample_fastqs(
            fastq_files=paths.filtered_files,
            output_files=paths.sampled_files,
            num_reads=settings.subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_polyT_exist = all(p.is_file() for p in paths.sampled_polyT_files)
    if not sampled_polyT_exist or settings.overwrite:
        # subsample polyT FastQs
        utils.subsample_fastqs(
            fastq_files=paths.polyT_files,
            output_files=paths.sampled_polyT_files,
            num_reads=settings.subsample_num,
            threads=settings.threads,
            logger=logger
        )

    sampled_randO_exist = all(p.is_file() for p in paths.sampled_randO_files)
    if not sampled_randO_exist or settings.overwrite:
        # subsample randO FastQs
        utils.subsample_fastqs(
            fastq_files=paths.randO_files,
            output_files=paths.sampled_randO_files,
            num_reads=settings.subsample_num,
            threads=settings.threads,
            logger=logger
        )

    utils.pseudoalignParse(
        paths=paths,
        fastq_files=paths.sampled_files,
        kb_out_dir=paths.kb_sub_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )

    utils.pseudoalignParse(
        paths=paths,
        fastq_files=paths.sampled_polyT_files,
        kb_out_dir=paths.kb_sub_polyT_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger
    )

    utils.pseudoalignParse(
        paths=paths,
        fastq_files=paths.sampled_randO_files,
        kb_out_dir=paths.kb_sub_randO_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger
    )
    
def subsample10X(settings: RunSettings, config_file:str, assay:str) -> None:
    '''Subsample the parse, polyT and randO FASTQ files for comparison.'''
    config = AnalysisConfig.from_yaml(config_file, assay)
    paths = PipelinePaths.build_10x(settings, config, assay)
    logger = helpers.setup_logger(paths.log_file)
    paths.ensure_dirs(logger)

    sampled_10x_exist = all(p.is_file() for p in paths.sampled_files)
    if not sampled_10x_exist or settings.overwrite:
        # subsample parse FASTQs 
        utils.subsample_fastqs(
            fastq_files=paths.multiplexed_files,
            output_files=paths.sampled_files,
            num_reads=settings.subsample_num,
            threads=settings.threads,
            logger=logger
        )

    utils.pseudoalign10X(
        paths=paths,
        fastq_files=paths.sampled_files,
        kb_out_dir=paths.kb_sub_dir,
        tech=config.technology,
        threads=settings.threads,
        logger=logger,
    )